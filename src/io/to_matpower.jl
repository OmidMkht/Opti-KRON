# --------------------------------------------------------------------------- #
# Write a reduced network back out as a MATPOWER case.
#
# A reduction is only useful to someone else if they can open it in the tool
# they already use. Opti-KRON produces a Kron-reduced admittance matrix and an
# assignment map; MATPOWER wants buses, branches and loads. This turns one into
# the other.
#
# ---- Recovering branches from an admittance matrix ------------------------ #
#
# Kron reduction returns Y_red, not a branch list, and the two are not the same
# kind of object: Y_red is what the kept buses *see*, with the eliminated ones
# folded in. But for a network without mutual coupling the map back is exact and
# elementary. Off the diagonal, an admittance matrix holds the negated branch
# admittance,
#
#     Y[i,j] = -y_ij        so    z_ij = -1 / Y_red[i,j]
#
# and whatever a row does not spend on its branches is a shunt to ground,
#
#     y_shunt_i = sum_j Y_red[i,j]      (the row sum)
#
# which is where the eliminated buses' own shunts, and the equivalent shunts
# created by eliminating them, end up. Together those reproduce Y_red exactly,
# so a power flow on the exported case and one on the reduced network are the
# same computation.
#
# Two things worth expecting in the output. Kron reduction can produce
# *negative* resistances -- an equivalent branch is not a physical conductor,
# and nothing requires it to look like one. MATPOWER solves such a case without
# complaint. And a reduced network is denser than the feeder it came from
# (Lemma 1), unless it was radialized first, which is normally what you want to
# export.
#
# ---- What cannot be exported ---------------------------------------------- #
#
# MATPOWER is a positive-sequence format, so a three-phase feeder has nowhere to
# go in it. Rather than silently collapsing phases into something that is not
# the network, this refuses. Three-phase reductions stay in the CSV format.
# --------------------------------------------------------------------------- #

"""
    write_matpower(path, net, A; scenario, name, base_mva, base_kv) -> path

Write the network reduced by assignment `A` as a MATPOWER case file.

`A` is an assignment matrix -- normally the radialized one, since exporting a
meshed equivalent is rarely what anyone wants. `scenario` picks which loading
column becomes the case's loads; MATPOWER holds one operating point per file, so
several scenarios mean several files.

The case is written in per unit on `base_mva`. The `Network` does not carry the
bases its CSVs were built with, so pass them if the exported file should be
readable in engineering units -- they are informational, and the power flow is
identical either way.

The header records which original buses each super-node represents, so the file
carries its own mapping back to the feeder it came from.
"""
function write_matpower(path::AbstractString, net::Network, A::AbstractMatrix;
    scenario::Int=1,
    name::AbstractString=_matpower_name(path),
    base_mva::Real=1.0,
    base_kv::Real=1.0,
    rtol::Real=1e-10)

    # 17 significant digits, not a prettier 10: these files are meant to be
    # reloaded and re-solved, and 10 digits loses enough on a large admittance
    # (case69 reaches 3e4) that the reassembled Ybus drifts from the reduction
    # it came from. 17 is what round-trips a Float64 exactly.

    is_three_phase(net) && error(
        "MATPOWER is a positive-sequence format and cannot represent a three-phase " *
        "feeder. Export $(net.name) through the CSV format instead.")
    1 <= scenario <= nscenarios(net) ||
        error("scenario $scenario is outside 1:$(nscenarios(net)).")

    kept = sort(super_nodes(A))
    Y_red, kept_rows = kron_reduce(net, kept)
    injections = aggregate_injections(A, net)[kept_rows, scenario]

    numbers = _bus_numbers(net, kept)
    slack_position = findfirst(==(net.slack), kept)
    slack_position === nothing && error("The slack bus is not a super-node in this assignment.")

    threshold = rtol * maximum(abs, Y_red)
    branches = _recover_branches(Y_red, numbers, threshold)
    shunts = vec(sum(Y_red, dims=2))          # whatever is not spent on branches

    open(path, "w") do io
        _write_matpower_header(io, name, net, A, kept, numbers, branches, scenario)
        println(io, "function mpc = ", name)
        println(io, "mpc.version = '2';")
        @printf(io, "mpc.baseMVA = %.17g;\n\n", base_mva)

        println(io, "%% bus data")
        println(io, "%\tbus_i\ttype\tPd\tQd\tGs\tBs\tarea\tVm\tVa\tbaseKV\tzone\tVmax\tVmin")
        println(io, "mpc.bus = [")
        for (position, bus) in enumerate(kept)
            # Injection convention flips: we carry generation-positive, MATPOWER
            # wants demand-positive.
            demand = -injections[position] * base_mva
            shunt = shunts[position] * base_mva
            @printf(io, "\t%d\t%d\t%.17g\t%.17g\t%.17g\t%.17g\t1\t1\t0\t%.17g\t1\t1.1\t0.9;\n",
                numbers[position], bus == net.slack ? 3 : 1,
                real(demand), imag(demand), real(shunt), imag(shunt), base_kv)
        end
        println(io, "];\n")

        println(io, "%% generator data")
        println(io, "%\tbus\tPg\tQg\tQmax\tQmin\tVg\tmBase\tstatus\tPmax\tPmin")
        println(io, "mpc.gen = [")
        @printf(io, "\t%d\t0\t0\t9999\t-9999\t1\t%.17g\t1\t9999\t-9999;\n",
            numbers[slack_position], base_mva)
        println(io, "];\n")

        println(io, "%% branch data")
        println(io, "%\tfbus\ttbus\tr\tx\tb\trateA\trateB\trateC\tratio\tangle\tstatus\tangmin\tangmax")
        println(io, "mpc.branch = [")
        for (from, to, impedance) in branches
            @printf(io, "\t%d\t%d\t%.17g\t%.17g\t0\t0\t0\t0\t0\t0\t1\t-360\t360;\n",
                from, to, real(impedance), imag(impedance))
        end
        println(io, "];")
    end
    return path
end

_matpower_name(path) = replace(splitext(basename(path))[1], r"[^A-Za-z0-9_]" => "_")

"""
Original bus numbers where they are integers, so the reduced case still refers
to the feeder's own labels; otherwise a fresh 1..n numbering.
"""
function _bus_numbers(net::Network, kept::AbstractVector{Int})
    parsed = [tryparse(Int, net.bus_ids[bus]) for bus in kept]
    return any(isnothing, parsed) || !allunique(parsed) ? collect(1:length(kept)) :
           Int[p for p in parsed]
end

"Series impedances between kept buses, from the off-diagonal admittances."
function _recover_branches(Y_red, numbers, threshold)
    n = length(numbers)
    branches = Tuple{Int,Int,ComplexF64}[]
    for i in 1:n, j in (i+1):n
        # Y is not exactly symmetric numerically; either direction describes the
        # same branch, so average them rather than privileging one.
        coupling = (Y_red[i, j] + Y_red[j, i]) / 2
        abs(coupling) > threshold || continue
        push!(branches, (numbers[i], numbers[j], -1 / coupling))
    end
    return branches
end

function _write_matpower_header(io, name, net, A, kept, numbers, branches, scenario)
    println(io, "%% ", name, " -- Opti-KRON reduction of ", net.name)
    println(io, "%")
    @printf(io, "%%   %d buses reduced to %d super-nodes (%.1f%% reduction), %d branches.\n",
        nnodes(net), length(kept), 100 * reduction_ratio(A), length(branches))
    @printf(io, "%%   Loads are scenario %d of %d.\n", scenario, nscenarios(net))
    println(io, "%")
    println(io, "%   Equivalent branches come from the Kron-reduced admittance matrix, so a")
    println(io, "%   negative resistance is possible and is not an error -- an equivalent")
    println(io, "%   branch stands in for eliminated network, not for a conductor.")
    println(io, "%")
    println(io, "%   Bus map -- each super-node and the original buses it represents:")
    for (position, bus) in enumerate(kept)
        represented = [net.bus_ids[j] for j in axes(A, 2) if A[bus, j] != 0]
        println(io, "%     ", numbers[position], " <- ", join(represented, " "))
    end
    println(io, "")
end
