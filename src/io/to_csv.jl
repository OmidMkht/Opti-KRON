# --------------------------------------------------------------------------- #
# Writing a reduced network back out in the schema it was read in.
#
# `write_matpower` refuses three-phase feeders, because a MATPOWER case has no
# way to say "this branch carries a and c but not b". That leaves the four
# published three-phase feeders and R100/R300 with nowhere to go, so reduced
# versions of those are written in the same CSV form the cases themselves use:
# `bus.csv` + `ybus.csv` + `load.csv` + `voltage.csv`, which `load_case` reads
# straight back through the Ybus fast path.
#
# `assignment.csv` is the extra one, and the only file here without an input
# counterpart. It is the reduction map -- which original bus each super-node
# stands for -- and it is what makes the reduced network interpretable rather
# than merely smaller.
#
# Everything is written at `%.17g`, full Float64 precision. Reduced admittances
# span orders of magnitude and a Kron reduction is not something anyone should
# have to re-run to recover exactly; see the precision note in to_matpower.jl.
# --------------------------------------------------------------------------- #

"""
    write_reduced_csv(dir, net, A; V=nothing, scenarios=1:nscenarios(net)) -> Vector{String}

Write the Kron reduction induced by assignment `A` into `dir`, in the CSV schema
[`read_network_ybus`](@ref) reads. Returns the paths written.

Five files:

- `bus.csv`        the kept buses, under their original ids and phases.
- `ybus.csv`       the reduced admittance, sparse triplets over the kept
                   node-phase rows -- the Schur complement, not a resampling.
- `load.csv`       injections after every eliminated bus hands its load to its
                   super-node.
- `voltage.csv`    the *full* network's voltage at the kept buses, which is the
                   operating point the reduction was certified against and what
                   the reduced network should reproduce. Omitted when `V` is not
                   given.
- `assignment.csv` the reduction map, one row per original bus.

Works for any feeder, single- or three-phase. Single-phase MATPOWER cases can
also go out through [`write_matpower`](@ref), which carries line parameters this
form does not.
"""
function write_reduced_csv(dir::AbstractString, net::Network, A::AbstractMatrix;
    V::Union{Nothing,AbstractMatrix}=nothing, scenarios=1:nscenarios(net))

    kept = sort(super_nodes(A))
    isempty(kept) && error("The assignment keeps no buses.")
    selected = collect(scenarios)
    all(s -> 1 <= s <= nscenarios(net), selected) ||
        error("`scenarios` indexes outside 1:$(nscenarios(net)).")

    Y_reduced, kept_rows = kron_reduce(net, kept)
    S_reduced = aggregate_injections(A, net)[kept_rows, selected]

    mkpath(dir)
    written = String[]

    push!(written, _write_reduced_bus(joinpath(dir, "bus.csv"), net, kept))
    push!(written, _write_reduced_ybus(joinpath(dir, "ybus.csv"), Y_reduced))
    push!(written, _write_reduced_injection(joinpath(dir, "load.csv"), net, kept,
        S_reduced, selected, ("p_pu", "q_pu")))
    if V !== nothing
        size(V, 1) == nphase_rows(net) ||
            error("V has $(size(V, 1)) rows but the network implies $(nphase_rows(net)).")
        push!(written, _write_reduced_injection(joinpath(dir, "voltage.csv"), net, kept,
            V[kept_rows, selected], selected, ("v_re_pu", "v_im_pu")))
    end
    push!(written, _write_assignment(joinpath(dir, "assignment.csv"), net, A))

    return written
end

function _write_reduced_bus(path::AbstractString, net::Network, kept::AbstractVector{Int})
    open(path, "w") do io
        println(io, "bus_id,phases,base_kv,type")
        for b in kept
            phases = join(String(p) for p in phases_of(net, b))
            println(io, net.bus_ids[b], ",", phases, ",1,", b == net.slack ? "slack" : "pq")
        end
    end
    return path
end

function _write_reduced_ybus(path::AbstractString, Y::AbstractMatrix)
    open(path, "w") do io
        println(io, "row,col,g,b")
        # The Schur complement is dense in general -- eliminating a bus makes its
        # neighbours a clique -- so this is written entry by entry rather than
        # from a sparsity pattern. Exact zeros are skipped; near-zeros are not,
        # because deciding what counts as structural here belongs to
        # `reduced_adjacency` and its relative tolerance, not to a file writer.
        for col in axes(Y, 2), row in axes(Y, 1)
            value = Y[row, col]
            iszero(value) && continue
            @printf(io, "%d,%d,%.17g,%.17g\n", row, col, real(value), imag(value))
        end
    end
    return path
end

"Shared writer for the two per-node-phase-per-scenario tables, which differ only in their value columns."
function _write_reduced_injection(path::AbstractString, net::Network,
    kept::AbstractVector{Int}, values::AbstractMatrix, selected::AbstractVector{Int},
    columns::Tuple{String,String})

    open(path, "w") do io
        println(io, "bus_id,phase,scenario,", columns[1], ",", columns[2])
        row = 0
        for b in kept, p in phases_of(net, b)
            row += 1
            for (k, s) in enumerate(selected)
                @printf(io, "%s,%s,h%03d,%.17g,%.17g\n", net.bus_ids[b], String(p), s,
                    real(values[row, k]), imag(values[row, k]))
            end
        end
    end
    return path
end

function _write_assignment(path::AbstractString, net::Network, A::AbstractMatrix)
    open(path, "w") do io
        println(io, "bus_id,super_node,kept")
        for j in 1:nnodes(net)
            i = findfirst(!iszero, view(A, :, j))
            i === nothing && error("Bus $(net.bus_ids[j]) is represented by nothing.")
            println(io, net.bus_ids[j], ",", net.bus_ids[i], ",", i == j ? 1 : 0)
        end
    end
    return path
end
