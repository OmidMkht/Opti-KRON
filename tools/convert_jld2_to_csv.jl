# --------------------------------------------------------------------------- #
# Convert a research .jld2 feeder into the public CSV format.
#
# Run from the repository root:
#     julia --project=. tools/convert_jld2_to_csv.jl reference/data/R100.jld2 data/R100
#
# Two things this does beyond a format change:
#
#   * Drops the source `node_phase_name` labels entirely and renames buses
#     B1..Bn in their original order, so the numbering still lines up with the
#     source for cross-checking while carrying none of its naming.
#   * Emits the Ybus fast-path format (bus.csv + ybus.csv), because these
#     feeders carry no per-line r/x -- only an assembled admittance matrix. See
#     src/io/from_ybus.jl for what that costs.
# --------------------------------------------------------------------------- #
using JLD2, SparseArrays, LinearAlgebra, Printf

const PHASE_CHARS = ('a', 'b', 'c')

function convert_case(jld2_path::AbstractString, out_dir::AbstractString; slack::Int=1)
    data = JLD2.load(jld2_path)
    Ybus = sparse(data["Ybus_3phase"])
    phase_node = Matrix{Bool}(data["phase_node"])
    B = size(phase_node, 2)
    nph = count(phase_node)

    size(Ybus, 1) == nph ||
        error("Ybus is $(size(Ybus,1)) rows but the phase mask implies $nph.")

    S = if haskey(data, "S_inj")
        Matrix{ComplexF64}(data["S_inj"])
    else
        V = Matrix{ComplexF64}(data["Voltage"])
        round.(V .* conj.(Ybus * V), digits=6)
    end

    # Rows of the slack bus carry the reference injection, not a load.
    slack_rows = node_row_range(phase_node, slack)
    S[slack_rows, :] .= 0.0 + 0.0im

    mkpath(out_dir)
    write_bus_csv(joinpath(out_dir, "bus.csv"), phase_node, slack)
    write_ybus_csv(joinpath(out_dir, "ybus.csv"), Ybus)
    write_load_csv(joinpath(out_dir, "load.csv"), S, phase_node)

    hl = get(data, "HL", nothing)
    ll = get(data, "LL", nothing)
    @printf("%-24s %4d buses, %4d node-phase rows, %3d scenarios, %5d Ybus nonzeros\n",
        basename(out_dir), B, nph, size(S, 2), nnz(Ybus))
    (hl === nothing || ll === nothing) ||
        println("    high-load scenario h", lpad(hl, 3, '0'), ", low-load h", lpad(ll, 3, '0'))
    return out_dir
end

"Rows of Ybus/S belonging to bus `i`."
function node_row_range(phase_node::AbstractMatrix{Bool}, i::Int)
    start = 1
    for k in 1:(i-1)
        start += count(view(phase_node, :, k))
    end
    return start:(start+count(view(phase_node, :, i))-1)
end

function write_bus_csv(path, phase_node, slack)
    open(path, "w") do io
        println(io, "bus_id,phases,base_kv,type")
        for i in axes(phase_node, 2)
            phases = join(PHASE_CHARS[p] for p in 1:3 if phase_node[p, i])
            # base_kv is not carried by the source data; the model is per-unit
            # throughout, so this column is informational only.
            println(io, "B$i,", phases, ",1.0,", i == slack ? "slack" : "pq")
        end
    end
end

function write_ybus_csv(path, Ybus)
    open(path, "w") do io
        println(io, "row,col,g,b")
        rows, vals = rowvals(Ybus), nonzeros(Ybus)
        for c in axes(Ybus, 2), k in nzrange(Ybus, c)
            v = vals[k]
            iszero(v) && continue
            println(io, rows[k], ",", c, ",", real(v), ",", imag(v))
        end
    end
end

function write_load_csv(path, S, phase_node)
    bus_of = Int[]
    phase_of = Char[]
    for i in axes(phase_node, 2), p in 1:3
        phase_node[p, i] || continue
        push!(bus_of, i)
        push!(phase_of, PHASE_CHARS[p])
    end

    open(path, "w") do io
        println(io, "bus_id,phase,scenario,p_pu,q_pu")
        for j in axes(S, 2)
            label = "h" * lpad(j, 3, '0')
            for r in axes(S, 1)
                iszero(S[r, j]) && continue      # absent rows draw zero on import
                println(io, "B", bus_of[r], ",", phase_of[r], ",", label, ",",
                    real(S[r, j]), ",", imag(S[r, j]))
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 2 || error("usage: julia tools/convert_jld2_to_csv.jl <in.jld2> <out_dir>")
    convert_case(ARGS[1], ARGS[2])
end
