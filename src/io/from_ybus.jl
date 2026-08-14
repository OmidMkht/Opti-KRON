# --------------------------------------------------------------------------- #
# The Ybus fast path, for feeders whose per-line parameters are not recoverable.
#
# Everything downstream works identically -- both solver families, radialization,
# every metric consume Ybus. Only exporting a reduced network as a feeder model
# is lost. Prefer from_csv.jl when branch data exists.
#
#   bus.csv    same schema as the branch-based path
#   ybus.csv   row, col, g, b -- sparse triplets over node-phase rows. Repeated
#              (row, col) pairs accumulate, and both triangles must be listed:
#              the loader does not assume symmetry, since transformers and
#              regulators break it legitimately.
#   load.csv   optional, same schema as the branch-based path
# --------------------------------------------------------------------------- #

"""
    read_network_ybus(dir; name, base_mva) -> Network

Load a feeder from `bus.csv` + `ybus.csv` (+ optional `load.csv`) in `dir`.
The resulting `Network` has `branches === nothing`; see `has_branch_data`.
"""
function read_network_ybus(dir::AbstractString;
    name::AbstractString=basename(rstrip(dir, ['/', '\\'])),
    base_mva::Real=1.0)

    bus_ids, phases, slack, index_of = _read_bus_table(dir)
    blocks = node_rows(phases)
    nph = count(phases)

    y_df = _require_csv(dir, "ybus.csv")
    Ybus = _assemble_ybus_triplets(y_df, nph)

    load_path = joinpath(dir, "load.csv")
    load_df = isfile(load_path) ? CSV.read(load_path, DataFrame) : nothing
    S = _read_injections(load_df, index_of, phases, blocks, base_mva)

    Lambda = adjacency_from_ybus(Ybus, blocks)
    return validate(Network(name, bus_ids, phases, Ybus, Lambda, slack,
        S, nothing))
end

function _assemble_ybus_triplets(df::DataFrame, nph::Int)
    rows = Int.(df.row)
    cols = Int.(df.col)

    bad = findfirst(i -> !(1 <= rows[i] <= nph) || !(1 <= cols[i] <= nph), eachindex(rows))
    bad === nothing ||
        error("ybus.csv row $bad indexes ($(rows[bad]), $(cols[bad])), outside 1:$nph. " *
              "Indices are node-phase rows, not bus numbers -- a three-phase bus occupies three rows.")

    vals = ComplexF64.(df.g) .+ im .* ComplexF64.(df.b)
    return sparse(rows, cols, vals, nph, nph, +)
end

"""
    network_from_matrices(Ybus, S; phases, slack, bus_ids, name) -> Network

Build a `Network` straight from in-memory matrices -- the seam for an existing
pipeline that already holds a `Ybus` and an injection matrix. `phases` defaults
to a balanced single-phase feeder, one row per bus.
"""
function network_from_matrices(Ybus::AbstractMatrix, S::AbstractMatrix;
    phases::Union{Nothing,AbstractMatrix{Bool}}=nothing,
    slack::Int=1,
    bus_ids::Union{Nothing,Vector{String}}=nothing,
    name::AbstractString="network")

    Y = sparse(ComplexF64.(Ybus))
    nph = size(Y, 1)

    ph = if phases === nothing
        mask = falses(3, nph)
        mask[1, :] .= true          # one phase per bus: node-phase rows == buses
        mask
    else
        Matrix{Bool}(phases)
    end

    count(ph) == nph ||
        error("The phase mask implies $(count(ph)) node-phase rows but Ybus is $(nph)x$(nph).")

    B = size(ph, 2)
    ids = bus_ids === nothing ? string.(1:B) : bus_ids
    blocks = node_rows(ph)
    Sc = Matrix{ComplexF64}(S)

    return validate(Network(name, ids, ph, Y, adjacency_from_ybus(Y, blocks),
        slack, Sc, nothing))
end
