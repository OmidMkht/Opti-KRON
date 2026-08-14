# --------------------------------------------------------------------------- #
# The one network type every Opti-KRON solver consumes.
#
# Both solver families -- the MILP in src/optimization and the exhaustive
# search in src/search -- take a `Network` and nothing else. Adding a new input
# format means writing one more constructor here, not touching either solver.
#
# Phase handling follows the layout the GPU search already uses: rows of Ybus,
# V and S are indexed by *node-phase pairs*, not by nodes. `phases` is the
# 3 x B mask saying which of {a,b,c} each node carries, and `node_rows` maps a
# node to its block of rows. A balanced single-phase feeder is just the case
# where `phases` has one true per column -- there is no separate single-phase
# code path, and adding one would be a mistake.
# --------------------------------------------------------------------------- #

"""
    Branch

One physical line/transformer, in per-unit. Retained separately from `Ybus`
because `Ybus` alone cannot be inverted back to per-line parameters once
shunts, transformers or regulators are present -- and released reduced
networks are only useful to other researchers if they carry line parameters.

`r`, `x`, `b` are `nph x nph` blocks, where `nph` is the number of phases the
branch carries (1 for a balanced single-phase model, up to 3 otherwise).
"""
struct Branch
    from::Int
    to::Int
    phases::Vector{Symbol}          # subset of [:a, :b, :c], in order
    r::Matrix{Float64}
    x::Matrix{Float64}
    b::Matrix{Float64}              # total shunt susceptance, split half per end
end

"""
    Network

A feeder ready to reduce.

- `bus_ids`     original bus labels, in node order. Preserved verbatim so a
                reduced network can be exported against the names the source
                data used rather than against internal indices.
- `phases`      3 x B `Bool` mask over {a,b,c}.
- `Ybus`        nph x nph admittance, complex, node-phase indexed.
- `Lambda`      B x B node-level adjacency (1 if any phase couples the nodes).
- `slack`       node index of the slack bus. Its voltage is not stored: it
                belongs to an operating point, not to the network, and
                `slack_voltage` derives a balanced default from `phases`.
- `S`           nph x nscenarios complex injections. One column per loading
                scenario; the papers reduce against 2-3 representative columns
                and validate against the rest.
- `branches`    per-line parameters, or `nothing` when the network came in
                through the Ybus fast path. When `nothing`, feeder export and
                the `:Ladder` MILP form are unavailable; everything else works.
"""
struct Network
    name::String
    bus_ids::Vector{String}
    phases::Matrix{Bool}
    Ybus::SparseMatrixCSC{ComplexF64,Int}
    Lambda::SparseMatrixCSC{Int,Int}
    slack::Int
    S::Matrix{ComplexF64}
    branches::Union{Nothing,Vector{Branch}}
end

"Number of buses (nodes), independent of how many phases each carries."
nnodes(net::Network) = size(net.phases, 2)

"Number of node-phase pairs -- the row dimension of `Ybus`, `S` and voltages."
nphase_rows(net::Network) = count(net.phases)

"Number of loading scenarios carried by the network."
nscenarios(net::Network) = size(net.S, 2)

"True when any bus carries more than one phase."
is_three_phase(net::Network) = any(>(1), vec(sum(net.phases, dims=1)))

"True when per-line parameters survived the import, so the reduced network can be exported as a feeder."
has_branch_data(net::Network) = net.branches !== nothing

"""
    node_rows(phases) -> Vector{Vector{Int}}

Rows of `Ybus` / `S` / `V` belonging to each node. Nodes are laid out in
order, each contributing one row per phase it carries.
"""
function node_rows(phases::AbstractMatrix{Bool})
    B = size(phases, 2)
    blocks = Vector{Vector{Int}}(undef, B)
    start = 1
    for i in 1:B
        len = count(view(phases, :, i))
        blocks[i] = collect(start:(start+len-1))
        start += len
    end
    return blocks
end

node_rows(net::Network) = node_rows(net.phases)

"""
    adjacency_from_ybus(Ybus, blocks) -> Lambda

Node-level adjacency from admittance sparsity: nodes i and j are adjacent when
their block of `Ybus` holds any nonzero. Used when the source data ships no
explicit adjacency -- which is the common case for Ybus-only inputs.

Note this recovers *topology*, not orientation; the radial orientation used by
radiality constraints is derived separately from the slack bus.
"""
function adjacency_from_ybus(Ybus::SparseMatrixCSC, blocks::Vector{Vector{Int}})
    B = length(blocks)
    row_node = zeros(Int, size(Ybus, 1))
    for i in 1:B, r in blocks[i]
        row_node[r] = i
    end

    I_idx, J_idx = Int[], Int[]
    rows, vals = rowvals(Ybus), nonzeros(Ybus)
    for c in 1:size(Ybus, 2), k in nzrange(Ybus, c)
        iszero(vals[k]) && continue
        i, j = row_node[rows[k]], row_node[c]
        (i == 0 || j == 0 || i == j) && continue
        push!(I_idx, i); push!(J_idx, j)
        push!(I_idx, j); push!(J_idx, i)
    end
    return sparse(I_idx, J_idx, ones(Int, length(I_idx)), B, B, (a, b) -> 1)
end

"""
    validate(net) -> Network

Fail loudly at import time rather than deep inside a solve. Checks the
invariants every solver assumes: consistent dimensions, a slack bus in range,
a connected graph, and no all-zero rows in `Ybus`.

An all-zero Ybus row means an ungrounded Laplacian, which makes the network
singular and every downstream solve meaningless -- worth catching here, since
the failure otherwise surfaces as an unhelpful factorization error.
"""
function validate(net::Network)
    B, nph = nnodes(net), nphase_rows(net)

    size(net.Ybus, 1) == size(net.Ybus, 2) == nph ||
        error("Ybus is $(size(net.Ybus)) but the phase mask implies $nph node-phase rows.")
    size(net.S, 1) == nph ||
        error("S has $(size(net.S, 1)) rows but the phase mask implies $nph.")
    length(net.bus_ids) == B ||
        error("Got $(length(net.bus_ids)) bus ids for $B buses.")
    1 <= net.slack <= B ||
        error("slack=$(net.slack) is outside 1:$B.")
    nscenarios(net) >= 1 ||
        error("Network carries no loading scenarios; S must have at least one column.")

    empties = [i for i in 1:nph if iszero(nnz(net.Ybus[i, :]))]
    isempty(empties) ||
        error("Ybus has all-zero row(s) at $(first(empties, 5))$(length(empties) > 5 ? " (+$(length(empties)-5) more)" : ""). " *
              "This is an ungrounded Laplacian and is singular -- check the source data's grounding/slack handling.")

    _assert_connected(net.Lambda, net.slack)
    return net
end

"Breadth-first reachability from the slack; a disconnected feeder cannot be reduced coherently."
function _assert_connected(Lambda::SparseMatrixCSC, slack::Int)
    B = size(Lambda, 1)
    seen = falses(B)
    seen[slack] = true
    queue = [slack]
    while !isempty(queue)
        i = pop!(queue)
        for j in findall(!iszero, view(Lambda, :, i))
            seen[j] && continue
            seen[j] = true
            push!(queue, j)
        end
    end
    all(seen) && return
    missing_nodes = findall(!, seen)
    error("$(length(missing_nodes)) bus(es) unreachable from the slack, e.g. $(first(missing_nodes, 5)). " *
          "Opti-KRON assumes one connected feeder.")
end

function Base.show(io::IO, net::Network)
    kind = is_three_phase(net) ? "three-phase" : "single-phase"
    print(io, "Network(\"", net.name, "\", ", nnodes(net), " buses, ", nphase_rows(net),
        " node-phase rows, ", kind, ", ", nscenarios(net), " scenario(s)",
        has_branch_data(net) ? "" : ", Ybus-only", ")")
end
