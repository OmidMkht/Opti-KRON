# --------------------------------------------------------------------------- #
# Radial orientation and the admissible-assignment set.
#
# Everything here works at the *node* level, matching the papers: a three-phase
# node is reduced with all of its phases at once, so orientation, hop distance
# and assignment eligibility are properties of the bus, never of one phase.
#
# Two things the solvers need from this file:
#
#   orient_radial     the parent/child structure rooted at the slack, which
#                     defines what "upstream" and "downstream" mean.
#   admissible_pairs  which (super-node, reduced-node) assignments are legal
#                     under a hop limit, a direction rule, and -- three-phase
#                     only -- phase availability.
# --------------------------------------------------------------------------- #

"""
    RadialTree

The feeder oriented as a tree rooted at the slack.

- `parents[i]`  upstream bus of `i`; `0` at the slack.
- `children[i]` downstream buses of `i`.
- `from`, `to`  edges in BFS discovery order, `from[e]` upstream of `to[e]`.
- `depth[i]`    hops from the slack; `0` at the slack.
"""
struct RadialTree
    parents::Vector{Int}
    children::Vector{Vector{Int}}
    from::Vector{Int}
    to::Vector{Int}
    depth::Vector{Int}
    slack::Int
end

nedges(t::RadialTree) = length(t.from)

"""
    orient_radial(Lambda, slack) -> RadialTree
    orient_radial(net)           -> RadialTree

Breadth-first orientation of a radial feeder from the slack bus.

Errors if the graph is disconnected or carries a cycle. The check `m == n - 1` is
exact for a connected graph, and rejecting here beats silently returning an
arbitrary spanning tree -- every downstream bound assumes the true radial paths.
"""
function orient_radial(Lambda::AbstractMatrix, slack::Int)
    n = size(Lambda, 1)
    size(Lambda, 2) == n || error("Adjacency must be square; got $(size(Lambda)).")
    1 <= slack <= n || error("slack=$slack is outside 1:$n.")

    visited = falses(n)
    parents = zeros(Int, n)
    children = [Int[] for _ in 1:n]
    depth = zeros(Int, n)
    from, to = Int[], Int[]

    queue = [slack]
    visited[slack] = true
    head = 1
    while head <= length(queue)
        i = queue[head]
        head += 1
        for j in findall(!iszero, view(Lambda, i, :))
            visited[j] && continue
            visited[j] = true
            parents[j] = i
            depth[j] = depth[i] + 1
            push!(children[i], j)
            push!(from, i)
            push!(to, j)
            push!(queue, j)
        end
    end

    if !all(visited)
        unreached = findall(!, visited)
        error("$(length(unreached)) bus(es) unreachable from the slack, e.g. $(first(unreached, 5)).")
    end
    # Count edges in the *graph*, not the BFS tree: the tree has n-1 edges either
    # way, so checking it would never detect a mesh.
    graph_edges = 0
    for i in 1:n, j in (i+1):n
        iszero(Lambda[i, j]) || (graph_edges += 1)
    end
    graph_edges == n - 1 ||
        error("Feeder is not radial: $n buses and $graph_edges branches means " *
              "$(graph_edges - n + 1) extra branch(es) closing a cycle. " *
              "Opti-KRON's error bounds assume a radial feeder.")

    return RadialTree(parents, children, from, to, depth, slack)
end

orient_radial(net::Network) = orient_radial(net.Lambda, net.slack)

"Buses from `i` up to and including the slack."
function path_to_root(tree::RadialTree, i::Int)
    path = [i]
    while tree.parents[path[end]] != 0
        push!(path, tree.parents[path[end]])
    end
    return path
end

"""
    interior_path(tree, i, j) -> Union{Nothing,Vector{Int}}

Buses strictly between `i` and `j` on the unique tree path, empty when they are
adjacent. A cluster is connected only if every bus here joins it too, which is
what the assignment constraints are built on.
"""
function interior_path(tree::RadialTree, i::Int, j::Int)
    i == j && return Int[]

    up_i = path_to_root(tree, i)
    up_j = path_to_root(tree, j)
    depth_of = Dict(b => k for (k, b) in enumerate(up_j))

    meet_i = findfirst(b -> haskey(depth_of, b), up_i)
    meet_i === nothing && return nothing
    meet = up_i[meet_i]

    # i -> meet (exclusive of i), then meet -> j (exclusive of j), deduped at the meet.
    left = up_i[2:meet_i]                       # excludes i, includes meet
    right = reverse(up_j[1:(depth_of[meet]-1)]) # excludes j and meet
    interior = vcat(left, right)
    return filter(!=(j), interior)
end

"Number of edges on the tree path between `i` and `j`."
function hop_distance(tree::RadialTree, i::Int, j::Int)
    i == j && return 0
    p = interior_path(tree, i, j)
    p === nothing && return typemax(Int)
    return length(p) + 1
end

"""
    admissible_pairs(net, tree; hops, direction) -> BitMatrix

`A[i, j] = true` when bus `j` may be assigned to super-node `i`.

Three rules, all from the papers:

- **hop limit** `hops` bounds the tree distance from a reduced bus to its
  super-node, keeping clusters local.
- **direction** `:any` allows either orientation; `:downstream` only into a bus
  nearer the slack, which radiality enforcement requires -- an upstream-absorbing
  assignment can orphan a subtree.
- **phase availability** a bus may only join a super-node carrying all of its
  phases (constraint (9h) of the three-phase paper). Vacuous single-phase.

The diagonal is always true: every bus may remain a super-node.
"""
function admissible_pairs(net::Network, tree::RadialTree=orient_radial(net);
    hops::Int=25,
    direction::Symbol=:any)

    direction in (:any, :downstream) ||
        error("direction must be :any or :downstream; got :$direction.")
    hops >= 1 || error("hops must be at least 1; got $hops.")

    B = nnodes(net)
    ok = falses(B, B)
    for i in 1:B
        ok[i, i] = true
    end

    for j in 1:B, i in 1:B
        i == j && continue
        interior = interior_path(tree, i, j)
        interior === nothing && continue
        length(interior) + 1 <= hops || continue

        # :downstream -- j is absorbed into i, so i must sit nearer the slack.
        direction === :downstream && tree.depth[i] >= tree.depth[j] && continue

        _phases_available(net, i, j) || continue
        ok[i, j] = true
    end
    return ok
end

"Super-node `i` must carry every phase that bus `j` does."
function _phases_available(net::Network, i::Int, j::Int)
    for p in 1:3
        net.phases[p, j] && !net.phases[p, i] && return false
    end
    return true
end
