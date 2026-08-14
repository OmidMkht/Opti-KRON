# --------------------------------------------------------------------------- #
# Radialization -- recovering a radial feeder from a meshed Kron-reduced one.
#
# Eliminating a bus of degree d makes its d neighbours mutually adjacent
# (Lemma 1 of the journal paper), so reducing a radial feeder generally yields a
# meshed one. Every mesh is a maximal clique of three or more buses, and those
# cliques are edge-disjoint -- which is what lets each one be repaired
# independently.
#
# Theorem 1: for each such clique, take the sub-tree of the *original* feeder
# spanning the clique's buses; the buses of degree >= 3 in that sub-tree are the
# minimal set whose reinsertion restores radiality.
#
# Reinserting them costs a little reduction and nothing in accuracy: the
# reinserted buses keep their own injections, so super-node voltages are
# unchanged. In exchange the reduced network regains sparsity, which in the
# paper's experiments made downstream OPF solve *faster* despite the extra buses.
# --------------------------------------------------------------------------- #

"""
    spanning_subtree(tree, nodes) -> Vector{Int}

Buses of the minimal sub-tree of `tree` connecting every bus in `nodes`
(the Steiner tree, which in a tree is just the union of the pairwise paths).

Taking paths between *consecutive* members is enough: the union of
`nodes[1]—nodes[2]`, `nodes[2]—nodes[3]`, ... already connects all of them and
never leaves the Steiner tree.
"""
function spanning_subtree(tree::RadialTree, nodes::AbstractVector{Int})
    length(nodes) <= 1 && return collect(nodes)

    members = Set{Int}(nodes)
    for k in 1:(length(nodes)-1)
        interior = interior_path(tree, nodes[k], nodes[k+1])
        interior === nothing &&
            error("Buses $(nodes[k]) and $(nodes[k+1]) are disconnected in the feeder.")
        union!(members, interior)
    end
    return sort(collect(members))
end

"""
    subtree_degrees(tree, members) -> Dict{Int,Int}

Degree of each bus *within* the sub-tree induced on `members` -- edges to buses
outside the set do not count.

The distinction matters: a bus can have degree 3 in the full feeder yet degree 2
in the sub-tree, in which case it is not critical and stays reduced. That is the
case Fig. 3 of the paper calls out with node 2.
"""
function subtree_degrees(tree::RadialTree, members::AbstractVector{Int})
    member_set = Set(members)
    degrees = Dict(b => 0 for b in members)
    for b in members
        p = tree.parents[b]
        if p != 0 && p in member_set
            degrees[b] += 1
            degrees[p] += 1
        end
    end
    return degrees
end

"""
    critical_nodes(net, A; tree, rtol) -> Vector{Int}

Buses that must be reinserted as super-nodes to make the Kron-reduced network
radial. Returns `Int[]` when the reduced network is already radial.
"""
function critical_nodes(net::Network, A::AbstractMatrix;
    tree::RadialTree=orient_radial(net), rtol::Real=1e-10)

    kept = super_nodes(A)
    length(kept) <= 2 && return Int[]

    Y_red, _ = kron_reduce(net, kept)
    adj = reduced_adjacency(Y_red, sort(kept), net; rtol=rtol)

    critical = Set{Int}()
    for clique_local in maximal_cliques(Graph(adj))
        length(clique_local) >= 3 || continue

        clique = sort(kept)[clique_local]              # local indices -> bus ids
        members = spanning_subtree(tree, clique)
        degrees = subtree_degrees(tree, members)

        for (bus, deg) in degrees
            # Only previously-reduced buses can be *re*inserted; clique members
            # are already super-nodes.
            deg >= 3 && iszero(A[bus, bus]) && push!(critical, bus)
        end
    end
    return sort(collect(critical))
end

"""
    radialize(net, A; tree, rtol) -> (A_radial, critical)

Reinsert the minimal set of buses that restores radiality.

Each reinserted bus becomes a super-node again and takes back its own injection.
Every other bus keeps its current super-node, so the voltages at the original
super-nodes are unchanged -- radialization trades reduction for structure, never
accuracy.
"""
function radialize(net::Network, A::AbstractMatrix;
    tree::RadialTree=orient_radial(net), rtol::Real=1e-10)

    critical = critical_nodes(net, A; tree=tree, rtol=rtol)
    isempty(critical) && return copy(A), critical

    A_radial = copy(A)
    for bus in critical
        # Detach from whichever super-node currently holds it, then restore it.
        for i in axes(A_radial, 1)
            A_radial[i, bus] = 0
        end
        A_radial[bus, bus] = 1
    end
    return A_radial, critical
end

"""
    is_radial(net, A; rtol) -> Bool

Whether the Kron-reduced network induced by `A` is a tree.
"""
function is_radial(net::Network, A::AbstractMatrix; rtol::Real=1e-10)
    kept = super_nodes(A)
    length(kept) <= 1 && return true

    Y_red, _ = kron_reduce(net, kept)
    adj = reduced_adjacency(Y_red, sort(kept), net; rtol=rtol)
    g = Graph(adj)
    return ne(g) == length(kept) - 1 && is_connected(g)
end
