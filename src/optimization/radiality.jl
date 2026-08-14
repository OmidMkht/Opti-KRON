# --------------------------------------------------------------------------- #
# Radiality as a constraint, rather than a repair.
#
# `radialize` in src/core fixes a meshed reduction afterwards, by reinserting
# buses. That is the journal paper's Theorem 1 and it is minimal *for the
# assignment it is given* -- but the assignment was chosen without knowing it
# would have to be repaired, so the pair (reduce, then repair) is not jointly
# optimal. Enforcing radiality inside the MILP instead lets the solver trade a
# merge it would have to undo for a different one it can keep.
#
# ---- Where the meshing comes from ----------------------------------------- #
#
# Eliminating a bus of degree d makes its d neighbours a clique (Lemma 1). On a
# tree the only buses with degree >= 3 are branching nodes, so a reduced network
# is meshed exactly when some branching node is eliminated while two or more of
# its branches still hold super-nodes -- the fill-in then ties those branches
# together, and a cycle appears.
#
# That is the whole condition, and it is expressible in the assignment diagonal
# alone. For a branching node j with m child branches, let d[b] indicate that
# branch b keeps at least one super-node:
#
#     sum_b d[b] <= 1 + (m - 1) * A[j,j]
#
# If j survives (A[j,j] = 1) the bound is m and nothing is restricted -- a kept
# branching node introduces no fill-in. If j is eliminated the bound is 1: at
# most one of its branches may keep anything, so no fill-in edge can form.
#
# `d[b]` needs no integrality. It is squeezed from both sides -- above every
# member's A[k,k], and below by their sum -- so at an integral A it is forced to
# the right value, and leaving it continuous keeps the branching set small.
# --------------------------------------------------------------------------- #

"""
    add_radiality_constraints!(model, Akk, tree; reach) -> Int

Constrain the assignment so the Kron-reduced network is radial. Returns the
number of branching nodes that needed constraining.

`Akk` is the diagonal of the assignment matrix and `tree` the feeder's radial
orientation. Only branching nodes matter; a degree-2 bus can be eliminated
freely, since collapsing a path leaves a path.

`reach[j,k]` says bus `k` is within reach of `j` geometrically -- the hop-limit
mask, and nothing else. Pass the mask **before** any screening. The presolve
below argues "`j` is the closest bus outside branch `b` to everything inside
it", which relies on reachability shrinking with distance; screening removes
pairs for reasons unrelated to distance, so a screened mask can wrongly report a
branch as unreachable and force a super-node that was not needed.
"""
function add_radiality_constraints!(model, Akk::AbstractVector, tree::RadialTree;
    reach::Union{Nothing,AbstractMatrix}=nothing)

    descendants = _descendants(tree)
    constrained = 0

    for j in eachindex(tree.children)
        branches = tree.children[j]
        length(branches) >= 2 || continue           # not a branching node
        subtrees = [vcat(c, descendants[c]) for c in branches]

        # A branch nobody outside it can represent must keep a super-node of its
        # own, whatever the solver would prefer. Every path from outside the
        # branch into it passes through j, so j is the nearest outside candidate
        # for all of it: if even j cannot reach some bus k, no outside bus can,
        # and k's representative has to live inside the branch.
        forced = reach === nothing ? falses(length(branches)) :
                 Bool[!all(reach[j, k] for k in subtree) for subtree in subtrees]

        if count(forced) >= 2
            # Two branches that must each keep something cannot both hang off an
            # eliminated j without meshing it. So j survives, and once it does
            # it introduces no fill-in at all -- no further constraints needed.
            @constraint(model, Akk[j] == 1)
            constrained += 1
            continue
        end

        alive = AffExpr(Float64(count(forced)))     # forced branches count as 1 each
        for (b, subtree) in enumerate(subtrees)
            forced[b] && continue
            keeps_any = @variable(model, lower_bound = 0, upper_bound = 1)
            for k in subtree
                @constraint(model, keeps_any >= Akk[k])
            end
            @constraint(model, keeps_any <= sum(Akk[k] for k in subtree))
            add_to_expression!(alive, keeps_any)
        end

        @constraint(model, alive <= 1 + (length(branches) - 1) * Akk[j])
        constrained += 1
    end
    return constrained
end

"Buses below each bus in the tree, computed once for all of them."
function _descendants(tree::RadialTree)
    below = [Int[] for _ in eachindex(tree.children)]
    # Deepest first, so a bus's children are already complete when it is reached.
    for j in sortperm(tree.depth, rev=true)
        for child in tree.children[j]
            append!(below[j], child)
            append!(below[j], below[child])
        end
    end
    return below
end
