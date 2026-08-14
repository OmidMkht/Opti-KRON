# --------------------------------------------------------------------------- #
# Kron reduction and the assignment algebra around it.
#
# The assignment matrix `A` is B x B at the *node* level: `A[i, j] = 1` means
# bus j is represented by super-node i, and `A[i, i] = 1` marks i as a
# super-node. Injections and voltages live at the node-phase level, so `A` is
# expanded through the phase mask before it touches them -- this is the
# three-phase paper's `A ⊗ I₃`, generalised to buses carrying different phases.
#
# The two directions:
#
#   aggregate_injections   I_agg = (A ⊗ I₃) I     push loads to super-nodes
#   lift_voltages          V_kron = (A ⊗ I₃)ᵀ V   give reduced buses their
#                                                 super-node's voltage
# --------------------------------------------------------------------------- #

"Phases carried by bus `i`, in a/b/c order."
phases_of(net::Network, i::Int) = [PHASE_SYMBOLS[p] for p in 1:3 if net.phases[p, i]]

"""
    expand_assignment(A, net) -> SparseMatrixCSC

Lift a node-level assignment matrix to node-phase rows.

This is `A ⊗ I₃` when every bus is three-phase. When buses carry different
phase subsets it maps each phase of a reduced bus onto the *same phase* of its
super-node, which is why `admissible_pairs` refuses assignments whose
super-node lacks a phase -- without that guarantee this mapping is undefined.
"""
function expand_assignment(A::AbstractMatrix, net::Network)
    B = nnodes(net)
    size(A) == (B, B) || error("Assignment must be $B x $B; got $(size(A)).")

    blocks = node_rows(net)
    rows, cols = Int[], Int[]

    for j in 1:B, i in 1:B
        iszero(A[i, j]) && continue
        phases_i, phases_j = phases_of(net, i), phases_of(net, j)
        for (local_j, p) in enumerate(phases_j)
            local_i = findfirst(==(p), phases_i)
            local_i === nothing &&
                error("Assignment sends bus $j (phase $p) to super-node $i, which lacks that phase.")
            push!(rows, blocks[i][local_i])
            push!(cols, blocks[j][local_j])
        end
    end

    nph = nphase_rows(net)
    return sparse(rows, cols, ones(Float64, length(rows)), nph, nph)
end

"""
    aggregate_injections(A, net) -> Matrix{ComplexF64}

Injections after every reduced bus hands its load to its super-node. Reduced
buses come out at exactly zero, which is the condition making them eligible for
Kron elimination in the first place.
"""
aggregate_injections(A::AbstractMatrix, net::Network) = expand_assignment(A, net) * net.S

"""
    lift_voltages(A, V, net) -> Matrix

Give every reduced bus the voltage of its super-node, mapping the full-size
reduced-network solution back onto the original bus set so it can be compared
against the full-network reference.
"""
lift_voltages(A::AbstractMatrix, V::AbstractMatrix, net::Network) =
    transpose(expand_assignment(A, net)) * V

"""
    super_nodes(A) -> Vector{Int}

Buses kept as super-nodes, i.e. the nonzero diagonal of `A`.
"""
super_nodes(A::AbstractMatrix) = findall(!iszero, diag(A))

"Fraction of buses eliminated, in [0, 1]."
reduction_ratio(A::AbstractMatrix) = 1.0 - count(!iszero, diag(A)) / size(A, 1)

"""
    kron_reduce(net, kept) -> (Y_reduced, kept_rows)

Schur complement of the admittance matrix onto the node-phase rows of `kept`.

    Y_kron = Y_KK - Y_KR pinv(Y_RR) Y_RK

The pseudo-inverse rather than a plain inverse is what makes this safe for
three-phase feeders: buses missing a phase contribute zero rows and columns to
`Y_RR`, so it is structurally singular and a direct solve would fail. On a
non-singular block `pinv` agrees with `inv`.

`kept` is a list of *buses*; all phases of a bus are kept or eliminated together.
"""
function kron_reduce(net::Network, kept::AbstractVector{Int})
    B = nnodes(net)
    all(b -> 1 <= b <= B, kept) || error("kept contains a bus outside 1:$B.")
    allunique(kept) || error("kept contains duplicate buses.")
    net.slack in kept || error("The slack bus ($(net.slack)) cannot be eliminated.")

    blocks = node_rows(net)
    kept_sorted = sort(kept)
    reduced = setdiff(1:B, kept_sorted)

    kept_rows = reduce(vcat, (blocks[b] for b in kept_sorted); init=Int[])
    isempty(reduced) && return Matrix{ComplexF64}(net.Ybus[kept_rows, kept_rows]), kept_rows

    red_rows = reduce(vcat, (blocks[b] for b in reduced); init=Int[])

    Y_kk = Matrix{ComplexF64}(net.Ybus[kept_rows, kept_rows])
    Y_kr = Matrix{ComplexF64}(net.Ybus[kept_rows, red_rows])
    Y_rk = Matrix{ComplexF64}(net.Ybus[red_rows, kept_rows])
    Y_rr = Matrix{ComplexF64}(net.Ybus[red_rows, red_rows])

    return Y_kk - Y_kr * pinv(Y_rr) * Y_rk, kept_rows
end

"""
    reduced_adjacency(Y_reduced, kept, net; rtol) -> Matrix{Int}

Node-level adjacency of a Kron-reduced network. Expect it denser than the
original: eliminating a bus of degree d makes its neighbours a d-clique (Lemma
1), which is what `radialization.jl` undoes.

`rtol` is *relative* to the largest magnitude in `Y_reduced`. The Schur
complement returns rounding residue rather than structural zeros, and comparing
against absolute zero counts that residue as branches -- enough to make a
correctly radialized network report as meshed. Keep it relative when overriding:
an absolute threshold turns the graph nearly complete, and
[`critical_nodes`](@ref) then runs `maximal_cliques` over it, which is
exponential. The default sits mid-plateau; on R100/R300 any `rtol` from 1e-14 to
1e-8 gives exactly `n-1` edges.

Coupling is `max(|Y[a,c]|, |Y[c,a]|)`. `Y_reduced` is asymmetric in general --
transformers and regulators break symmetry legitimately -- and testing one
direction alone yields an adjacency matrix `Graphs.Graph` rejects.
"""
function reduced_adjacency(Y_reduced::AbstractMatrix, kept::AbstractVector{Int}, net::Network;
    rtol::Real=1e-10)
    kept_sorted = sort(kept)
    n = length(kept_sorted)

    lengths = [count(view(net.phases, :, b)) for b in kept_sorted]
    offsets = cumsum(vcat(0, lengths[1:end-1]))
    threshold = isempty(Y_reduced) ? 0.0 : rtol * maximum(abs, Y_reduced)

    block(a, c) = view(Y_reduced, (offsets[a]+1):(offsets[a]+lengths[a]),
        (offsets[c]+1):(offsets[c]+lengths[c]))

    adj = zeros(Int, n, n)
    for a in 1:n, c in (a+1):n
        coupled = any(v -> abs(v) > threshold, block(a, c)) ||
                  any(v -> abs(v) > threshold, block(c, a))
        coupled && (adj[a, c] = adj[c, a] = 1)
    end
    return adj
end

"""
    identity_assignment(net) -> Matrix{Float64}

The unreduced starting point: every bus is its own super-node. This is
`A_prev` at the first iteration of Algorithm 2.
"""
identity_assignment(net::Network) = Matrix{Float64}(I, nnodes(net), nnodes(net))

"""
    assign!(A, super, reduced) -> A

Move `reduced` into `super`'s cluster, carrying along anything already assigned
to `reduced`.

That second part is what makes iteration work: by the time a bus is absorbed it
may already represent a cluster of its own, and those buses have to follow it
rather than being orphaned at a bus that is no longer a super-node.
"""
function assign!(A::AbstractMatrix, super::Int, reduced::Int)
    super == reduced && error("Cannot assign bus $super to itself.")
    iszero(A[super, super]) && error("Bus $super is not a super-node.")
    iszero(A[reduced, reduced]) && error("Bus $reduced has already been reduced.")

    for j in axes(A, 2)
        if !iszero(A[reduced, j])
            A[reduced, j] = 0
            A[super, j] = 1
        end
    end
    A[reduced, reduced] = 0
    return A
end
