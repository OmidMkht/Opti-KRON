# --------------------------------------------------------------------------- #
# Zero-injection warm start.
#
# Absorbing a bus that injects nothing perturbs nothing: C = conj(S ./ V) = 0,
# so (A - I)C is untouched and e = Z (A - I) C is *exactly* zero. What remains
# is only the nominal gap between the two voltages, which does not depend on the
# assignment -- so the error model collapses to one linear constraint per pair,
#
#     c_ij * A[i,j] <= Ē,   c_ij = max over phases/scenarios of | |V_i| - |V_j| |
#
# with every structural constraint unchanged. That little MILP needs no error
# machinery, no indicators and no dense Z, and returns a genuinely feasible
# point of the full problem.
#
# Sequencing rule: the pairs used here must be a subset of those the *target*
# solve has binaries for. Pick a representative from outside that set and the
# target has no variable to hold it, so the start reads back all-zero and is
# rejected as infeasible for a reason that looks nothing like its cause. Hence
# `admissible` is an argument, never derived -- `solve_milp` passes the screened
# mask, so the rule holds by construction.
# --------------------------------------------------------------------------- #

"""
    zero_injection_buses(net; scenarios, tol) -> Vector{Int}

Buses drawing no net power in any of `scenarios`, excluding the slack.

Judged on the injection matrix, which already folds in shunts -- so a bus
carrying only a shunt is correctly *not* zero-injection. All phases of a bus
must be quiet for it to count, since buses are reduced whole.
"""
function zero_injection_buses(net::Network;
    scenarios=1:nscenarios(net), tol::Real=0.0)

    blocks = node_rows(net)
    quiet = Int[]
    for bus in 1:nnodes(net)
        bus == net.slack && continue
        all(abs(net.S[r, s]) <= tol for r in blocks[bus], s in scenarios) &&
            push!(quiet, bus)
    end
    return quiet
end

"""
    zero_injection_warmstart(net, V, Ē, admissible, tree, sel; kwargs...) -> Matrix{Float64}

A feasible assignment that merges only zero-injection buses, for use as a MIP
start.

`admissible` must be the pair set the target solve will use -- see the header.
Returns the identity assignment if no zero-injection bus can be merged, which is
feasible too, just less useful.
"""
function zero_injection_warmstart(net::Network, V::AbstractMatrix, Ē::Real,
    admissible::BitMatrix, tree::RadialTree, sel::AbstractVector{Int};
    max_reduction::Real=1.0,
    enforce_radiality::Bool=false,
    reach::Union{Nothing,BitMatrix}=nothing,
    tol::Real=0.0,
    prefer::Symbol=:auto,
    time_limit::Union{Nothing,Real}=nothing)

    B = nnodes(net)
    quiet = zero_injection_buses(net; scenarios=sel, tol=tol)
    isempty(quiet) && return identity_assignment(net)

    # "Merging a zero-injection bus is allowed, nothing else." Contiguity then
    # rules out reaching past a loaded bus on its own, with no extra work: the
    # interior pair was already removed here, so the path constraint kills it.
    allowed = copy(admissible)
    is_quiet = falses(B)
    is_quiet[quiet] .= true
    for ij in findall(allowed)
        i, j = Tuple(ij)
        i == j && continue
        is_quiet[j] || (allowed[ij] = false)
    end
    paths = _cluster_paths!(allowed, tree)

    absV = abs.(V)
    row_of = _phase_rows(net)
    pairs = findall(allowed)

    factory, _ = select_optimizer(prefer=prefer, verbose=false,
        time_limit=time_limit, mip_gap=1e-3)
    model = Model(factory)

    @variable(model, a[1:length(pairs)], Bin)
    A = Matrix{AffExpr}(undef, B, B)
    A .= AffExpr(0.0)
    for (k, ij) in enumerate(pairs)
        A[ij] = a[k]
    end

    @objective(model, Min, sum(A[i, i] for i in 1:B))
    _assignment_constraints!(model, A, allowed, paths, net.slack, B, max_reduction)
    enforce_radiality &&
        add_radiality_constraints!(model, [A[i, i] for i in 1:B], tree; reach=reach)

    # The whole error model, for this subproblem: a constant against the budget.
    for (k, ij) in enumerate(pairs)
        i, j = Tuple(ij)
        i == j && continue
        @constraint(model, _voltage_gap(net, absV, row_of, i, j, sel) * a[k] <= Ē)
    end

    optimize!(model)
    has_values(model) || return identity_assignment(net)

    start = zeros(Float64, B, B)
    for (k, ij) in enumerate(pairs)
        start[ij] = round(value(a[k]))
    end
    return start
end

"Worst gap between two buses' voltage magnitudes, over shared phases and scenarios."
function _voltage_gap(net::Network, absV, row_of, i::Int, j::Int, sel)
    worst = 0.0
    for p in 1:3
        net.phases[p, j] || continue
        ri, rj = row_of[p, i], row_of[p, j]
        for s in sel
            worst = max(worst, abs(absV[ri, s] - absV[rj, s]))
        end
    end
    return worst
end
