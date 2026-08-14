# --------------------------------------------------------------------------- #
# The Opti-KRON MILP -- one phase-aware formulation.
#
#   minimise   sum_i A[i,i]                    (super-nodes kept)
#   subject to sum_i A[i,j] = 1     for all j  (every bus represented once)
#              A[i,j] <= A[i,i]                (only super-nodes represent)
#              A[i,j] <= A[i,k]                (clusters are connected)
#              A[i,j] = 1  =>  annulus_ij      (merge respects the budget)
#
# A balanced feeder is the degenerate case where every bus carries one phase.
#
# Error model: moving bus j's injection to super-node i perturbs every voltage by
# e = Z (A - I) C, with C = conj(S ./ V) fixed at the operating point. Linear in
# A, which is what makes this a MILP rather than a bilinear program. Z is the
# slack-referenced bus impedance, because the slack holds its voltage and absorbs
# whatever the merge pushes at it; plain inv(Ybus) needs a shunt to ground to
# exist at all, and differs by a uniform shift when it does.
#
# The annulus | |V_i + e_i| - |V_j| | <= E is on voltage *magnitude*, so it
# excludes a disk and is nonconvex. Both linearisations bound the aligned error
# Re(conj(V_i) e_i): `:R` against a first-order bound on |V_i + e_i|, `:S`
# against the secant of its square. Neither is robust to simultaneous merges, so
# `annulus_violation` re-checks the answer against the nonconvex constraint.
#
# On/off is by JuMP indicator constraints rather than big-M: no constant to size
# wrong, and nothing dilutes the LP relaxation.
# --------------------------------------------------------------------------- #

"""
    MilpSolution

What the MILP returned, and what it cost.

- `A`          node-level assignment, `B x B`. `A[i,j] = 1` means bus `j` is
               represented by super-node `i`; `A[i,i] = 1` marks a super-node.
- `status`     JuMP termination status. Check it: a time-limited run can return
               a feasible but non-optimal assignment.
- `scenarios`  the loading scenarios the budget was enforced over -- the
               assignment is certified against these only.
- `screening`  what preprocessing removed, or `nothing` when `screen=false`.
"""
struct MilpSolution <: ReductionSolution
    A::Matrix{Float64}
    kept::Vector{Int}
    reduction::Float64
    objective::Float64
    bound::Float64
    gap::Float64
    status::MOI.TerminationStatusCode
    solver::Symbol
    scenarios::Vector{Int}
    nbinaries::Int
    nindicators::Int
    build_time::Float64
    solve_time::Float64
    screening::Union{Nothing,ScreeningReport}
end

function Base.show(io::IO, s::MilpSolution)
    print(io, "MilpSolution(", length(s.kept), " super-nodes, ",
        round(100 * s.reduction, digits=1), "% reduction, ", s.status,
        _gap_text(s.gap), ", $(s.solver) in $(round(s.solve_time, digits=1))s)")
end

"""
Render a MIP gap for people. A solver that stopped before finding any bound
reports a sentinel around 1e97, where "unknown" is the honest answer.
"""
function _gap_text(gap)
    iszero(gap) && return ""
    isfinite(gap) && gap < 1e6 || return ", gap unknown"
    return ", gap $(round(100 * gap, digits=2))%"
end

"""
    solve_milp(net, V, Ē; kwargs...) -> MilpSolution

Reduce `net` as far as possible while holding every merged bus within `Ē` per
unit of its original voltage magnitude.

`V` is the power-flow solution of the *unreduced* feeder, the operating point the
constant-current model linearises around. It is an argument because it is a
modelling choice -- delta connections, ZIP loads and taps belong in a tool that
models them. Check a `V` from elsewhere with [`powerflow_residual`](@ref): one
inconsistent with `Ybus` and `S` invalidates the bound while still solving.

- `scenarios`      columns of `S`/`V` to enforce over. Defaults to all, which
                   scales the model linearly; the papers use 2-3 representative
                   scenarios and validate against the rest.
- `hops`           hop limit between a reduced bus and its super-node.
- `direction`      `:any`, or `:downstream` to absorb only towards the slack.
- `approach`       `:R` or `:S`, the two annulus linearisations (see header).
- `max_reduction`  refuse to eliminate more than this fraction of buses.
- `pin`            buses that must survive. Pin both terminals of a device and
                   its admittance survives exactly -- see [`read_devices`](@ref).
- `warm_start`     `:zero_injection` (default) first merges only buses that
                   inject nothing, perturbing no voltage at all. `:identity` or
                   an assignment matrix also work.
- `prefer`, `verbose`, `time_limit`, `mip_gap` go to [`select_optimizer`](@ref).

The result is *not* radialised -- Kron reduction meshes a radial feeder. Follow
with [`radialize`](@ref).
"""
function solve_milp(net::Network, V::AbstractMatrix, Ē::Real;
    scenarios=1:nscenarios(net),
    hops::Int=25,
    direction::Symbol=:any,
    approach::Symbol=:R,
    max_reduction::Real=1.0,
    pin::AbstractVector{Int}=Int[],
    screen::Bool=true,
    enforce_radiality::Bool=false,
    tree::RadialTree=orient_radial(net),
    Z::Union{Nothing,AbstractMatrix}=nothing,
    warm_start::Union{Symbol,AbstractMatrix}=:zero_injection,
    prefer::Symbol=:auto,
    verbose::Bool=false,
    time_limit::Union{Nothing,Real}=nothing,
    mip_gap::Union{Nothing,Real}=1e-3)

    approach in (:R, :S) || error("approach must be :R or :S; got :$approach.")
    Ē > 0 || error("Ē must be positive; got $Ē.")
    0 <= max_reduction <= 1 || error("max_reduction must be in [0, 1]; got $max_reduction.")

    B, nph = nnodes(net), nphase_rows(net)
    size(V) == (nph, nscenarios(net)) ||
        error("V is $(size(V)) but the network implies $((nph, nscenarios(net))). " *
              "V must be the full-network power-flow solution, one column per scenario.")

    sel = collect(scenarios)
    isempty(sel) && error("`scenarios` selected no columns; the error budget would be vacuous.")
    all(s -> 1 <= s <= nscenarios(net), sel) ||
        error("`scenarios` indexes outside 1:$(nscenarios(net)).")

    t_build = time()

    # ---- geometry ---------------------------------------------------------- #
    admissible = admissible_pairs(net, tree; hops=hops, direction=direction)
    # Pinning is geometry, not an extra constraint. Applying it here rather than
    # afterwards means screening, the warm start and the binary count all see a
    # smaller problem instead of a same-sized one with a fixed variable.
    _apply_pins!(admissible, pin, B)
    # Radiality reasons about who can reach whom geometrically, so it needs the
    # hop mask before screening removes pairs for unrelated reasons.
    reach = enforce_radiality ? copy(admissible) : nothing
    paths = _cluster_paths!(admissible, tree)

    # ---- operating point --------------------------------------------------- #
    # Constant current: a bus keeps the current it drew in the base case, so
    # moving the bus moves that current intact.
    C = conj.(Matrix{ComplexF64}(net.S) ./ Matrix{ComplexF64}(V))
    Zfull = Z === nothing ? bus_impedance(net) : Matrix{ComplexF64}(Z)
    size(Zfull) == (nph, nph) || error("Z must be $nph x $nph; got $(size(Zfull)).")
    absV = abs.(V)

    # ---- screening --------------------------------------------------------- #
    # Kills pairs that cannot meet the budget under any completion, and marks
    # which annulus rows can bind. Everything it removes was provably irrelevant.
    needed, report = if screen
        result, stats = screen!(admissible, net, tree, V, Ē, Zfull, C, sel)
        paths = _cluster_paths!(admissible, tree)
        result, stats
    else
        nothing, nothing
    end
    pairs = findall(admissible)

    # ---- model ------------------------------------------------------------- #
    factory, solver_name = select_optimizer(prefer=prefer, verbose=verbose,
        time_limit=time_limit, mip_gap=mip_gap)
    model = Model(factory)

    # Resolved *after* screening, so the start can only pick pairs this model has
    # a binary for -- see warmstart.jl for what goes wrong otherwise.
    start_map = _resolve_warm_start(warm_start, net, V, Ē, admissible, tree, sel,
        max_reduction, enforce_radiality, reach, prefer, time_limit)

    @variable(model, a[1:length(pairs)], Bin)
    A = Matrix{AffExpr}(undef, B, B)
    A .= AffExpr(0.0)
    for (k, ij) in enumerate(pairs)
        A[ij] = a[k]
        set_start_value(a[k], start_map[ij] > 0.5 ? 1.0 : 0.0)
    end

    @objective(model, Min, sum(A[i, i] for i in 1:B))

    _assignment_constraints!(model, A, admissible, paths, net.slack, B, max_reduction)
    enforce_radiality &&
        add_radiality_constraints!(model, [A[i, i] for i in 1:B], tree; reach=reach)
    nind = _annulus_constraints!(model, a, pairs, net, C, Zfull, V, absV, sel, Ē,
        approach, needed)

    build_time = time() - t_build

    optimize!(model)
    status = termination_status(model)
    has_values(model) ||
        error("The MILP finished with no feasible assignment (status: $status). " *
              "Ē=$Ē may be too tight for hops=$hops on this feeder.")

    A_val = zeros(Float64, B, B)
    for (k, ij) in enumerate(pairs)
        A_val[ij] = round(value(a[k]))
    end

    return MilpSolution(A_val, super_nodes(A_val), reduction_ratio(A_val),
        objective_value(model), objective_bound(model), relative_gap(model),
        status, solver_name, sel, length(pairs), nind, build_time, solve_time(model),
        report)
end

"""
    bus_impedance(net) -> Matrix{ComplexF64}

The slack-referenced bus impedance: `inv(Ybus)` on the non-slack node-phase rows,
zero on the slack's.

Those zeros are the model, not a truncation: the slack holds its voltage against
any current asked of it, so a perturbation there propagates nowhere (zero column)
and the slack's own voltage never moves (zero row).

Inverting the full `Ybus` instead only works when the feeder carries a shunt to
ground. Without one it is a singular Laplacian, and `inv` quietly returns entries
around 1e12 that cancel to noise in the constraint rows rather than failing --
exactly what feeders assembled from branch data do.
"""
function bus_impedance(net::Network)
    nph = nphase_rows(net)
    non_slack = setdiff(1:nph, node_rows(net)[net.slack])
    Y = Matrix{ComplexF64}(net.Ybus[non_slack, non_slack])

    Z = zeros(ComplexF64, nph, nph)
    Z[non_slack, non_slack] = inv(Y)

    # inv() on a numerically singular block returns huge entries rather than
    # throwing, and the damage surfaces far downstream. One matvec catches it.
    probe = randn(ComplexF64, length(non_slack))
    residual = norm(Y * (Z[non_slack, non_slack] * probe) - probe) / norm(probe)
    residual < 1e-6 ||
        error("The admittance matrix is numerically singular (relative residual " *
              "$(round(residual, sigdigits=3)) on inversion). Check the feeder for " *
              "an isolated island or a branch with near-zero impedance.")
    return Z
end

"""
Interior paths of every admissible pair, dropping pairs that could only form a
disconnected cluster.

A cluster must be connected, so every bus `k` strictly between `i` and `j` must
join `i` too. On an unbalanced feeder `k` may be *ineligible*: two phase-`a`
laterals off a three-phase backbone bus are admissible to each other while the
bus between them lacks phases `i` has. Skipping the constraint would permit a
cluster in two pieces, so the pair is dropped -- on R100 at `hops=5`, 70 of 1198
off-diagonal pairs.

One pass suffices: paths in a tree compose, so if `(i,k)` is dropped for some `m`
on its interior path, `m` lies on `i`-to-`j` too.
"""
function _cluster_paths!(admissible::BitMatrix, tree::RadialTree)
    _prune_disconnected!(admissible, tree)

    paths = Dict{CartesianIndex{2},Vector{Int}}()
    for ij in findall(admissible)
        i, j = Tuple(ij)
        i == j && continue
        interior = interior_path(tree, i, j)
        isempty(interior) || (paths[ij] = interior)
    end
    return paths
end

"""
Drop pairs whose cluster could only be disconnected, returning how many went.
Screening calls this too: removing a pair can orphan the path of another.
"""
function _prune_disconnected!(admissible::BitMatrix, tree::RadialTree)
    removed = 0
    for ij in findall(admissible)
        i, j = Tuple(ij)
        i == j && continue
        interior = interior_path(tree, i, j)
        if interior === nothing || any(k -> !admissible[i, k], interior)
            admissible[ij] = false
            removed += 1
        end
    end
    return removed
end

"""
Force each bus in `pin` to survive, by clearing every off-diagonal way of
representing it. With "represented exactly once" that fixes `A[i,i] = 1` without
adding a row. Pinned buses may still represent *others* -- only their own
elimination is removed.
"""
function _apply_pins!(admissible::AbstractMatrix{Bool}, pin::AbstractVector{Int}, B::Int)
    isempty(pin) && return admissible
    for i in pin
        1 <= i <= B || error("pin names bus $i, which is outside 1:$B.")
        for k in 1:B
            k == i || (admissible[k, i] = false)
        end
        admissible[i, i] = true
    end
    return admissible
end

"Structural constraints: who may represent whom, and what a cluster may look like."
function _assignment_constraints!(model, A, admissible, paths, slack, B, max_reduction)
    # Every bus is represented exactly once -- also what keeps e = Z(A-I)C valid.
    for j in 1:B
        @constraint(model, sum(A[i, j] for i in 1:B if admissible[i, j]) == 1)
    end

    # Only a super-node may represent anything, and the slack always is one.
    for ij in findall(admissible)
        i, j = Tuple(ij)
        i == j || @constraint(model, A[ij] <= A[i, i])
    end
    @constraint(model, A[slack, slack] == 1)

    # Clusters are connected: everything on the path joins too.
    for (ij, interior) in paths
        i, _ = Tuple(ij)
        for k in interior
            @constraint(model, A[ij] <= A[i, k])
        end
    end

    # Optional floor on how much network is allowed to survive.
    max_reduction < 1 &&
        @constraint(model, sum(A[i, i] for i in 1:B) >= B * (1 - max_reduction))
    return model
end

"""
The voltage-error budget, as one indicator-gated pair of inequalities per
(merge, phase, scenario). Returns how many were added.

The phase loop is what makes this the three-phase model: a merge is checked on
every phase the *reduced* bus carries, against the same phase of the super-node,
which `admissible_pairs` has already guaranteed carries them all. On a
single-phase feeder the body runs once, giving the single-phase paper's
constraint.
"""
function _annulus_constraints!(model, a, pairs, net::Network, C, Z, V, absV,
    sel, Ē, approach, needed)

    nph = nphase_rows(net)
    row_of = _phase_rows(net)
    count_added = 0

    # Without screening every (phase, scenario) of every pair gets a row.
    binding(ij, p, k) = needed === nothing || (p, k) in needed[ij]

    # Re(conj(V_i) e_i) at every row/scenario a constraint will ask for. This is
    # the only place the dense Z is touched, and building the aligned scalar
    # directly rather than e_real and e_imag separately halves the work.
    #
    # Z is dense, so each expression spans every binary in the model. They are
    # therefore held as *variables* pinned by one equality each, not substituted
    # into the rows that use them: a super-node is shared by many merges, and
    # substituting would copy a several-thousand-term expression into each --
    # on case533mt at hops=5, 27282 dense rows where 1066 suffice.
    Dr, Di = _injection_shift(a, pairs, net, C, sel, row_of, nph)
    aligned = Dict{Tuple{Int,Int},VariableRef}()

    for ij in pairs, p in 1:3
        i, j = Tuple(ij)
        net.phases[p, j] || continue
        ri = row_of[p, i]
        for k in eachindex(sel)
            binding(ij, p, k) || continue
            haskey(aligned, (ri, k)) && continue
            g = conj(V[ri, sel[k]]) .* view(Z, ri, :)
            w = AffExpr(0.0)
            for m in 1:nph
                # Re(g_m * D_m) = Re(g_m) Re(D_m) - Im(g_m) Im(D_m)
                add_to_expression!(w, real(g[m]), Dr[m, k])
                add_to_expression!(w, -imag(g[m]), Di[m, k])
            end
            shared = @variable(model)
            @constraint(model, shared == w)
            aligned[(ri, k)] = shared
        end
    end

    for (idx, ij) in enumerate(pairs), p in 1:3
        i, j = Tuple(ij)
        net.phases[p, j] || continue
        ri, rj = row_of[p, i], row_of[p, j]

        for k in eachindex(sel)
            binding(ij, p, k) || continue
            s = sel[k]
            mag_i, mag_j = absV[ri, s], absV[rj, s]

            # S's aligned term is exactly twice R's, so halve the bounds rather
            # than scale the expression -- one fewer AffExpr per row.
            lower, upper = if approach === :R
                mag_i * ((mag_j - Ē) - mag_i), mag_i * ((mag_j + Ē) - mag_i)
            else
                ((mag_j - Ē)^2 - mag_i^2) / 2, ((mag_j + Ē)^2 - mag_i^2) / 2
            end

            w = aligned[(ri, k)]
            @constraint(model, a[idx] => {w >= lower})
            @constraint(model, a[idx] => {w <= upper})
            count_added += 2
        end
    end
    return count_added
end

"""
`(A - I)C` at the node-phase level, split into real and imaginary parts: the
current each row gains or loses under the assignment, as an affine expression.

This is [`expand_assignment`](@ref) written for JuMP expressions rather than
numbers. Row `(i,p)` collects the current of every bus `j` that could be assigned
to `i` and carries phase `p`, and gives up its own.
"""
function _injection_shift(a, pairs, net::Network, C, sel, row_of, nph)
    Dr = [AffExpr(0.0) for _ in 1:nph, _ in eachindex(sel)]
    Di = [AffExpr(0.0) for _ in 1:nph, _ in eachindex(sel)]

    for (idx, ij) in enumerate(pairs), p in 1:3
        i, j = Tuple(ij)
        net.phases[p, j] || continue          # subset rule: p in phases(j) => p in phases(i)
        ri, rj = row_of[p, i], row_of[p, j]
        for k in eachindex(sel)
            c = C[rj, sel[k]]
            add_to_expression!(Dr[ri, k], real(c), a[idx])
            add_to_expression!(Di[ri, k], imag(c), a[idx])
        end
    end

    for r in 1:nph, k in eachindex(sel)     # the -I term
        c = C[r, sel[k]]
        add_to_expression!(Dr[r, k], -real(c))
        add_to_expression!(Di[r, k], -imag(c))
    end
    return Dr, Di
end

"`row_of[p, i]` is the Ybus row of phase `p` at bus `i`, or 0 when absent."
function _phase_rows(net::Network)
    blocks = node_rows(net)
    row_of = zeros(Int, 3, nnodes(net))
    for i in 1:nnodes(net)
        local_idx = 0
        for p in 1:3
            net.phases[p, i] || continue
            local_idx += 1
            row_of[p, i] = blocks[i][local_idx]
        end
    end
    return row_of
end

"""
    annulus_violation(net, A, V, Ē; scenarios) -> Float64

Worst violation of the *original nonconvex* annulus, over every active merge and
every scenario. `<= 0` means the assignment genuinely respects the budget.

This is the correctness bar, and it deliberately shares no code with the model
above: it solves `Ybus[ns,ns] e = ((A-I)C)[ns,:]` on the untouched sparse `Ybus`
and evaluates `| |V_i + e_i| - |V_j| |` directly. The MILP has to produce an
assignment that survives *this* test, not one that matches a linearisation.

Pass the scenarios the model did *not* see to find out whether a reduction fitted
to representative scenarios holds up across the rest.
"""
function annulus_violation(net::Network, A::AbstractMatrix, V::AbstractMatrix,
    Ē::Real; scenarios=1:nscenarios(net))

    sel = collect(scenarios)
    nph = nphase_rows(net)
    blocks = node_rows(net)
    row_of = _phase_rows(net)

    C = conj.(Matrix{ComplexF64}(net.S) ./ Matrix{ComplexF64}(V))
    rhs = (expand_assignment(A, net) - I) * C[:, sel]

    non_slack = setdiff(1:nph, blocks[net.slack])
    e = zeros(ComplexF64, nph, length(sel))
    e[non_slack, :] = Matrix(net.Ybus[non_slack, non_slack]) \ rhs[non_slack, :]

    worst = -Inf
    for j in axes(A, 2), i in axes(A, 1)
        iszero(A[i, j]) && continue
        for p in 1:3
            net.phases[p, j] || continue
            ri, rj = row_of[p, i], row_of[p, j]
            for k in eachindex(sel)
                deviation = abs(abs(e[ri, k] + V[ri, sel[k]]) - abs(V[rj, sel[k]]))
                worst = max(worst, deviation - Ē)
            end
        end
    end
    return worst
end

"Turn the `warm_start` option into an assignment the model can start from."
function _resolve_warm_start(spec, net, V, Ē, admissible, tree, sel,
    max_reduction, enforce_radiality, reach, prefer, time_limit)

    spec isa AbstractMatrix && return Matrix{Float64}(spec)
    spec === :identity && return identity_assignment(net)
    spec === :zero_injection || error(
        "warm_start must be :zero_injection, :identity, or an assignment matrix; got :$spec.")

    return zero_injection_warmstart(net, V, Ē, admissible, tree, sel;
        max_reduction=max_reduction, enforce_radiality=enforce_radiality,
        reach=reach, prefer=prefer, time_limit=time_limit)
end
