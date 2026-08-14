# --------------------------------------------------------------------------- #
# The Opti-KRON MILP -- one phase-aware formulation.
#
# There is no separate single-phase model. A balanced feeder is the degenerate
# case where every bus carries one phase, so the phase loops below run once and
# the model collapses to the journal paper's formulation on its own. Writing two
# models would mean maintaining two, and they would drift.
#
#   minimise   sum_i A[i,i]                    (number of super-nodes kept)
#   subject to sum_i A[i,j] = 1     for all j  (every bus is represented once)
#              A[i,j] <= A[i,i]                (only super-nodes represent)
#              A[i,j] <= A[i,k]                (clusters are connected)
#              A[i,j] = 1  =>  annulus_ij      (merge respects the error budget)
#
# ---- The error model ------------------------------------------------------ #
#
# Moving bus j's injection to super-node i perturbs every voltage in the feeder.
# Under the constant-current model (C = conj(S ./ V) held fixed at the operating
# point) that perturbation is exactly
#
#     e = Z (A - I) C
#
# which is *linear in A* -- the whole reason the problem is a MILP rather than a
# bilinear program. This is single-shot: one model, one solve. No iterative
# decomposition, no cutting planes, no big-M.
#
# Two things about Z are worth stating, because both look wrong at first glance.
#
# 1. Z is the *slack-referenced* bus impedance matrix: zero on the slack's rows
#    and columns, and inv(Ybus[ns,ns]) on the rest. That is the exact model --
#    the slack holds its voltage, so it absorbs whatever current the merge
#    pushes at it and passes none of it on -- and it is what `annulus_violation`
#    re-derives independently to check the answer.
#
#    Plain inv(Ybus) is a near miss rather than an equivalent, and it is worth
#    being precise about the gap. It needs Ybus to be non-singular at all, which
#    means a shunt to ground somewhere; a feeder assembled from branch data
#    alone has none, every row of Ybus sums to zero, and inv returns entries
#    around 1e12 that cancel to millivolt-scale garbage in the annulus rows
#    instead of failing outright. Even when a grounding exists, solving
#    (Ybus + G) x = r and the slack-referenced system give answers differing by
#    a uniform shift c*1, and c vanishes only for particular r. On R100 -- whose
#    grounding sits on the slack rows -- the two agreed to 1.2e-16 over a real
#    assignment, with the slack row of e at 7e-17, so passing a precomputed
#    inv(Ybus) as `Z` is reasonable there. It is not safe in general, which is
#    why it is not the default.
#
#    What *is* structural, and what the model does rely on, is that the
#    column-sum constraint makes A column-stochastic, so (A-I)C sums to zero
#    over the node-phase rows of each phase: aggregation moves current between
#    buses, it never creates any.
#
# 2. Z is dense and Ybus is sparse, so using Z looks like it throws away
#    sparsity for nothing. It buys conditioning. A near-zero-impedance branch (a
#    jumper or a closed switch) puts a near-*infinite* entry in Ybus -- on R100,
#    |Ybus| spans 0.83 to 4.7e6, a condition number of 2.9e9 -- and that single
#    outlier would set the conditioning of the whole constraint matrix. The same
#    branch puts a near-*zero* entry in Z, whose magnitudes span only 4.5e-10 to
#    7.1. Z is the better-behaved object precisely because it is the inverse.
#
# ---- The annulus ---------------------------------------------------------- #
#
# The accuracy requirement is on voltage *magnitude*: after the merge, bus i's
# voltage must sit within Ē of the magnitude bus j had before it,
#
#     | |V_i + e_i| - |V_j| |  <=  Ē
#
# an annulus in the complex plane. It excludes a disk, so it is nonconvex and
# cannot go to a MILP as written. Both linearisations below come from the
# report "Linear and Robust Reformulations of a Complex-Magnitude (Annulus)
# Constraint" and both work on the aligned error term Re(conj(V_i) e_i), which
# is the component of the perturbation that actually moves the magnitude:
#
#   :R  tracks that term against a first-order bound on |V_i + e_i|
#   :S  uses the full secant of |V_i + e_i|^2, so its terms are squared
#
# Neither carries the report's robustness margin (Q_i for R, E_i for S). That
# margin assumes A activates one merge at a time; here A aggregates every
# simultaneous merge, and the only bound available then is a worst-case sum over
# every admissible column, which on realistic feeders exceeds |V_j| + Ē and
# collapses the upper bound below the lower one -- making even the do-nothing
# assignment infeasible. So the margins are set to zero. Approach R's lower
# bound is exact regardless; the rest is a linearisation, not a certified bound,
# which is why `annulus_violation` re-checks the answer against the untouched
# nonconvex constraint.
#
# ---- On/off ---------------------------------------------------------------- #
#
# The annulus only applies to merges the solver actually selects, so each pair's
# two inequalities hang off A[i,j] as JuMP indicator constraints. The research
# code used big-M for this and had to size M from a worst-case deviation bound.
# Indicators need no such constant: the solver enforces the implication directly,
# there is no bound to get wrong, and nothing dilutes the LP relaxation. Gurobi
# and HiGHS both support them.
# --------------------------------------------------------------------------- #

"""
    MilpSolution

What the MILP returned, and what it cost.

- `A`            node-level assignment, `B x B`. `A[i,j] = 1` means bus `j` is
                 represented by super-node `i`; `A[i,i] = 1` marks a super-node.
- `kept`         the super-nodes, i.e. `super_nodes(A)`.
- `reduction`    fraction of buses eliminated, in `[0, 1]`.
- `status`       JuMP termination status. Check it: a time-limited run can
                 return a feasible but non-optimal assignment.
- `gap`          relative MIP gap at termination, `0.0` when proven optimal.
- `scenarios`    the loading scenarios the error budget was enforced over --
                 the assignment is only certified against these.
- `screening`    what preprocessing removed, or `nothing` when `screen=false`.
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
reports a sentinel around 1e97 rather than a gap; printing that as a percentage
produces a hundred digits of nonsense where "unknown" is the honest answer.
"""
function _gap_text(gap)
    iszero(gap) && return ""
    isfinite(gap) && gap < 1e6 || return ", gap unknown"
    return ", gap $(round(100 * gap, digits=2))%"
end

"""
    solve_milp(net, V, Ē; kwargs...) -> MilpSolution

Reduce `net` as far as possible while holding every merged bus within `Ē`
per unit of its original voltage magnitude.

`V` is the AC power-flow solution of the *unreduced* feeder, `nphase_rows(net)`
x `nscenarios(net)`, and is the operating point the constant-current model
linearises around. [`powerflow`](@ref) will produce one; it is passed in rather
than computed here because the operating point is a modelling choice. A feeder
with delta connections, ZIP loads or regulator taps should be solved in a tool
that models them, and those voltages handed over directly. Check any `V` from
elsewhere with [`powerflow_residual`](@ref) -- an operating point inconsistent
with `Ybus` and `S` invalidates the error bound while still solving happily.

Keyword arguments:

- `scenarios`      which columns of `S`/`V` the budget is enforced over.
                   Defaults to all of them, which is the honest default but
                   scales the model linearly -- the papers enforce 2-3
                   representative scenarios and validate against the rest.
- `hops`   hop limit between a reduced bus and its super-node.
- `direction`      `:any`, or `:downstream` to absorb only towards the slack.
- `approach`       `:R` or `:S`, the two annulus linearisations (see header).
- `max_reduction`  refuse to eliminate more than this fraction of buses.
- `pin`            bus indices that must survive as super-nodes. Used to keep
                   transformers, regulators and switches intact: pin both
                   terminals and the device's admittance comes through the Schur
                   complement unchanged, exactly. See `src/io/devices.jl` for
                   why terminal-pinning is sufficient, and `read_devices` for
                   getting the list out of a case directory.
- `warm_start`     `:zero_injection` (default) solves a small exact subproblem
                   first -- merging only buses that inject nothing, which
                   perturbs no voltage at all -- and starts from that feasible
                   point. `:identity` starts from no reduction. An assignment
                   matrix is used as given.
- `prefer`, `verbose`, `time_limit`, `mip_gap` go to [`select_optimizer`](@ref).

The returned assignment is *not* radialised -- Kron reduction meshes a radial
feeder, and repairing that is a separate step. Follow with [`radialize`](@ref).
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
    # Pinning is expressed as geometry, not as an extra constraint: take away
    # every off-diagonal way of representing a pinned bus and "represented
    # exactly once" forces A[i,i] = 1 on its own. Doing it here rather than
    # after the fact means screening, the warm start and the binary count all
    # see a smaller problem instead of a same-sized one with a fixed variable.
    _apply_pins!(admissible, pin, B)
    # Radiality reasons about who can reach whom geometrically, so it needs the
    # hop mask before screening removes pairs for unrelated reasons.
    reach = enforce_radiality ? copy(admissible) : nothing
    paths = _cluster_paths!(admissible, tree)

    # ---- operating point --------------------------------------------------- #
    # Constant current: the injection a bus draws is fixed at whatever current it
    # was drawing in the base case, so moving the bus moves that current intact.
    C = conj.(Matrix{ComplexF64}(net.S) ./ Matrix{ComplexF64}(V))
    Zfull = Z === nothing ? bus_impedance(net) : Matrix{ComplexF64}(Z)
    size(Zfull) == (nph, nph) || error("Z must be $nph x $nph; got $(size(Zfull)).")
    absV = abs.(V)

    # ---- screening --------------------------------------------------------- #
    # Kills pairs that cannot satisfy the budget under any completion, and marks
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

    # Resolved *after* screening, so the start can only ever pick pairs this
    # model has a binary for -- see the header of warmstart.jl for what goes
    # wrong otherwise.
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

The slack-referenced bus impedance matrix: `inv(Ybus)` on the non-slack
node-phase rows, zero on the slack's.

Those zero rows and columns are the model, not a truncation. The slack holds its
voltage against any current the network asks of it, so a perturbation there
propagates nowhere (zero column) and the slack's own voltage never moves (zero
row).

Inverting the full `Ybus` instead only works when the feeder carries a shunt to
ground somewhere -- otherwise `Ybus` is a singular Laplacian, every row summing
to zero, and `inv` quietly returns entries around 1e12 that cancel to noise in
the constraint rows rather than failing outright. Feeders assembled from branch
data with no shunts are exactly that case.
"""
function bus_impedance(net::Network)
    nph = nphase_rows(net)
    non_slack = setdiff(1:nph, node_rows(net)[net.slack])
    Y = Matrix{ComplexF64}(net.Ybus[non_slack, non_slack])

    Z = zeros(ComplexF64, nph, nph)
    Z[non_slack, non_slack] = inv(Y)

    # inv() on a numerically singular block returns huge entries rather than
    # throwing, and the damage only shows up much later as constraint rows that
    # cancel to noise. One matvec is enough to catch it here instead.
    probe = randn(ComplexF64, length(non_slack))
    residual = norm(Y * (Z[non_slack, non_slack] * probe) - probe) / norm(probe)
    residual < 1e-6 ||
        error("The admittance matrix is numerically singular (relative residual " *
              "$(round(residual, sigdigits=3)) on inversion). Check the feeder for " *
              "an isolated island or a branch with near-zero impedance.")
    return Z
end

"""
    _cluster_paths!(admissible, tree) -> Dict{CartesianIndex,Vector{Int}}

Interior paths of every admissible pair, dropping the pairs that could only
form a disconnected cluster.

A cluster has to be a connected piece of the feeder, which means every bus `k`
strictly between `i` and `j` must join `i` too. That is the `A[i,j] <= A[i,k]`
constraint below. But on an unbalanced feeder `k` may be *ineligible* to join
`i` -- two single-phase laterals on phase `a`, hanging off a common three-phase
backbone bus, are admissible to each other under the subset rule while the
backbone bus between them carries phases `i` does not have. Merging across it
would leave a cluster in two disconnected pieces.

The research code skipped those constraints (`admissibles[i,k] || continue`),
which silently permits exactly that. Here the pair is dropped instead, which
costs a little reduction and keeps every cluster connected. On R100 at
`hops=5` this removes 70 of 1198 off-diagonal pairs.

One pass suffices. If `(i,k)` is itself dropped for some bus `m` on its own
interior path, then `m` also lies on the path from `i` to `j` -- paths in a tree
compose -- so `(i,j)` is dropped in the same sweep.
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
    _prune_disconnected!(admissible, tree) -> Int

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
    _apply_pins!(admissible, pin, B) -> admissible

Force each bus in `pin` to survive, by clearing every off-diagonal way of
representing it. Combined with the "represented exactly once" constraint this
fixes `A[i,i] = 1` without adding a row, and shrinks the model rather than
merely constraining it. Pinned buses may still represent *other* buses -- only
their own elimination is removed.
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
    # Every bus is represented exactly once. This is also what keeps
    # e = Z(A-I)C valid with Z = inv(Ybus) -- see the header.
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
    _annulus_constraints!(...) -> Int

The voltage-error budget, as one indicator-gated pair of inequalities per
(merge, phase, scenario). Returns how many were added.

The phase loop is what makes this the three-phase model: a merge is checked on
every phase the *reduced* bus carries, comparing each against the same phase of
the super-node. `admissible_pairs` has already guaranteed the super-node carries
them all. On a single-phase feeder the loop body runs once and this is the
journal paper's constraint.
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
    # directly -- rather than e_real and e_imag separately -- halves the work,
    # since every constraint uses only this combination.
    #
    # Each of these expressions is dense: Z is dense, so it spans every binary
    # in the model. They are therefore held as *variables*, pinned by one
    # equality each, rather than as raw expressions substituted into the
    # constraints that use them. A super-node is shared by many merges, so
    # substituting would copy a several-thousand-term expression into every row
    # referencing it -- on case533mt at hops=5 that is 27282 dense rows where
    # 1066 suffice. Through a variable each later use is a single coefficient,
    # and the dense part is paid once.
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

            # Approach S's aligned term is exactly twice Approach R's, so rather
            # than scaling the expression we halve the bounds -- same constraint,
            # one fewer AffExpr per row.
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
    _injection_shift(...) -> (Dr, Di)

`(A - I)C` at the node-phase level, split into real and imaginary parts: the
current each row gains or loses under the assignment, as an affine expression.

This is `expand_assignment` written for JuMP expressions rather than numbers --
the same `A ⊗ I₃` generalised to mixed-phase buses, but carrying the variables
through instead of a sparsity pattern. Row `(i,p)` collects the current of every
bus `j` that could be assigned to `i` and that carries phase `p`, and gives up
its own.
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
with the slack held fixed, and evaluates `| |V_i + e_i| - |V_j| |` directly. The
MILP only has to produce an assignment that survives *this* test -- it does not
have to match any particular linearisation.

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

"""
    _resolve_warm_start(spec, ...) -> Matrix{Float64}

Turn the `warm_start` option into an assignment the model can start from.
"""
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
