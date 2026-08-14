# --------------------------------------------------------------------------- #
# Constraint screening: shrink the MILP before the solver sees it.
#
#   kill  a pair whose merge violates the budget under *every* completion of the
#         assignment -- the binary is dropped entirely
#   skip  a row that holds under every completion -- no indicators are built
#
# Both are conservative, so the screen removes only what could not have changed
# the answer -- and it screens the *true* annulus, not its linearisation.
#
# The error at row i is e[i,s] = sum_m (Z[i,rep(m)] - Z[i,m]) C[m,s], where
# rep(m) is unknown but must be one of m's admissible representatives. Bounding
# term by term gives a radius, and the reverse triangle inequality turns that
# into a window |e_i + V_i| in [|V_i| - dev, |V_i| + dev]. The merge (i,j) needs
# that inside [|V_j| - Ē, |V_j| + Ē]: disjoint means kill, contained means skip,
# overlapping means build it and let the solver decide.
#
# Two radii, neither dominant. `_magnitude_radius` is a pure triangle inequality;
# `_rectangular_radius` brackets Re and Im separately and takes the far corner of
# the box, which treats them as independent when they come from one shared choice
# of representative -- with Z[i,m]=0, Z[i,r1]=1, Z[i,r2]=j, C=1 the truth is 1,
# magnitude gives 1, the box gives sqrt(2). Both valid, so is their minimum.
# --------------------------------------------------------------------------- #

"""
    ScreeningReport

What the screen removed, and what it left behind.

- `pairs_removed`     binaries eliminated, split by the rule that caught them.
- `rows_skipped`      annulus rows proved redundant.
- `rows_kept`         annulus rows actually built.
- `deviation`         the radius bound on `|e|`, per node-phase row and scenario.
"""
struct ScreeningReport
    pairs_before::Int
    pairs_after::Int
    voltage_pairs_removed::Int
    bound_pairs_removed::Int
    rows_kept::Int
    rows_skipped::Int
    passes::Int
    deviation::Matrix{Float64}
end

function Base.show(io::IO, r::ScreeningReport)
    removed = r.pairs_before - r.pairs_after
    print(io, "ScreeningReport(", r.pairs_after, "/", r.pairs_before, " pairs kept, ",
        removed, " removed [", r.voltage_pairs_removed, " voltage-spread, ",
        r.bound_pairs_removed, " bound]; ", r.rows_kept, " rows kept, ",
        r.rows_skipped, " skipped; ", r.passes, " pass(es))")
end

"""
    screen!(admissible, net, tree, V, Ē, Z, C, sel; max_passes) -> (needed, report)

Prune `admissible` in place and decide which annulus rows to build.

`needed` maps a surviving pair to the `(phase, scenario-slot)` combinations that
can actually bind. Runs to a fixed point: killing pairs shrinks the
representative sets, which tightens the radius, which can kill more.
"""
function screen!(admissible::BitMatrix, net::Network, tree::RadialTree,
    V::AbstractMatrix, Ē::Real, Z::AbstractMatrix, C::AbstractMatrix,
    sel::AbstractVector{Int}; max_passes::Int=3)

    pairs_before = count(admissible)
    absV = abs.(V)

    # Cheap and exact, so it runs first and shrinks what the costly screen sees.
    voltage_removed = screen_voltage_spread!(admissible, net, tree, absV, Ē, sel)

    row_of = _phase_rows(net)
    needed = Dict{CartesianIndex{2},Vector{Tuple{Int,Int}}}()
    deviation = zeros(Float64, nphase_rows(net), length(sel))
    bound_removed = 0
    passes = 0

    for pass in 1:max_passes
        passes = pass
        deviation = deviation_bound(Z, C, admissible, net, sel)
        empty!(needed)
        killed_this_pass = 0

        for ij in findall(admissible)
            i, j = Tuple(ij)
            binding = Tuple{Int,Int}[]
            doomed = false

            for p in 1:3
                net.phases[p, j] || continue
                ri, rj = row_of[p, i], row_of[p, j]
                for k in eachindex(sel)
                    s = sel[k]
                    dev = deviation[ri, k]
                    reachable = (max(0.0, absV[ri, s] - dev), absV[ri, s] + dev)
                    allowed = (absV[rj, s] - Ē, absV[rj, s] + Ē)

                    if reachable[1] > allowed[2] || reachable[2] < allowed[1]
                        doomed = true                       # no completion can satisfy it
                        break
                    elseif reachable[1] < allowed[1] || reachable[2] > allowed[2]
                        push!(binding, (p, k))              # might bind; build it
                    end
                end
                doomed && break
            end

            # The diagonal is never killable -- e = 0 is always reachable -- so a
            # doomed diagonal means the bound is wrong, not the pair.
            if doomed && i != j
                admissible[ij] = false
                killed_this_pass += 1
            else
                needed[ij] = binding
            end
        end

        bound_removed += killed_this_pass
        killed_this_pass == 0 && break
        # Killing pairs can orphan a path; restore contiguity before re-bounding.
        _prune_disconnected!(admissible, tree)
    end

    rows_kept = sum(length, values(needed); init=0)
    total_rows = sum(_rows_for_pair(net, ij, length(sel)) for ij in findall(admissible); init=0)

    report = ScreeningReport(pairs_before, count(admissible), voltage_removed,
        bound_removed, rows_kept, total_rows - rows_kept, passes, deviation)
    return needed, report
end

"Annulus rows a pair would need without screening: one per phase of j, per scenario."
_rows_for_pair(net::Network, ij::CartesianIndex, nscen::Int) =
    count(view(net.phases, :, ij[2])) * nscen

"""
    screen_voltage_spread!(admissible, net, tree, absV, Ē, sel) -> Int

Kill pairs whose cluster could never fit the budget, from voltages alone.

Exact, and it sees what the per-pair bound cannot: contiguity forces every bus
between `i` and `j` into one cluster, so all take one post-merge voltage, so
their original magnitudes must fit a window of width `2*Ē`. If they do not, no
assignment rescues the pair, for any injections -- which bites hardest on
long-range pairs, where the injection bound is weakest.

Only scenarios in `sel` count: a spread appearing solely in an unmodelled
scenario is not grounds to remove a pair the model was free to use.
"""
function screen_voltage_spread!(admissible::BitMatrix, net::Network, tree::RadialTree,
    absV::AbstractMatrix, Ē::Real, sel::AbstractVector{Int})

    row_of = _phase_rows(net)
    removed = 0

    for ij in findall(admissible)
        i, j = Tuple(ij)
        i == j && continue
        interior = interior_path(tree, i, j)
        interior === nothing && continue
        members = vcat(interior, j)         # i itself is the centre, not a member

        for p in 1:3
            net.phases[p, j] || continue
            for s in sel
                lo, hi = Inf, -Inf
                for k in members
                    net.phases[p, k] || continue
                    v = absV[row_of[p, k], s]
                    lo, hi = min(lo, v), max(hi, v)
                end
                if hi - lo > 2 * Ē
                    admissible[ij] = false
                    removed += 1
                    @goto next_pair
                end
            end
        end
        @label next_pair
    end
    return removed
end

"""
    deviation_bound(Z, C, admissible, net, sel) -> Matrix{Float64}

Radius bounding `|e[i,s]|` over every assignment the current mask allows, per
node-phase row and scenario. Pointwise minimum of the magnitude and rectangular
bounds; see the header for why neither dominates.
"""
function deviation_bound(Z::AbstractMatrix, C::AbstractMatrix, admissible::BitMatrix,
    net::Network, sel::AbstractVector{Int})

    reps = _representative_rows(admissible, net)
    magnitude = _magnitude_radius(Z, C, reps, sel)
    rectangular = _rectangular_radius(Z, C, reps, sel)
    return min.(magnitude, rectangular)
end

"""
For each node-phase row, the rows it could be represented by. Phase `p` of bus
`j` only ever moves to phase `p` of a bus admitting `j`, which keeps this a
per-row question rather than a per-bus one.
"""
function _representative_rows(admissible::BitMatrix, net::Network)
    row_of = _phase_rows(net)
    reps = [Int[] for _ in 1:nphase_rows(net)]
    for j in 1:nnodes(net), p in 1:3
        net.phases[p, j] || continue
        rj = row_of[p, j]
        for i in 1:nnodes(net)
            admissible[i, j] || continue
            push!(reps[rj], row_of[p, i])
        end
    end
    return reps
end

"Triangle-inequality radius: sum over members of the widest term they can contribute."
function _magnitude_radius(Z, C, reps, sel)
    nph = length(reps)
    spread = zeros(Float64, nph, nph)          # spread[i,m] = max_r |Z[i,r] - Z[i,m]|
    for m in 1:nph, i in 1:nph
        widest = 0.0
        base = Z[i, m]
        for r in reps[m]
            widest = max(widest, abs(Z[i, r] - base))
        end
        spread[i, m] = widest
    end
    return spread * abs.(view(C, :, sel))
end

"""
Bracket `Re(e)` and `Im(e)` separately, then take the far corner of the box.

Each term is `(Z[i,r] - Z[i,m]) C[m,s]`, with real part `dRe*Cr - dIm*Ci` and
imaginary part `dRe*Ci + dIm*Cr`. Both are linear in `(dRe, dIm)` with
coefficients independent of `r`, so intervals for `Re(Z[i,r])` and `Im(Z[i,r])`
over the admissible reps -- computed once, for all scenarios -- combine with the
sign of each `C` entry to give a per-scenario bound.
"""
function _rectangular_radius(Z, C, reps, sel)
    nph = length(reps)
    re_lo = zeros(Float64, nph, nph)
    re_hi = zeros(Float64, nph, nph)
    im_lo = zeros(Float64, nph, nph)
    im_hi = zeros(Float64, nph, nph)

    for m in 1:nph, i in 1:nph
        base = Z[i, m]
        rlo = rhi = ilo = ihi = 0.0            # r == m is always available
        for r in reps[m]
            dre, dim = real(Z[i, r] - base), imag(Z[i, r] - base)
            rlo, rhi = min(rlo, dre), max(rhi, dre)
            ilo, ihi = min(ilo, dim), max(ihi, dim)
        end
        re_lo[i, m], re_hi[i, m] = rlo, rhi
        im_lo[i, m], im_hi[i, m] = ilo, ihi
    end

    radius = zeros(Float64, nph, length(sel))
    for (k, s) in enumerate(sel), i in 1:nph
        re_min = re_max = im_min = im_max = 0.0
        for m in 1:nph
            cr, ci = real(C[m, s]), imag(C[m, s])
            # Re = dRe*cr - dIm*ci ; Im = dRe*ci + dIm*cr
            re_min += _term_min(re_lo[i, m], re_hi[i, m], cr) + _term_min(im_lo[i, m], im_hi[i, m], -ci)
            re_max += _term_max(re_lo[i, m], re_hi[i, m], cr) + _term_max(im_lo[i, m], im_hi[i, m], -ci)
            im_min += _term_min(re_lo[i, m], re_hi[i, m], ci) + _term_min(im_lo[i, m], im_hi[i, m], cr)
            im_max += _term_max(re_lo[i, m], re_hi[i, m], ci) + _term_max(im_lo[i, m], im_hi[i, m], cr)
        end
        radius[i, k] = hypot(max(abs(re_min), abs(re_max)), max(abs(im_min), abs(im_max)))
    end
    return radius
end

_term_min(lo, hi, coefficient) = coefficient >= 0 ? coefficient * lo : coefficient * hi
_term_max(lo, hi, coefficient) = coefficient >= 0 ? coefficient * hi : coefficient * lo
