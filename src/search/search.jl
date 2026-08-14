# --------------------------------------------------------------------------- #
# Exhaustive search at q = 1 -- Algorithm 1 of the PSCC/EPSR three-phase paper.
#
# Where the MILP answers "what is the best assignment?" exactly, the search asks
# "what is the best single merge available now?" and repeats. One bus is
# eliminated per iteration; every feasible pair (s, r) is scored, the best taken,
# the graph updated. Greedy, so not optimal -- but each step is embarrassingly
# parallel and needs no solver, which keeps it going where the MILP runs out.
#
# The budget is the MILP's, exactly, so "within E" means one thing whichever
# backend produced the answer: C = conj(S ./ V), e = Z (A - I) C with Z
# slack-referenced, deviation | |V_i + e_i| - |V_j| |. Scoring against Y directly
# -- what the research code does, on a feeder it has already grounded -- lets the
# slack float, understating the error until the search appears to beat the MILP.
#
# Per-candidate cost is an axpy, not a solve: candidates within an iteration
# differ in rank one, so with Z formed once the update is a scaled column
# difference. See gpu.jl, where that arithmetic lives.
# --------------------------------------------------------------------------- #

"""
    SearchSolution

What the exhaustive search returned, and what it cost.

- `A`          node-level assignment, `B x B`, same convention as the MILP's.
- `kept`       the super-nodes.
- `reduction`  fraction of buses eliminated.
- `objective`  the worst voltage-magnitude error of the accepted assignment.
- `status`     `:exhausted` (no feasible merge left -- Algorithm 1's `Δ = 0`),
               `:max_reduction`, or `:time_limit`.
- `iterations` merges accepted, which is also the number of buses eliminated.
- `evaluated`  candidates scored; `screened` is how many were skipped.

There is no optimality gap here: the search does not know what it missed. A
`SearchSolution` says "this reduction holds the budget", never "no better one
exists" -- that is what the MILP is for.
"""
struct SearchSolution <: ReductionSolution
    A::Matrix{Float64}
    kept::Vector{Int}
    reduction::Float64
    objective::Float64
    status::Symbol
    solver::Symbol
    scenarios::Vector{Int}
    iterations::Int
    evaluated::Int
    screened::Int
    build_time::Float64
    solve_time::Float64
end

function Base.show(io::IO, s::SearchSolution)
    print(io, "SearchSolution(", length(s.kept), " super-nodes, ",
        round(100 * s.reduction, digits=1), "% reduction, ", s.status,
        ", worst err ", @sprintf("%.1e", s.objective),
        ", $(s.solver) in $(round(s.solve_time, digits=1))s)")
end

"""
    search_reduce(net, V, Ē; kwargs...) -> SearchSolution

Reduce `net` by greedy exhaustive search, holding every bus within `Ē` of its
original voltage magnitude on the selected scenarios.

- `scenarios`      columns of `V` to enforce the budget over. Default: all.
- `objective`      `:max` scores by worst error (Eq. 15), `:sum` by the total.
- `radial`         only eliminate buses of current degree <= 2, which keeps the
                   reduced network radial by construction. A bus becomes
                   eligible later as its neighbours are absorbed.
- `upstream`       only merge towards the slack, by BFS depth.
- `max_reduction`  stop once this fraction of buses is gone.
- `pin`            buses that must survive; never offered for elimination.
                   Keeps equipment intact -- see `src/io/devices.jl`.
- `screen`         drop provably-infeasible candidates first (rigorous).
- `backend`        `:cpu` / `:zbus` (host) or `:gpu`; same arithmetic.
- `batch`          candidates scored per pass.
- `time_limit`     seconds, checked between iterations.
"""
function search_reduce(net::Network, V::AbstractMatrix, Ē::Real;
    scenarios=1:nscenarios(net),
    objective::Symbol=:max,
    radial::Bool=true,
    upstream::Bool=false,
    max_reduction::Real=1.0,
    pin::AbstractVector{Int}=Int[],
    screen::Bool=true,
    backend::Symbol=:cpu,
    batch::Int=256,
    Zbus::Union{Nothing,AbstractMatrix}=nothing,
    time_limit::Union{Nothing,Real}=nothing,
    verbose::Bool=false)

    objective in (:max, :sum) || error("objective must be :max or :sum; got :$objective.")
    backend in (:cpu, :zbus, :gpu) ||
        error("backend must be :cpu, :zbus or :gpu; got :$backend.")
    Ē > 0 || error("Ē must be positive; got $Ē.")
    0 <= max_reduction <= 1 || error("max_reduction must be in [0, 1]; got $max_reduction.")

    B, nph = nnodes(net), nphase_rows(net)
    size(V) == (nph, nscenarios(net)) ||
        error("V is $(size(V)) but the network implies $((nph, nscenarios(net))).")

    sel = collect(scenarios)
    isempty(sel) && error("`scenarios` selected no columns; the error budget would be vacuous.")
    all(s -> 1 <= s <= nscenarios(net), sel) ||
        error("`scenarios` indexes outside 1:$(nscenarios(net)).")

    # A missing GPU is a normal state, not a failure: CUDA is optional, and most
    # machines running this will not have one. Fall back, but say so -- a run
    # that silently took ten times longer than expected is worse than a warning.
    if backend === :gpu && !gpu_available()
        @warn "CUDA is unavailable, so the search is running on the CPU. " *
              "Same result, same error model, slower on large feeders. " *
              "Install CUDA.jl and check `CUDA.functional()` to use the GPU."
        backend = :cpu
    end

    t_build = time()

    V̂ = Matrix{ComplexF64}(V[:, sel])
    absV = abs.(V̂)
    # The same error model the MILP is stated in, and the same one
    # `annulus_violation` re-checks against: constant current C = conj(S ./ V),
    # perturbation e = Z (A - I) C, with Z slack-referenced so the slack holds
    # its voltage. Using Y directly instead -- as the research code does, on a
    # feeder it has already grounded -- lets the slack float, which understates
    # the error and lets the search merge past the budget.
    C = conj.(Matrix{ComplexF64}(net.S[:, sel]) ./ V̂)
    Z = Zbus === nothing ? bus_impedance(net) : Matrix{ComplexF64}(Zbus)
    size(Z) == (nph, nph) || error("Zbus must be $nph x $nph; got $(size(Z)).")

    blocks = node_rows(net)
    neighbours = _neighbour_sets(net.Lambda)
    depth = upstream ? _bus_depths(neighbours, net.slack) : nothing

    protected = falses(B)
    for i in pin
        1 <= i <= B || error("pin names bus $i, which is outside 1:$B.")
        protected[i] = true
    end

    A = Matrix{Float64}(LinearAlgebra.I, B, B)
    alive = trues(B)                       # currently a super-node
    Cagg = copy(C)                         # aggregated current, per phase row
    Vcur = copy(V̂)                         # V + e, the voltage each row now holds
    owner = collect(1:nph)                 # phase row -> row holding its voltage
    members = [[i] for i in 1:B]           # buses each super-node stands for

    solver = backend === :gpu ? :search_gpu : :search_cpu
    context = _search_context(Val(backend), Z, Vcur, absV, batch)

    build_time = time() - t_build
    t_solve = time()

    floor_kept = max(1, ceil(Int, B * (1 - max_reduction)))
    evaluated = screened = iterations = 0
    worst = 0.0
    status = :exhausted

    while count(alive) > floor_kept
        if time_limit !== nothing && time() - t_solve > time_limit
            status = :time_limit
            break
        end

        pairs = _candidate_pairs(net, neighbours, alive, depth, radial, protected)
        isempty(pairs) && break
        if screen
            kept_pairs = _screen_pairs(pairs, net, blocks, absV, Ē)
            screened += length(pairs) - length(kept_pairs)
            pairs = kept_pairs
            isempty(pairs) && break
        end

        best, score, err = _best_candidate(context, pairs, net, blocks, Cagg, owner,
            absV, Ē, objective)
        evaluated += length(pairs)
        best === nothing && break

        s, r = pairs[best]
        # The context carries state derived from `Cagg` (the Zbus path caches the
        # current voltages), so it is updated from the pre-merge currents before
        # `_accept!` overwrites them.
        _accept_context!(context, net, blocks, Cagg, s, r)
        _accept!(A, alive, Cagg, owner, members, neighbours, net, blocks, s, r)
        iterations += 1
        worst = err
        verbose && iterations % 25 == 0 &&
            @printf("  %4d merges, %4d kept, worst err %.2e, %d candidates\n",
                iterations, count(alive), err, length(pairs))
    end

    count(alive) <= floor_kept && (status = :max_reduction)
    solve_time = time() - t_solve
    kept = findall(alive)

    return SearchSolution(A, kept, 1 - length(kept) / B, worst, status, solver,
        sel, iterations, evaluated, screened, build_time, solve_time)
end

# ---- graph -------------------------------------------------------------- #

"Adjacency as neighbour sets, so merging is O(degree) rather than O(B)."
function _neighbour_sets(Lambda::SparseMatrixCSC)
    B = size(Lambda, 1)
    sets = [Set{Int}() for _ in 1:B]
    rows = rowvals(Lambda)
    for c in 1:B, k in nzrange(Lambda, c)
        i = rows[k]
        i == c && continue
        push!(sets[c], i)
        push!(sets[i], c)
    end
    return sets
end

"BFS depth from the slack, for the `upstream` rule."
function _bus_depths(neighbours::Vector{Set{Int}}, slack::Int)
    depth = fill(typemax(Int), length(neighbours))
    depth[slack] = 0
    queue = [slack]
    head = 1
    while head <= length(queue)
        i = queue[head]
        head += 1
        for j in neighbours[i]
            depth[j] == typemax(Int) || continue
            depth[j] = depth[i] + 1
            push!(queue, j)
        end
    end
    return depth
end

"""
    _candidate_pairs(net, neighbours, alive, depth, radial) -> Vector{Tuple{Int,Int}}

`𝒜 = {(s, r) | Λ_{s,r} >= 1, φ_r ⊆ φ_s}`, restricted by the `radial` and
`upstream` rules. `r` is the bus about to be eliminated, `s` the one absorbing
it; the slack and any pinned bus are never eliminated.
"""
function _candidate_pairs(net::Network, neighbours::Vector{Set{Int}},
    alive::BitVector, depth, radial::Bool, protected::BitVector)

    pairs = Tuple{Int,Int}[]
    for r in 1:nnodes(net)
        (alive[r] && r != net.slack && !protected[r]) || continue
        radial && length(neighbours[r]) > 2 && continue
        for s in neighbours[r]
            alive[s] || continue
            depth === nothing || depth[s] < depth[r] || continue
            # φ_r ⊆ φ_s, per phase and never a phase count: a super-node can
            # only stand for a bus whose phases it actually carries.
            all(p -> !net.phases[p, r] || net.phases[p, s], 1:3) || continue
            push!(pairs, (s, r))
        end
    end
    return pairs
end

"""
    _screen_pairs(pairs, net, blocks, absV, Ē) -> Vector{Tuple{Int,Int}}

Drop candidates that cannot meet the bound under any redistribution.

When `r` takes `s`'s voltage its error is `| |V̂_r| − |V^c_s| |`, and `s`'s own
error stays within `Ē`, so by the triangle inequality the candidate's error is
at least `Vdif(s,r) − Ē`. Anything with `Vdif > 2Ē` is therefore infeasible no
matter what the network does elsewhere -- this prunes without ever being wrong.
`Vdif` uses only the reference voltages, so it is static and costs one pass.
"""
function _screen_pairs(pairs::Vector{Tuple{Int,Int}}, net::Network,
    blocks::Vector{Vector{Int}}, absV::AbstractMatrix, Ē::Real)

    threshold = 2 * Ē
    keep = Tuple{Int,Int}[]
    for (s, r) in pairs
        spread = 0.0
        for (row_r, row_s) in _phase_pairs(net, blocks, s, r),
            k in axes(absV, 2)

            spread = max(spread, abs(absV[row_r, k] - absV[row_s, k]))
            spread > threshold && break
        end
        spread <= threshold && push!(keep, (s, r))
    end
    return keep
end

"""
    _phase_pairs(net, blocks, s, r) -> Vector{Tuple{Int,Int}}

Rows of `r` paired with the rows of `s` carrying the same phase. Assumes
`φ_r ⊆ φ_s`, which `_candidate_pairs` has already checked.
"""
function _phase_pairs(net::Network, blocks::Vector{Vector{Int}}, s::Int, r::Int)
    phases_r = phases_of(net, r)
    phases_s = phases_of(net, s)
    out = Tuple{Int,Int}[]
    for (local_r, p) in enumerate(phases_r)
        local_s = findfirst(==(p), phases_s)
        local_s === nothing && continue
        push!(out, (blocks[r][local_r], blocks[s][local_s]))
    end
    return out
end

# ---- accepting a merge --------------------------------------------------- #

"Fold `r`'s cluster into `s`: assignment, adjacency, current and voltage ownership."
function _accept!(A, alive, Cagg, owner, members, neighbours, net, blocks, s::Int, r::Int)
    for j in members[r]
        A[s, j] = 1.0
        A[r, j] = 0.0
    end
    append!(members[s], members[r])
    empty!(members[r])
    alive[r] = false

    for (row_r, row_s) in _phase_pairs(net, blocks, s, r)
        for k in axes(Cagg, 2)
            Cagg[row_s, k] += Cagg[row_r, k]
            Cagg[row_r, k] = 0
        end
        # Everything that was following r's voltage now follows s's.
        for i in eachindex(owner)
            owner[i] == row_r && (owner[i] = row_s)
        end
    end

    for j in neighbours[r]
        j == s && continue
        push!(neighbours[s], j)
        push!(neighbours[j], s)
        delete!(neighbours[j], r)
    end
    delete!(neighbours[s], r)
    empty!(neighbours[r])
    return A
end
