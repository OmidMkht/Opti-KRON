# --------------------------------------------------------------------------- #
# The Zbus reformulation, and the GPU backend built on it.
#
# The candidates of one iteration differ from each other only in *rank one*.
# Moving bus r's aggregated current onto s changes `(A − I)C` by
# `Cagg[r_p]·(e_{s_p} − e_{r_p})` per phase, so with the slack-referenced
# `Z` in hand
#
#     V^c = V_cur + Σ_p  Cagg[r_p] · ( Z[:, s_p] − Z[:, r_p] )
#
# and a candidate costs a few scaled column differences instead of a solve. `Z`
# is formed once; per-candidate work is an axpy.
#
# It is also what makes the GPU worthwhile. Scoring is then a wide, shallow,
# perfectly parallel sweep over (row, candidate, scenario) -- exactly the shape a
# GPU wants, and exactly the shape a sparse solve is not. Memory does not grow
# with the scenario count either, which is what the reference implementation was
# reaching for.
#
# The cost is `Z`: dense, `nphase_rows^2`. 4 MB on R300, 229 MB on ieee8500.
#
# The context below is generic over its array type, so the identical arithmetic
# runs on `Matrix` or `CuMatrix`. That is deliberate: it means the GPU path can
# be validated against the CPU path running the same formulas, and any
# disagreement is hardware or precision, never a second implementation drifting.
# --------------------------------------------------------------------------- #

const CUDA_MODULE = Ref{Union{Nothing,Module}}(nothing)

"""
    load_cuda!()

Load CUDA.jl if it is installed, without touching a device. Optional in the same
way Gurobi is: declared nowhere, loaded here, and absent is not an error until
someone actually asks for `:gpu`.
"""
function load_cuda!()
    CUDA_MODULE[] === nothing || return CUDA_MODULE[]
    try
        CUDA_MODULE[] = Base.require(Base.PkgId(
            Base.UUID("052768ef-5323-5732-b1bb-66c8b64840ba"), "CUDA"))
    catch
        CUDA_MODULE[] = nothing
    end
    return CUDA_MODULE[]
end

"""
    gpu_available() -> Bool

True when CUDA.jl is installed *and* a device is actually usable. Both halves
matter: an installed package on a machine with no driver is the common case on
a laptop, and it should route to the CPU rather than fail late.
"""
function gpu_available()
    cuda = load_cuda!()
    cuda === nothing && return false
    try
        # `invokelatest`, because `load_cuda!` may have introduced these methods
        # after this function was compiled. Calling straight through would raise
        # a world-age MethodError, the `catch` would swallow it, and a perfectly
        # good GPU would be reported as absent -- exactly the failure
        # `load_gurobi!` exists to prevent on the solver side.
        return Base.invokelatest(cuda.functional)
    catch
        return false
    end
end

"""
Buffers for the Zbus formulation. `Z` is the dense inverse, `Vcur` the voltages
of the assignment accepted so far, and `absV` the reference magnitudes. All
three live wherever the backend puts them -- host or device.
"""
struct Rank1SearchContext{M,R}
    Z::M
    Vcur::M
    absV::R
    batch::Int
    nscenarios::Int
    device::Symbol
end

function _search_context(::Val{:gpu}, Z::Matrix{ComplexF64}, Vcur::Matrix{ComplexF64},
    absV::Matrix{Float64}, batch::Int)

    # `search_reduce` has already downgraded to :cpu if there is no device, so
    # reaching here without one means someone built the context directly.
    gpu_available() || error("""
        backend :gpu needs CUDA.jl and a working device. `gpu_available()` says no.
        Install CUDA.jl and check `CUDA.functional()`, or use backend=:search_cpu,
        which computes the same reduction from the same formulas.""")
    return _rank1_context(Z, Vcur, absV, batch, :gpu)
end

"""
`:zbus` names the host path explicitly, for when a caller wants it regardless of
what `backend` resolution would pick. It is the same context `:cpu` builds.
"""
_search_context(::Val{:zbus}, Z::Matrix{ComplexF64}, Vcur::Matrix{ComplexF64},
    absV::Matrix{Float64}, batch::Int) = _rank1_context(Z, Vcur, absV, batch, :cpu)

"""
    _rank1_context(Z, Vcur, absV, batch, device) -> Rank1SearchContext

Place the slack-referenced impedance and the current voltages on `device`.
"""
function _rank1_context(Z::Matrix{ComplexF64}, Vcur::Matrix{ComplexF64},
    absV::Matrix{Float64}, batch::Int, device::Symbol)

    device === :cpu &&
        return Rank1SearchContext(Z, copy(Vcur), absV, batch, size(Vcur, 2), :cpu)

    CuArray = load_cuda!().CuArray
    return Rank1SearchContext(Base.invokelatest(CuArray, Z),
        Base.invokelatest(CuArray, Vcur), Base.invokelatest(CuArray, absV),
        batch, size(Vcur, 2), :gpu)
end

"""
    _accept_context!(context::Rank1SearchContext, net, blocks, Cagg, s, r)

Carry the accepted merge into the cached voltages, by the same rank-1 update
used to score it. Recomputing `Z * Cagg` here instead would cost a dense
matrix product per iteration and undo the point of the reformulation.
"""
function _accept_context!(context::Rank1SearchContext, net::Network,
    blocks::Vector{Vector{Int}}, Cagg::Matrix{ComplexF64}, s::Int, r::Int)

    for (row_r, row_s) in _phase_pairs(net, blocks, s, r), k in axes(Cagg, 2)
        coefficient = Cagg[row_r, k]
        iszero(coefficient) && continue
        @views context.Vcur[:, k] .+= coefficient .* (context.Z[:, row_s] .- context.Z[:, row_r])
    end
    return context
end

"""
    _best_candidate(context::Rank1SearchContext, ...) -> (index, score, error)

Score every candidate through the rank-1 update. The arithmetic is written in
whole-array operations so it runs unchanged on host and device.
"""
function _best_candidate(context::Rank1SearchContext, pairs::Vector{Tuple{Int,Int}},
    net::Network, blocks::Vector{Vector{Int}}, Cagg::Matrix{ComplexF64},
    owner::Vector{Int}, absV::Matrix{Float64}, Ē::Real, objective::Symbol)

    nph, nscen = size(context.Vcur, 1), context.nscenarios
    best, best_score, best_error = nothing, Inf, Inf

    for start in 1:context.batch:length(pairs)
        stop = min(start + context.batch - 1, length(pairs))
        width = stop - start + 1

        # Up to three phase pairs per candidate; unused slots carry a zero
        # coefficient and a dummy column, so every candidate is the same shape.
        source = ones(Int, 3, width)
        sink = ones(Int, 3, width)
        coefficient = zeros(ComplexF64, 3, width, nscen)
        ownermap = repeat(owner, 1, width)

        for (slot, index) in enumerate(start:stop)
            s, r = pairs[index]
            for (p, (row_r, row_s)) in enumerate(_phase_pairs(net, blocks, s, r))
                sink[p, slot] = row_s
                source[p, slot] = row_r
                for k in 1:nscen
                    coefficient[p, slot, k] = Cagg[row_r, k]
                end
                for i in 1:nph
                    ownermap[i, slot] == row_r && (ownermap[i, slot] = row_s)
                end
            end
        end

        score, error = _score_batch(context, source, sink, coefficient, ownermap,
            width, nph, nscen, Ē, objective)

        for slot in 1:width
            if score[slot] < best_score
                best, best_score, best_error = start + slot - 1, score[slot], error[slot]
            end
        end
    end

    return (isfinite(best_score) ? best : nothing), best_score, best_error
end

"""
    _score_batch(context, source, sink, coefficient, ownermap, ...) -> (scores, worsts)

The rank-1 update and its error, for a whole batch at once.

`Z[:, sink] - Z[:, source]` is the direction each phase's transport moves every
voltage in; scaling by the current being moved and summing over the (at most
three) phases gives the perturbation. Adding `Vcur` gives the candidate's
voltages, and gathering through `ownermap` gives what each bus actually sees.
"""
function _score_batch(context::Rank1SearchContext, source::Matrix{Int},
    sink::Matrix{Int}, coefficient::Array{ComplexF64,3}, ownermap::Matrix{Int},
    width::Int, nph::Int, nscen::Int, Ē::Real, objective::Symbol)

    to_device = context.device === :gpu ? load_cuda!().CuArray : identity

    # nph x (3*width): the column difference for every phase slot in the batch.
    difference = context.Z[:, vec(to_device(sink))] .- context.Z[:, vec(to_device(source))]
    difference = reshape(difference, nph, 3, width)

    scores = fill(Inf, width)
    worsts = fill(Inf, width)
    accumulated = zeros(Float64, width)
    running = zeros(Float64, width)

    for k in 1:nscen
        weights = reshape(to_device(coefficient[:, :, k]), 1, 3, width)
        perturbation = dropdims(sum(difference .* weights, dims=2), dims=2)  # nph x width
        candidate = context.Vcur[:, k] .+ perturbation

        # Each bus reads its super-node's voltage: a flat gather, so that the
        # per-candidate owner map costs one indexing pass rather than a permute.
        offsets = to_device(collect(0:(width-1))' .* nph)
        assigned = reshape(vec(candidate)[vec(to_device(ownermap) .+ offsets)], nph, width)

        deviation = abs.(context.absV[:, k] .- abs.(assigned))
        running .= Array(vec(maximum(deviation, dims=1)))
        accumulated .+= objective === :sum ? Array(vec(sum(deviation, dims=1))) : 0.0

        for slot in 1:width
            worsts[slot] = k == 1 ? running[slot] : max(worsts[slot], running[slot])
        end
    end

    for slot in 1:width
        scores[slot] = worsts[slot] > Ē ? Inf :
                       (objective === :max ? worsts[slot] : accumulated[slot])
    end
    return scores, worsts
end
