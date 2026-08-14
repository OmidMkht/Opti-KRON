# --------------------------------------------------------------------------- #
# Solver selection: Gurobi when it is usable, HiGHS otherwise.
#
# Gurobi is markedly faster on these MILPs and is what the papers used, but it
# is commercial. Requiring it would put a licence between an outside researcher
# and their first successful run -- which is the barrier this package exists to
# remove. So HiGHS is a hard dependency and always works, and Gurobi is picked
# up automatically when it is both installed and licensed.
#
# "Installed" and "licensed" are different failures and both are handled: a
# missing package throws on load, an absent or expired licence throws when the
# environment is constructed. Either way we fall back, once, with a notice.
# --------------------------------------------------------------------------- #

const GUROBI_UUID = Base.UUID("2e9cd046-0924-5485-92f1-d5272153d98b")

# Cached across calls: probing Gurobi means constructing an environment, which
# is slow and prints a banner. `:unprobed` -> not yet tried.
const _GUROBI_STATE = Ref{Symbol}(:unprobed)
const _GUROBI_MODULE = Ref{Any}(nothing)

"""
    load_gurobi!() -> Bool

Bring Gurobi's methods into scope without touching a licence. Whether it is
*usable* is [`gurobi_available`](@ref); this answers only whether it is
installed.

Must run from `__init__`, for pure world-age reasons: `Base.require` defines new
methods, and those are invisible to any call chain already running. Load and use
Gurobi inside one dynamic call and MOI's `is_empty` dispatch fails with
`MethodError: no method matching is_empty(::Gurobi.Optimizer)` -- which is
exactly what a first `solve_milp` used to do on a licensed machine. The licence
probe stays lazy; that is the slow part.
"""
function load_gurobi!()
    _GUROBI_MODULE[] === nothing || return true
    try
        _GUROBI_MODULE[] = Base.require(Base.PkgId(GUROBI_UUID, "Gurobi"))
        return true
    catch
        return false            # not installed; HiGHS covers us
    end
end

"""
    gurobi_available() -> Bool

Whether Gurobi can actually solve here: the package loads *and* a licence is
present. Probed once and cached.
"""
function gurobi_available()
    _GUROBI_STATE[] === :unprobed || return _GUROBI_STATE[] === :available

    _GUROBI_STATE[] = :unavailable
    if load_gurobi!()
        try
            # Constructing an Env is what actually validates the licence; loading
            # the package alone proves nothing.
            Base.invokelatest(getfield(_GUROBI_MODULE[], :Env))
            _GUROBI_STATE[] = :available
        catch
            # Licence missing, expired, or size-limited -- all the same to us.
        end
    end
    return _GUROBI_STATE[] === :available
end

"""
    select_optimizer(; prefer=:auto, verbose=false, time_limit=nothing, mip_gap=nothing)
        -> (optimizer_factory, name::Symbol)

Build a JuMP optimizer factory.

- `:auto` (default) uses Gurobi when available and HiGHS otherwise.
- `:gurobi` errors rather than silently falling back — use it when a benchmark
  result would be misleading if it quietly came from a different solver.
- `:highs` forces HiGHS even where Gurobi is available, which is what makes a
  Gurobi-vs-HiGHS comparison possible on the same machine.

Attributes are set through the solver-independent MOI interface, so the same
`time_limit` and `mip_gap` mean the same thing to both solvers.
"""
function select_optimizer(; prefer::Symbol=:auto,
    verbose::Bool=false,
    time_limit::Union{Nothing,Real}=nothing,
    mip_gap::Union{Nothing,Real}=nothing)

    prefer in (:auto, :gurobi, :highs) ||
        error("prefer must be :auto, :gurobi or :highs; got :$prefer.")

    name = if prefer === :gurobi
        gurobi_available() ||
            error("prefer=:gurobi but Gurobi is not usable here (package missing, or no valid licence). " *
                  "Use prefer=:auto to fall back to HiGHS.")
        :gurobi
    elseif prefer === :highs
        :highs
    else
        gurobi_available() ? :gurobi : :highs
    end

    base = if name === :gurobi
        # A zero-argument closure rather than the type itself. `gurobi_available`
        # may have loaded the package during *this* call, in which case its
        # bindings are not visible in the caller's world age and handing the
        # type straight to JuMP fails with "must be callable with zero
        # arguments". Deferring through invokelatest sidesteps that.
        let gurobi = _GUROBI_MODULE[]
            () -> Base.invokelatest(getfield(gurobi, :Optimizer))
        end
    else
        _warn_highs_once()
        HiGHS.Optimizer
    end

    return optimizer_with_attributes(base, _attributes(verbose, time_limit, mip_gap)...), name
end

"Solver-independent MOI attributes, so the same knobs mean the same thing to both solvers."
function _attributes(verbose, time_limit, mip_gap)
    # Must be `Pair[]`, not `Pair{Any,Any}[]`: MOI dispatches on the pair's
    # key type, and an Any-typed key fails its "optimizer attribute or string"
    # check even when the value it holds is a genuine attribute.
    attrs = Pair[]
    verbose || push!(attrs, MOI.Silent() => true)
    time_limit === nothing || push!(attrs, MOI.TimeLimitSec() => Float64(time_limit))
    mip_gap === nothing || push!(attrs, MOI.RelativeGapTolerance() => Float64(mip_gap))
    return attrs
end

const _HIGHS_NOTICE_SHOWN = Ref(false)

function _warn_highs_once()
    _HIGHS_NOTICE_SHOWN[] && return
    _HIGHS_NOTICE_SHOWN[] = true
    gurobi_available() && return          # deliberate :highs choice, not a fallback
    @info "Gurobi is unavailable (not installed, or no valid licence) -- solving with HiGHS. " *
          "Results are equivalent; expect longer solve times on large feeders."
end

"Reset the cached solver probe. Testing aid; not needed in normal use."
function reset_solver_cache!()
    _GUROBI_STATE[] = :unprobed
    _HIGHS_NOTICE_SHOWN[] = false
    # `_GUROBI_MODULE` is deliberately left alone: it is the loaded module, not
    # probe state, and dropping it would send the next call back through
    # `Base.require` from whatever call chain is running -- the world-age
    # failure `load_gurobi!` exists to prevent.
    return nothing
end
