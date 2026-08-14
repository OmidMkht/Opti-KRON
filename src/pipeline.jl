# --------------------------------------------------------------------------- #
# The whole method, end to end, as one call.
#
# Feeder -> operating point -> optimal reduction -> radial equivalent -> the
# accuracy actually achieved, including on scenarios the optimizer never saw.
# This is what `run_optikron.jl` drives, which leaves the runner as its options.
# --------------------------------------------------------------------------- #

"""
    Reduction

The result of a full run. `assignment` is the one to use downstream -- it is the
radialized map when radialization ran, and the raw one otherwise.
"""
struct Reduction
    network::Network
    V::Matrix{ComplexF64}
    solution::ReductionSolution
    assignment::Matrix{Float64}
    reinserted::Vector{Int}
    Ē::Float64
    scenarios::Vector{Int}
    radiality::Symbol
    devices::Vector{Device}
end

"""
    optikron(case; Ē, scenarios, hops, radiality, backend, preserve, V) -> Reduction

Reduce a feeder, start to finish.

`case` is a directory under `data/`, a path to one, or a `Network` in hand. `V`
defaults to [`powerflow`](@ref). `directory` says where a pre-loaded `Network`
came from, so its shipped `voltage.csv` and device tables are still found.

`radiality` is `:in_model` (default), `:post` (reinsert the minimal repair set,
Theorem 1) or `:none`. `:post` is not the default because it can leave the budget
behind: reinserting a bus changes `(A-I)C` and so `e`, and on case533mt at Ē=0.001
the worst violation moved from -1.4e-06 to +3.4e-05.

`backend` is `:milp` (exact, certified optimum) or `:search_cpu` / `:search_gpu`
(greedy, no hop limit, scales further). `solver` forces `:gurobi` / `:highs`
rather than `:auto`, worth doing for a benchmark.

`preserve` keeps equipment intact where the case ships `transformer.csv` /
`switch.csv`: `:devices` (default) pins transformers, regulators and switches;
`:transformers` merges across closed switches; `:none` pins nothing. The device
*model* is preserved exactly; the loading through it still moves.
"""
function optikron(case;
    Ē::Real=0.01,
    scenarios=Int[],
    hops::Int=5,
    radiality::Symbol=:in_model,
    backend::Symbol=:milp,
    solver::Symbol=:auto,
    preserve::Symbol=:devices,
    directory::Union{Nothing,AbstractString}=nothing,
    time_limit::Union{Nothing,Real}=3600,
    V::Union{Nothing,AbstractMatrix}=nothing,
    kwargs...)

    radiality in (:post, :in_model, :none) ||
        error("radiality must be :post, :in_model or :none; got :$radiality.")
    preserve in (:devices, :transformers, :none) ||
        error("preserve must be :devices, :transformers or :none; got :$preserve.")
    backend in (:milp, :search_cpu, :search_gpu) ||
        error("backend must be :milp, :search_cpu or :search_gpu; got :$backend.")

    network = case isa Network ? case : load_case(case)
    selected = collect(isempty(scenarios) ? (1:nscenarios(network)) : scenarios)

    # Shipped voltages and device tables live in the case directory, not in the
    # `Network`, so a caller that pre-loaded the feeder needs `directory` to keep
    # both.
    source = directory !== nothing ? case_directory(directory) :
             case isa AbstractString ? case_directory(case) : nothing

    # A case shipping its own voltages is stating its operating point, possibly
    # from a richer model than this package's power flow. Prefer it.
    shipped = (V === nothing && source !== nothing) ? read_voltage(source, network) : nothing
    operating_point = V !== nothing ? Matrix{ComplexF64}(V) :
                      shipped !== nothing ? shipped : powerflow(network)

    devices = (preserve !== :none && source !== nothing) ?
              read_devices(source, network; switches=(preserve === :devices)) : Device[]

    solution = if backend === :milp
        solve_milp(network, operating_point, Ē;
            scenarios=selected,
            hops=hops,
            pin=preserved_buses(devices),
            enforce_radiality=(radiality === :in_model),
            # Radiality can only be enforced when merges move towards the slack:
            # absorbing upstream can orphan a subtree the constraint cannot see.
            direction=(radiality === :in_model ? :downstream : :any),
            prefer=solver,
            time_limit=time_limit,
            kwargs...)
    else
        # No hop limit here -- the search merges neighbours and lets distance
        # accumulate through the chain -- and radiality comes from eliminating
        # only degree <= 2 buses rather than from a constraint.
        search_reduce(network, operating_point, Ē;
            scenarios=selected,
            radial=(radiality !== :none),
            upstream=(radiality === :in_model),
            pin=preserved_buses(devices),
            backend=(backend === :search_gpu ? :gpu : :cpu),
            time_limit=time_limit,
            kwargs...)
    end

    assignment, reinserted = (radiality === :post && backend === :milp) ?
                             radialize(network, solution.A) : (solution.A, Int[])

    return Reduction(network, operating_point, solution, assignment, reinserted,
        Float64(Ē), selected, radiality, devices)
end

"""
    load_case(name) -> Network

Load a case by directory name (resolved under `data/` when it is not a path),
picking the branch-based or Ybus reader by what the directory holds.
"""
function load_case(name::AbstractString)
    directory = case_directory(name)
    return isfile(joinpath(directory, "ybus.csv")) ?
           read_network_ybus(directory) : read_network_csv(directory)
end

"""
    list_cases([io]; root)

Print the installed cases and their sizes. Unit-test fixtures are skipped.
"""
function list_cases(io::IO=stdout; root::AbstractString=joinpath(dirname(@__DIR__), "data"))
    println(io, "Cases in ", root, ":\n")
    for name in sort(filter(d -> isdir(joinpath(root, d)), readdir(root)))
        startswith(name, "test") && continue
        net = load_case(joinpath(root, name))
        @printf(io, "  %-13s %5d buses  %5d node-phase rows  %3d scenario(s)  %s\n",
            name, nnodes(net), nphase_rows(net), nscenarios(net),
            is_three_phase(net) ? "three-phase" : "single-phase")
    end
end

"Resolve a case name to its directory, under `data/` when it is not a path."
function case_directory(name::AbstractString)
    directory = isdir(name) ? name : joinpath(dirname(@__DIR__), "data", name)
    isdir(directory) || error("No such case: $name (looked in $directory)")
    return directory
end

"Backend-specific line in a `Reduction` summary: each reports where its effort went."
function _solution_detail(io::IO, sol::MilpSolution)
    sol.screening === nothing && return
    s = sol.screening
    @printf(io, "Screening         %d of %d pairs kept, %d of %d error rows built\n",
        s.pairs_after, s.pairs_before, s.rows_kept, s.rows_kept + s.rows_skipped)
end

_solution_detail(io::IO, sol::SearchSolution) =
    @printf(io, "Search            %d merges, %d candidates scored, %d screened out\n",
        sol.iterations, sol.evaluated, sol.screened)

"""
    export_reduced(r::Reduction, dir; scenarios) -> Vector{String}

Write the reduction out in whichever form the feeder supports, and say so.

Single-phase feeders that came in with line parameters go out as MATPOWER, which
carries per-branch r/x. Everything else goes out as CSV, in the schema
`load_case` reads, plus `assignment.csv`. Three-phase feeders can only take the
second: a MATPOWER branch cannot say it carries a and c but not b.
"""
function export_reduced(r::Reduction, dir::AbstractString; scenarios=r.scenarios)
    if !is_three_phase(r.network) && has_branch_data(r.network)
        return export_matpower(r, dir; scenarios=scenarios)
    end
    return write_reduced_csv(dir, r.network, r.assignment; V=r.V, scenarios=scenarios)
end

"""
    export_matpower(r::Reduction, dir; scenarios, base_mva, base_kv) -> Vector{String}

Write the reduced network as MATPOWER case files, one per scenario.

Defaults to the scenarios the reduction was certified against: a file per loading
of a 168-scenario case is rarely what anyone means, and those are the ones the
equivalent is answerable for. Three-phase feeders are refused, see
[`write_matpower`](@ref).
"""
function export_matpower(r::Reduction, dir::AbstractString;
    scenarios=r.scenarios, base_mva::Real=1.0, base_kv::Real=1.0)

    mkpath(dir)
    selected = collect(scenarios)
    single = length(selected) == 1
    written = String[]

    for s in selected
        label = single ? "" : "_s$(s)"
        path = joinpath(dir, "$(r.network.name)_reduced$(label).m")
        write_matpower(path, r.network, r.assignment;
            scenario=s, base_mva=base_mva, base_kv=base_kv)
        push!(written, path)
    end
    return written
end

"Worst violation of the true annulus on the scenarios the model enforced."
enforced_violation(r::Reduction) =
    annulus_violation(r.network, r.assignment, r.V, r.Ē; scenarios=r.scenarios)

"Worst violation on the scenarios the model never saw; `NaN` when it saw all of them."
function held_out_violation(r::Reduction)
    rest = setdiff(1:nscenarios(r.network), r.scenarios)
    return isempty(rest) ? NaN :
           annulus_violation(r.network, r.assignment, r.V, r.Ē; scenarios=rest)
end

function Base.show(io::IO, ::MIME"text/plain", r::Reduction)
    net, sol = r.network, r.solution
    kept = length(super_nodes(r.assignment))
    println(io, net)
    @printf(io, "\nOperating point   |V| in [%.4f, %.4f] pu   (residual %.1e)\n",
        extrema(abs.(r.V))..., powerflow_residual(net, r.V))
    @printf(io, "Reduction         %d buses -> %d super-nodes (%.1f%%)  [%s, %s, %.2fs]\n",
        nnodes(net), length(sol.kept), 100 * sol.reduction, sol.status, sol.solver,
        sol.solve_time)

    _solution_detail(io, sol)
    if !isempty(r.devices)
        transformers = count(d -> d.kind === :transformer, r.devices)
        @printf(io, "Preserved         %d device(s) kept exactly (%d transformer, %d switch), %d buses pinned\n",
            length(r.devices), transformers, length(r.devices) - transformers,
            length(preserved_buses(r.devices)))
    end
    if r.radiality === :post
        @printf(io, "Radialization     +%d buses -> %d kept (%.1f%%)\n",
            length(r.reinserted), kept, 100 * reduction_ratio(r.assignment))
    end
    r.radiality === :none ||
        @printf(io, "                  radial: %s\n", is_radial(net, r.assignment))

    @printf(io, "\nAccuracy          budget Ē = %.4f pu\n", r.Ē)
    @printf(io, "                  %3d enforced scenarios: %+.2e\n",
        length(r.scenarios), enforced_violation(r))
    held_out = held_out_violation(r)
    isnan(held_out) ||
        @printf(io, "                  %3d held-out scenarios: %+.2e  (not certified)\n",
            nscenarios(net) - length(r.scenarios), held_out)
    print(io, "\nA[i,j] = 1 means bus j is represented by bus i.")
end
