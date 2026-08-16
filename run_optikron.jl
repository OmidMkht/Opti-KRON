#!/usr/bin/env julia
# --------------------------------------------------------------------------- #
# Opti-KRON runner.
#
#     julia --project=. run_optikron.jl              # the options below
#     julia --project=. run_optikron.jl ieee123      # another case
#     julia --project=. run_optikron.jl ieee123 0.01 # ...and another Ē
#     julia --project=. run_optikron.jl --list       # what cases are installed
#
# The method itself is `optikron` in src/pipeline.jl; this file is only its
# options.
# --------------------------------------------------------------------------- #

# ------------------------------- OPTIONS ----------------------------------- #
case      = "ieee8500"   # any case in data/ (see --list), or a path to a folder
Ē         = 0.001       # allowed voltage-magnitude error, per unit
scenarios = []          # loadings to enforce it on ([] = all)
hops      = 5          # how far a bus may travel to its super-node
radiality = :in_model   # :in_model (a model constraint) | :post (repair after) | :none
preserve  = :both       # :both | :switches | :regulators | :none
contract  = false       # short-circuit closed switches before reducing, then expand back
backend   = :milp       # :milp | :search_cpu | :search_gpu
solver    = :auto       # :auto (Gurobi if licensed, else HiGHS) | :gurobi | :highs
export_to = ""          # folder to write the reduced network to ("" = none)
time_limit = 3600       # seconds the solver may spend (nothing for no limit)
verbose   = true        # stage timing, the search's merge trace, and the
                        # Gurobi/HiGHS log; false runs silent
# --------------------------------------------------------------------------- #

using Printf

# A large case spends minutes inside a single call with nothing on stdout, which
# is indistinguishable from a hang. Every step long enough to worry about is
# announced before it starts, with the clock since launch.
const T0 = time()
stage(msg) = verbose && (@printf("[%6.1fs] %s\n", time() - T0, msg); flush(stdout))

# Loaded as the project's own package rather than `include`d. Including the
# source builds a second copy of the module, and then a session that has also
# done `using OptiKRON` carries two of every name and two of every type -- which
# surfaces as ambiguous bindings, or worse as `Network` not matching `Network`.
using OptiKRON

"--list" in ARGS && (list_cases(); exit(0))
length(ARGS) >= 1 && (case = ARGS[1])
length(ARGS) >= 2 && (Ē = parse(Float64, ARGS[2]))

# Loaded here so `scenarios` can be checked against the case: naming scenario 67
# of a case that has one is a confusing way to fail. `directory` is then what
# lets the loaded network still find its voltage.csv and device tables.
stage("Loading $case ...")
network = load_case(case)
scenarios = filter(s -> s <= nscenarios(network), scenarios)
isempty(scenarios) && (scenarios = collect(1:nscenarios(network)))
stage("  $(nnodes(network)) buses, $(nphase_rows(network)) node-phase rows, " *
      "$(nscenarios(network)) scenario(s)")

stage("Reducing: backend=$backend, Ē=$Ē, hops=$hops, preserve=$preserve ...")
elapsed = @elapsed result = optikron(network; directory=case, Ē, scenarios, hops,
    radiality, preserve, contract, backend, solver, time_limit, verbose)

stage("Reduced in $(round(elapsed, digits=1))s; checking accuracy " *
      "(a dense solve -- seconds to a minute on a large case) ...")
show(stdout, MIME"text/plain"(), result)

# `enforced_violation` reports the margin, `deviation - Ē`, so it is negative
# whenever the reduction is inside budget. The error itself is that plus Ē --
# printed here because how much of the budget was actually spent is what says
# whether a bigger Ē would buy more reduction.
worst = enforced_violation(result) + result.Ē
@printf("\nMax error         %.3e pu of Ē = %.3e  (%.1f%% of budget)\n",
    worst, result.Ē, 100 * worst / result.Ē)
held = held_out_violation(result)
isnan(held) || @printf("Max error         %.3e pu on held-out scenarios\n", held + result.Ē)

@printf("\nTotal             %.1fs (power flow, reduction and radiality)\n", elapsed)

if !isempty(export_to)
    stage("Writing to $export_to ...")
    foreach(f -> println("Wrote ", f), export_reduced(result, export_to))
end

# `result.assignment` is the reduction map and `result.network` the feeder it
# came from; see src/pipeline.jl for the other fields.
