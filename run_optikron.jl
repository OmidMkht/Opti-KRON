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
case      = "ieee123"   # any case in data/ (see --list), or a path to a folder
Ē         = 0.001       # allowed voltage-magnitude error, per unit
scenarios = []          # loadings to enforce it on ([] = all)
hops      = 10          # how far a bus may travel to its super-node
radiality = :in_model   # :in_model (a model constraint) | :post (repair after) | :none
preserve  = :devices    # keep equipment intact: :devices | :transformers | :none
backend   = :milp       # :milp (exact) | :search_cpu | :search_gpu (greedy, scales)
solver    = :auto       # :auto (Gurobi if licensed, else HiGHS) | :gurobi | :highs
export_to = ""          # folder to write the reduced network to ("" = none)
time_limit = 3600       # seconds the solver may spend (nothing for no limit)
# --------------------------------------------------------------------------- #

using Printf
include(joinpath(@__DIR__, "src", "OptiKRON.jl"))
using .OptiKRON

"--list" in ARGS && (list_cases(); exit(0))
length(ARGS) >= 1 && (case = ARGS[1])
length(ARGS) >= 2 && (Ē = parse(Float64, ARGS[2]))

# Loaded here so `scenarios` can be checked against the case: naming scenario 67
# of a case that has one is a confusing way to fail. `directory` is then what
# lets the loaded network still find its voltage.csv and device tables.
network = load_case(case)
scenarios = filter(s -> s <= nscenarios(network), scenarios)
isempty(scenarios) && (scenarios = collect(1:nscenarios(network)))

elapsed = @elapsed result = optikron(network; directory=case, Ē, scenarios, hops,
    radiality, preserve, backend, solver, time_limit)
show(stdout, MIME"text/plain"(), result)
@printf("\n\nTotal             %.1fs (power flow, reduction and radiality)\n", elapsed)

isempty(export_to) || foreach(f -> println("Wrote ", f), export_reduced(result, export_to))

# `result.assignment` is the reduction map and `result.network` the feeder it
# came from; see src/pipeline.jl for the other fields.
