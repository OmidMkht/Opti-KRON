#!/usr/bin/env julia
# --------------------------------------------------------------------------- #
# Opti-KRON runner: edit the options, run it.
#
#     julia --project=. run_optikron.jl              # the options below
#     julia --project=. run_optikron.jl ieee123 0.01 # override case and Ē
#     julia --project=. run_optikron.jl --list       # what cases are installed
#
# The method itself is `optikron` in src/pipeline.jl.
# --------------------------------------------------------------------------- #

# ------------------------------- OPTIONS ----------------------------------- #
case      = "R100"      # any case in data/, or a path to a folder
Ē         = 0.001       # allowed voltage-magnitude error, per unit
scenarios = []          # loadings to enforce it on ([] = all)
hops      = 5           # how far a bus may travel to its super-node
radiality = :in_model   # :in_model | :post | :none
preserve  = :required   # :required | :all | :switches | :regulators | :none | tuple
backend   = :milp       # :milp | :search_cpu | :search_gpu
contract  = false       # short-circuit closed switches, then expand back
solver    = :auto       # :auto | :gurobi | :highs
time_limit = 3600       # seconds the solver may spend (nothing for no limit)

export_to = ""          # folder for the reduced network ("" = none)
export_dss = false      # also write an OpenDSS circuit (needs an OpenDSS origin)
verbose   = true        # stage timing and the solver log
# --------------------------------------------------------------------------- #

using Printf, OptiKRON

const T0 = time()
stage(msg) = verbose && (@printf("[%6.1fs] %s\n", time() - T0, msg); flush(stdout))

"--list" in ARGS && (list_cases(); exit(0))
length(ARGS) >= 1 && (case = ARGS[1])
length(ARGS) >= 2 && (Ē = parse(Float64, ARGS[2]))

# Loaded up front so `scenarios` can be checked against the case: naming scenario
# 67 of a case that has one is a confusing way to fail.
stage("Loading $case ...")
network = load_case(case)
scenarios = filter(s -> s <= nscenarios(network), scenarios)
isempty(scenarios) && (scenarios = collect(1:nscenarios(network)))
stage("  $(nnodes(network)) buses, $(nphase_rows(network)) rows, " *
      "$(nscenarios(network)) scenario(s)")

# Dropping a required kind is allowed -- the reduction is still valid and still
# certified -- but the OpenDSS writer cannot honour equipment it was told to fold
# away, so say so before spending the solve rather than after.
absent = missing_for_dss(preserve)
export_dss && !isempty(absent) && @warn """
    preserve=$preserve leaves $(join(absent, ", ")) unpinned.
    The reduction stays valid and inside its budget, but the OpenDSS circuit will
    have absorbed that equipment: a center-tapped transformer or phase shifter
    cannot be rebuilt from a Schur complement at all, and a folded regulator or
    switch freezes at its present tap or position. Use preserve=:required for a
    circuit whose equipment still means what it did."""

stage("Reducing: backend=$backend, Ē=$Ē, hops=$hops, preserve=$preserve ...")
elapsed = @elapsed result = optikron(network; directory=case, Ē, scenarios, hops,
    radiality, preserve, contract, backend, solver, time_limit, verbose)

stage("Reduced in $(round(elapsed, digits=1))s; checking accuracy " *
      "(a dense solve -- up to a minute on a large case) ...")
show(stdout, MIME"text/plain"(), result)

# `enforced_violation` is the margin `deviation - Ē`, negative inside budget. The
# error itself is that plus Ē -- how much of the budget was spent is what says
# whether a larger Ē would buy more reduction.
worst = enforced_violation(result) + result.Ē
@printf("\nMax error         %.3e pu of Ē = %.3e  (%.1f%% of budget)\n",
    worst, result.Ē, 100 * worst / result.Ē)
held = held_out_violation(result)
isnan(held) || @printf("Max error         %.3e pu on held-out scenarios\n", held + result.Ē)
@printf("\nTotal             %.1fs\n", elapsed)

if !isempty(export_to)
    stage("Writing to $export_to ...")
    foreach(f -> println("Wrote ", f), export_reduced(result, export_to))

    # OpenDSS is driven from Python, so this shells out to converter/ rather than
    # linking it in. Only feeders that came from OpenDSS can take it; the script
    # says so and skips the rest.
    if export_dss
        stage("Writing OpenDSS ...")
        script = joinpath(@__DIR__, "converter", "build_reduced_dss.py")
        run(`python $script --reduced $export_to --case $case`)
    end
end
