# Opti-KRON

Optimal Kron-based network reduction for distribution feeders — single-phase and unbalanced three-phase.

> **Status: pre-release.** The input layer is in place and tested. The solvers are
> being ported from the research code. See [Roadmap](#roadmap) for what works today.

## What problem does Opti-KRON solve?

Detailed feeder models with thousands of nodes make optimal power flow, hosting-capacity
and control studies expensive. Network reduction shrinks them — but classical Kron reduction
can only eliminate zero-injection nodes, and most other methods require you to decide up
front which buses to keep.

Opti-KRON chooses the kept buses **for** you, by optimization, subject to an explicit bound
on voltage error. You give it a feeder and a voltage-error budget `E̅`; it returns the most
reduced network that respects that budget, plus the map saying which original bus each load
and device now lives at.

## Why use it?

- **Explicit error bound.** Reduction is constrained by a user-set per-unit voltage error, not assessed after the fact.
- **Interpretable aggregation.** A binary assignment matrix maps every reduced bus to a specific super-node, so loads, DERs and controllable devices keep known locations.
- **Radiality preserved**, two ways: enforced inside the optimization, or repaired afterwards by reinserting a minimal set of nodes (Theorem 1). A radial equivalent also runs *faster* downstream, because of sparsity.
- **Unbalanced three-phase supported**, not just balanced single-phase equivalents.

## Two solver families

| | What it is | When to use it |
|---|---|---|
| `src/optimization/` | MILP (JuMP), exact | Any feeder, single- or three-phase; when you want a certified optimum |
| `src/search/` | Exhaustive search at q=1, CPU and GPU | Large networks; nonconvex voltage-magnitude objectives |

The exhaustive search operates directly on complex values, so it evaluates the exact
voltage-magnitude error the MILP can only approximate through rectangular decomposition.
On a 1000-node feeder the GPU version runs ~16× faster than CPU.

## Solvers

The MILP runs on **Gurobi** when a licence is present and falls back to **HiGHS**
otherwise — no configuration needed. HiGHS is a hard dependency, so the package
works out of the box without a commercial licence; Gurobi is optional and picked
up automatically.

```julia
factory, name = select_optimizer()                  # :gurobi if licensed, else :highs
factory, name = select_optimizer(prefer = :highs)   # force the open solver
factory, name = select_optimizer(prefer = :gurobi)  # error rather than fall back
```

## Input format

One format, both solvers. Three CSVs in a directory:

```
bus.csv      bus_id, phases, base_kv, type          # phases: "abc", "a", "ac"; exactly one type=slack
branch.csv   from_bus, to_bus, phases, r, x, b      # or r_aa..r_cc, x_aa..x_cc, b_aa..b_cc for 3-phase
load.csv     bus_id, phase, scenario, p_pu, q_pu    # injections: generation positive, load negative
```

```julia
include("src/OptiKRON.jl"); using .OptiKRON
net = read_network_csv("data/test4")
# Network("test4", 4 buses, 4 node-phase rows, single-phase, 2 scenario(s))
```

Branch parameters are carried through — rather than only a `Ybus` — because a `Ybus` cannot
be inverted back to per-line `r`/`x` once shunts, transformers or regulators are present.
Without them a reduced network can only be emitted as another `Ybus`, which is far less
useful to someone who wants to open it in their own tool.

When branch data genuinely isn't available, there's a fast path that gives up feeder export
but is otherwise fully supported:

```julia
net = read_network_ybus("path/to/case")             # bus.csv + ybus.csv (sparse triplets)
net = network_from_matrices(Ybus, S; slack=1)       # straight from memory
```

## Quick start

```
julia --project=. run_optikron.jl              # reduce R100 with the defaults
julia --project=. run_optikron.jl R300 0.002   # another case, another Ē
```

[`run_optikron.jl`](run_optikron.jl) is *only* an options block — case, `Ē`,
scenarios, hops, radiality mode, backend. The method itself is `optikron` in
[`src/pipeline.jl`](src/pipeline.jl), which loads the feeder, solves its power
flow, reduces it, recovers radiality and reports the accuracy actually achieved,
including on the scenarios the optimizer never saw.

```julia
result = optikron("R100"; Ē = 0.001, scenarios = [67, 91], hops = 5,
                  radiality = :in_model)
result.assignment      # the reduction map, radialized
```

## Reducing a feeder

```julia
net = read_network_ybus("data/R100")
V   = powerflow(net)                    # operating point: one column per scenario

sol = solve_milp(net, V, 0.001;
                 scenarios = [67, 91],  # enforce the budget on representative loadings
                 hops = 5)      # how far a bus may travel to its super-node

sol.kept          # the super-nodes
sol.reduction     # 0.88  -- fraction of buses eliminated
sol.A             # assignment matrix: A[i,j] = 1 means bus j is represented by i

annulus_violation(net, sol.A, V, 0.001) # <= 0 means the budget really holds
A_radial, added = radialize(net, sol.A) # recover a radial equivalent
```

The error budget is a *hard constraint*, not a penalty term: the solver returns
the smallest network that respects it, or reports infeasibility.

It is a hard constraint on the *linearized* annulus, though, and that distinction
is why `annulus_violation` exists. The annulus is nonconvex, both linearizations
bound the aligned error rather than the magnitude itself, and neither is robust to
several merges landing at once — so an assignment can be optimal, at zero gap, and
still miss the true constraint by a little. Measured across the shipped feeders it
is rare and small: three runs in two hundred, all on `case141`, all under 0.4% of
the budget. `annulus_violation` re-derives the original nonconvex constraint from
the untouched sparse `Ybus` and is the check that actually certifies a reduction —
run it. Passing it the scenarios the model did *not* see also tells you how a
reduction fitted to representative loadings holds up across the rest.

`V` is passed in rather than computed inside the solver, because the operating
point is a modelling choice — a feeder with delta connections, ZIP loads or
regulator taps belongs in a tool that models them, with those voltages handed
over directly. `powerflow_residual(net, V)` checks that any `V` you bring is
actually consistent with `Ybus` and `S`; one that isn't will still solve, but
the error bound it produces means nothing.

## Included test cases

Two radial three-phase feeders ship in [`data/`](data/), ready to reduce:

| Case | Buses | Node-phase rows | Scenarios | Phasing |
|---|---|---|---|---|
| `R100` | 100 | 166 | 168 hourly | unbalanced three-phase, radial |
| `R300` | 300 | 506 | 168 hourly | unbalanced three-phase, radial |

Both are genuinely unbalanced — a three-phase backbone with single-phase
laterals — so they exercise the phase-availability rules rather than behaving
like single-phase equivalents. See [`data/README.md`](data/README.md) for the
full description.

Ten further benchmark networks ship alongside them: five single-phase MATPOWER
feeders (`case69`, `case85`, `case141`, `case533mt`, `case1197`) and five
published three-phase feeders (`ieee34`, `ieee37`, `ieee123`, `european_lv`,
`ieee8500`) carrying the operating point they were published with, along with
their transformer, regulator, phase-shift, switch and capacitor tables. See
[`data/README.md`](data/README.md) for all of them.

Worked reductions live in
[`reduced_networks/`](reduced_networks) — two per feeder, each recording the
error threshold, hop limit and scenarios it was produced under, and each verified
against the exact nonconvex annulus.

## Roadmap

- [x] `Network` type, phase-aware, shared by both solver families
- [x] CSV import (branch-based), Ybus fast path, in-memory constructor, import validation
- [x] Core: radial topology, Kron reduction, radialization
- [x] MILP solver — one phase-aware model, single- and three-phase
- [x] Constraint screening — provably-redundant pairs and rows removed up front
- [x] Zero-injection warm start — the exact part of the reduction, solved first
- [x] Radiality both ways — repaired after (Theorem 1), or enforced in the model
- [x] Core: AC power flow, operating-point consistency check, error metrics
- [x] `run_optikron.jl`: feeder → reduction → radial equivalent → accuracy
- [x] MATPOWER import and export — reduce a `.m` case, get a `.m` case back
- [x] Twelve benchmark feeders, single- and three-phase, up to 4876 buses
- [x] Exhaustive search — CPU and GPU backends, same error model as the MILP
      (GPU is ~25x on a 907-bus feeder; CUDA optional, absent falls back to CPU)
- [x] Preserve transformers, regulators, phase shifters and switches through a
      reduction, selectable per kind
- [x] Re-export the IEEE cases at ~14 decimals, then regenerate their reduced
      examples

## Citation

If you use the Opti-KRON methodology, please cite the journal paper:

> O. Mokhtari, S. Chevalier, and M. Almassalkhi, "Optimal Kron-based Reduction of Networks
> (Opti-KRON) for Three-phase Distribution Feeders," *Electric Power Systems Research*,
> vol. 263, art. 113615, 2027.
> [doi:10.1016/j.epsr.2026.113615](https://doi.org/10.1016/j.epsr.2026.113615) ·
> [arXiv:2510.19608](https://arxiv.org/abs/2510.19608)

For the single-phase MILP and the radialization step:

> O. Mokhtari, S. Chevalier, and M. Almassalkhi, "Structure-preserving Optimal Kron-based
> Reduction of Radial Distribution Networks," arXiv preprint, 2025.
> [arXiv:2508.15006](https://arxiv.org/abs/2508.15006)

## License

MIT — see [LICENSE](LICENSE). You may use, modify and redistribute this freely,
including commercially, provided the copyright notice is retained.

The license does not require citation; that is a scholarly norm, and
[CITATION.cff](CITATION.cff) records how to do it.

## Contact

Omid Mokhtari — omid.mokhtari@uvm.edu
Department of Electrical and Biomedical Engineering, University of Vermont
