# Reduced networks

Worked reductions of every feeder in [`data/`](../data), two per case: one
moderate, one aggressive. Each is a complete network in its own right — load it
and solve it like any other case.

```julia
net = load_case("reduced_networks/ieee123/57pct")
```

One folder per feeder, one subfolder per reduction level, named for the
reduction it achieves. Every level carries a `README.md` recording the error
threshold, hop limit and scenarios it was produced under, and the violation it
was re-checked at.

## What is here

| Case | Level | Reduction | Ē | Hops | Scenarios | PQ error | OpenDSS error |
|---|---|---|---|---|---|---|---|
| `ieee37` | [`67pct`](ieee37/67pct) | 39 → 13 | 0.003 | 10 | all | 2.9e-03 | 5.6e-04 |
| `ieee37` | [`77pct`](ieee37/77pct) | 39 → 9 | 0.008 | 10 | all | 7.1e-03 | 1.4e-03 |
| `ieee34` | [`51pct`](ieee34/51pct) | 37 → 18 | 0.003 | 10 | all | 2.6e-03 | 1.7e-04 |
| `ieee34` | [`62pct`](ieee34/62pct) | 37 → 14 | 0.01 | 10 | all | 6.6e-03 | 5.0e-04 |
| `ieee123` | [`57pct`](ieee123/57pct) | 130 → 56 | 0.002 | 10 | all | 2.0e-03 | 9.1e-04 |
| `ieee123` | [`78pct`](ieee123/78pct) | 130 → 29 | 0.006 | 10 | all | 5.9e-03 | 3.9e-03 |
| `european_lv` | [`91pct`](european_lv/91pct) | 907 → 82 | 0.01 | 10 | all | 5.2e-04 | 1.5e-04 |
| `european_lv` | [`98pct`](european_lv/98pct) | 907 → 21 | 0.01 | 25 | all | 1.2e-03 | 8.8e-04 |
| `case69` | [`59pct`](case69/59pct) | 69 → 28 | 0.001 | 10 | all | 9.9e-04 | — |
| `case69` | [`77pct`](case69/77pct) | 69 → 16 | 0.006 | 10 | all | 5.8e-03 | — |
| `case85` | [`59pct`](case85/59pct) | 85 → 35 | 0.002 | 10 | all | 1.9e-03 | — |
| `case85` | [`78pct`](case85/78pct) | 85 → 19 | 0.006 | 10 | all | 5.9e-03 | — |
| `case141` | [`74pct`](case141/74pct) | 141 → 37 | 0.001 | 10 | all | 9.9e-04 | — |
| `case141` | [`85pct`](case141/85pct) | 141 → 21 | 0.003 | 10 | all | 3.0e-03 | — |
| `case533mt` | [`66pct`](case533mt/66pct) | 533 → 180 | 0.001 | 10 | all | 1.0e-03 | — |
| `case533mt` | [`83pct`](case533mt/83pct) | 533 → 88 | 0.003 | 10 | all | 3.0e-03 | — |
| `case1197` | [`50pct`](case1197/50pct) | 1197 → 593 | 0.006 | 10 | all | 6.0e-03 | — |
| `case1197` | [`65pct`](case1197/65pct) | 1197 → 417 | 0.009 | 10 | all | 9.0e-03 | — |
| `R100` | [`82pct`](R100/82pct) | 100 → 18 | 0.001 | 10 | 67, 91 | 9.5e-04 | — |
| `R100` | [`91pct`](R100/91pct) | 100 → 9 | 0.005 | 10 | 67, 91 | 4.8e-03 | — |
| `R300` | [`84pct`](R300/84pct) | 300 → 48 | 0.001 | 10 | 67, 91 | 9.5e-04 | — |
| `R300` | [`91pct`](R300/91pct) | 300 → 26 | 0.003 | 10 | 67, 91 | 3.0e-03 | — |
| `ieee8500` | [`31pct`](ieee8500/31pct) | 4876 → 3346 | 0.005 | none | all | 5.0e-03 | 3.0e-03 |
| `ieee8500` | [`41pct`](ieee8500/41pct) | 4876 → 2856 | 0.02 | none | all | 1.6e-02 | 2.5e-03 |

All but `ieee8500` are proven optimal by the MILP. `ieee8500` does not close at
4876 buses, so its two levels come from the GPU exhaustive search — greedy, no
hop limit, same error budget, 269s and 236s. `dss_build.json` records the
OpenDSS build.

Every level is reduced with `preserve = :required` — **center-tapped
transformers, phase shifters, regulators and switches** kept intact. Those four
are what an OpenDSS form needs to mean anything: a center-tapped transformer or
phase shifter cannot be rebuilt from a Schur complement at all, and a folded
regulator or switch freezes at its present tap or position. Plain transformers
and capacitor banks are *not* pinned; an ordinary two-winding crossing is
recoverable by rebasing each winding's kV, which is what lets `ieee37` reduce to
13 buses rather than 15.

Nothing here is mandatory. `preserve` accepts `:all`, `:none` or any explicit
tuple, and the reduction stays valid and certified under every choice — only the
OpenDSS form suffers, and [`run_optikron.jl`](../run_optikron.jl) warns about
that when you ask for one.

**Regenerating any of these** takes one edit: each level's note carries the exact
option block that reproduces it, checked byte-for-byte against the shipped files.
`ieee8500` is the one case whose `ybus.csv` is not shipped (~1.2 MB a level);
everything else is complete in the repository.

## Every one is verified

All twenty-four are radial and inside their stated budget on the *exact
nonconvex* annulus — re-derived from the files against the untouched sparse
`Ybus`, not against the linearisation the solver used. Worst margins run from
`-1.9e-07` on `case533mt/66pct` to `-9.5e-03` on `european_lv/91pct`; each
level's own note carries its number.

Three qualifications worth reading before using one:

- **`R100` and `R300` are certified on two hours out of 168.** Across the other
  166 the error reaches 3.6× to 4.6× `Ē`. That is the expected cost of fitting to
  representative loadings, and it is why the notes report enforced and held-out
  separately.
- **`european_lv` does not respond to the budget.** It reduces identically
  anywhere in Ē = 0.001–0.01, because its annulus never binds — geometry and the
  hop limit decide the answer. Its two examples are therefore separated by `hops`,
  not by error threshold.
- **All ten OpenDSS round trips now solve inside budget**, after two converter
  corrections (delta loads staying delta, and per-feeder load representation);
  see [`converter/README.md`](../converter/README.md#the-round-trip-and-how-it-is-solved).

## OpenDSS form

Where a feeder came from OpenDSS, its reduction also ships as a self-contained,
solved OpenDSS circuit in `dss/` — `Master.dss`, the element files, and
`validation.csv` giving the solved-versus-original voltage at every surviving
super-node. All ten land inside their budget; see
[`converter/README.md`](../converter/README.md#the-round-trip-and-how-it-is-solved).

## Formats

Single-phase feeders that arrived with line parameters go out as MATPOWER `.m`,
with the bus map in the header comment. Three-phase feeders have no MATPOWER form
— a branch there cannot carry `a` and `c` but not `b` — so they go out as CSV:
`bus.csv`, `ybus.csv` (the Schur complement), `load.csv` (injections after each
eliminated bus hands its load over), `voltage.csv`, and `assignment.csv`, the
reduction map that makes the result interpretable rather than merely smaller.
