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

| Case | Level | Reduction | Ē | Hops | Scenarios | Format | OpenDSS |
|---|---|---|---|---|---|---|---|
| `ieee37` | [`62pct`](ieee37/62pct) | 39 → 15 | 0.003 | 10 | 1 of 1 | CSV | yes |
| `ieee37` | [`74pct`](ieee37/74pct) | 39 → 10 | 0.008 | 10 | 1 of 1 | CSV | yes |
| `ieee34` | [`51pct`](ieee34/51pct) | 37 → 18 | 0.003 | 10 | 1 of 1 | CSV | yes |
| `ieee34` | [`59pct`](ieee34/59pct) | 37 → 15 | 0.010 | 10 | 1 of 1 | CSV | yes |
| `ieee123` | [`57pct`](ieee123/57pct) | 130 → 56 | 0.002 | 10 | 1 of 1 | CSV | yes |
| `ieee123` | [`78pct`](ieee123/78pct) | 130 → 29 | 0.006 | 10 | 1 of 1 | CSV | yes |
| `european_lv` | [`91pct`](european_lv/91pct) | 907 → 82 | not binding | 10 | 1 of 1 | CSV | yes |
| `european_lv` | [`98pct`](european_lv/98pct) | 907 → 21 | not binding | 25 | 1 of 1 | CSV | yes |
| `ieee8500` | [`42pct`](ieee8500/42pct) | 4876 → 2821 | 0.005 | none | 1 of 1 | map only | — |
| `ieee8500` | [`54pct`](ieee8500/54pct) | 4876 → 2224 | 0.020 | none | 1 of 1 | map only | — |
| `case69` | [`59pct`](case69/59pct) | 69 → 28 | 0.001 | 10 | 1 of 1 | MATPOWER | — |
| `case69` | [`77pct`](case69/77pct) | 69 → 16 | 0.006 | 10 | 1 of 1 | MATPOWER | — |
| `case85` | [`59pct`](case85/59pct) | 85 → 35 | 0.002 | 10 | 1 of 1 | MATPOWER | — |
| `case85` | [`78pct`](case85/78pct) | 85 → 19 | 0.006 | 10 | 1 of 1 | MATPOWER | — |
| `case141` | [`74pct`](case141/74pct) | 141 → 37 | 0.001 | 10 | 1 of 1 | MATPOWER | — |
| `case141` | [`85pct`](case141/85pct) | 141 → 21 | 0.003 | 10 | 1 of 1 | MATPOWER | — |
| `case533mt` | [`66pct`](case533mt/66pct) | 533 → 180 | 0.001 | 10 | 2 of 2 | MATPOWER | — |
| `case533mt` | [`83pct`](case533mt/83pct) | 533 → 88 | 0.003 | 10 | 2 of 2 | MATPOWER | — |
| `case1197` | [`50pct`](case1197/50pct) | 1197 → 593 | 0.006 | 10 | 1 of 1 | MATPOWER | — |
| `case1197` | [`65pct`](case1197/65pct) | 1197 → 417 | 0.009 | 10 | 1 of 1 | MATPOWER | — |
| `R100` | [`82pct`](R100/82pct) | 100 → 18 | 0.001 | 10 | 67, 91 of 168 | CSV | — |
| `R100` | [`91pct`](R100/91pct) | 100 → 9 | 0.005 | 10 | 67, 91 of 168 | CSV | — |
| `R300` | [`84pct`](R300/84pct) | 300 → 48 | 0.001 | 10 | 67, 91 of 168 | CSV | — |
| `R300` | [`91pct`](R300/91pct) | 300 → 26 | 0.003 | 10 | 67, 91 of 168 | CSV | — |

All but `ieee8500` are proven optimal by the MILP. `ieee8500` does not close at
4876 buses, so its two levels come from the GPU exhaustive search — greedy, no
hop limit, same error budget, 263s and 303s. `dss_build.json` records the
OpenDSS build.

`ieee8500` ships as **map only** — no `ybus.csv`. Kron fill-in runs to 289–631 MB
there, 99.5% round-off; each level's note carries the one-liner that rebuilds it
with `kron_reduce`.

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
- **A solved OpenDSS round trip is a harder test than the budget** — on `ieee37`
  it runs 28× over. Not a conversion defect; see
  [`converter/README.md`](../converter/README.md#the-round-trip-is-a-harder-test-than-the-budget).

## OpenDSS form

Where a feeder came from OpenDSS, its reduction also ships as a self-contained,
solved OpenDSS circuit in `dss/` — `Master.dss`, the element files, and
`validation.csv` giving the solved-versus-original voltage at every surviving
super-node. `ieee8500` has none; see
[`converter/README.md`](../converter/README.md#the-ieee8500-limitation) for why
and what it costs to get one anyway.

## Formats

Single-phase feeders that arrived with line parameters go out as MATPOWER `.m`,
with the bus map in the header comment. Three-phase feeders have no MATPOWER form
— a branch there cannot carry `a` and `c` but not `b` — so they go out as CSV:
`bus.csv`, `ybus.csv` (the Schur complement), `load.csv` (injections after each
eliminated bus hands its load over), `voltage.csv`, and `assignment.csv`, the
reduction map that makes the result interpretable rather than merely smaller.
