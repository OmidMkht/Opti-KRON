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

`ieee8500` ships as **map only** — `assignment.csv` plus the bus, load and
voltage tables, without `ybus.csv`. Kron reduction fills in, and at this size the
Schur complement runs to 289 MB and 631 MB. Measured, 99.5% of those entries are
round-off: only 28312 of 5390534 exceed 1e-12 of the largest, and the count is
unchanged at 1e-10 and 1e-8. Each level's note carries the one-liner that rebuilds
it with `kron_reduce`.

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
- **A solved OpenDSS round trip is a harder test than the budget**, and on
  `ieee37` it is much harder — see below.

## OpenDSS form

Where a feeder came from OpenDSS, its reduction also ships as a self-contained,
solved OpenDSS circuit in `dss/`, rebuilt by
[`converter/build_reduced_dss.py`](../converter/build_reduced_dss.py). Each one
carries `Master.dss`, the element files, `bus_mapping.csv`, `report.json`, and
`validation.csv` giving the solved-versus-original voltage at every surviving
super-node.

That solved error is **not** the certified budget. The budget holds a
constant-current linearisation; OpenDSS re-solves with constant-power loads and
the real winding connections:

| Case | Level | Ē | Solved DSS error | |
|---|---|---|---|---|
| `european_lv` | 91pct | not binding | 1.6e-04 | well inside |
| `european_lv` | 98pct | not binding | 5.1e-04 | well inside |
| `ieee123` | 57pct | 0.002 | 9.4e-04 | inside |
| `ieee123` | 78pct | 0.006 | 5.2e-03 | inside |
| `ieee34` | 51pct | 0.003 | 9.4e-03 | **3× over** |
| `ieee34` | 59pct | 0.010 | 1.0e-02 | at the edge |
| `ieee37` | 62pct | 0.003 | 8.4e-02 | **28× over** |
| `ieee37` | 74pct | 0.008 | 8.4e-02 | **11× over** |

The conversion itself is exact — power aggregation error `0.0` and a synthesized
Kron `Ybus` matching to `2.9e-16` on `ieee37`, the *worst* row of the table. The
gradient tracks how far each feeder sits from the constant-current assumption:
`european_lv` is 98% unloaded, `ieee123` is wye-connected, and `ieee37` is the
one all-delta IEEE feeder, where load responds to line-to-line voltage.

### Why `ieee8500` has no `dss/`

A known limitation, not a failure to try. The reduction absorbs `ieee8500`'s 2354
service transformers rather than pinning them — pinning would fix 2400 buses, 49%
of the feeder, as a hard ceiling. The OpenDSS converter needs both terminals of a
preserved transformer to survive, so a reduction that has eliminated them cannot
be re-instantiated as a circuit.

Preserving them makes the conversion possible and costs real reduction, measured
on this feeder with the GPU search:

| Ē | absorbed (shipped) | transformers pinned |
|---|---|---|
| 0.005 | 42% | 31% |
| 0.020 | 54% | 41% |

Eleven to thirteen points. The shipped levels keep the reduction; run
`optikron("ieee8500"; preserve = (:phase_shift, :regulator, :switch, :transformer),
backend = :search_gpu)` if you would rather have the circuit.

## Formats

Single-phase feeders that arrived with line parameters go out as MATPOWER `.m`,
with the bus map in the header comment. Three-phase feeders have no MATPOWER form
— a branch there cannot carry `a` and `c` but not `b` — so they go out as CSV:
`bus.csv`, `ybus.csv` (the Schur complement), `load.csv` (injections after each
eliminated bus hands its load over), `voltage.csv`, and `assignment.csv`, the
reduction map that makes the result interpretable rather than merely smaller.
