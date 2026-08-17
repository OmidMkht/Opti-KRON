# R300 — 84pct

300 → 48 buses (84.0%), radial, CSV/MATPOWER, MILP.

| | |
|---|---|
| Ē | 0.001 |
| Backend | MILP, proven optimal |
| Hops | 10 |
| Scenarios | 67, 91 |
| Radiality | `:in_model` |
| Preserved | `:required` — center-tap transformers, phase shifters, regulators, switches |

## Rebuild it

Set these in [`run_optikron.jl`](../../../run_optikron.jl), then run
`julia --project=. run_optikron.jl`:

```julia
case      = "R300"
Ē         = 0.001
scenarios = [67, 91]
hops      = 10
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/R300/84pct"
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `9.505e-04` pu, 95% of `Ē`.
- Held-out scenarios: `4.166e-03` pu — **not certified**, 4.2× `Ē`.

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/R300/84pct")`
