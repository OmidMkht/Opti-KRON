# ieee123 — 57pct

130 → 56 buses (56.9%), radial, CSV/MATPOWER, MILP.

| | |
|---|---|
| Ē | 0.002 |
| Backend | MILP, proven optimal |
| Hops | 10 |
| Scenarios | all |
| Radiality | `:in_model` |
| Preserved | `:required` — center-tap transformers, phase shifters, regulators, switches |

## Rebuild it

Set these in [`run_optikron.jl`](../../../run_optikron.jl), then run
`julia --project=. run_optikron.jl`:

```julia
case      = "ieee123"
Ē         = 0.002
hops      = 10
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/ieee123/57pct"
export_dss = true
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `1.994e-03` pu, 100% of `Ē`.
- **OpenDSS re-solve** (constant-power, real windings): `9.149e-04` pu, 46% of `Ē`.
- Equipment audit: location-only.

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `dss` — solved OpenDSS circuit
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/ieee123/57pct")`
