# ieee123 — 78pct

130 → 29 buses (77.7%), radial, CSV/MATPOWER, MILP.

| | |
|---|---|
| Ē | 0.006 |
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
Ē         = 0.006
hops      = 10
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/ieee123/78pct"
export_dss = true
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `5.880e-03` pu, 98% of `Ē`.
- **OpenDSS re-solve** (constant-power, real windings): `3.949e-03` pu, 66% of `Ē`.
- Equipment audit: location-only.

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `dss` — solved OpenDSS circuit
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/ieee123/78pct")`
