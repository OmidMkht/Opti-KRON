# ieee34 — 51pct

37 → 18 buses (51.4%), radial, CSV/MATPOWER, MILP.

| | |
|---|---|
| Ē | 0.003 |
| Backend | MILP, proven optimal |
| Hops | 10 |
| Scenarios | all |
| Radiality | `:in_model` |
| Preserved | `:required` — center-tap transformers, phase shifters, regulators, switches |

## Rebuild it

Set these in [`run_optikron.jl`](../../../run_optikron.jl), then run
`julia --project=. run_optikron.jl`:

```julia
case      = "ieee34"
Ē         = 0.003
hops      = 10
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/ieee34/51pct"
export_dss = true
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `2.551e-03` pu, 85% of `Ē`.
- **OpenDSS re-solve** (constant-power, real windings): `1.706e-04` pu, 6% of `Ē`.
- Equipment audit: location-only.

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `dss` — solved OpenDSS circuit
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/ieee34/51pct")`
