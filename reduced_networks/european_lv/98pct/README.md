# european_lv — 98pct

907 → 21 buses (97.7%), radial, CSV/MATPOWER, MILP.

| | |
|---|---|
| Ē | 0.01 |
| Backend | MILP, proven optimal |
| Hops | 25 |
| Scenarios | all |
| Radiality | `:in_model` |
| Preserved | `:required` — center-tap transformers, phase shifters, regulators, switches |

## Rebuild it

Set these in [`run_optikron.jl`](../../../run_optikron.jl), then run
`julia --project=. run_optikron.jl`:

```julia
case      = "european_lv"
Ē         = 0.01
hops      = 25
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/european_lv/98pct"
export_dss = true
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `1.169e-03` pu, 12% of `Ē`.
- **OpenDSS re-solve** (constant-power, real windings): `8.828e-04` pu, 9% of `Ē`.
- Equipment audit: location-only.

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `dss` — solved OpenDSS circuit
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/european_lv/98pct")`
