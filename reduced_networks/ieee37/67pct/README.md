# ieee37 — 67pct

39 → 13 buses (66.7%), radial, CSV/MATPOWER, MILP.

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
case      = "ieee37"
Ē         = 0.003
hops      = 10
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/ieee37/67pct"
export_dss = true
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `2.909e-03` pu, 97% of `Ē`.
- **OpenDSS re-solve** (constant-power, real windings): `5.605e-04` pu, 19% of `Ē`.
- Equipment audit: 1 winding(s) rebased.

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `dss` — solved OpenDSS circuit
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/ieee37/67pct")`
