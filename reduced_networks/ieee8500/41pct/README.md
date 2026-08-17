# ieee8500 — 41pct

4876 → 2856 buses (41.4%), radial, CSV, search.

| | |
|---|---|
| Ē | 0.02 |
| Backend | exhaustive search, GPU |
| Hops | none |
| Scenarios | all |
| Radiality | `:in_model` |
| Preserved | `:required` — center-tap transformers, phase shifters, regulators, switches |

## Rebuild it

Set these in [`run_optikron.jl`](../../../run_optikron.jl), then run
`julia --project=. run_optikron.jl`:

```julia
case      = "ieee8500"
Ē         = 0.02
backend   = :search_gpu   # the MILP does not close at this size
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/ieee8500/41pct"
export_dss = true
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `1.611e-02` pu, 81% of `Ē`.
- **OpenDSS re-solve** (constant-power, real windings): `2.540e-03` pu, 13% of `Ē`.
- Equipment audit: location-only.

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `dss` — solved OpenDSS circuit
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/ieee8500/41pct")`
