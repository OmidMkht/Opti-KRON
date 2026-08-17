# case141 — 74pct

141 → 37 buses (73.8%), radial, CSV/MATPOWER, MILP.

| | |
|---|---|
| Ē | 0.001 |
| Backend | MILP, proven optimal |
| Hops | 10 |
| Scenarios | all |
| Radiality | `:in_model` |
| Preserved | `:required` — center-tap transformers, phase shifters, regulators, switches |

## Rebuild it

Set these in [`run_optikron.jl`](../../../run_optikron.jl), then run
`julia --project=. run_optikron.jl`:

```julia
case      = "case141"
Ē         = 0.001
hops      = 10
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/case141/74pct"
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `9.920e-04` pu, 99% of `Ē`.

## Files

- `case141_reduced.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case141/74pct")`
