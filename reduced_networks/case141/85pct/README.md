# case141 — 85pct

141 → 21 buses (85.1%), radial, CSV/MATPOWER, MILP.

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
case      = "case141"
Ē         = 0.003
hops      = 10
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/case141/85pct"
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `2.990e-03` pu, 100% of `Ē`.

## Files

- `case141_reduced.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case141/85pct")`
