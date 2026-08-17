# case85 — 78pct

85 → 19 buses (77.6%), radial, CSV/MATPOWER, MILP.

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
case      = "case85"
Ē         = 0.006
hops      = 10
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/case85/78pct"
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `5.919e-03` pu, 99% of `Ē`.

## Files

- `case85_reduced.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case85/78pct")`
