# case85 — 59pct

85 → 35 buses (58.8%), radial, CSV/MATPOWER, MILP.

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
case      = "case85"
Ē         = 0.002
hops      = 10
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/case85/59pct"
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `1.944e-03` pu, 97% of `Ē`.

## Files

- `case85_reduced.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case85/59pct")`
