# case1197 — 50pct

1197 → 593 buses (50.5%), radial, CSV/MATPOWER, MILP.

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
case      = "case1197"
Ē         = 0.006
hops      = 10
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/case1197/50pct"
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `5.990e-03` pu, 100% of `Ē`.

## Files

- `case1197_reduced.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case1197/50pct")`
