# case1197 — 65pct

1197 → 417 buses (65.2%), radial, CSV/MATPOWER, MILP.

| | |
|---|---|
| Ē | 0.009 |
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
Ē         = 0.009
hops      = 10
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/case1197/65pct"
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `8.962e-03` pu, 100% of `Ē`.

## Files

- `case1197_reduced.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case1197/65pct")`
