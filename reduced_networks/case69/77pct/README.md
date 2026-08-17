# case69 — 77pct

69 → 16 buses (76.8%), radial, CSV/MATPOWER, MILP.

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
case      = "case69"
Ē         = 0.006
hops      = 10
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/case69/77pct"
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `5.779e-03` pu, 96% of `Ē`.

## Files

- `case69_reduced.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case69/77pct")`
