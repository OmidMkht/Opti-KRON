# case533mt — 66pct

533 → 180 buses (66.2%), radial, CSV/MATPOWER, MILP.

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
case      = "case533mt"
Ē         = 0.001
hops      = 10
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/case533mt/66pct"
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `9.998e-04` pu, 100% of `Ē`.

## Files

- `case533mt_reduced_s1.m` — MATPOWER, bus map in header
- `case533mt_reduced_s2.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case533mt/66pct")`
