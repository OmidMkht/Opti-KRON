# R100 — 91pct

100 → 9 buses (91.0%), radial, CSV/MATPOWER, MILP.

| | |
|---|---|
| Ē | 0.005 |
| Backend | MILP, proven optimal |
| Hops | 10 |
| Scenarios | 67, 91 |
| Radiality | `:in_model` |
| Preserved | `:required` — center-tap transformers, phase shifters, regulators, switches |

## Rebuild it

Set these in [`run_optikron.jl`](../../../run_optikron.jl), then run
`julia --project=. run_optikron.jl`:

```julia
case      = "R100"
Ē         = 0.005
scenarios = [67, 91]
hops      = 10
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/R100/91pct"
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `4.764e-03` pu, 95% of `Ē`.
- Held-out scenarios: `1.159e-02` pu — **not certified**, 2.3× `Ē`.

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/R100/91pct")`
