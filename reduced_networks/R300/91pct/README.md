# R300 — 91pct

300 → 26 buses (91.3%), radial, CSV/MATPOWER, MILP.

| | |
|---|---|
| Ē | 0.003 |
| Backend | MILP, proven optimal |
| Hops | 10 |
| Scenarios | 67, 91 |
| Radiality | `:in_model` |
| Preserved | `:required` — center-tap transformers, phase shifters, regulators, switches |

## Rebuild it

Set these in [`run_optikron.jl`](../../../run_optikron.jl), then run
`julia --project=. run_optikron.jl`:

```julia
case      = "R300"
Ē         = 0.003
scenarios = [67, 91]
hops      = 10
preserve  = :required
radiality = :in_model
export_to = "reduced_networks/R300/91pct"
```

## Accuracy

- **Constant-PQ** (what the budget certifies): `2.986e-03` pu, 100% of `Ē`.
- Held-out scenarios: `1.393e-02` pu — **not certified**, 4.6× `Ē`.

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/R300/91pct")`
