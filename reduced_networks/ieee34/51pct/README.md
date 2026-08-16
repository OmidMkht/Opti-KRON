# ieee34 — 51pct

`37` buses reduced to `18` super-nodes (51.4%), radial, CSV format.

## How it was produced

| | |
|---|---|
| Error threshold | `Ē = 0.003` per unit. |
| Hops | `10` |
| Scenarios enforced | all 1 — the case carries a single loading |
| Radiality | enforced in the model (`:in_model`) |
| Devices | switches, regulators and phase shifters preserved |

```julia
result = optikron("ieee34"; Ē = 0.003, hops = 10,
                  radiality = :in_model)
```

## Verified

Re-checked against the exact nonconvex annulus on the untouched sparse
`Ybus`, not against the linearisation the solver used:

- Enforced scenarios: worst violation `-3.66e-06` — inside budget.

## Files

- `assignment.csv` — the reduction map, `bus_id → super_node`
- `bus.csv` — surviving buses and their phasing
- `load.csv` — injections after each eliminated bus hands its load over
- `voltage.csv` — the operating point at the surviving buses
- `ybus.csv` — the Schur complement

Load it with `load_case("reduced_networks/ieee34/51pct")`.
