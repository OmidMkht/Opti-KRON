# R100 — 91pct

`100` buses reduced to `9` super-nodes (91.0%), radial, CSV format.

## How it was produced

| | |
|---|---|
| Error threshold | `Ē = 0.005` per unit. |
| Hops | `10` |
| Scenarios enforced | 67 and 91 of 168 — peak and minimum hour |
| Radiality | enforced in the model (`:in_model`) |
| Devices | switches, regulators and phase shifters preserved |

```julia
result = optikron("R100"; Ē = 0.005, hops = 10, scenarios = [67, 91],
                  radiality = :in_model)
```

## Verified

Re-checked against the exact nonconvex annulus on the untouched sparse
`Ybus`, not against the linearisation the solver used:

- Enforced scenarios: worst violation `-2.36e-04` — inside budget.
- Held-out scenarios (166 of 168): `+6.59e-03`. **Not certified** — the budget
  was never enforced here, and the error reaches 2.3× `Ē` across them.

## Files

- `assignment.csv` — the reduction map, `bus_id → super_node`
- `bus.csv` — surviving buses and their phasing
- `load.csv` — injections after each eliminated bus hands its load over
- `voltage.csv` — the operating point at the surviving buses
- `ybus.csv` — the Schur complement

Load it with `load_case("reduced_networks/R100/91pct")`.
