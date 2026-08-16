# R300 — 84pct

`300` buses reduced to `48` super-nodes (84.0%), radial, CSV format.

## How it was produced

| | |
|---|---|
| Error threshold | `Ē = 0.001` per unit. |
| Backend | MILP, proven optimal (`:milp`) |
| Hops | `10` |
| Scenarios enforced | 67 and 91 of 168 — peak and minimum hour |
| Radiality | enforced in the model (`:in_model`) |
| Devices | switches, regulators and phase shifters preserved |

```julia
result = optikron("R300"; Ē = 0.001, hops = 10, scenarios = [67, 91],
                  radiality = :in_model)
```

## Verified

Re-checked against the exact nonconvex annulus on the untouched sparse
`Ybus`, not against the linearisation the solver used:

- Enforced scenarios: worst violation `-4.95e-05` — inside budget.
- Held-out scenarios (166 of 168): `+3.17e-03`. **Not certified** — the budget
  was never enforced here, and the error reaches 4.2× `Ē` across them.

## Files

- `assignment.csv` — the reduction map, `bus_id → super_node`
- `bus.csv` — surviving buses and their phasing
- `load.csv` — injections after each eliminated bus hands its load over
- `voltage.csv` — the operating point at the surviving buses
- `ybus.csv` — the Schur complement

Load it with `load_case("reduced_networks/R300/84pct")`.
