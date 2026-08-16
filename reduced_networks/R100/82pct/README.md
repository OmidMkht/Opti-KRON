# R100 — 82pct

100 → 18 buses (82.0%), radial, CSV.

| | |
|---|---|
| Ē | 0.001 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 67, 91 of 168 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("R100"; Ē=0.001, hops=10, scenarios=[67, 91], radiality=:in_model)
```

**Violation** (enforced): `-5.33e-05` — inside budget.  
**Violation** (held-out, 166 of 168): `+2.63e-03` — 3.6× `Ē`, not certified.

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/R100/82pct")`
