# R100 — 91pct

100 → 9 buses (91.0%), radial, CSV.

| | |
|---|---|
| Ē | 0.005 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 67, 91 of 168 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("R100"; Ē=0.005, hops=10, scenarios=[67, 91], radiality=:in_model)
```

**Violation** (enforced): `-2.36e-04` — inside budget.  
**Violation** (held-out, 166 of 168): `+6.59e-03` — 2.3× `Ē`, not certified.

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/R100/91pct")`
