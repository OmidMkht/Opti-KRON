# R300 — 91pct

300 → 26 buses (91.3%), radial, CSV.

| | |
|---|---|
| Ē | 0.003 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 67, 91 of 168 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("R300"; Ē=0.003, hops=10, scenarios=[67, 91], radiality=:in_model)
```

**Violation** (enforced): `-3.92e-04` — inside budget.  
**Violation** (held-out, 166 of 168): `+1.09e-02` — 4.6× `Ē`, not certified.

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/R300/91pct")`
