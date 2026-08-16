# ieee123 — 57pct

130 → 56 buses (56.9%), radial, CSV.

| | |
|---|---|
| Ē | 0.002 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 1 of 1 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("ieee123"; Ē=0.002, hops=10, radiality=:in_model)
```

**Violation** (enforced): `-6.47e-06` — inside budget.  

## OpenDSS

Solved error `9.41e-04` pu — inside `Ē`.

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `dss` — solved OpenDSS circuit
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/ieee123/57pct")`
