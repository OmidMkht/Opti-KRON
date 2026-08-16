# ieee123 — 78pct

130 → 29 buses (77.7%), radial, CSV.

| | |
|---|---|
| Ē | 0.006 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 1 of 1 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("ieee123"; Ē=0.006, hops=10, radiality=:in_model)
```

**Violation** (enforced): `-1.11e-04` — inside budget.  

## OpenDSS

Solved error `5.20e-03` pu — inside `Ē`.

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `dss` — solved OpenDSS circuit
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/ieee123/78pct")`
