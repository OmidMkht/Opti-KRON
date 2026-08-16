# european_lv — 91pct

907 → 82 buses (91.0%), radial, CSV.

| | |
|---|---|
| Ē | 0.01 (not binding) |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 1 of 1 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("european_lv"; Ē=0.01, hops=10, radiality=:in_model)
```

**Violation** (enforced): `-9.52e-03` — inside budget.  

## OpenDSS

Solved error `1.55e-04` pu — inside `Ē`.

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `dss` — solved OpenDSS circuit
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/european_lv/91pct")`
