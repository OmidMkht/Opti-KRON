# european_lv — 98pct

907 → 21 buses (97.7%), radial, CSV.

| | |
|---|---|
| Ē | 0.01 (not binding) |
| Backend | MILP, optimal |
| Hops | 25 |
| Scenarios | 1 of 1 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("european_lv"; Ē=0.01, hops=25, radiality=:in_model)
```

**Violation** (enforced): `-9.05e-03` — inside budget.  

## OpenDSS

Solved error `5.12e-04` pu — inside `Ē`.

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `dss` — solved OpenDSS circuit
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/european_lv/98pct")`
