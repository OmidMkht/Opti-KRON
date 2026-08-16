# case141 — 85pct

141 → 21 buses (85.1%), radial, MATPOWER.

| | |
|---|---|
| Ē | 0.003 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 1 of 1 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("case141"; Ē=0.003, hops=10, radiality=:in_model)
```

**Violation** (enforced): `-9.90e-06` — inside budget.  

## Files

- `case141_reduced.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case141/85pct")`
