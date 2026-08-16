# case69 — 77pct

69 → 16 buses (76.8%), radial, MATPOWER.

| | |
|---|---|
| Ē | 0.006 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 1 of 1 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("case69"; Ē=0.006, hops=10, radiality=:in_model)
```

**Violation** (enforced): `-2.21e-04` — inside budget.  

## Files

- `case69_reduced.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case69/77pct")`
