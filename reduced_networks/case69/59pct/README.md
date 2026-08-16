# case69 — 59pct

69 → 28 buses (59.4%), radial, MATPOWER.

| | |
|---|---|
| Ē | 0.001 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 1 of 1 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("case69"; Ē=0.001, hops=10, radiality=:in_model)
```

**Violation** (enforced): `-1.29e-05` — inside budget.  

## Files

- `case69_reduced.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case69/59pct")`
