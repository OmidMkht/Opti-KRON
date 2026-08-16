# case1197 — 50pct

1197 → 593 buses (50.5%), radial, MATPOWER.

| | |
|---|---|
| Ē | 0.006 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 1 of 1 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("case1197"; Ē=0.006, hops=10, radiality=:in_model)
```

**Violation** (enforced): `-2.82e-06` — inside budget.  

## Files

- `case1197_reduced.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case1197/50pct")`
