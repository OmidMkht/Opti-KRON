# case85 — 59pct

85 → 35 buses (58.8%), radial, MATPOWER.

| | |
|---|---|
| Ē | 0.002 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 1 of 1 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("case85"; Ē=0.002, hops=10, radiality=:in_model)
```

**Violation** (enforced): `-5.58e-05` — inside budget.  

## Files

- `case85_reduced.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case85/59pct")`
