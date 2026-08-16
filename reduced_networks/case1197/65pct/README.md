# case1197 — 65pct

1197 → 417 buses (65.2%), radial, MATPOWER.

| | |
|---|---|
| Ē | 0.009 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 1 of 1 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("case1197"; Ē=0.009, hops=10, radiality=:in_model)
```

**Violation** (enforced): `-3.82e-05` — inside budget.  

## Files

- `case1197_reduced.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case1197/65pct")`
