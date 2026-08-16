# case85 — 78pct

85 → 19 buses (77.6%), radial, MATPOWER.

| | |
|---|---|
| Ē | 0.006 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 1 of 1 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("case85"; Ē=0.006, hops=10, radiality=:in_model)
```

**Violation** (enforced): `-8.12e-05` — inside budget.  

## Files

- `case85_reduced.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case85/78pct")`
