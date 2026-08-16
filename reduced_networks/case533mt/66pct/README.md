# case533mt — 66pct

533 → 180 buses (66.2%), radial, MATPOWER.

| | |
|---|---|
| Ē | 0.001 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 2 of 2 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("case533mt"; Ē=0.001, hops=10, radiality=:in_model)
```

**Violation** (enforced): `-1.91e-07` — inside budget.  

## Files

- `case533mt_reduced_s1.m` — MATPOWER, bus map in header
- `case533mt_reduced_s2.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case533mt/66pct")`
