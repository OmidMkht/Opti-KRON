# case533mt — 83pct

533 → 88 buses (83.5%), radial, MATPOWER.

| | |
|---|---|
| Ē | 0.003 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 2 of 2 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("case533mt"; Ē=0.003, hops=10, radiality=:in_model)
```

**Violation** (enforced): `-5.17e-05` — inside budget.  

## Files

- `case533mt_reduced_s1.m` — MATPOWER, bus map in header
- `case533mt_reduced_s2.m` — MATPOWER, bus map in header

`load_case("reduced_networks/case533mt/83pct")`
