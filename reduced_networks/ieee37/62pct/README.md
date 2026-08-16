# ieee37 — 62pct

39 → 15 buses (61.5%), radial, CSV.

| | |
|---|---|
| Ē | 0.003 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 1 of 1 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("ieee37"; Ē=0.003, hops=10, radiality=:in_model)
```

**Violation** (enforced): `-1.84e-04` — inside budget.  

## OpenDSS

Solved error `8.39e-02` pu — **28× `Ē`**. Not a conversion fault; see
[converter/README.md](../../../converter/README.md).

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `dss` — solved OpenDSS circuit
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/ieee37/62pct")`
