# ieee34 — 59pct

37 → 15 buses (59.5%), radial, CSV.

| | |
|---|---|
| Ē | 0.01 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 1 of 1 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("ieee34"; Ē=0.01, hops=10, radiality=:in_model)
```

**Violation** (enforced): `-3.42e-03` — inside budget.  

## OpenDSS

Solved error `1.03e-02` pu — **1× `Ē`**. Not a conversion fault; see
[converter/README.md](../../../converter/README.md).

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `dss` — solved OpenDSS circuit
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/ieee34/59pct")`
