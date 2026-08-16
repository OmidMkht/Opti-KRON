# ieee34 — 51pct

37 → 18 buses (51.4%), radial, CSV.

| | |
|---|---|
| Ē | 0.003 |
| Backend | MILP, optimal |
| Hops | 10 |
| Scenarios | 1 of 1 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("ieee34"; Ē=0.003, hops=10, radiality=:in_model)
```

**Violation** (enforced): `-3.66e-06` — inside budget.  

## OpenDSS

Solved error `9.38e-03` pu — **3× `Ē`**. Not a conversion fault; see
[converter/README.md](../../../converter/README.md).

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `dss` — solved OpenDSS circuit
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point
- `ybus.csv` — Schur complement

`load_case("reduced_networks/ieee34/51pct")`
