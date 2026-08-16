# ieee123 — 57pct

`130` buses reduced to `56` super-nodes (56.9%), radial, CSV format.

## How it was produced

| | |
|---|---|
| Error threshold | `Ē = 0.002` per unit. |
| Backend | MILP, proven optimal (`:milp`) |
| Hops | `10` |
| Scenarios enforced | all 1 — the case carries a single loading |
| Radiality | enforced in the model (`:in_model`) |
| Devices | switches, regulators and phase shifters preserved |

```julia
result = optikron("ieee123"; Ē = 0.002, hops = 10,
                  radiality = :in_model)
```

## Verified

Re-checked against the exact nonconvex annulus on the untouched sparse
`Ybus`, not against the linearisation the solver used:

- Enforced scenarios: worst violation `-6.47e-06` — inside budget.

## Solved as OpenDSS

`dss/` is this reduction as a self-contained OpenDSS circuit,
rebuilt by `converter/build_reduced_dss.py` and solved.

- Power aggregation error: `0.0` — exact.
- Synthesized Kron `Ybus` vs this reduction's: `4.216173281961903e-14` relative.
- **Solved voltage error against the original: `9.41e-04` pu.**

  Inside the certified budget even under a full constant-power re-solve.

  Per-node detail is in `dss/validation.csv`.

## Files

- `assignment.csv` — the reduction map, `bus_id → super_node`
- `bus.csv` — surviving buses and their phasing
- `dss` — the same reduction as a solved OpenDSS circuit
- `load.csv` — injections after each eliminated bus hands its load over
- `voltage.csv` — the operating point at the surviving buses
- `ybus.csv` — the Schur complement

Load it with `load_case("reduced_networks/ieee123/57pct")`.
