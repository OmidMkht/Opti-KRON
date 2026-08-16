# ieee37 — 62pct

`39` buses reduced to `15` super-nodes (61.5%), radial, CSV format.

## How it was produced

| | |
|---|---|
| Error threshold | `Ē = 0.003` per unit. |
| Backend | MILP, proven optimal (`:milp`) |
| Hops | `10` |
| Scenarios enforced | all 1 — the case carries a single loading |
| Radiality | enforced in the model (`:in_model`) |
| Devices | switches, regulators and phase shifters preserved |

```julia
result = optikron("ieee37"; Ē = 0.003, hops = 10,
                  radiality = :in_model)
```

## Verified

Re-checked against the exact nonconvex annulus on the untouched sparse
`Ybus`, not against the linearisation the solver used:

- Enforced scenarios: worst violation `-1.84e-04` — inside budget.

## Solved as OpenDSS

`dss/` is this reduction as a self-contained OpenDSS circuit,
rebuilt by `converter/build_reduced_dss.py` and solved.

- Power aggregation error: `0.0` — exact.
- Synthesized Kron `Ybus` vs this reduction's: `2.92728359381024e-16` relative.
- **Solved voltage error against the original: `8.39e-02` pu.**

  That is `28×` the certified `Ē`. Not a conversion fault — the
  budget is certified against a constant-current linearisation, and
  OpenDSS re-solves with constant-power loads and the real winding
  connections. See [converter/README.md](../../../converter/README.md).

  Per-node detail is in `dss/validation.csv`.

## Files

- `assignment.csv` — the reduction map, `bus_id → super_node`
- `bus.csv` — surviving buses and their phasing
- `dss` — the same reduction as a solved OpenDSS circuit
- `load.csv` — injections after each eliminated bus hands its load over
- `voltage.csv` — the operating point at the surviving buses
- `ybus.csv` — the Schur complement

Load it with `load_case("reduced_networks/ieee37/62pct")`.
