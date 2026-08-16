# european_lv — 98pct

`907` buses reduced to `21` super-nodes (97.7%), radial, CSV format.

## How it was produced

| | |
|---|---|
| Error threshold | `Ē = 0.01` — **not binding**: this feeder reduces identically anywhere in 0.001–0.01, so the pair of examples is separated by `hops`, not by budget. |
| Backend | MILP, proven optimal (`:milp`) |
| Hops | `25` |
| Scenarios enforced | all 1 — the case carries a single loading |
| Radiality | enforced in the model (`:in_model`) |
| Devices | switches, regulators and phase shifters preserved |

```julia
result = optikron("european_lv"; Ē = 0.01, hops = 25,
                  radiality = :in_model)
```

## Verified

Re-checked against the exact nonconvex annulus on the untouched sparse
`Ybus`, not against the linearisation the solver used:

- Enforced scenarios: worst violation `-9.05e-03` — inside budget.

## Solved as OpenDSS

`dss/` is this reduction as a self-contained OpenDSS circuit,
rebuilt by `converter/build_reduced_dss.py` and solved.

- Power aggregation error: `0.0` — exact.
- Synthesized Kron `Ybus` vs this reduction's: `1.7839580851209356e-13` relative.
- **Solved voltage error against the original: `5.12e-04` pu.**

  Inside the certified budget even under a full constant-power re-solve.

  Per-node detail is in `dss/validation.csv`.

## Files

- `assignment.csv` — the reduction map, `bus_id → super_node`
- `bus.csv` — surviving buses and their phasing
- `dss` — the same reduction as a solved OpenDSS circuit
- `load.csv` — injections after each eliminated bus hands its load over
- `voltage.csv` — the operating point at the surviving buses
- `ybus.csv` — the Schur complement

Load it with `load_case("reduced_networks/european_lv/98pct")`.
