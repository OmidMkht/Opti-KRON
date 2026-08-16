# european_lv — 91pct

`907` buses reduced to `82` super-nodes (91.0%), radial, CSV format.

## How it was produced

| | |
|---|---|
| Error threshold | `Ē = 0.01` — **not binding**: this feeder reduces identically anywhere in 0.001–0.01, so the pair of examples is separated by `hops`, not by budget. |
| Hops | `10` |
| Scenarios enforced | all 1 — the case carries a single loading |
| Radiality | enforced in the model (`:in_model`) |
| Devices | switches, regulators and phase shifters preserved |

```julia
result = optikron("european_lv"; Ē = 0.01, hops = 10,
                  radiality = :in_model)
```

## Verified

Re-checked against the exact nonconvex annulus on the untouched sparse
`Ybus`, not against the linearisation the solver used:

- Enforced scenarios: worst violation `-9.52e-03` — inside budget.

## Files

- `assignment.csv` — the reduction map, `bus_id → super_node`
- `bus.csv` — surviving buses and their phasing
- `load.csv` — injections after each eliminated bus hands its load over
- `voltage.csv` — the operating point at the surviving buses
- `ybus.csv` — the Schur complement

Load it with `load_case("reduced_networks/european_lv/91pct")`.
