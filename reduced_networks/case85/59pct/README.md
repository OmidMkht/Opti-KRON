# case85 — 59pct

`85` buses reduced to `35` super-nodes (58.8%), radial, MATPOWER format.

## How it was produced

| | |
|---|---|
| Error threshold | `Ē = 0.002` per unit. |
| Hops | `10` |
| Scenarios enforced | all 1 — the case carries a single loading |
| Radiality | enforced in the model (`:in_model`) |
| Devices | switches, regulators and phase shifters preserved |

```julia
result = optikron("case85"; Ē = 0.002, hops = 10,
                  radiality = :in_model)
```

## Verified

Re-checked against the exact nonconvex annulus on the untouched sparse
`Ybus`, not against the linearisation the solver used:

- Enforced scenarios: worst violation `-5.58e-05` — inside budget.

## Files

- `case85_reduced.m` — MATPOWER case, bus map in the header comment

Load it with `load_case("reduced_networks/case85/59pct")`.
