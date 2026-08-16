# case533mt — 83pct

`533` buses reduced to `88` super-nodes (83.5%), radial, MATPOWER format.

## How it was produced

| | |
|---|---|
| Error threshold | `Ē = 0.003` per unit. |
| Hops | `10` |
| Scenarios enforced | all 2 — the case carries 2 loadings |
| Radiality | enforced in the model (`:in_model`) |
| Devices | switches, regulators and phase shifters preserved |

```julia
result = optikron("case533mt"; Ē = 0.003, hops = 10,
                  radiality = :in_model)
```

## Verified

Re-checked against the exact nonconvex annulus on the untouched sparse
`Ybus`, not against the linearisation the solver used:

- Enforced scenarios: worst violation `-5.17e-05` — inside budget.

## Files

- `case533mt_reduced_s1.m` — MATPOWER case, bus map in the header comment
- `case533mt_reduced_s2.m` — MATPOWER case, bus map in the header comment

Load it with `load_case("reduced_networks/case533mt/83pct")`.
