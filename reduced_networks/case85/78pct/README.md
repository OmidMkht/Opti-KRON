# case85 — 78pct

`85` buses reduced to `19` super-nodes (77.6%), radial, MATPOWER format.

## How it was produced

| | |
|---|---|
| Error threshold | `Ē = 0.006` per unit. |
| Hops | `10` |
| Scenarios enforced | all 1 — the case carries a single loading |
| Radiality | enforced in the model (`:in_model`) |
| Devices | switches, regulators and phase shifters preserved |

```julia
result = optikron("case85"; Ē = 0.006, hops = 10,
                  radiality = :in_model)
```

## Verified

Re-checked against the exact nonconvex annulus on the untouched sparse
`Ybus`, not against the linearisation the solver used:

- Enforced scenarios: worst violation `-8.12e-05` — inside budget.

## Files

- `case85_reduced.m` — MATPOWER case, bus map in the header comment

Load it with `load_case("reduced_networks/case85/78pct")`.
