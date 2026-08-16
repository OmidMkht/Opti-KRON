# case69 — 77pct

`69` buses reduced to `16` super-nodes (76.8%), radial, MATPOWER format.

## How it was produced

| | |
|---|---|
| Error threshold | `Ē = 0.006` per unit. |
| Hops | `10` |
| Scenarios enforced | all 1 — the case carries a single loading |
| Radiality | enforced in the model (`:in_model`) |
| Devices | switches, regulators and phase shifters preserved |

```julia
result = optikron("case69"; Ē = 0.006, hops = 10,
                  radiality = :in_model)
```

## Verified

Re-checked against the exact nonconvex annulus on the untouched sparse
`Ybus`, not against the linearisation the solver used:

- Enforced scenarios: worst violation `-2.21e-04` — inside budget.

## Files

- `case69_reduced.m` — MATPOWER case, bus map in the header comment

Load it with `load_case("reduced_networks/case69/77pct")`.
