# case141 — 74pct

`141` buses reduced to `37` super-nodes (73.8%), radial, MATPOWER format.

## How it was produced

| | |
|---|---|
| Error threshold | `Ē = 0.001` per unit. |
| Hops | `10` |
| Scenarios enforced | all 1 — the case carries a single loading |
| Radiality | enforced in the model (`:in_model`) |
| Devices | switches, regulators and phase shifters preserved |

```julia
result = optikron("case141"; Ē = 0.001, hops = 10,
                  radiality = :in_model)
```

## Verified

Re-checked against the exact nonconvex annulus on the untouched sparse
`Ybus`, not against the linearisation the solver used:

- Enforced scenarios: worst violation `-6.55e-06` — inside budget.

## Files

- `case141_reduced.m` — MATPOWER case, bus map in the header comment

Load it with `load_case("reduced_networks/case141/74pct")`.
