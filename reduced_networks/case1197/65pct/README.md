# case1197 — 65pct

`1197` buses reduced to `417` super-nodes (65.2%), radial, MATPOWER format.

## How it was produced

| | |
|---|---|
| Error threshold | `Ē = 0.009` per unit. |
| Hops | `10` |
| Scenarios enforced | all 1 — the case carries a single loading |
| Radiality | enforced in the model (`:in_model`) |
| Devices | switches, regulators and phase shifters preserved |

```julia
result = optikron("case1197"; Ē = 0.009, hops = 10,
                  radiality = :in_model)
```

## Verified

Re-checked against the exact nonconvex annulus on the untouched sparse
`Ybus`, not against the linearisation the solver used:

- Enforced scenarios: worst violation `-3.82e-05` — inside budget.

## Files

- `case1197_reduced.m` — MATPOWER case, bus map in the header comment

Load it with `load_case("reduced_networks/case1197/65pct")`.
