# case85 — 59pct

`85` buses reduced to `35` super-nodes (58.8%), radial, MATPOWER format.

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
result = optikron("case85"; Ē = 0.002, hops = 10,
                  radiality = :in_model)
```

## Verified

Re-checked against the exact nonconvex annulus on the untouched sparse
`Ybus`, not against the linearisation the solver used:

- Enforced scenarios: worst violation `-5.58e-05` — inside budget.

## No `ybus.csv` here

Kron reduction fills in. The Schur complement of this level runs to
millions of nonzeros — 289 MB at `42pct` and 631 MB at `54pct` — because
the writer keeps every entry that is not exactly zero, on the principle
that deciding what counts as structural is not a file writer's job.

Measured, **99.5% of those entries are round-off**: only 28312 of
5390534 at `42pct`, and 20593 of 11743569 at `54pct`, exceed 1e-12 of the
largest, and that count is unchanged at 1e-10 and 1e-8. Too large to ship,
and one line to rebuild:

```julia
using OptiKRON, CSV, DataFrames
net = load_case("case85")
df  = CSV.read("reduced_networks/case85/59pct/assignment.csv", DataFrame)
idx = Dict(id => i for (i, id) in enumerate(net.bus_ids))
A   = zeros(nnodes(net), nnodes(net))
for r in eachrow(df); A[idx[string(r.super_node)], idx[string(r.bus_id)]] = 1; end
Y_reduced = kron_reduce(net, A)
```

`assignment.csv` is the reduction; the rest is derived from it.

## Files

- `case85_reduced.m` — MATPOWER case, bus map in the header comment

Load it with `load_case("reduced_networks/case85/59pct")`.
