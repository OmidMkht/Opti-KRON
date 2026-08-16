# ieee8500 — 54pct

`4876` buses reduced to `2224` super-nodes (54.4%), radial, CSV format.

## How it was produced

| | |
|---|---|
| Error threshold | `Ē = 0.02` per unit. |
| Backend | exhaustive search at q=1, GPU (`:search_gpu`) — the MILP does not close at this size |
| Hops | no limit — the search merges neighbours and lets distance accumulate |
| Scenarios enforced | all 1 — the case carries a single loading |
| Radiality | enforced in the model (`:in_model`) |
| Devices | switches, regulators and phase shifters preserved |

```julia
result = optikron("ieee8500"; Ē = 0.02, backend = :search_gpu,
                  radiality = :in_model)
```

## Verified

Re-checked against the exact nonconvex annulus on the untouched sparse
`Ybus`, not against the linearisation the solver used:

- Enforced scenarios: worst violation `-8.23e-05` — inside budget.

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
net = load_case("ieee8500")
df  = CSV.read("reduced_networks/ieee8500/54pct/assignment.csv", DataFrame)
idx = Dict(id => i for (i, id) in enumerate(net.bus_ids))
A   = zeros(nnodes(net), nnodes(net))
for r in eachrow(df); A[idx[string(r.super_node)], idx[string(r.bus_id)]] = 1; end
Y_reduced = kron_reduce(net, A)
```

`assignment.csv` is the reduction; the rest is derived from it.

## Files

- `assignment.csv` — the reduction map, `bus_id → super_node`
- `bus.csv` — surviving buses and their phasing
- `load.csv` — injections after each eliminated bus hands its load over
- `voltage.csv` — the operating point at the surviving buses

Load it with `load_case("reduced_networks/ieee8500/54pct")`.
