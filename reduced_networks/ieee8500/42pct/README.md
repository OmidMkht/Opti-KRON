# ieee8500 — 42pct

4876 → 2821 buses (42.1%), radial, CSV, search.

| | |
|---|---|
| Ē | 0.005 |
| Backend | search (`:search_gpu`) |
| Hops | none |
| Scenarios | 1 of 1 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("ieee8500"; Ē=0.005, backend=:search_gpu, radiality=:in_model)
```

**Violation** (enforced): `-9.30e-07` — inside budget.  

## No `ybus.csv`

Kron fill-in is 289 MB; only 28312 of 5390534 entries (0.5%) exceed 1e-12 of the
largest, rest is round-off. Too large to ship. Rebuild:

```julia
net = load_case("ieee8500")
df  = CSV.read("reduced_networks/ieee8500/42pct/assignment.csv", DataFrame)
idx = Dict(id => i for (i, id) in enumerate(net.bus_ids))
A   = zeros(nnodes(net), nnodes(net))
for r in eachrow(df); A[idx[string(r.super_node)], idx[string(r.bus_id)]] = 1; end
Y_reduced = kron_reduce(net, A)
```

## Files

- `assignment.csv` — reduction map
- `bus.csv` — surviving buses
- `load.csv` — injections post-reduction
- `voltage.csv` — operating point

`load_case("reduced_networks/ieee8500/42pct")`
