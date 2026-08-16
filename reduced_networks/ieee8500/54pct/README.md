# ieee8500 — 54pct

4876 → 2224 buses (54.4%), radial, CSV, search.

| | |
|---|---|
| Ē | 0.02 |
| Backend | search (`:search_gpu`) |
| Hops | none |
| Scenarios | 1 of 1 |
| Radiality | `:in_model` |
| Devices | switches, regulators, phase shifters preserved |

```julia
optikron("ieee8500"; Ē=0.02, backend=:search_gpu, radiality=:in_model)
```

**Violation** (enforced): `-8.23e-05` — inside budget.  

## No `ybus.csv`

Kron fill-in is 631 MB; only 20593 of 11743569 entries (0.2%) exceed 1e-12 of the
largest, rest is round-off. Too large to ship. Rebuild:

```julia
net = load_case("ieee8500")
df  = CSV.read("reduced_networks/ieee8500/54pct/assignment.csv", DataFrame)
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

`load_case("reduced_networks/ieee8500/54pct")`
