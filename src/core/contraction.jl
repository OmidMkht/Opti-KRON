# --------------------------------------------------------------------------- #
# Switch contraction -- reduce with switches short-circuited, then put them back.
#
# A closed switch is a jumper of near-zero impedance, so its two terminals are
# electrically one node. Modelling it as an edge costs on both sides: the huge
# admittance (1.15e6 pu on ieee123, against 3e4 for the next largest branch) sets
# sigma_max and so the condition number, and preserving the switch means pinning
# both terminals, which on ieee123 freezes 19 buses.
#
# Contracting removes both costs. Reduce the contracted feeder, then expand the
# assignment back to the original bus set: every original bus is represented by
# whichever member of its super-node's group it sits nearest, which is well posed
# because deleting the switch edges splits a radial feeder into one component per
# group member. Both terminals survive, so the switch comes back through
# `kron_reduce` exactly -- preserved for free rather than by pinning.
#
# Ordering is the thing to be careful with. Contraction renumbers buses, so every
# index crossing the boundary goes through `SwitchContraction`: `parent` maps an
# original bus to its contracted one, `members` maps back. Groups are ordered by
# their smallest original index, so the numbering is deterministic and as close
# to the original as it can be.
# --------------------------------------------------------------------------- #

"""
    SwitchContraction

The index map between an original feeder and its switch-contracted form.

- `parent[j]`    contracted bus holding original bus `j`.
- `members[m]`   original buses held by contracted bus `m`, ascending.
- `nbuses`       original bus count, so an expansion can size its output.

`length(members)` is the contracted bus count. A group of one is a bus no switch
touched, so an untouched feeder gives `parent == 1:nbuses`.
"""
struct SwitchContraction
    parent::Vector{Int}
    members::Vector{Vector{Int}}
    nbuses::Int
end

ncontracted(c::SwitchContraction) = length(c.members)

"""
    contract_switches(net, switches) -> (Network, SwitchContraction)

Short-circuit each switch in `switches`, returning the contracted feeder and the
index map back.

Contraction is exact, not an approximation: it states `V_a = V_b` and adds the
two buses' injected currents, which for the admittance is "sum the rows, sum the
columns, drop one". The result is what the feeder already behaves like, since the
jumper between them carries no meaningful voltage drop.

`switches` is normally `read_devices(dir, net; kinds=(:switch,))`. Open switches
are not edges and must not be passed; `read_devices` already drops them.
"""
function contract_switches(net::Network, switches::AbstractVector{Device})
    B = nnodes(net)
    parent = collect(1:B)
    find(x) = (while parent[x] != x
        parent[x] = parent[parent[x]]
        x = parent[x]
    end; x)

    for d in switches
        a, b = find(d.from), find(d.to)
        a == b || (parent[max(a, b)] = min(a, b))
    end

    # Groups ordered by their smallest member, so the numbering is deterministic
    # and stays as close to the original as contraction allows.
    roots = sort!(unique(find(j) for j in 1:B))
    index_of_root = Dict(r => m for (m, r) in enumerate(roots))
    group = [index_of_root[find(j)] for j in 1:B]
    members = [Int[] for _ in roots]
    for j in 1:B
        push!(members[group[j]], j)
    end

    contraction = SwitchContraction(group, members, B)
    return _contracted_network(net, contraction), contraction
end

contract_switches(net::Network, dir::AbstractString) =
    contract_switches(net, read_devices(dir, net; kinds=(:switch,)))

"Build the contracted `Network`: phases unioned, rows and columns summed."
function _contracted_network(net::Network, c::SwitchContraction)
    n = ncontracted(c)
    blocks = node_rows(net)

    phases = falses(3, n)
    for m in 1:n, j in c.members[m], p in 1:3
        net.phases[p, j] && (phases[p, m] = true)
    end
    new_blocks = node_rows(phases)

    # Original node-phase row -> contracted node-phase row, the only index map
    # the numerical work needs. Both sides list phases in a/b/c order, so the
    # match is by phase, never by position -- a single-phase switch can put a
    # phase-`a` bus into a group whose other member carries `abc`.
    rowmap = zeros(Int, nphase_rows(net))
    for j in 1:nnodes(net)
        m = c.parent[j]
        group_phases = [p for p in 1:3 if phases[p, m]]
        for (local_j, p) in enumerate(p for p in 1:3 if net.phases[p, j])
            local_m = findfirst(==(p), group_phases)
            rowmap[blocks[j][local_j]] = new_blocks[m][local_m]
        end
    end

    nph = count(phases)
    rows, cols, vals = Int[], Int[], ComplexF64[]
    Y = net.Ybus
    yrows, yvals = rowvals(Y), nonzeros(Y)
    for col in 1:size(Y, 2), k in nzrange(Y, col)
        push!(rows, rowmap[yrows[k]])
        push!(cols, rowmap[col])
        push!(vals, yvals[k])
    end
    Ybus = sparse(rows, cols, vals, nph, nph, +)

    S = zeros(ComplexF64, nph, nscenarios(net))
    for r in 1:nphase_rows(net)
        S[rowmap[r], :] .+= net.S[r, :]
    end

    bus_ids = [net.bus_ids[first(c.members[m])] for m in 1:n]
    Lambda = adjacency_from_ybus(Ybus, new_blocks)
    return validate(Network(net.name * "_contracted", bus_ids, phases, Ybus,
        Lambda, c.parent[net.slack], S, nothing))
end

"""
    uncontract(A_contracted, contraction, tree) -> Matrix{Float64}

Expand an assignment over the contracted feeder back to the original bus set.

A contracted super-node stands for a group of original buses joined by switches.
When the group survives, every member survives with it -- the switch is an edge
again, so both terminals are kept and [`kron_reduce`](@ref) reproduces its
admittance exactly. When the group is absorbed, each original bus goes to the
member of the absorbing group it lies nearest in `tree`.

That nearest-member rule is the primary/secondary choice made precise. Deleting
the switch edges splits a radial feeder into one component per group member, so
every bus has exactly one member on its own side and there is no tie to break.
"""
function uncontract(A_contracted::AbstractMatrix, c::SwitchContraction, tree::RadialTree)
    n = ncontracted(c)
    size(A_contracted) == (n, n) ||
        error("Assignment is $(size(A_contracted)) but the contracted feeder has $n buses.")

    A = zeros(Float64, c.nbuses, c.nbuses)
    for j in 1:c.nbuses
        m = c.parent[j]
        host = findfirst(!iszero, view(A_contracted, :, m))
        host === nothing && error("Contracted bus $m is represented by nothing.")
        if host == m
            A[j, j] = 1.0                       # its own group survived
        else
            A[_nearest(c.members[host], j, tree), j] = 1.0
        end
    end
    return A
end

"The candidate lying nearest `j` in the original feeder -- which side of the switch it came from."
function _nearest(candidates::AbstractVector{Int}, j::Int, tree::RadialTree)
    length(candidates) == 1 && return only(candidates)
    best, best_hops = first(candidates), typemax(Int)
    for g in candidates
        h = hop_distance(tree, j, g)
        h < best_hops && ((best, best_hops) = (g, h))
    end
    return best
end
