# --------------------------------------------------------------------------- #
# Equipment that must survive a reduction intact.
#
# A transformer, regulator or switch is a *modelled device*, not just impedance.
# Its tap ratio, winding connection and phase shift are meaningful to whatever
# reads the reduced network, and none of them can be recovered once the device
# has been folded into a Schur complement. So the reduction is told to leave
# both of its terminals standing.
#
# Why pinning the terminals is enough -- and exactly enough:
#
#   `kron_reduce` eliminates buses by Schur complement,
#   Y_red = Y_kk - Y_kr * pinv(Y_rr) * Y_rk. Take a device on edge (i,j) of a
#   radial feeder with both i and j kept. Deleting that edge splits the tree
#   into T_i and T_j, and no eliminated bus is adjacent to both sides -- such a
#   bus would be a second i-j path, which a tree does not have. So Y_rr is
#   block-diagonal across the two sides and so is its pseudo-inverse, while
#   Y_kr[i,:] is supported only on T_i and Y_rk[:,j] only on T_j. The
#   correction term at entry (i,j) is therefore identically zero:
#
#       Y_red[i,j] == Y[i,j]
#
#   exactly, to the last bit. The same argument rules out fill-in bridging the
#   two sides, so the device also stays the unique cut edge -- the reduced
#   network cannot grow a parallel path around it.
#
# What this does *not* preserve is loading. The assignment matrix relocates
# injections, and a bus on the far side of a device may be represented by one
# on the near side, so the current through the device changes. That is intended:
# the error budget is what bounds it. The physics model is exact; the loading
# moves.
#
# Two filters matter when reading these files, and skipping either is costly:
#
#   - Only devices whose `model_scope` is `explicit_ybus` exist as an edge in
#     our Ybus at all. Feeders exported with service transformers already
#     collapsed mark those `collapsed_primary_load`; their secondary bus is not
#     in `bus.csv`, so there is no edge to preserve and pinning would freeze the
#     primary for nothing. On `ieee8500` that distinction is 2 pinned buses
#     against 1140 -- 0.08% of the feeder against 45.6% of it.
#   - Both terminals must actually appear in `bus.csv`. Source data routinely
#     names buses that its own bus list does not carry.
# --------------------------------------------------------------------------- #

"""
    Device

One piece of equipment the reduction must keep intact, as a bus-index pair.

`kind` is `:transformer` (covering regulators and tap changers, which the
source data ships in the same table) or `:switch`. `note` carries the source
file's own qualifier verbatim -- `model_scope` for transformers, `status` for
switches -- so a reduction can report *why* a bus was pinned.
"""
struct Device
    id::String
    kind::Symbol
    from::Int
    to::Int
    note::String
end

"""
    read_devices(dir, net; switches=true, transformers=true) -> Vector{Device}

Devices in `dir` whose terminals both exist in `net` and which are genuine edges
of `net.Ybus`. Returns an empty vector when the case ships no device files,
which is the normal case for the MATPOWER feeders.

Reads `transformer.csv` (`transformer_id, from_bus, from_phases, to_bus,
to_phases, model_scope`) and `switch.csv` (`switch_id, from_bus, from_phases,
to_bus, to_phases, status`). Phase columns are ignored -- pinning is per bus,
and a pinned bus keeps all of its phases.

Only `model_scope == "explicit_ybus"` transformers are returned; see the header
for why the distinction is worth 45% of `ieee8500`. Switches marked `open` are
dropped too: an open switch is not an edge in `Ybus`, so there is nothing to
preserve exactly. Pass `switches=false` to merge across closed switches, which
is nearly free -- a closed switch is a near-zero-impedance jumper whose two ends
are electrically one node -- and worth doing unless you need the switch to stay
operable for reconfiguration studies.
"""
function read_devices(dir::AbstractString, net::Network;
    switches::Bool=true, transformers::Bool=true)

    index_of = Dict(id => k for (k, id) in enumerate(net.bus_ids))
    found = Device[]

    transformers && _read_device_file!(found, joinpath(dir, "transformer.csv"),
        :transformer, index_of, :model_scope, scope -> scope == "explicit_ybus")
    switches && _read_device_file!(found, joinpath(dir, "switch.csv"),
        :switch, index_of, :status, status -> status != "open")

    return found
end

"Shared reader for the two device tables, which differ only in their qualifier column."
function _read_device_file!(found::Vector{Device}, path::AbstractString, kind::Symbol,
    index_of::Dict{String,Int}, qualifier::Symbol, admit)

    isfile(path) || return found
    df = CSV.read(path, DataFrame)
    isempty(df) && return found

    id_column = Symbol(kind == :transformer ? "transformer_id" : "switch_id")
    for column in (id_column, :from_bus, :to_bus)
        String(column) in names(df) ||
            error("$(basename(path)) has no `$column` column.")
    end
    has_qualifier = String(qualifier) in names(df)

    for row in eachrow(df)
        # A blank qualifier is unknown, not a value: read it as "" and let the
        # admission rule decide. Source data leaves switch status blank often.
        raw = has_qualifier ? row[qualifier] : missing
        note = ismissing(raw) ? "" : _as_string(raw)
        admit(note) || continue

        from = get(index_of, _as_string(row.from_bus), nothing)
        to = get(index_of, _as_string(row.to_bus), nothing)
        # A terminal missing from bus.csv means the device was already collapsed
        # upstream of us; there is no edge here to keep.
        (from === nothing || to === nothing || from == to) && continue

        push!(found, Device(_as_string(row[id_column]), kind, from, to, note))
    end
    return found
end

"""
    preserved_buses(devices) -> Vector{Int}

The bus indices to pin, sorted and deduplicated. Feed to `solve_milp(...;
pin=...)` or let [`optikron`](@ref) do it.
"""
preserved_buses(devices::AbstractVector{Device}) =
    sort!(unique!(reduce(vcat, ([d.from, d.to] for d in devices); init=Int[])))

function Base.show(io::IO, d::Device)
    print(io, "Device(", d.id, ", ", d.kind, ", ", d.from, "->", d.to,
        isempty(d.note) ? "" : ", \"$(d.note)\"", ")")
end
