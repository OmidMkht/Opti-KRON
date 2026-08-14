# --------------------------------------------------------------------------- #
# Equipment that must survive a reduction intact.
#
# A transformer, regulator or switch is a modelled device, not just impedance,
# and its tap ratio, winding connection and phase shift cannot be recovered once
# folded into a Schur complement. Pinning both terminals preserves it exactly:
# deleting the device edge splits a radial feeder into T_i and T_j, and no
# eliminated bus touches both sides -- such a bus would be a second i-j path,
# which a tree does not have. So Y_rr is block-diagonal across the split, the
# correction at entry (i,j) is identically zero, and Y_red[i,j] == Y[i,j] to the
# last bit. The same argument rules out fill-in bridging the sides, so the device
# also stays the unique cut edge.
#
# Loading is not preserved: the assignment relocates injections, so the current
# through the device changes. The error budget bounds that.
#
# Two filters matter. Only `model_scope = explicit_ybus` devices are an edge in
# our Ybus at all -- feeders exported with service transformers already collapsed
# mark those `collapsed_primary_load`, worth 2 pinned buses against 1140 on
# ieee8500. Both terminals must also appear in bus.csv.
# --------------------------------------------------------------------------- #

"""
    Device

One piece of equipment the reduction must keep intact, as a bus-index pair.

`kind` is `:transformer` (covering regulators and tap changers, which ship in the
same table) or `:switch`. `note` carries the source file's qualifier verbatim --
`model_scope` or `status` -- so a reduction can report *why* a bus was pinned.
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
of `net.Ybus`. Empty when the case ships no device files, as the MATPOWER feeders
do not.

Reads `transformer.csv` (`transformer_id, from_bus, from_phases, to_bus,
to_phases, model_scope`) and `switch.csv` (`switch_id, ..., status`). Phase
columns are ignored -- pinning is per bus, and a pinned bus keeps every phase.

Only `explicit_ybus` transformers are returned (see the header, worth 45% of
`ieee8500`), and `open` switches are dropped -- an open switch is not an edge in
`Ybus`. Pass `switches=false` to merge across closed ones, which is nearly free
since a closed switch is a jumper whose ends are electrically one node, unless
the switch must stay operable for reconfiguration studies.
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
        # A blank qualifier is unknown, not a value -- switch status often is.
        raw = has_qualifier ? row[qualifier] : missing
        note = ismissing(raw) ? "" : _as_string(raw)
        admit(note) || continue

        from = get(index_of, _as_string(row.from_bus), nothing)
        to = get(index_of, _as_string(row.to_bus), nothing)
        # A missing terminal means the device was collapsed before it reached us.
        (from === nothing || to === nothing || from == to) && continue

        push!(found, Device(_as_string(row[id_column]), kind, from, to, note))
    end
    return found
end

"Bus indices to pin, sorted and deduplicated. Goes to `solve_milp(...; pin=...)`."
preserved_buses(devices::AbstractVector{Device}) =
    sort!(unique!(reduce(vcat, ([d.from, d.to] for d in devices); init=Int[])))

function Base.show(io::IO, d::Device)
    print(io, "Device(", d.id, ", ", d.kind, ", ", d.from, "->", d.to,
        isempty(d.note) ? "" : ", \"$(d.note)\"", ")")
end
