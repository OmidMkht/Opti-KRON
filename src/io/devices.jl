# --------------------------------------------------------------------------- #
# Equipment that must survive a reduction intact.
#
# A transformer, regulator or switch is a modelled device, not just impedance,
# and its tap ratio, winding connection and phase shift cannot be recovered once
# folded into a Schur complement. Pinning both terminals preserves it exactly:
# deleting the device edge splits a radial feeder into T_i and T_j, and no
# eliminated bus touches both sides -- such a bus would be a second i-j path,
# which a tree does not have. So Y_rr is block-diagonal across the split and the
# correction at entry (i,j) is identically zero. The same argument rules out
# fill-in bridging the sides, so the device also stays the unique cut edge.
#
# Exact in arithmetic; in Float64 the residue is under one ULP of the block's own
# scale, because `pinv` reaches the inverse through an SVD that does not preserve
# block-diagonality bit-for-bit.
#
# Loading is not preserved: the assignment relocates injections, so the current
# through the device changes. The error budget bounds that.
#
# The stronger reason to pin is *state*. A regulator taps and a switch opens; a
# reduction that folded either in would silently encode one tap position or one
# switch state and be wrong the moment the device moved. A delta-wye phase shift
# never changes, so pinning it is about keeping the network interpretable.
#
# `transformer.csv` is the edge list; `regulator.csv`,
# `phase_shift_equipment.csv` and `center_tapped_transformer.csv` classify rows
# of it rather than adding edges. A transformer named by the phase-shift table is
# `:phase_shift`, one named by the center-tap table is `:center_tap`, one named by
# a regulator's `transformer_id` is `:regulator`, and the rest are
# `:transformer`. Rows with `enabled = 0`, open switches, and terminals absent
# from bus.csv are all dropped -- none of them is an edge of this Ybus.
#
# A center-tapped transformer is the one kind that cannot be crossed even with
# rebasing. Its two secondaries sit on one bus in opposite orientation (`.1.0`
# and `.0.2` about the centre conductor), so the primary-to-secondary map is not
# a diagonal scaling and no per-unit change of base reproduces it. Everything
# downstream of one is a separate phase reference domain.
#
# `capacitor_bank.csv` is a shunt, not an edge, so it pins one bus. It is off by
# default: the Schur complement carries its admittance correctly, and pinning
# only matters when the bank's switching state has to stay operable.
# --------------------------------------------------------------------------- #

const DEVICE_KINDS = (:transformer, :center_tap, :regulator, :phase_shift, :switch, :capacitor)

"""
The kinds that are *edges* of this Ybus, and so the ones the exactness argument
above applies to. `:capacitor` is not among them: a bank is a shunt, so pinning
its bus keeps the bank in place and relocatable but does not leave the diagonal
block untouched. It is opt-in for that reason, and `read_devices` defaults here.
"""
const EDGE_KINDS = (:transformer, :center_tap, :regulator, :phase_shift, :switch)

"""
    Device

One piece of equipment the reduction must keep intact, as a bus-index pair.

`kind` is one of `$(DEVICE_KINDS)`. `note` carries the source row's own
qualifier -- winding connections for a transformer, `status` for a switch -- so a
reduction can report *why* a bus was pinned.
"""
struct Device
    id::String
    kind::Symbol
    from::Int
    to::Int
    note::String
end

"""
    read_devices(dir, net; kinds=EDGE_KINDS) -> Vector{Device}

Devices in `dir` of the requested `kinds` whose terminals both exist in `net`.
Empty when the case ships no device files, as the MATPOWER feeders do not.

Reads `transformer.csv`, `regulator.csv`, `phase_shift_equipment.csv`,
`center_tapped_transformer.csv`, `switch.csv` and `capacitor_bank.csv`. Phase
columns are ignored -- pinning is per bus, and a pinned bus keeps every phase.

`kinds` accepts any subset of `$(DEVICE_KINDS)`, or `:all`. Dropping `:switch`
merges across closed switches, which is nearly free since a closed switch is a
jumper whose ends are electrically one node -- worth it unless the switch must
stay operable for reconfiguration studies. Dropping `:regulator` lets the
reduction absorb a regulator at its present tap, which is only safe if the taps
are frozen. Dropping `:center_tap` is not safe on a feeder that has any: see the
header.
"""
function read_devices(dir::AbstractString, net::Network; kinds=EDGE_KINDS)
    wanted = _device_kinds(kinds)
    index_of = Dict(id => k for (k, id) in enumerate(net.bus_ids))
    found = Device[]

    shifters = _ids_in(joinpath(dir, "phase_shift_equipment.csv"), :equipment_id)
    regulated = _ids_in(joinpath(dir, "regulator.csv"), :transformer_id)
    centre = _ids_in(joinpath(dir, "center_tapped_transformer.csv"), :transformer_id)

    _each_row(joinpath(dir, "transformer.csv")) do row
        id = _as_string(row.transformer_id)
        kind = id in shifters ? :phase_shift :
               id in centre ? :center_tap :
               id in regulated ? :regulator : :transformer
        kind in wanted || return
        note = string(_column(row, :from_connection, "?"), "-",
            _column(row, :to_connection, "?"))
        _push_device!(found, index_of, id, kind, row, note)
    end

    if :switch in wanted
        _each_row(joinpath(dir, "switch.csv")) do row
            status = _column(row, :status, "")
            status == "open" && return          # not an edge of this Ybus
            _push_device!(found, index_of, _as_string(row.switch_id), :switch, row, status)
        end
    end

    if :capacitor in wanted
        _each_row(joinpath(dir, "capacitor_bank.csv")) do row
            hasproperty(row, :bus_id) || return
            bus = get(index_of, _as_string(row.bus_id), nothing)
            bus === nothing && return
            push!(found, Device(_as_string(row.capacitor_id), :capacitor, bus, bus,
                _column(row, :connection, "")))
        end
    end
    return found
end

"Normalise the `kinds` argument, rejecting anything that is not a device kind."
function _device_kinds(kinds)
    kinds === :all && return DEVICE_KINDS
    kinds isa Symbol && (kinds = (kinds,))
    for k in kinds
        k in DEVICE_KINDS || error("$k is not a device kind; expected some of $DEVICE_KINDS.")
    end
    return Tuple(kinds)
end

"Values of `column` in `path`, as a set. Empty when the file is absent or has no such column."
function _ids_in(path::AbstractString, column::Symbol)
    ids = Set{String}()
    _each_row(path) do row
        hasproperty(row, column) && push!(ids, _as_string(getproperty(row, column)))
    end
    return ids
end

"Apply `f` to every enabled row of `path`, skipping the file when it is absent or empty."
function _each_row(f, path::AbstractString)
    isfile(path) || return
    df = CSV.read(path, DataFrame)
    isempty(df) && return
    for row in eachrow(df)
        _column(row, :enabled, "1") == "0" && continue
        f(row)
    end
end

"A row's `column` as a lowercase string, or `default` when absent or blank."
function _column(row, column::Symbol, default::AbstractString)
    hasproperty(row, column) || return default
    value = getproperty(row, column)
    ismissing(value) && return default
    text = lowercase(strip(_as_string(value)))
    return isempty(text) ? default : text
end

function _push_device!(found::Vector{Device}, index_of::Dict{String,Int},
    id::AbstractString, kind::Symbol, row, note::AbstractString)

    for column in (:from_bus, :to_bus)
        hasproperty(row, column) || error("A device row has no `$column` column.")
    end
    from = get(index_of, _as_string(row.from_bus), nothing)
    to = get(index_of, _as_string(row.to_bus), nothing)
    # A missing terminal means the device was collapsed before it reached us.
    (from === nothing || to === nothing || from == to) && return
    push!(found, Device(String(id), kind, from, to, String(note)))
end

"Bus indices to pin, sorted and deduplicated. Goes to `solve_milp(...; pin=...)`."
preserved_buses(devices::AbstractVector{Device}) =
    sort!(unique!(reduce(vcat, ([d.from, d.to] for d in devices); init=Int[])))

function Base.show(io::IO, d::Device)
    print(io, "Device(", d.id, ", ", d.kind, ", ", d.from, "->", d.to,
        isempty(d.note) ? "" : ", \"$(d.note)\"", ")")
end
