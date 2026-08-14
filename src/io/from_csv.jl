# --------------------------------------------------------------------------- #
# The canonical public input format: three CSVs describing buses, lines and
# loads. Line parameters rather than a raw Ybus, because a Ybus cannot be
# inverted back to per-line r/x once shunts or transformers are present -- so a
# Ybus-only pipeline can only ever emit a Ybus. See from_ybus.jl for the fast
# path when branch data genuinely is not available.
#
#   bus.csv     bus_id, phases, base_kv, type
#                 phases over {a,b,c}, e.g. "abc"; exactly one type = "slack"
#
#   branch.csv  from_bus, to_bus, phases, <impedance columns>
#                 single-phase: r, x, b
#                 three-phase:  r_aa..r_cc, x_aa..x_cc, b_aa..b_cc, upper
#                 triangle only (the matrices are symmetric), per unit
#
#   load.csv    bus_id, phase, scenario, p_pu, q_pu
#                 absent buses/phases draw zero; `scenario` is any label, and
#                 columns of S come out sorted by first appearance
#
# p_pu/q_pu are *injections* -- generation positive, load negative.
# --------------------------------------------------------------------------- #

const PHASE_SYMBOLS = (:a, :b, :c)
const PHASE_INDEX = Dict(:a => 1, :b => 2, :c => 3)

# CSV.jl hands back InlineStrings (String3, String7, ...) for short text columns
# and raw Ints for numeric ids. Everything downstream keys on `String`, so
# normalise once at the boundary rather than sprinkling conversions around.
_as_string(x) = String(string(x))

"Parse a phase string like \"abc\" or \"ac\" into ordered symbols."
function parse_phases(s::AbstractString)
    cleaned = lowercase(strip(s))
    isempty(cleaned) && error("Empty phase string.")
    out = Symbol[]
    for ch in cleaned
        sym = Symbol(ch)
        sym in PHASE_SYMBOLS || error("Unknown phase '$ch' in \"$s\"; expected some ordering of a, b, c.")
        sym in out && error("Phase '$ch' listed twice in \"$s\".")
        push!(out, sym)
    end
    return sort(out, by=p -> PHASE_INDEX[p])
end

"""
    read_network_csv(dir; name, base_mva) -> Network

Load a feeder from `bus.csv`, `branch.csv` and `load.csv` in `dir`.
"""
function read_network_csv(dir::AbstractString;
    name::AbstractString=basename(rstrip(dir, ['/', '\\'])),
    base_mva::Real=1.0)

    bus_ids, phases, slack, index_of = _read_bus_table(dir)
    branch_df = _require_csv(dir, "branch.csv")
    load_path = joinpath(dir, "load.csv")
    load_df = isfile(load_path) ? CSV.read(load_path, DataFrame) : nothing

    branches = _read_branches(branch_df, index_of)
    blocks = node_rows(phases)
    Ybus = _build_ybus(branches, phases, blocks)
    Lambda = adjacency_from_ybus(Ybus, blocks)
    S = _read_injections(load_df, index_of, phases, blocks, base_mva)

    return validate(Network(name, bus_ids, phases, Ybus, Lambda, slack,
        S, branches))
end

function _require_csv(dir::AbstractString, file::AbstractString)
    path = joinpath(dir, file)
    isfile(path) || error("Expected $file in $dir. See src/io/from_csv.jl for the schema.")
    return CSV.read(path, DataFrame)
end

"""
    _read_bus_table(dir) -> (bus_ids, phases, slack, index_of)

Read `bus.csv`, which both the branch-based and the Ybus fast path depend on:
it is what fixes bus ordering, the phase mask and the slack bus.
"""
function _read_bus_table(dir::AbstractString)
    bus_df = _require_csv(dir, "bus.csv")

    bus_ids = _as_string.(bus_df.bus_id)
    allunique(bus_ids) || error("bus.csv contains duplicate bus_id values.")
    index_of = Dict(id => i for (i, id) in enumerate(bus_ids))

    phases = falses(3, length(bus_ids))
    for (i, spec) in enumerate(bus_df.phases)
        for p in parse_phases(_as_string(spec))
            phases[PHASE_INDEX[p], i] = true
        end
    end

    slack_rows = findall(t -> lowercase(strip(_as_string(t))) == "slack", bus_df.type)
    length(slack_rows) == 1 ||
        error("bus.csv must mark exactly one bus as type=slack; found $(length(slack_rows)).")

    return bus_ids, phases, only(slack_rows), index_of
end

"Read branch rows, accepting either the single-phase (r, x, b) or three-phase (r_aa ...) column set."
function _read_branches(df::DataFrame, index_of::Dict{String,Int})
    cols = Set(string.(names(df)))
    scalar_form = "r" in cols
    scalar_form || "r_aa" in cols ||
        error("branch.csv needs either an `r` column (single-phase) or `r_aa` (three-phase).")

    branches = Branch[]
    for (row_number, row) in enumerate(eachrow(df))
        from_id, to_id = _as_string(row.from_bus), _as_string(row.to_bus)
        haskey(index_of, from_id) || error("branch.csv row $row_number references unknown bus \"$from_id\".")
        haskey(index_of, to_id) || error("branch.csv row $row_number references unknown bus \"$to_id\".")
        ph = parse_phases(_as_string(row.phases))

        if scalar_form
            length(ph) == 1 ||
                error("branch.csv row $row_number lists phases $(row.phases) but uses the scalar " *
                      "r/x/b columns. Use the r_aa ... r_cc form for multi-phase branches.")
            r = fill(Float64(row.r), 1, 1)
            x = fill(Float64(row.x), 1, 1)
            b = fill("b" in cols ? Float64(row.b) : 0.0, 1, 1)
        else
            r = _symmetric_block(row, "r", ph, row_number)
            x = _symmetric_block(row, "x", ph, row_number)
            b = "b_aa" in cols ? _symmetric_block(row, "b", ph, row_number) : zeros(length(ph), length(ph))
        end
        push!(branches, Branch(index_of[from_id], index_of[to_id], ph, r, x, b))
    end
    return branches
end

"Assemble the symmetric `prefix`_pq block for the phases this branch carries."
function _symmetric_block(row, prefix::String, ph::Vector{Symbol}, row_number::Int)
    n = length(ph)
    M = zeros(n, n)
    for a in 1:n, c in a:n
        pa, pc = ph[a], ph[c]
        key = "$(prefix)_$(pa)$(pc)"
        alt = "$(prefix)_$(pc)$(pa)"
        col = hasproperty(row, Symbol(key)) ? Symbol(key) :
              hasproperty(row, Symbol(alt)) ? Symbol(alt) :
              error("branch.csv row $row_number is missing column $key.")
        value = getproperty(row, col)
        ismissing(value) && error("branch.csv row $row_number has an empty $col.")
        M[a, c] = M[c, a] = Float64(value)
    end
    return M
end

"""
    _build_ybus(branches, phases, blocks) -> SparseMatrixCSC

Standard pi-model assembly over node-phase rows. Series admittance is the
matrix inverse of (r + jx) restricted to the phases the branch carries; the
total shunt `b` is split half to each end.
"""
function _build_ybus(branches::Vector{Branch}, phases::AbstractMatrix{Bool}, blocks::Vector{Vector{Int}})
    nph = count(phases)
    I_idx, J_idx, vals = Int[], Int[], ComplexF64[]

    push_block!(rows, cols, M) = for (a, r) in enumerate(rows), (c, col) in enumerate(cols)
        iszero(M[a, c]) && continue
        push!(I_idx, r); push!(J_idx, col); push!(vals, M[a, c])
    end

    for (k, br) in enumerate(branches)
        z = br.r .+ im .* br.x
        det_guard = length(br.phases) == 1 ? abs(z[1, 1]) : abs(det(z))
        det_guard > eps() ||
            error("Branch $k ($(br.from) -> $(br.to)) has a singular series impedance; " *
                  "a zero-impedance jumper must be merged before import, not passed through.")
        y = inv(z)
        shunt = (br.b ./ 2) .* im

        rows_f = _phase_rows(br, br.from, phases, blocks, k)
        rows_t = _phase_rows(br, br.to, phases, blocks, k)

        push_block!(rows_f, rows_f, y .+ shunt)
        push_block!(rows_t, rows_t, y .+ shunt)
        push_block!(rows_f, rows_t, -y)
        push_block!(rows_t, rows_f, -y)
    end

    return sparse(I_idx, J_idx, vals, nph, nph, +)
end

"Map a branch's phases to the global rows of one of its end buses."
function _phase_rows(br::Branch, bus::Int, phases::AbstractMatrix{Bool}, blocks::Vector{Vector{Int}}, k::Int)
    bus_phases = [PHASE_SYMBOLS[p] for p in 1:3 if phases[p, bus]]
    map(br.phases) do p
        local_idx = findfirst(==(p), bus_phases)
        if local_idx === nothing
            available = join(bus_phases, ", ")
            error("Branch $k carries phase $p but bus $bus only has $available. " *
                  "Every branch phase must exist at both end buses.")
        end
        blocks[bus][local_idx]
    end
end

"Build the nph x nscenarios injection matrix; absent entries draw zero."
function _read_injections(df, index_of::Dict{String,Int}, phases::AbstractMatrix{Bool},
    blocks::Vector{Vector{Int}}, base_mva::Real)
    nph = count(phases)
    df === nothing && return zeros(ComplexF64, nph, 1)

    labels = unique(_as_string.(df.scenario))
    column_of = Dict(label => j for (j, label) in enumerate(labels))
    S = zeros(ComplexF64, nph, length(labels))

    for (row_number, row) in enumerate(eachrow(df))
        bus_id = _as_string(row.bus_id)
        haskey(index_of, bus_id) || error("load.csv row $row_number references unknown bus \"$bus_id\".")
        bus = index_of[bus_id]
        p = only(parse_phases(_as_string(row.phase)))
        phases[PHASE_INDEX[p], bus] ||
            error("load.csv row $row_number puts load on phase $p of bus \"$bus_id\", which does not carry it.")

        bus_phases = [PHASE_SYMBOLS[q] for q in 1:3 if phases[q, bus]]
        global_row = blocks[bus][findfirst(==(p), bus_phases)]
        S[global_row, column_of[_as_string(row.scenario)]] += (Float64(row.p_pu) + im * Float64(row.q_pu)) / base_mva
    end
    return S
end

"""
    read_voltage(dir, net) -> Union{Nothing,Matrix{ComplexF64}}

Read `voltage.csv` from a case directory, or `nothing` when the case does not
ship one.

The operating point is problem data -- the reduction is linearised around it --
so a case that carries its own voltages saves guessing at them. Columns are
`bus_id, phase, v_re_pu, v_im_pu`, with an optional `scenario` column when the
case holds more than one; without it every scenario shares the one profile.

Check what you load with [`powerflow_residual`](@ref) before trusting it: a
voltage profile that does not solve the `Ybus` and `S` beside it will still
reduce, and the resulting error bound will mean nothing.
"""
function read_voltage(dir::AbstractString, net::Network)
    path = joinpath(dir, "voltage.csv")
    isfile(path) || return nothing

    df = CSV.read(path, DataFrame)
    index_of = Dict(id => k for (k, id) in enumerate(net.bus_ids))
    blocks = node_rows(net)

    has_scenarios = "scenario" in names(df)
    labels = has_scenarios ? unique(_as_string.(df.scenario)) : ["_"]
    column_of = Dict(label => k for (k, label) in enumerate(labels))
    length(labels) in (1, nscenarios(net)) ||
        error("voltage.csv holds $(length(labels)) scenarios but the case has " *
              "$(nscenarios(net)).")

    V = zeros(ComplexF64, nphase_rows(net), length(labels))
    for row in eachrow(df)
        bus = get(index_of, _as_string(row.bus_id), nothing)
        bus === nothing && error("voltage.csv names bus $(row.bus_id), which bus.csv does not.")
        phase = parse_phases(_as_string(row.phase))[1]
        local_index = findfirst(==(phase), phases_of(net, bus))
        local_index === nothing &&
            error("voltage.csv gives phase $phase at bus $(row.bus_id), which does not carry it.")
        column = has_scenarios ? column_of[_as_string(row.scenario)] : 1
        V[blocks[bus][local_index], column] = Float64(row.v_re_pu) + im * Float64(row.v_im_pu)
    end

    any(iszero, V) && error("voltage.csv leaves some node-phase rows unset.")
    # One profile stands for every scenario; several must match the case exactly.
    return size(V, 2) == nscenarios(net) ? V : repeat(V, 1, nscenarios(net))
end
