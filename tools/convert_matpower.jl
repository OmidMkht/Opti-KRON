# --------------------------------------------------------------------------- #
# Convert MATPOWER `.m` cases into the Opti-KRON CSV format.
#
# Run from the repository root:
#     julia --project=. tools/convert_matpower.jl data/MATPOWER/case69.m data/case69
#
# Several files sharing one topology become one case with one scenario each,
# which is how a high/low loading pair belongs together:
#     julia --project=. tools/convert_matpower.jl \
#         data/MATPOWER/case533mt_hi.m,data/MATPOWER/case533mt_lo.m data/case533mt hi,lo
#
# ---- Why this is more than reading three matrices ------------------------- #
#
# A MATPOWER case is a *program*, not a data file, and the distribution cases
# use that. Three things bite:
#
#   1. Units. case69, case85 and case141 list branch impedances in ohms and
#      loads in kW, then convert to per-unit in executable MATLAB at the bottom
#      of the file. A parser that reads only the matrices gets ohms, silently,
#      and every per-unit result downstream is wrong by a factor of Vbase^2/Sbase.
#      case141 additionally splits a kVA load into P and Q through a power
#      factor. So the trailing statements are parsed and applied, and anything
#      executable that is *not* recognised is a hard error -- emitting plausible
#      numbers in the wrong units is the one failure this converter must never
#      have.
#
#   2. Expressions inside the matrices. case533mt writes base voltages as
#      `135/sqrt(3)` and its MVA base as `50/3`. Entries are therefore evaluated
#      as arithmetic, not parsed as floats.
#
#   3. Open switches. case533mt carries 577 branches of which 45 are out of
#      service -- tie switches held open. Keeping them would close 45 loops and
#      the feeder would be rejected as meshed. Dropping them leaves 532 branches
#      across 533 buses: a tree, as intended.
# --------------------------------------------------------------------------- #
using Printf

const BUS_I, BUS_TYPE, PD, QD, GS, BS = 1, 2, 3, 4, 5, 6
const BASE_KV = 10
const GEN_BUS, PG, QG = 1, 2, 3
const F_BUS, T_BUS, BR_R, BR_X, BR_B = 1, 2, 3, 4, 5
const TAP, SHIFT, BR_STATUS = 9, 10, 11

"""
    convert_matpower(paths, out_dir; scenarios, name) -> out_dir

Convert one or more MATPOWER files into `out_dir` as `bus.csv`, `branch.csv`
and `load.csv`, and copy the sources in beside them so the case stays traceable
to what it came from.

Every file must describe the same network; only their loads may differ. Each
becomes one loading scenario.
"""
function convert_matpower(paths::Vector{String}, out_dir::AbstractString;
    scenarios::Vector{String}=String[],
    name::AbstractString=basename(rstrip(out_dir, ['/', '\\'])))

    isempty(paths) && error("No input files given.")
    labels = isempty(scenarios) ? [_case_label(p) for p in paths] : scenarios
    length(labels) == length(paths) ||
        error("Got $(length(paths)) files but $(length(labels)) scenario names.")

    cases = [read_matpower(p) for p in paths]
    reference = first(cases)
    for (path, case) in zip(paths[2:end], cases[2:end])
        _assert_same_network(reference, case, first(paths), path)
    end

    mkpath(out_dir)
    _write_bus(joinpath(out_dir, "bus.csv"), reference)
    _write_branch(joinpath(out_dir, "branch.csv"), reference)
    _write_load(joinpath(out_dir, "load.csv"), cases, labels)
    for path in paths
        destination = joinpath(out_dir, basename(path))
        # Re-running in place is a normal thing to do when checking a conversion;
        # copying a file onto itself is not.
        abspath(path) == abspath(destination) || cp(path, destination; force=true)
    end

    live = count(row -> row[BR_STATUS] != 0, reference.branch)
    @printf("%-14s %5d buses, %4d branches (%d open switch(es) dropped), %d scenario(s)%s\n",
        name, length(reference.bus), live, length(reference.branch) - live,
        length(cases), live == length(reference.bus) - 1 ? ", radial" : " -- NOT radial")
    return out_dir
end

convert_matpower(path::AbstractString, out_dir::AbstractString; kwargs...) =
    convert_matpower([String(path)], out_dir; kwargs...)

"Strip the leading `case` from a filename, so case533mt_hi -> hi where it helps."
function _case_label(path)
    stem = splitext(basename(path))[1]
    parts = split(stem, '_')
    return length(parts) > 1 ? String(parts[end]) : "base"
end

# ---- reading -------------------------------------------------------------- #

"""
    read_matpower(path) -> NamedTuple

Parse a MATPOWER case: `baseMVA`, `bus`, `gen`, `branch`, with all unit
conversions in the file already applied.
"""
function read_matpower(path::AbstractString)
    text = read(path, String)
    base_mva = _scalar(text, "mpc.baseMVA")
    bus = _matrix(text, "bus")
    gen = _matrix(text, "gen")
    branch = _matrix(text, "branch")

    isempty(bus) && error("$path: no mpc.bus matrix found.")
    isempty(branch) && error("$path: no mpc.branch matrix found.")

    _apply_conversions!(bus, branch, base_mva, text, path)
    return (base_mva=base_mva, bus=bus, gen=gen, branch=branch, path=path)
end

"Body of `mpc.<field> = [ ... ];`, one vector per row, entries evaluated."
function _matrix(text::AbstractString, field::AbstractString)
    m = match(Regex("mpc\\.$field\\s*=\\s*\\[(.*?)\\n\\s*\\];", "s"), text)
    m === nothing && return Vector{Float64}[]

    rows = Vector{Float64}[]
    for line in split(m.captures[1], '\n')
        cleaned = strip(rstrip(strip(split(line, '%')[1]), ';'))
        isempty(cleaned) && continue
        push!(rows, [_evaluate(token) for token in split(cleaned)])
    end
    return rows
end

function _scalar(text::AbstractString, lhs::AbstractString)
    m = match(Regex("$(replace(lhs, "." => "\\."))\\s*=\\s*([^;%\\n]+)"), text)
    m === nothing && error("Could not find $lhs.")
    return _evaluate(m.captures[1])
end

"""
    _evaluate(token) -> Float64

Evaluate a numeric entry, which may be arithmetic rather than a literal --
case533mt writes base voltages as `135/sqrt(3)`.

Deliberately not `eval` of arbitrary parsed code: only numbers, the four
operations, powers and `sqrt` are permitted, so a malformed or unexpected file
fails loudly instead of running.
"""
function _evaluate(token::AbstractString)
    expression = Meta.parse(strip(token))
    return Float64(_evaluate(expression))
end

_evaluate(x::Real) = x
function _evaluate(expression::Expr)
    expression.head === :call ||
        error("Unsupported expression in case file: $expression")
    operator = expression.args[1]
    operands = [_evaluate(a) for a in expression.args[2:end]]

    operator === :+ && return sum(operands)
    operator === :- && return length(operands) == 1 ? -operands[1] : operands[1] - operands[2]
    operator === :* && return prod(operands)
    operator === :/ && return operands[1] / operands[2]
    operator === :^ && return operands[1]^operands[2]
    operator === :sqrt && return sqrt(operands[1])
    error("Unsupported function `$operator` in case file.")
end

# ---- the executable tail --------------------------------------------------- #

"""
    _apply_conversions!(bus, branch, base_mva, text, path)

Apply the unit conversions the case file performs in MATLAB after its matrices,
and refuse to continue if it does anything else.

The refusal is the point. These cases state impedances in ohms and loads in kW
and fix that in code; a converter that ignored the code would produce a case
that loads cleanly, solves cleanly, and is wrong by orders of magnitude.
"""
function _apply_conversions!(bus, branch, base_mva, text, path)
    tail = _executable_tail(text)

    if occursin(r"mpc\.branch\(:,\s*\[BR_R\s+BR_X\]\)\s*=\s*mpc\.branch\(:,\s*\[BR_R\s+BR_X\]\)\s*/", tail)
        v_base = bus[1][BASE_KV] * 1e3
        s_base = base_mva * 1e6
        z_base = v_base^2 / s_base
        for row in branch
            row[BR_R] /= z_base
            row[BR_X] /= z_base
        end
    end

    # A kVA magnitude split into P and Q by power factor (case141).
    pf_match = match(r"pf\s*=\s*([0-9.]+)\s*;", tail)
    if pf_match !== nothing && occursin("sin(acos(pf))", tail)
        pf = parse(Float64, pf_match.captures[1])
        for row in bus
            row[QD] = row[PD] * sin(acos(pf))
            row[PD] = row[PD] * pf
        end
    end

    if occursin(r"mpc\.bus\(:,\s*\[PD,?\s*QD\]\)\s*=\s*mpc\.bus\(:,\s*\[PD,?\s*QD\]\)\s*/\s*1e3", tail)
        for row in bus
            row[PD] /= 1e3
            row[QD] /= 1e3
        end
    end

    unknown = _unrecognised_statements(tail)
    isempty(unknown) || error("""
        $path performs conversions this converter does not recognise:
            $(join(unknown, "\n    "))
        Refusing to continue: the matrices alone may be in ohms or kW, and
        emitting them as per-unit would be silently wrong. Teach
        _apply_conversions! the statement, or convert the file in MATPOWER first.""")
    return bus, branch
end

"Everything after the last matrix, which is where MATPOWER cases put their code."
function _executable_tail(text::AbstractString)
    last_close = findlast("];", text)
    return last_close === nothing ? "" : text[last(last_close)+1:end]
end

"""
Statements in the tail that are not among the conversions handled above.

MATLAB continues a statement across lines with `...`, and these cases use it for
the long `idx_bus` / `idx_brch` constant lists. Those are joined back into single
logical lines first, or each fragment would be reported as its own mystery.
"""
function _unrecognised_statements(tail::AbstractString)
    joined = replace(tail, r"\.\.\.[^\n]*\n" => " ")
    known = [
        r"^\s*$", r"^\s*%", r"^\s*end\s*;?\s*$", r"^\s*function\b",
        r"^\s*define_constants\s*;?\s*$",
        r"^\s*\[[^\]]*\]\s*=\s*idx_\w+\s*;",     # column-index constants
        r"^\s*(V|S)base\s*=", r"^\s*pf\s*=",
        r"^\s*mpc\.(bus|branch)\(:,",
        r"^\s*mpc\.(version|baseMVA)\b",
    ]
    return [strip(line) for line in split(joined, '\n')
            if !any(pattern -> occursin(pattern, line), known)]
end

# ---- writing --------------------------------------------------------------- #

function _assert_same_network(a, b, path_a, path_b)
    length(a.bus) == length(b.bus) ||
        error("$path_a has $(length(a.bus)) buses, $path_b has $(length(b.bus)).")
    for (row_a, row_b) in zip(a.bus, b.bus)
        row_a[BUS_I] == row_b[BUS_I] ||
            error("$path_a and $path_b disagree on bus numbering.")
    end
    length(a.branch) == length(b.branch) ||
        error("$path_a and $path_b have different branch counts.")
    for (row_a, row_b) in zip(a.branch, b.branch)
        (row_a[F_BUS], row_a[T_BUS], row_a[BR_STATUS]) ==
        (row_b[F_BUS], row_b[T_BUS], row_b[BR_STATUS]) ||
            error("$path_a and $path_b disagree on topology; they cannot be scenarios of one case.")
    end
    return true
end

_bus_id(row) = string(Int(row[BUS_I]))

function _write_bus(path, case)
    slack_rows = [row for row in case.bus if Int(row[BUS_TYPE]) == 3]
    length(slack_rows) == 1 ||
        error("Expected exactly one slack bus (type 3); found $(length(slack_rows)).")

    open(path, "w") do io
        println(io, "bus_id,phases,base_kv,type")
        for row in case.bus
            shunt = row[GS] != 0 || row[BS] != 0
            shunt && error("Bus $(_bus_id(row)) carries a shunt (Gs/Bs), which " *
                           "branch.csv cannot express. Use the Ybus path for this case.")
            println(io, _bus_id(row), ",a,", row[BASE_KV], ",",
                Int(row[BUS_TYPE]) == 3 ? "slack" : "pq")
        end
    end
end

function _write_branch(path, case)
    open(path, "w") do io
        println(io, "from_bus,to_bus,phases,r,x,b")
        for row in case.branch
            row[BR_STATUS] == 0 && continue           # open tie switch
            tap = length(row) >= TAP ? row[TAP] : 0.0
            shift = length(row) >= SHIFT ? row[SHIFT] : 0.0
            (tap in (0.0, 1.0) && shift == 0.0) ||
                error("Branch $(Int(row[F_BUS]))-$(Int(row[T_BUS])) has tap $tap and " *
                      "shift $shift; branch.csv models neither. Use the Ybus path.")
            println(io, Int(row[F_BUS]), ",", Int(row[T_BUS]), ",a,",
                row[BR_R], ",", row[BR_X], ",", row[BR_B])
        end
    end
end

"""
Injections, per scenario: generation positive, load negative, per unit on the
case's own MVA base. The slack carries none -- it is the reference, and whatever
it supplies is a result rather than an input.
"""
function _write_load(path, cases, labels)
    open(path, "w") do io
        println(io, "bus_id,phase,scenario,p_pu,q_pu")
        for (case, label) in zip(cases, labels)
            generation = Dict{Int,Tuple{Float64,Float64}}()
            for row in case.gen
                bus = Int(row[GEN_BUS])
                p, q = get(generation, bus, (0.0, 0.0))
                generation[bus] = (p + row[PG], q + row[QG])
            end

            for row in case.bus
                Int(row[BUS_TYPE]) == 3 && continue
                bus = Int(row[BUS_I])
                gen_p, gen_q = get(generation, bus, (0.0, 0.0))
                p = (gen_p - row[PD]) / case.base_mva
                q = (gen_q - row[QD]) / case.base_mva
                (p == 0 && q == 0) && continue        # absent rows import as zero
                println(io, bus, ",a,", label, ",", p, ",", q)
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 ||
        error("usage: julia tools/convert_matpower.jl <in1.m[,in2.m,...]> <out_dir> [scen1,scen2,...]")
    inputs = String.(split(ARGS[1], ','))
    names = length(ARGS) >= 3 ? String.(split(ARGS[3], ',')) : String[]
    convert_matpower(inputs, ARGS[2]; scenarios=names)
end
