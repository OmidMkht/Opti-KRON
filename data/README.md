# Benchmark cases

| Case | Buses | Node-phase rows | Scenarios | Phasing |
|---|---|---|---|---|
| `R100` | 100 | 166 | 168 | Unbalanced three-phase, radial |
| `R300` | 300 | 506 | 168 | Unbalanced three-phase, radial |
| `case69` | 69 | 69 | 1 | Single-phase, radial |
| `case85` | 85 | 85 | 1 | Single-phase, radial |
| `case141` | 141 | 141 | 1 | Single-phase, radial |
| `case533mt` | 533 | 533 | 2 (`hi`, `lo`) | Single-phase, radial |
| `case1197` | 1197 | 1197 | 1 | Single-phase, radial |
| `ieee37` | 38 | 114 | 1 | Three-phase throughout, radial |
| `ieee123` | 126 | 265 | 1 | Unbalanced three-phase, radial |
| `european906` | 907 | 2721 | 1 | Three-phase throughout, radial |
| `ieee8500` | 2501 | 3786 | 1 | Unbalanced three-phase, radial |
| `test4` | 4 | 4 | 2 | Single-phase; unit tests |
| `test3ph` | 3 | 7 | 1 | Mixed `abc`/`a`; unit tests |

```julia
net = load_case("case69")                  # picks the right reader for the case
net = read_network_ybus("data/R100")       # R100/R300 use the Ybus fast path
net = read_network_csv("data/case69")      # the rest carry branch parameters
```

## Reduced networks

Every case ships two worked reductions under `reduced/`, named for the reduction
they achieve — a moderate one and an aggressive one, so the accuracy/size
trade-off is visible rather than asserted:

| Case | Moderate | Aggressive | Format |
|---|---|---|---|
| `ieee37` | `63pct` — Ē=0.003, 38 → 14 | `76pct` — Ē=0.009, 38 → 9 | CSV |
| `ieee123` | `59pct` — Ē=0.002, 126 → 52 | `77pct` — Ē=0.005, 126 → 29 | CSV |
| `european906` | `79pct` — Ē=0.001, 907 → 191 | `89pct` — Ē=0.004, 907 → 100 | CSV |
| `case69` | `59pct` — Ē=0.001, 69 → 28 | `77pct` — Ē=0.006, 69 → 16 | MATPOWER |
| `case85` | `59pct` — Ē=0.002, 85 → 35 | `78pct` — Ē=0.006, 85 → 19 | MATPOWER |
| `case141` | `74pct` — Ē=0.001, 141 → 37 | `85pct` — Ē=0.003, 141 → 21 | MATPOWER |
| `case533mt` | `66pct` — Ē=0.001, 533 → 180 | `83pct` — Ē=0.003, 533 → 88 | MATPOWER |
| `case1197` | `50pct` — Ē=0.006, 1197 → 593 | `65pct` — Ē=0.009, 1197 → 417 | MATPOWER |
| `R100` | `82pct` — Ē=0.001, 100 → 18 | `91pct` — Ē=0.005, 100 → 9 | CSV |
| `R300` | `84pct` — Ē=0.001, 300 → 48 | `91pct` — Ē=0.003, 300 → 26 | CSV |

All at `hops = 10`, radiality enforced inside the optimization, every one proven
optimal and inside its budget on the exact nonconvex check. `R100`/`R300` are
enforced on the peak and minimum hours (scenarios 67 and 91); the rest carry a
single loading. `ieee8500` ships no reduced example — at 2501 buses the tight
budgets take tens of minutes, so it stays a live benchmark rather than a
precomputed one.

The aggressive example keeps **at least 30% fewer buses** than the moderate one.
That is the pairing rule, rather than a fixed gap in reduction percentage: once
a feeder is past 80% reduced another ten points is nearly unreachable, while
`R300` going from 48 buses to 26 is obviously a meaningful further step.
Counting buses measures that; counting percentage points does not.

### Two output formats

Single-phase feeders that arrived with line parameters go out as **MATPOWER**
`.m` files, ready to open and solve, with the bus map in the header comment.

Three-phase feeders have no MATPOWER form — a MATPOWER branch cannot say it
carries `a` and `c` but not `b` — so they go out as **CSV**, in the same schema
the cases themselves use, which `load_case` reads straight back:

| File | What it holds |
|---|---|
| `bus.csv` | the kept buses, under their original ids and phases |
| `ybus.csv` | the reduced admittance — the Schur complement, at full precision |
| `load.csv` | injections after each eliminated bus hands its load to its super-node |
| `voltage.csv` | the full network's voltage at the kept buses: the operating point the budget was certified against |
| `assignment.csv` | the reduction map, one row per original bus |

`assignment.csv` is the one file with no input counterpart, and the one that
makes a reduced network interpretable rather than merely smaller.

```julia
reduced = load_case("data/R300/reduced/91pct")   # 26 buses, 2 scenarios
```

`case1197` is the hardest of these: its low-voltage sections drop enough voltage
that even Ē = 0.009 only reaches 65%. See the note on its impedance spread below.

## MATPOWER cases

Five public distribution feeders converted from MATPOWER with
[`tools/convert_matpower.jl`](../tools/convert_matpower.jl). Each folder keeps
the original `.m` file beside the CSVs, so any conversion can be checked against
its source.

| Case | Source |
|---|---|
| `case69` | Baran & Wu, "Optimal capacitor placement on radial distribution systems," *IEEE Trans. Power Delivery* 4(1), 1989. Derived from a portion of the PG&E system. |
| `case85` | Das et al., radial distribution test system. |
| `case141` | Khodr et al., 141-bus feeder. |
| `case533mt` | Malmer & Thorin, *Network reconfiguration for renewable generation maximization*, 2023. The 533-bus system of the journal paper. |
| `case1197` | Moses et al., hybrid MV/LV system: a 30-bus 22 kV network with 22 copies of a 415 V residential LV network (Western Power, Australia) attached. |

Three things the converter handles that a plain matrix reader would not, and
which are worth knowing when reading the `.m` files:

- **`case69`, `case85`, `case141` state impedances in ohms and loads in kW**,
  converting to per unit in executable MATLAB at the foot of the file. The CSVs
  here are per unit. `case141` additionally splits a kVA magnitude into P and Q
  through a power factor.
- **`case533mt` carries 45 open tie switches** (`status = 0`), dropped on
  conversion: 577 branches become 532 across 533 buses, which is the radial
  configuration the case intends. It also writes base voltages as expressions
  (`135/sqrt(3)`), and its `hi`/`lo` files share one topology, so they are two
  loading scenarios of a single case rather than two cases. `lo` is net
  *generating* (−1.61 MW against `hi`'s +14.87 MW).
- **`case1197` spans two voltage levels** on one 100 MVA base, so its low-voltage
  cables carry per-unit impedances up to 1007. That is arithmetically correct —
  `Z_base` at 415 V is 1.7 mΩ — but it makes the case a genuine conditioning
  stress test rather than a routine feeder.

## Published test feeders

Four three-phase feeders from the IEEE PES Distribution Test Feeder Working
Group and the European LV test case, in the same `bus.csv` + `ybus.csv` form as
`R100`/`R300`. Each also ships `voltage.csv`, the operating point the feeder was
published with; `optikron` reduces against that rather than recomputing it.

| Case | Buses | Node-phase rows | Ē = 0.001 | Ē = 0.01 |
|---|---|---|---|---|
| `ieee37` | 38 | 114 | 23 buses (39.5%) | 8 buses (78.9%) |
| `ieee123` | 126 | 265 | 73 buses (42.1%) | 17 buses (86.5%) |
| `european906` | 907 | 2721 | 235 buses (74.1%) | 143 buses (84.2%) |
| `ieee8500` | 2501 | 3786 | 1672 buses (33.1%)\* | 326 buses (87.0%) |

At `hops = 5`, radiality enforced inside the optimization, all within budget on
the exact nonconvex check. Every entry is proven optimal except the one marked
\*, where the solver stopped at a 1800 s limit holding a feasible reduction
rather than a proven-best one — see below.

Three things to know about these:

- **`ieee8500` does not hold 8500 buses.** The name counts *nodes* the way the
  source case does, including neutral and secondary service conductors. Our
  phase model carries `a`/`b`/`c` only, so the same feeder arrives as 2501 buses
  over 3786 node-phase rows.
- **`european906` and `ieee8500` run voltage regulators.** Their shipped
  `voltage.csv` reflects tap positions that a constant-power wye power flow
  cannot reproduce, which is exactly why the operating point is an input to the
  reduction rather than something it derives. `ieee123`'s shipped voltages,
  which have no such boost, match `powerflow(net)` to 9.1e-08 pu.
- **`ieee8500` at `Ē = 0.001` is where tightening the budget stops being free.**
  Screening is what decides the model's size, and it is far less effective at
  the tighter budget: 363 constraint rows survive at `Ē = 0.01` against 27510 at
  `Ē = 0.001`, a 75x larger model from the same feeder. The loose budget then
  solves to proven optimality in 2.3 s; the tight one is still searching at
  1800 s. Both builds cost about 25 s. If you need `Ē = 0.001` here, budget
  solver time in the tens of minutes, or accept the feasible answer — it is
  inside the budget either way.

`ieee37` is delta-connected and three-phase at every bus, so phase availability
never binds there; `ieee123` and `ieee8500` are heavily single-phase-laterals,
where it does.

### Equipment kept through the reduction

These four also ship `transformer.csv` and `switch.csv`. A transformer,
regulator or switch is a modelled device, not just impedance, and none of its
tap ratio, winding connection or phase shift survives being folded into a Schur
complement. So `optikron` pins both terminals of each, and the device then comes
through **exactly** — `Y_red[i,j] == Y[i,j]` to the last bit, not approximately.
The reason is in [`src/io/devices.jl`](../src/io/devices.jl): with both terminals
kept, no eliminated bus touches both sides of the device edge, so the Schur
correction at that entry is identically zero.

What is *not* preserved is loading. The assignment relocates injections, so the
current through the device changes; the error budget is what bounds that.

| Case | Devices kept | Buses pinned | Cost at `Ē = 0.01`, `hops = 5` |
|---|---|---|---|
| `ieee37` | 2 transformers | 4 | 8 → 9 buses |
| `ieee123` | 1 transformer, 6 switches | 13 | 17 → 24 buses |
| `european906` | 1 transformer | 2 | none measured |
| `ieee8500` | 1 transformer, 25 switches | 52 | — |

Control it with `preserve`: `:devices` (default), `:transformers` (merge across
closed switches — a closed switch is a near-zero-impedance jumper whose ends are
electrically one node, so merging it is nearly free and removes a real numerical
hazard), or `:none`. On `ieee123` that ladder is 24 / 19 / 17 buses kept.

Two filters matter when reading these tables. Only transformers marked
`model_scope = explicit_ybus` are edges in our Ybus at all — `ieee8500` ships
1178 transformers of which **1177 are `collapsed_primary_load`**, service
transformers already folded into their primary at export, with no secondary bus
here to pin. Preserving on the name alone would freeze 1140 of 2501 buses
(45.6%) instead of 2 (0.08%), capping reduction near 54% and destroying the 87%
result above.

## R100 and R300

Radial three-phase distribution feeders at 100 and 300 buses. Both are genuinely
unbalanced: buses carry either all three phases or a single phase, so a reduction
has to respect phase availability — a bus can only be absorbed into one carrying
all of its phases. This is what makes them useful test cases rather than
single-phase equivalents in disguise.

| | `R100` | `R300` |
|---|---|---|
| Buses | 100 | 300 |
| Node-phase rows | 166 | 506 |
| Three-phase buses | 33 | 103 |
| Single-phase buses | 67 | 197 |
| Ybus nonzeros | 1,023 | 3,279 |
| Maximum depth from slack | 17 | 42 |

Note the phase mix: most buses are single-phase laterals hanging off a
three-phase backbone. That asymmetry is the point — it constrains which
assignments are legal and is invisible in a balanced single-phase model.

Each case carries **168 hourly loading scenarios** spanning one week. The
reduction methods are designed to run against a small representative subset and
then generalize; shipping the full series is what makes that test reproducible.

| Scenario | Meaning |
|---|---|
| `h067` | peak loading of the week |
| `h091` | minimum loading of the week |
| `h001`…`h168` | full hourly series |

Buses are numbered `B1..Bn`, with `B1` the slack.

### Format

These carry **no per-line parameters** — the underlying model is an assembled
admittance matrix rather than a branch table — so they ship as `bus.csv` +
`ybus.csv` and load through the Ybus fast path. Everything works on that path
except exporting a reduced network as a feeder model and the `:Ladder` MILP
form. See [`src/io/from_ybus.jl`](../src/io/from_ybus.jl).

`ybus.csv` holds sparse triplets over node-phase rows (`row, col, g, b`), not
bus indices — a three-phase bus occupies three consecutive rows. `bus.csv` is
what maps rows back to buses.
