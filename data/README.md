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
| `ieee34` | 37 | 95 | 1 | Unbalanced three-phase, radial |
| `ieee37` | 39 | 117 | 1 | Three-phase throughout, radial |
| `ieee123` | 130 | 274 | 1 | Unbalanced three-phase, radial |
| `european_lv` | 907 | 2721 | 1 | Three-phase throughout, radial |
| `ieee8500` | 4876 | 8541 | 1 | Unbalanced three-phase, radial |
| `test4` | 4 | 4 | 2 | Single-phase; unit tests |
| `test3ph` | 3 | 7 | 1 | Mixed `abc`/`a`; unit tests |

```julia
net = load_case("case69")                          # picks the right reader for the case
net = load_case("reduced_networks/R300/91pct")     # reduced networks load the same way
```

`R100`/`R300` carry 168 hourly scenarios spanning a week — `h067` is the peak
hour and `h091` the minimum. Their buses are `B1..Bn` with `B1` the slack, and
most are single-phase laterals off a three-phase backbone, which is what makes
phase availability bind. The IEEE feeders each ship `voltage.csv` (the published
operating point, which `optikron` reduces against rather than recomputing),
device tables, and `validation.json` recording the export's own `max |Y·V − I|`.

## Reduced networks

Two worked reductions per case, moved out to
[`reduced_networks/`](../reduced_networks) so the originals and the reductions
stay separate. Each level records the error threshold, hop limit and scenarios it
was produced under, and the violation it was re-checked at.

## Sources

Five MATPOWER feeders, converted with
[`tools/convert_matpower.jl`](../tools/convert_matpower.jl). Each folder keeps
its original `.m` beside the CSVs so the conversion can be checked.

| Case | Source |
|---|---|
| `case69` | Baran & Wu, "Optimal capacitor placement on radial distribution systems," *IEEE Trans. Power Delivery* 4(1), 1989. From a portion of the PG&E system. |
| `case85` | Das et al., radial distribution test system. |
| `case141` | Khodr et al., 141-bus feeder. |
| `case533mt` | Malmer & Thorin, *Network reconfiguration for renewable generation maximization*, 2023. |
| `case1197` | Moses et al., hybrid MV/LV: a 30-bus 22 kV network with 22 copies of a 415 V residential LV network (Western Power, Australia). |

Five three-phase feeders from the IEEE PES Distribution Test Feeder Working Group
and the European LV test case. Note that **`ieee8500` does not hold 8500 buses**:
the name counts nodes including neutral and secondary service conductors, while
this package carries `a`/`b`/`c` only.

Three conversion details worth knowing:

- `case69`, `case85` and `case141` state impedances in ohms and loads in kW,
  converting to per unit in executable MATLAB at the foot of the file. The CSVs
  here are per unit.
- `case533mt` carries 45 open tie switches dropped on conversion (577 branches
  become 532 across 533 buses, the radial configuration intended), and its
  `hi`/`lo` files are two loadings of one topology, `lo` net *generating*.
- `case1197` spans two voltage levels on one 100 MVA base, so its LV cables carry
  per-unit impedances up to 1007. Arithmetically correct — `Z_base` at 415 V is
  1.7 mΩ — but it makes the case a conditioning stress test.
