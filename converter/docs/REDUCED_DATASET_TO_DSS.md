# Reduced CSV Dataset to OpenDSS

## Purpose

`convert_reduced_dataset.py` converts a supplied radial reduced CSV dataset
into a self-contained OpenDSS snapshot, runs its power flow, and compares the
solved voltage magnitude with the original voltage at every surviving
super-node.

This workflow is distinct from `reduce_dss.py`: the assignment and Kron Y-bus
have already been computed by the reduction method and are treated as inputs.

## Inputs and assignment semantics

```powershell
python convert_reduced_dataset.py FULL_MASTER FULL_DATASET REDUCED_DATASET OUTPUT `
  --sbase-mva 5
```

The reduced directory must contain:

- `assignment.csv`: `bus_id,super_node,kept`
- `bus.csv`
- `ybus.csv` (optional when the full dataset Y-bus is available)
- `load.csv`
- `voltage.csv`

For every original bus `i`, `super_node` is the surviving bus `j` to which it
is assigned. A kept row must itself appear in the reduced bus list. Child
phases must be a subset of the super-node phases.

## Relocation and Kron-network rules

The converter verifies, phase by phase, that reduced power is exactly

```text
S_reduced[j,phase] = sum(S_full[i,phase] for i assigned to j)
```

This is an input-consistency check used by the reduction method; it does not
replace the models in the generated circuit. The default `original-models`
mode recompiles the full circuit, relocates each regular device, and preserves
its model in per unit. Device name, load model, nominal power, connection,
phase terminals, per-unit voltage limits, and supporting curves remain the
same. For example, a model-2 load on `701.1.2` remains a model-2 phase-to-phase
delta load after it moves to its super-node.

If the original and mapped buses have different physical voltage bases, define

```text
r = Vbase_new / Vbase_old
```

The converter rebases physical parameters while preserving per-unit behavior:

```text
V_rated,new = r V_rated,old
I_rated,new = I_rated,old / r
Z_new       = r^2 Z_old
Y_new       = Y_old / r^2
```

Power ratings, ZIP fractions, load/generator model numbers, connections, and
per-unit voltage thresholds are unchanged. This covers standard OpenDSS
`Load`, `Generator`, `PVSystem`, `Storage`, `Isource`, `VCCS`, `Capacitor`, and
one-terminal `Reactor` models. In particular, changing a load's rated `kV`
lets OpenDSS reproduce its original constant-PQ, constant-current,
constant-impedance, or ZIP law at the new physical voltage base.

Supporting load shapes, growth shapes, temperature shapes, XY curves, and
spectra are copied when they exist.

The scalar/diagonal transformation is deliberately limited to regular
devices. Before conversion, the code divides the original circuit into phase
reference domains. Multi-winding transformers, connection-changing
transformers, and equipment listed in `phase_shift_equipment.csv` form domain
boundaries. An assignment crossing one of those detectable boundaries is
rejected with the offending bus mappings. The input contract additionally
forbids bank-level open-delta or zigzag conversions that are not represented
as a connection-changing OpenDSS transformer object or listed in
`phase_shift_equipment.csv`.

The constant-PQ assumption therefore applies to the ordered `load.csv` used by
the reduction calculations, not to the final DSS device definitions. Two
diagnostic compatibility modes remain available: `mapped-elements` freezes
each solved load as model-1 PQ while retaining its connection, and
`nodal-wye` realizes the phase powers in `load.csv` as independent wye loads.

When `ybus.csv` is supplied, it is used directly. When it is omitted, the
converter reconstructs the conventional Schur-complement Y-bus from the full
ordered matrix and retained phase nodes. The calculation uses sparse LU and
column blocks, applies a relative `1e-12` sparsification threshold, always
retains the diagonal, verifies sparse-LU rank, and writes the compact result as
`ybus.csv` with `ybus_rebuild.json` metadata.

Passive shunts and capacitor-bank effects are already present in the supplied
or rebuilt Kron Y-bus. In `original-models` mode, the converter
subtracts each relocated physical shunt's YPrim from the synthetic equivalent
before writing the original shunt object. The total generated Y-bus therefore
remains equal to the supplied reduced Y-bus without counting these devices
twice. `equipment_mapping.csv` records each capacitor, switch, transformer,
regulator, phase-shift terminal, and center-tapped winding's original and
assigned bus for audit.

Physical transformers are a special case. Each terminal is moved through the
assignment when necessary. For an allowed ordinary two-winding crossing, each
winding `kV` is rebased using its terminal's voltage-base ratio. The original
YPrim is transformed as `Ynew = D Yold D`, stamped at the mapped nodes, and
subtracted from the synthetic Kron network. Solved taps and per-unit winding
impedances are retained. Crossings that require a non-diagonal voltage
transformation are rejected by the phase-domain check.

All remaining radial off-diagonal blocks are written as matrix-valued Reactor
branches; diagonal residuals are written as grounded matrix shunts. Tiny
cross-phase terms below `1e-7` of an edge's dominant coupling can be omitted
when needed to represent a one-phase edge with compatible OpenDSS terminals.

## Power flow and voltage error

The original full circuit provides physical bus voltage bases, Vsource
properties, equipment definitions, and the snapshot load multiplier. The
reduced per-unit Y-bus is converted back to SI, relocated devices are written,
controls are disabled, and OpenDSS solves a snapshot.

The target is the full-network voltage stored for each surviving super-node.
For every retained phase:

```text
error_mag_pu = abs(abs(V_reduced_pu) - abs(V_supernode_target_pu))
```

`validation.csv` contains every target, solution, magnitude error, and complex
error. `report.json` reports the maximum magnitude error with its bus and phase.
The converter also reassembles the generated passive OpenDSS Y-bus and compares
it with the supplied Kron matrix, ensuring the reported power-flow error is not
silently caused by a bad DSS serialization.

The report's `equipment_per_unit_transformation` section lists transformed
regular-device counts and voltage-base-ratio examples. `voltage_level_audit`
records cross-level terminal counts, transformed transformer windings, and the
phase-domain validation result.

## Output files

- `Master.dss`
- `Source.dss`
- `Transformers.dss`
- `Branches.dss`
- `Shunts.dss`
- `PhysicalShunts.dss`
- `Loads.dss`
- `Injections.dss`
- `Models.dss`
- `bus_mapping.csv`
- `equipment_mapping.csv`
- `validation.csv`
- `report.json`

## IEEE 123 examples

```powershell
python convert_reduced_dataset.py data\test_cases\123Bus\IEEE123Master.dss `
  data\datasets\ieee123 data\datasets\ieee123\reduced\57pct `
  data\datasets\ieee123\reduced\57pct\dss --sbase-mva 5

python convert_reduced_dataset.py data\test_cases\123Bus\IEEE123Master.dss `
  data\datasets\ieee123 data\datasets\ieee123\reduced\78pct `
  data\datasets\ieee123\reduced\78pct\dss --sbase-mva 5
```

## Verified IEEE 34 and IEEE 37 results

All values below are the maximum retained-node error
`abs(abs(V_reduced) - abs(V_full))`; the base is 2.5 MVA. The generated
passive OpenDSS Y-bus agrees with the supplied reduced Y-bus to machine
precision.

| Case | Full -> reduced buses | Phase nodes | Max error (pu) | Percent | Location |
|---|---:|---:|---:|---:|---|
| IEEE 34, 51pct | 37 -> 18 | 50 | 0.000155104 | 0.015510% | 834.b |
| IEEE 34, 59pct | 37 -> 15 | 43 | 0.000493878 | 0.049388% | 832.b |
| IEEE 37, 62pct | 39 -> 15 | 45 | 0.000411213 | 0.041121% | 720.a |
| IEEE 37, 74pct | 39 -> 10 | 30 | 0.001700091 | 0.170009% | 734.a |

The corresponding self-contained circuits and detailed per-node validation
tables are under each reduced dataset's `dss` directory.

## IEEE 8500 assignment status

The two supplied IEEE 8500 reductions omit their dense Y-buses. Compact,
full-rank Kron matrices were previously reconstructed for them, but both
assignments cross center-tapped three-winding service-transformer boundaries.
The current converter therefore rejects both before writing a DSS circuit.

| Case | Buses | Phase nodes | Sparse Y nnz | Rejected bus assignments | Reason |
|---|---:|---:|---:|---:|---|
| 42pct | 2821 | 5324 | 28324 | 4 | Crosses center-tapped three-winding service transformers |
| 54pct | 2224 | 4073 | 20613 | 170 | Crosses center-tapped three-winding service transformers |

The `dss` folders currently stored below these reductions are historical
location-only diagnostic outputs and are not outputs accepted by the present
converter. New assignments must keep each center-tapped service-secondary
domain intact. Ordinary scalar voltage-ratio crossings through two-winding,
same-connection transformers remain permitted.
