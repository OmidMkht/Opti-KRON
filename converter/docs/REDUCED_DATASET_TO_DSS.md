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
- `ybus.csv`
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

Thus solved constant-PQ load, generation, storage, and other net injections at
an eliminated bus move to its super-node. Positive CSV power is network
injection; the OpenDSS `Load` sign is reversed when writing `Loads.dss`.

The supplied `ybus.csv` is used for the radial network. It is not rebuilt from
the bus assignment. Passive shunts and capacitor-bank effects already present
in the full passive Y-bus are therefore embedded in the supplied Kron Y-bus;
they are not written a second time as capacitors. `equipment_mapping.csv`
records each capacitor, switch, transformer, regulator, and phase-shift
terminal's original and assigned bus for audit.

Physical transformers are a special case. Their terminal buses must survive.
The converter subtracts their original YPrim contribution from the synthetic
Kron network and writes the solved physical transformer definitions and taps
to `Transformers.dss`. This preserves winding connection, voltage ratio, and
tap behavior without double counting their admittance.

All remaining radial off-diagonal blocks are written as matrix-valued Reactor
branches; diagonal residuals are written as grounded matrix shunts. Tiny
cross-phase terms below `1e-7` of an edge's dominant coupling can be omitted
when needed to represent a one-phase edge with compatible OpenDSS terminals.

## Power flow and voltage error

The original full circuit provides physical bus voltage bases and Vsource
properties. The reduced per-unit Y-bus is converted back to SI, aggregated
constant-PQ devices are written, controls are disabled, and OpenDSS solves a
snapshot.

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

## Output files

- `Master.dss`
- `Source.dss`
- `Transformers.dss`
- `Branches.dss`
- `Shunts.dss`
- `Loads.dss`
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

