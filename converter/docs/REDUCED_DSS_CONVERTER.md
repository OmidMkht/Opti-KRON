# Full-to-Reduced OpenDSS Converter

## Purpose and exactness contract

`reduce_dss.py` accepts a full solved OpenDSS circuit and a bus-level binary
assignment matrix `A`, performs a phase-domain Kron reduction, and writes a
self-contained radial reduced OpenDSS circuit.

The assignment convention is:

```text
A[i,j] = 1  means original bus j is represented by super-bus i
A[i,i] = 1  means bus i survives
```

The generated constant-PQ equivalents are fitted so the original converged
complex voltages at every surviving phase node are a solution of the reduced
circuit. This is an exact **snapshot** target, subject to OpenDSS and numerical
tolerances. It does not promise exact voltages after loads, controls, topology,
or time-series conditions change.

## Installation and implementation

```powershell
python -m pip install -r requirements.txt
```

Maintained files:

- CLI: `reduce_dss.py`
- library: `opendss_tools/converter.py`
- integration test: `tests/test_converter.py`

## Inputs

### Full OpenDSS circuit

The first argument is the entry-point `.dss` file. It must compile, establish
valid voltage bases, solve in snapshot mode, and currently contain exactly one
OpenDSS `Vsource`.

### `assignment.csv` or `.npy`

`A` must be a dense square binary matrix with one row and column per physical
bus. CSV files may contain a header and optional row-label column; `.npy` files
are also accepted.

The converter enforces:

- shape `N x N`, where `N` is the number of actual OpenDSS buses;
- binary entries within `--matrix-tolerance`;
- exactly one `1` in every column;
- every nonzero assignment points to a surviving diagonal row;
- a child's active phases are a subset of its super-bus phases;
- the source bus survives and maps to itself;
- every transformer/regulator terminal bus survives and maps to itself.

The supplied assignment must lead to a connected radial reduced admittance
graph. The converter validates this; it does not search for or optimize `A`.

### `bus_order.csv`

This file contains one `bus_id` per row in exactly the row/column order of `A`.
A header named `bus`, `bus_id`, or `name` is optional. Names are normalized to
lower case and stripped of phase suffixes. The set must exactly match the
compiled OpenDSS bus set—missing, extra, or duplicate names are rejected.

Example:

```csv
bus_id
sourcebus
1
2
3
```

## Usage

```powershell
python reduce_dss.py <full-master.dss> <assignment.csv> `
  <bus_order.csv> <output-directory>
```

Example:

```powershell
python reduce_dss.py data\test_cases\123Bus\IEEE123Master.dss `
  assignment.csv bus_order.csv data\reduced\ieee123
```

Numerical controls:

```powershell
python reduce_dss.py FULL\Master.dss assignment.csv bus_order.csv reduced `
  --matrix-tolerance 1e-9 --max-voltage-error-pu 1e-6
```

| Argument | Default | Meaning |
|---|---:|---|
| `--matrix-tolerance` | `1e-9` | Binary, graph-block, reciprocity, and decomposition tolerance |
| `--max-voltage-error-pu` | `1e-6` | Maximum accepted solved voltage-magnitude error |

If final validation exceeds the voltage limit, the command fails but retains
the output files for diagnosis.

## Mathematical procedure

Let the original ordered phase nodes be partitioned into surviving nodes `K`
and eliminated nodes `R`. The passive network is reduced as

```text
Yred = YKK - YKR pinv(YRR) YRK
```

where `pinv` is the numerical pseudoinverse. Transformers and autotransformers
are not absorbed into a same-voltage series equivalent. Their terminal buses
must survive, their solved definitions/taps are copied explicitly, and their
admittance contribution is added back at the retained nodes.

For the original surviving voltage vector `VK`, the converter fits

```text
Ieq = Yred VK
Seq = VK conj(Ieq)
```

and writes a one-phase constant-PQ `Load` equivalent for each nonzero phase
injection. OpenDSS loads use the opposite sign from network injection, so the
written kW/kvar values are `-Seq`. Source-bus fitted injection is set to zero
because the retained voltage source supplies the slack current.

The assignment lift is used to validate bus/phase aggregation, but direct
`A I` injection aggregation is intentionally not used to construct the final
snapshot. In general, direct aggregation is an approximation and does not make
the original `VK` an exact solution of the Kron network.

## Turning the reduced Y-bus into OpenDSS elements

The reduced off-diagonal bus blocks define graph edges. Each reciprocal branch
admittance block is inverted with a pseudoinverse and written as a matrix-valued
OpenDSS `Reactor` using `rmatrix` and `xmatrix`. Remaining diagonal admittance
is written as a grounded matrix-valued shunt reactor.

This representation preserves general phase coupling. The equivalents should
not be interpreted as reconstructed physical line geometry, conductor length,
or equipment nameplate data.

The converter rejects:

- a disconnected or meshed reduced graph;
- incompatible phase dimensions across an edge;
- non-reciprocal branch blocks;
- an admittance block that cannot be reconstructed from its pseudoinverse
  within tolerance.

## Output circuit

The output directory contains:

| File | Purpose |
|---|---|
| `Master.dss` | Self-contained entry point; disables controls and solves a snapshot |
| `Source.dss` | Retained Vsource properties |
| `Transformers.dss` | Explicit preserved transformers and fixed solved taps |
| `Branches.dss` | Matrix-valued Kron-equivalent series reactors |
| `Shunts.dss` | Matrix-valued grounded shunt reactors |
| `Loads.dss` | Operating-point constant-PQ phase equivalents |
| `bus_mapping.csv` | Original bus, assigned super-bus, and survivor flag |
| `validation.csv` | Target and solved complex voltages and per-node errors |
| `report.json` | Full/reduced sizes and maximum voltage errors |

`bus_mapping.csv` is the audit trail for `A`. `validation.csv` reports target
and reduced rectangular voltage values, absolute complex error in volts, and
relative magnitude error. `report.json` summarizes the maximum values.

## Validation sequence

The converter does not stop after writing algebraic equivalents. It:

1. writes all reduced DSS files;
2. compiles the new `Master.dss` in a fresh OpenDSS context;
3. solves the reduced snapshot;
4. confirms that every surviving phase node exists;
5. compares solved voltages to the original target `VK`;
6. enforces `--max-voltage-error-pu`.

This catches DSS serialization, phase-terminal, sign, and element-decomposition
errors that an in-memory Kron check would miss.

## Controls and operating-state limitations

The output deliberately uses `controlmode=off`. Regulator transformers are
preserved at their solved taps, but `RegControl` logic is not transferred.
Capacitor controls, switch automation, relays, fuses, reclosers, load shapes,
ZIP dependence, generator controls, storage dispatch, and protection behavior
are not dynamically preserved. Capacitors and other passive shunts present in
the solved network are embedded in the reduced admittance.

For a different scenario, solve the full circuit at that scenario and generate
a new reduced snapshot. A multi-scenario or control-preserving equivalent is a
different modeling problem from this converter.

## Current supported scope

- Unbalanced one-, two-, and three-phase nodes numbered 1, 2, and 3
- One OpenDSS Vsource
- Bus-level `A` with phase-compatible assignments
- Explicitly retained transformer/regulator terminal buses
- A supplied assignment whose reduced graph is radial
- Snapshot-equivalent constant-PQ injections

## Common failures

- **Bus-order mismatch:** generate `bus_order.csv` from the actual compiled
  OpenDSS circuit, not the feeder's nominal bus count.
- **Transformer terminal eliminated:** change `A` so every transformer and
  regulator endpoint survives and maps to itself.
- **Phase incompatibility:** never map a child phase to a super-bus that lacks
  that phase.
- **Reduced graph is not a tree:** revise `A`; the converter does not radialize
  a meshed result automatically.
- **Voltage error too large:** inspect `validation.csv`, confirm the full model
  converged at the intended operating point, and consider a tighter/cleaner
  assignment or numerical tolerance.

## Running tests

```powershell
python -m unittest discover -s tests -v
```
