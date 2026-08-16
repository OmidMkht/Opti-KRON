# Ordered OpenDSS Dataset Exporter

## Purpose

`export_dss_dataset.py` converts a solved, unbalanced OpenDSS distribution
circuit into a phase-consistent dataset. Every electrical table uses one
canonical ordering and one common power base. The exporter also verifies the
network identity

```text
Ypu Vpu = Ipu
```

after serializing and reading the CSV files back from disk.

This is an operating-point export. Native ZIP, constant-current, motor,
generator, storage, and other power-conversion models are first solved by
OpenDSS. Their solved terminal currents are then represented by equivalent
complex power/current injections at that solution. Consequently, `load.csv`
is suitable as a constant-PQ snapshot even when the original elements were not
constant PQ.

## Installation

From the project directory:

```powershell
python -m pip install -r requirements.txt
```

The maintained implementation is:

- CLI: `export_dss_dataset.py`
- library: `opendss_tools/exporter.py`
- integration test: `tests/test_dataset_exporter.py`

## Basic use

```powershell
python export_dss_dataset.py <Master.dss> <output-directory> `
  --sbase-mva <MVA>
```

Example:

```powershell
python export_dss_dataset.py data\test_cases\123Bus\IEEE123Master.dss `
  data\datasets\ieee123 --sbase-mva 5
```

Optional arguments:

| Argument | Default | Meaning |
|---|---:|---|
| `--scenario` | `h001` | Label written into load and voltage rows |
| `--max-iterations` | `100` | OpenDSS power-flow iteration limit |
| `--max-control-iterations` | `100` | Regulator/capacitor control iteration limit |
| `--serialized-tolerance` | `1e-12` | Maximum permitted CSV round-trip `|YV-I|` in pu |
| `--independent-tolerance` | `1e-5` | Maximum permitted independently assembled KCL residual in pu |

The command fails if the OpenDSS solution does not converge or either
validation tolerance is exceeded.

## Canonical ordering

OpenDSS supplies its solved Y-node order. The exporter preserves physical-bus
first-occurrence order and, within each bus, sorts active node numbers as 1, 2,
3 (`a`, `b`, `c`). Ground node 0 is not an exported phase node.

`bus.csv` is the authoritative physical-bus order. Expanding each bus's
`phases` string in `a`, `b`, `c` order produces the authoritative
`node_index`. That node ordering is shared by:

- the row and column indices in `ybus.csv`;
- every row in `voltage.csv`;
- every row in `load.csv`.

Equipment tables reference buses by `bus_id`; lookup uses the ordered
`bus.csv`. Bus names are normalized to lower case and stripped of OpenDSS node
suffixes.

## Per-unit convention

The user supplies one common three-phase/phase-domain power base `Sbase`.
OpenDSS supplies the line-to-neutral voltage base for each bus. For node `i`:

```text
Vpu_i = Vsi_i / Vbase_i
Ibase_i = Sbase / Vbase_i
Ipu_i = Isi_i / Ibase_i
Ypu = diag(Vbase/Sbase) Ysi diag(Vbase)
Spu_i = Vpu_i conj(Ipu_i)
```

This scaling deliberately uses the same `Sbase` for every phase. Therefore,
the sum of phase powers equals the total complex power divided by `Sbase`; no
extra factor of three is required. To match the required sample format,
`bus.csv` reports `base_kv=1.0`: all exported electrical values are already in
per unit.

The signs in `load.csv` are net-injection signs:

- positive `p_pu`/`q_pu`: injection into the passive network;
- negative `p_pu`/`q_pu`: consumption from the passive network;
- slack rows include the source-bus network injection;
- zero-injection nodes are retained.

## Output files

### `bus.csv`

One row per ordered physical bus.

| Column | Meaning |
|---|---|
| `bus_id` | Normalized OpenDSS bus name |
| `phases` | Active phases in canonical `abc` order |
| `base_kv` | `1.0`, because the dataset is fully per unit |
| `type` | `slack` for a Vsource terminal bus, otherwise `pq` |

### `ybus.csv`

Sparse coordinate representation of the passive phase-domain Y-bus.
`row` and `col` are one-based `node_index` values; the entry is
`g + j b`, in per unit. Power-conversion elements and the ideal source are excluded;
passive lines, transformers, reactors, capacitive/shunt primitives, and other
network admittances are included.

### `voltage.csv`

One row per ordered phase node with exactly the identifying fields
`bus_id,phase,scenario` and rectangular values `v_re_pu,v_im_pu`. Magnitude
and angle are deliberately omitted because they are directly derivable.

### `load.csv`

Despite the filename, this is a complete ordered **net-injection** table, not
only a list of OpenDSS `Load` objects. It uses the original five-column schema
`bus_id,phase,scenario,p_pu,q_pu` and retains one row for every ordered node,
including zero and slack injections. Current is reconstructed as
`Ipu = conj(Spu/Vpu)` when validating `Ypu Vpu = Ipu`.

### `transformer.csv`

One row per transformer secondary winding. It retains transformer/winding ID,
endpoint buses and phases, both winding connections, solved `tap_pu`, and
enabled state. The dimensionless `voltage_ratio` is the shifted-winding rated
kV divided by winding-1 rated kV.

### `switch.csv`

One row per detected switch, including endpoint buses/phases, enabled state,
and `open`/`closed` status. Detection uses the OpenDSS switch flag and common
switch naming patterns.

### `bus_coordinates.csv`

One row for every bus in canonical order using only `bus_id,x,y`. Missing `x`
and `y` values are blank; missing coordinates never remove a bus. If no
coordinates were loaded by the master file, the exporter looks in the master
file's directory for a colocated `BusCoords`-named `.csv`, `.dat`, `.txt`, or
`.dss` file and loads it automatically. Coordinates are geometric metadata,
not electrical per-unit quantities.

### `regulator.csv`

One row per OpenDSS `RegControl`. It links the control to its transformer and
terminal buses and includes controlled winding, solved `tap_pu`, `vreg_pu`,
`band_pu`, and enabled state. The per-unit setpoint is calculated against the
regulated bus's line-to-neutral base.

### `capacitor_bank.csv`

One row per capacitor bank/control association. It includes bus, phases,
connection, `rated_q_pu`, step count, solved states, control mode, and enabled
state. Uncontrolled banks have a blank control mode.

### `phase_shift_equipment.csv`

One row for every enabled or disabled three-phase transformer winding pair
whose wye/delta connection differs from winding 1. Rows retain equipment ID,
endpoint buses/phases, both connections, enabled state, and signed
`phase_shift_deg`. The sign convention is the to-bus positive-sequence voltage
angle minus the from-bus angle.
For the OpenDSS default `Lag`/`ANSI` setting, the low-voltage side lags the
high-voltage side by 30 degrees; `Lead`/`Euro` reverses the sign. A header-only
file means no directly identifiable three-phase mixed-connection transformer
exists. General voltage-angle drops caused by impedance are intentionally not
classified as fixed phase-shift equipment.

### `validation.json`

Summary data include convergence, bus/node counts, Y-bus nonzeros, base MVA,
voltage range, nodes below 0.5 pu, and both KCL checks.

## Validation

The exporter performs two distinct checks:

1. **Independent in-memory KCL check.** It assembles the passive Y-bus from
   each enabled element's `YPrim`, reads currents independently from enabled
   power-conversion elements, and compares the two at non-source nodes.
2. **Serialized round-trip check.** It reads `ybus.csv`, `voltage.csv`, and
   `load.csv` back from disk, reconstructs `Ipu=conj(Spu/Vpu)`, and evaluates
   `Ypu Vpu - Ipu`. This detects order, sign, scaling, and serialization mistakes.

The serialized residual is bounded by the eight-decimal quantization tolerance.
The independent residual is the stronger comparison against OpenDSS element
currents and is normally small but nonzero due to solver and floating-point
tolerances.

Coordinates and equipment values are written with at most eight digits after
the decimal point. `ybus.csv`, `voltage.csv`, and `load.csv` retain up to 15
decimal places so small structural admittances are not discarded and the
operating-point identity remains near machine precision. Values are written in
ordinary decimal notation rather than scientific notation. The exporter
serializes Y and V first and then refits the constant-PQ `load.csv` values to
those serialized quantities.

Open switches are never stamped into the exported passive Y-bus. This includes
both conductors opened through the OpenDSS API and the IEEE test-case convention
of parking a switch terminal at a bus whose name ends in `_open`. Those parked
dummy buses are excluded from the canonical bus/node order; `switch.csv` maps
the endpoint back to its logical bus name and records `status=open`.

## Generated test-case datasets

The maintained outputs are under `data/datasets/`:

| Dataset | Buses | Phase nodes | Coordinates | Regulators | Capacitors | Phase shifters | Sbase (MVA) |
|---|---:|---:|---:|---:|---:|---:|---:|
| IEEE 34 | 37 | 95 | 37 | 6 | 2 | 1 | 2.5 |
| IEEE 37 | 39 | 117 | 39 | 2 | 0 | 0 | 2.5 |
| IEEE 123 | 130 | 274 | 130 | 7 | 4 | 0 | 5.0 |
| European LV | 907 | 2721 | 906 | 0 | 0 | 1 | 0.8 |
| IEEE 8500 | 4876 | 8541 | 2459 | 12 | 10 | 1 | 27.5 |

Counts reflect the actual OpenDSS circuit graph after compilation, including
source, regulator-created, secondary, auxiliary, and internal buses. They are
not expected to equal the marketing name of a feeder.

## Important reduction cautions

Before reducing a dataset, preserve or explicitly model source terminals,
transformer ratios/connections, regulator terminal buses and solved taps,
switch status, capacitor states, voltage-base changes, and phase availability.
The snapshot injections reproduce one solved operating point; they do not
preserve future ZIP behavior, control actions, time-series behavior, protection
logic, or regulator/capacitor switching thresholds.

## Running tests

```powershell
python -m unittest discover -s tests -v
```
