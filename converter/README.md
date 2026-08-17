# OpenDSS converter

The bridge between OpenDSS and this package, in both directions. It is what
produced the three-phase feeders in [`data/`](../data), and what turns a
reduction back into a circuit you can open in OpenDSS.

Python rather than Julia because OpenDSS is driven through `opendssdirect.py`;
there is no Julia binding that exposes the element-level detail this needs.

```bash
python -m pip install -r converter/requirements.txt
python -m unittest discover -s converter/tests -v
```

## What it does

| Entry point | Direction |
|---|---|
| [`export_dss_dataset.py`](export_dss_dataset.py) | OpenDSS feeder → ordered, phase-consistent per-unit CSV, verifying `YV = I` |
| [`build_reduced_dss.py`](build_reduced_dss.py) | every reduction in `reduced_networks/` → a solved OpenDSS snapshot |
| [`convert_reduced_dataset.py`](convert_reduced_dataset.py) | one reduced CSV dataset → OpenDSS, solved, with super-node voltage errors |
| [`reduce_dss.py`](reduce_dss.py) | OpenDSS feeder + an assignment → a radial reduced OpenDSS snapshot |

```bash
python converter/build_reduced_dss.py            # rebuild every dss/ folder
python converter/build_reduced_dss.py ieee123    # just one case
```

`opendss_cases/` holds the five upstream feeders the shipped datasets were
exported from, so the whole pipeline reproduces from this repository alone.

## Why `load.csv` is not the loads

The exporter does not read load objects out of OpenDSS. It rounds `Y` and `V`,
then back-computes

```text
S = V ⊙ conj(Y V)
```

so the serialized dataset satisfies `YV = I` *exactly* — that is the
`serialized_max_abs_yv_minus_i_pu ≈ 7e-16` in each `validation.json`, and it is
what lets a consumer trust the three files against each other.

The cost is that an unloaded node carries round-off rather than a hard zero,
around `1e-15`. This is not removable here: zeroing an injection of magnitude
*f* leaves a residual of *f*/|V| at that node, so the round-trip residual tracks
any floor one-for-one and the exporter's own `1e-12` tolerance fails above about
`1e-13` — well below where the dust ends. Consumers should treat "no load" as a
relative threshold; `zero_injection_buses` in the Julia package uses `1e-9 ×`
the feeder's largest injection.

## Which reductions get an OpenDSS form

Only the feeders that came from OpenDSS. The reduction map and a Kron `Ybus` are
not enough on their own — the converter needs the original master to recover the
line geometry, transformer windings and regulator settings that an admittance
matrix cannot carry.

- **`ieee34`, `ieee37`, `ieee123`, `european_lv`, `ieee8500`** — full `dss/`
  output beside the CSVs. Ten reductions in all.
- **`R100`, `R300`** — no OpenDSS origin, CSV only.
- **`case69`, `case85`, `case141`, `case533mt`, `case1197`** — single-phase
  MATPOWER, exported as `.m` instead.

The CSV files are the primary artefact throughout; `dss/` sits beside them and
never replaces them.

## The round trip, and how it is solved

Each generated snapshot is solved and every surviving super-node compared with
its original voltage, in `dss/validation.csv` and `dss/report.json`. All ten
reductions land inside their budget:

| Case | Level | Ē | Solved DSS error | Audit |
|---|---|---|---|---|
| `ieee37` | 67pct | 0.003 | 5.6e-04 | 1 winding(s) rebased |
| `ieee37` | 77pct | 0.008 | 1.4e-03 | 1 winding(s) rebased |
| `ieee34` | 51pct | 0.003 | 1.7e-04 | location-only |
| `ieee34` | 62pct | 0.01 | 5.0e-04 | location-only |
| `ieee123` | 57pct | 0.002 | 9.1e-04 | location-only |
| `ieee123` | 78pct | 0.006 | 3.9e-03 | location-only |
| `european_lv` | 91pct | 0.01 | 1.5e-04 | location-only |
| `european_lv` | 98pct | 0.01 | 8.8e-04 | location-only |
| `ieee8500` | 31pct | 0.005 | 3.0e-03 | location-only |
| `ieee8500` | 41pct | 0.02 | 2.5e-03 | location-only |

Getting there took two corrections worth recording, because both look like
reduction error and neither is.

**Delta loads must stay delta.** An earlier converter derived loads from the
per-phase injections in `load.csv` and emitted them as wye. That is exact at the
operating point but wrong once you re-solve: a delta load holds constant P,Q
against `V₁−V₂`, a wye load against `V₁` alone. On `ieee37` — the one all-delta
IEEE feeder — the solved error was **8.4e-02, twenty-eight times the budget**,
and the size of the error tracked each feeder's delta fraction almost exactly
(0% on `european_lv`, 9% on `ieee123`, 54% on `ieee34`, all of `ieee37`).
`original-models` fixes it by relocating each original device untouched, changing
only its bus, so a `701.1.2 Conn=Delta Model=2` load stays exactly that. That is
worth 204× on `ieee37` and 60× on `ieee34`.

**But relocation is only right when the original device states the operating
point.** `european_lv`'s loads are `kW=1 PF=0.95 Yearly=Shape_N` — nominal
placeholders scaled by a profile, with the master running a 1440-step yearly
simulation — so relocating them reproduces no particular snapshot and the error
went to 1.9e-02. That feeder takes `nodal-wye`, which reads the powers the
dataset actually pins down. The mode is therefore chosen per feeder in
`build_reduced_dss.py`, with the reasoning beside it.

## Moving equipment across a voltage level

Relocating a device to its super-node is only sound while both sit at the same
nominal voltage. A **constant-impedance or shunt element, a capacitor bank, a ZIP
load, or an inverter** does not carry its behaviour with it: those models are
written against a voltage base, and moving them to a bus at a different base
changes what they do unless every parameter is re-referred through the
transformer ratio.

The converter now does that re-referring for an **ordinary two-winding
crossing**: each winding's kV is rebased by its terminals' voltage ratio, the
original YPrim is transformed as `Ynew = D Yold D` and stamped at the mapped
nodes. `report.json` records it under `equipment_per_unit_transformation`, and
`voltage_level_audit` reports three separate things:

| field | meaning |
|---|---|
| `location_only_voltage_level_valid` | nothing changed level at all |
| `per_unit_transformed_equivalent_valid` | crossings were rebased correctly |
| `incompatible_phase_domain_assignment_count` | a boundary was crossed that **cannot** be rebased |

Only the third is a failure. `ieee37` reduces across its 0.277/2.771 kV winding
and reports `1 winding rebased` — valid, and worth two extra buses of reduction.

### What cannot be crossed

A **center-tapped transformer** is the case no change of base repairs. Its two
secondaries sit on one bus in opposite orientation about the centre conductor
(`.1.0` and `.0.2`), so the primary-to-secondary map is not a diagonal scaling.
Multi-winding transformers, connection-changing transformers and anything in
`phase_shift_equipment.csv` therefore define **phase reference domains**, and an
assignment crossing one is rejected with the offending bus mappings.

`ieee8500` carries 1177 of them across 2315 buses, which is why its reductions
run to 31% and 41% where the smaller feeders reach 60–90%. That is the honest
ceiling for a feeder whose service transformers are split-phase, not a defect.
`center_tapped_transformer.csv` ships beside each dataset so the Julia side can
pin them; `preserve = :required` does so by default.

### Two error figures on `ieee8500`

`report.json` reports both `max_voltage_magnitude_error_pu` and
`max_energized_voltage_magnitude_error_pu`. On `ieee8500` the raw maximum is
6.3e-02 at `e182723.b`, one of ten neutral-like conductors the source dataset
carries as phase rows at ~0.06 pu; the energized maximum, 3.0e-03, is the figure
that describes the reduction. `build_reduced_dss.py` prints the energized one and
notes the raw one in brackets.

## Documentation

- [Dataset exporter](docs/DATASET_EXPORTER.md)
- [Reduced CSV dataset to DSS](docs/REDUCED_DATASET_TO_DSS.md)
- [Reduced DSS converter](docs/REDUCED_DSS_CONVERTER.md)
- [CSV format](docs/CSV_FORMAT.md)
