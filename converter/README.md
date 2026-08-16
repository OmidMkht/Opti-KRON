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

- **`ieee34`, `ieee37`, `ieee123`, `european_lv`** — full `dss/` output beside
  the CSVs. Eight reductions in all.
- **`ieee8500`** — OpenDSS origin, but its shipped reductions cannot be converted;
  see below.
- **`R100`, `R300`** — no OpenDSS origin, CSV only.
- **`case69`, `case85`, `case141`, `case533mt`, `case1197`** — single-phase
  MATPOWER, exported as `.m` instead.

The CSV files are the primary artefact throughout; `dss/` sits beside them and
never replaces them.

### The `ieee8500` limitation

The converter requires both terminals of a preserved transformer to survive the
reduction, so it can re-instantiate it as a real OpenDSS `Transformer`. Opti-KRON
deliberately does **not** pin plain transformers — on `ieee8500` the 2354 service
transformers would fix 2400 buses, 49% of the feeder, as a hard ceiling on
reduction — so the shipped reductions have absorbed them into the Kron equivalent
and there is nothing left to re-instantiate. The conversion stops with

```text
Transformer Transformer.t21396254a has non-surviving terminal buses ['l2804253']
```

which is the converter refusing to invent a circuit it cannot justify, not a bug.
The two are individually correct and jointly incompatible on this feeder.

Preserving the transformers resolves it and costs 11–13 points of reduction
(42%→31% at Ē=0.005, 54%→41% at Ē=0.02, measured with the GPU search). The
shipped levels keep the reduction; pass
`preserve = (:phase_shift, :regulator, :switch, :transformer)` if you want the
circuit instead.

A deeper fix — emitting the absorbed transformer's equivalent as a branch — is
not obviously sound: eliminating a service transformer merges an MV bus with an
LV one, and the resulting super-node has no single base kV for OpenDSS to carry.
That is left open rather than guessed at.

## The round trip is a harder test than the budget

Each generated snapshot is solved and every surviving super-node compared with
its original voltage, in `dss/validation.csv` and `dss/report.json`. That number
is *not* the same as the reduction's certified budget, and on some feeders it is
much larger:

| Case | Level | Ē | Solved DSS error | |
|---|---|---|---|---|
| `european_lv` | 91pct | not binding | 1.6e-04 | well inside |
| `european_lv` | 98pct | not binding | 5.1e-04 | well inside |
| `ieee123` | 57pct | 0.002 | 9.4e-04 | inside |
| `ieee123` | 78pct | 0.006 | 5.2e-03 | inside |
| `ieee34` | 51pct | 0.003 | 9.4e-03 | **3× over** |
| `ieee34` | 59pct | 0.010 | 1.0e-02 | at the edge |
| `ieee37` | 62pct | 0.003 | 8.4e-02 | **28× over** |
| `ieee37` | 74pct | 0.008 | 8.4e-02 | **11× over** |

This is not a conversion defect, and the reports rule that out directly: power
aggregation is exact to `0.0` on every case, and the synthesized Kron `Ybus`
reproduces the reduction's own to `2.9e-16` on `ieee37` — the *best* of the six,
against the worst voltage error.

It is the modelling gap. Opti-KRON certifies its budget against a **constant
current** linearisation at a fixed operating point; OpenDSS then re-solves the
reduced circuit with **constant power** loads and the real winding connections.
The gradient across the table tracks exactly how far each feeder sits from that
assumption: `european_lv` is 98% unloaded so the distinction barely exists,
`ieee123` is wye-connected and holds, and `ieee37` — the one all-delta IEEE
feeder, where load responds to line-to-line voltage — is worst by an order of
magnitude.

The README already says the operating point is a modelling choice and that
"a feeder with delta connections, ZIP loads or regulator taps belongs in a tool
that models them". These numbers are what that sentence costs in practice, and
they are published rather than hidden: use `ieee37`'s reduced snapshot knowing
its solved error is 8%, not 0.3%.

## Documentation

- [Dataset exporter](docs/DATASET_EXPORTER.md)
- [Reduced CSV dataset to DSS](docs/REDUCED_DATASET_TO_DSS.md)
- [Reduced DSS converter](docs/REDUCED_DSS_CONVERTER.md)
- [CSV format](docs/CSV_FORMAT.md)
