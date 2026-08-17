"""Emit an OpenDSS snapshot beside every reduced dataset that has an OpenDSS origin.

Writes `reduced_networks/<case>/<level>/dss/`, solves each one, and reports the
voltage error at every surviving super-node against the original feeder.

Only the feeders that came from OpenDSS can be done: the reduction map alone is
not enough, the converter needs the original master to recover line geometry,
transformer windings and regulator settings that a Ybus cannot carry.

    python converter/build_reduced_dss.py            # all of them
    python converter/build_reduced_dss.py ieee123    # one case
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

CONVERTER = Path(__file__).resolve().parent
REPO = CONVERTER.parent
CASES_DIR = CONVERTER / "opendss_cases"

# case -> (master relative to opendss_cases, sbase_mva, load representation)
#
# `original-models` relocates each original device untouched, so a delta load
# stays delta -- worth 200x on ieee37, which is all-delta. It is only right where
# the original load objects already state the operating point. european_lv's do
# not: every load there is a nominal `kW=1 PF=0.95` scaled by its own
# `Yearly=Shape_N` profile, and its master runs a 1440-step yearly simulation, so
# relocating those originals reproduces no particular snapshot. That feeder takes
# the powers from the dataset instead, which is what its CSV actually pins down.
SOURCES = {
    # ieee34Mod1.dss is the circuit itself; Run_IEEE34Mod1.dss is an interactive
    # demo that compiles it twice and issues plot commands.
    "ieee34": ("34Bus/ieee34Mod1.dss", 2.5, "original-models"),
    "ieee37": ("37Bus/ieee37.dss", 2.5, "original-models"),
    "ieee123": ("123Bus/IEEE123Master.dss", 5.0, "original-models"),
    "european_lv": ("LVTestCase/Master.dss", 0.8, "nodal-wye"),
    "ieee8500": ("8500-Node/Master.dss", 27.5, "original-models"),
}


def main() -> int:
    # `--reduced DIR --case NAME` converts one directory, which is what
    # run_optikron.jl calls after exporting. Bare names convert every level of
    # those cases under reduced_networks/.
    single = None
    argv = sys.argv[1:]
    if "--reduced" in argv:
        i = argv.index("--reduced")
        single = Path(argv[i + 1])
        del argv[i:i + 2]
        i = argv.index("--case")
        argv_case = argv[i + 1]
        del argv[i:i + 2]
        if argv_case not in SOURCES:
            print(f"{argv_case}: no OpenDSS origin, nothing to write")
            return 0
        argv = [argv_case]

    wanted = argv or sorted(SOURCES)
    rows, failures = [], []

    for case in wanted:
        if case not in SOURCES:
            print(f"{case}: no OpenDSS origin, skipping")
            continue
        master_rel, sbase, representation = SOURCES[case]
        master = CASES_DIR / master_rel
        full = REPO / "data" / case
        case_dir = REPO / "reduced_networks" / case
        if single is not None:
            levels = [single]
        elif case_dir.is_dir():
            levels = sorted((p for p in case_dir.iterdir() if p.is_dir()),
                            key=lambda p: p.name)
        else:
            print(f"{case}: no reductions yet, skipping")
            continue

        for reduced in levels:
            level = reduced.name
            out = reduced / "dss"
            # One conversion per process. OpenDSS is a global singleton behind
            # opendssdirect and `ClearAll` does not fully reset it -- line codes
            # and wire data from a previous feeder survive and collide, which
            # surfaces as "Y matrix build aborted" on the *second* case in a run
            # even though it converts cleanly on its own.
            proc = subprocess.run(
                [sys.executable, str(CONVERTER / "convert_reduced_dataset.py"),
                 str(master), str(full), str(reduced), str(out),
                 "--sbase-mva", str(sbase),
                 "--load-representation", representation],
                capture_output=True, text=True,
            )
            if proc.returncode != 0:
                detail = (proc.stderr.strip().splitlines() or ["unknown error"])[-1]
                failures.append((case, level, detail[:110]))
                print(f"  {case:<12} {level:<7} FAILED  {detail[:80]}")
                continue
            report = json.loads(proc.stdout)
            rows.append((case, level, report))
            # Report the error over *energized* nodes. The raw maximum is set by
            # near-zero-voltage neutral-like conductors that the source dataset
            # carries as phase rows -- on ieee8500 a 0.06 pu node dominates it and
            # says nothing about the reduction.
            err = report.get("max_energized_voltage_magnitude_error_pu",
                             report["max_voltage_magnitude_error_pu"])
            raw = report["max_voltage_magnitude_error_pu"]
            tail = "" if abs(raw - err) < 1e-12 else f"  (raw {raw:.2e})"
            # Three different things, only two of which are problems.
            # `location_only` false just means a device changed voltage level;
            # that is fine when the converter rebased it, which
            # `per_unit_transformed_equivalent_valid` reports. A crossed phase
            # domain -- a center-tapped transformer or phase shifter -- is not
            # recoverable by rebasing and is a real failure.
            audit = report.get("voltage_level_audit", {})
            domains = audit.get("incompatible_phase_domain_assignment_count", 0)
            if domains:
                flag = f"  CROSSED {domains} PHASE DOMAIN(S)"
            elif not audit.get("per_unit_transformed_equivalent_valid", True):
                flag = "  INVALID PER-UNIT EQUIVALENT"
            elif not audit.get("location_only_voltage_level_valid", True):
                n = audit.get("transformed_transformer_winding_count", 0)
                flag = f"  ({n} winding(s) rebased)"
            else:
                flag = ""
            print(
                f"  {case:<12} {level:<7} {report['full_bus_count']:>5} -> "
                f"{report['reduced_bus_count']:<5} "
                f"max |V| err {err:.3e} pu{tail}{flag}"
            )

    # A one-directory run drops its summary beside that directory; a full sweep
    # writes the index for the whole collection.
    summary = (single / "dss_build.json") if single is not None \
        else (REPO / "reduced_networks" / "dss_build.json")
    summary.write_text(
        json.dumps(
            {
                "converted": [
                    {"case": c, "level": l, **r} for c, l, r in rows
                ],
                "failed": [
                    {"case": c, "level": l, "error": e} for c, l, e in failures
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"\n{len(rows)} converted, {len(failures)} failed -> {summary}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
