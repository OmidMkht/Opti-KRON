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

# case -> (master relative to opendss_cases, sbase_mva)
SOURCES = {
    "ieee34": ("34Bus/Run_IEEE34Mod1.dss", 2.5),
    "ieee37": ("37Bus/ieee37.dss", 2.5),
    "ieee123": ("123Bus/IEEE123Master.dss", 5.0),
    "european_lv": ("LVTestCase/Master.dss", 0.8),
    "ieee8500": ("8500-Node/Master.dss", 27.5),
}


def main() -> int:
    wanted = sys.argv[1:] or sorted(SOURCES)
    rows, failures = [], []

    for case in wanted:
        if case not in SOURCES:
            print(f"{case}: no OpenDSS origin, skipping")
            continue
        master_rel, sbase = SOURCES[case]
        master = CASES_DIR / master_rel
        full = REPO / "data" / case
        case_dir = REPO / "reduced_networks" / case
        if not case_dir.is_dir():
            print(f"{case}: no reductions yet, skipping")
            continue

        for level in sorted(p.name for p in case_dir.iterdir() if p.is_dir()):
            reduced = case_dir / level
            out = reduced / "dss"
            # One conversion per process. OpenDSS is a global singleton behind
            # opendssdirect and `ClearAll` does not fully reset it -- line codes
            # and wire data from a previous feeder survive and collide, which
            # surfaces as "Y matrix build aborted" on the *second* case in a run
            # even though it converts cleanly on its own.
            proc = subprocess.run(
                [sys.executable, str(CONVERTER / "convert_reduced_dataset.py"),
                 str(master), str(full), str(reduced), str(out),
                 "--sbase-mva", str(sbase)],
                capture_output=True, text=True,
            )
            if proc.returncode != 0:
                detail = (proc.stderr.strip().splitlines() or ["unknown error"])[-1]
                failures.append((case, level, detail[:110]))
                print(f"  {case:<12} {level:<7} FAILED  {detail[:80]}")
                continue
            report = json.loads(proc.stdout)
            rows.append((case, level, report))
            print(
                f"  {case:<12} {level:<7} {report['full_bus_count']:>5} -> "
                f"{report['reduced_bus_count']:<5} "
                f"max |V| err {report['max_voltage_magnitude_error_pu']:.3e} pu  "
                f"agg err {report['aggregation_max_abs_power_error_pu']:.1e}"
            )

    summary = REPO / "reduced_networks" / "dss_build.json"
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
