from __future__ import annotations

import argparse
import json
from pathlib import Path

from opendss_tools import ConversionError, convert


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert a full OpenDSS circuit and assignment matrix into a radial Kron-equivalent DSS circuit.")
    parser.add_argument("master", type=Path, help="Full circuit Master.dss")
    parser.add_argument("assignment", type=Path, help="Dense binary assignment CSV or .npy")
    parser.add_argument("bus_order", type=Path, help="CSV containing bus_id in A row/column order")
    parser.add_argument("output", type=Path, help="Output directory")
    parser.add_argument("--matrix-tolerance", type=float, default=1e-9)
    parser.add_argument("--max-voltage-error-pu", type=float, default=1e-6)
    args = parser.parse_args()
    try:
        report = convert(args.master, args.assignment, args.bus_order, args.output, args.matrix_tolerance, args.max_voltage_error_pu)
    except ConversionError as exc:
        parser.error(str(exc))
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
