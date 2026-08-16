from __future__ import annotations

import argparse
import json
from pathlib import Path

from opendss_tools import DatasetError, export_dataset


def main() -> int:
    p=argparse.ArgumentParser(description='Export an ordered phase-domain per-unit dataset from OpenDSS.')
    p.add_argument('master',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--sbase-mva',type=float,required=True,help='Common phase-domain power base in MVA')
    p.add_argument('--scenario',default='h001')
    p.add_argument('--max-iterations',type=int,default=100)
    p.add_argument('--max-control-iterations',type=int,default=100)
    p.add_argument('--serialized-tolerance',type=float,default=1e-12)
    p.add_argument('--independent-tolerance',type=float,default=1e-5)
    args=p.parse_args()
    try:
        report=export_dataset(args.master,args.output,args.sbase_mva,args.scenario,args.max_iterations,args.max_control_iterations,args.serialized_tolerance,args.independent_tolerance)
    except DatasetError as exc:
        p.error(str(exc))
    print(json.dumps(report,indent=2)); return 0


if __name__=='__main__':
    raise SystemExit(main())
