from __future__ import annotations

import argparse
import json
from pathlib import Path

from opendss_tools import ConversionError, convert_reduced_dataset


def main() -> int:
    parser=argparse.ArgumentParser(description='Convert a radial reduced CSV dataset into a solved OpenDSS circuit.')
    parser.add_argument('full_master',type=Path)
    parser.add_argument('full_dataset',type=Path)
    parser.add_argument('reduced_dataset',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--sbase-mva',type=float,required=True)
    parser.add_argument('--matrix-tolerance',type=float,default=1e-9)
    parser.add_argument('--aggregation-tolerance',type=float,default=1e-12)
    parser.add_argument(
        '--kron-relative-tolerance',type=float,default=1e-12,
        help='Relative threshold for rebuilding an omitted reduced Y-bus while retaining its diagonal and rank.',
    )
    parser.add_argument(
        '--load-representation',choices=('original-models','mapped-elements','nodal-wye'),
        default='original-models',
        help='Relocate original models with per-unit voltage-base rebasing (default), freeze solved loads as PQ, or use CSV nodal-wye powers.',
    )
    args=parser.parse_args()
    try:
        report=convert_reduced_dataset(
            args.full_master,args.full_dataset,args.reduced_dataset,args.output,args.sbase_mva,
            args.matrix_tolerance,args.aggregation_tolerance,args.load_representation,
            args.kron_relative_tolerance,
        )
    except ConversionError as exc:
        parser.error(str(exc))
    print(json.dumps(report,indent=2));return 0


if __name__=='__main__':raise SystemExit(main())
