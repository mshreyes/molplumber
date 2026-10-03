"""Command-line interface for molplumber."""
from __future__ import annotations

import argparse
import json
import sys

import pandas as pd

from .core import DEFAULT_THRESHOLDS, assess_leakage, summarize


def _read(path: str) -> pd.DataFrame:
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    return pd.read_csv(path, sep=sep)


def _get_col(df: pd.DataFrame, col: str, path: str) -> pd.Series:
    if col not in df.columns:
        sys.exit(f"Column '{col}' not found in {path}. Available: {list(df.columns)}")
    return df[col]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="molplumber",
        description="Report how similar each test molecule is to the training set.",
    )
    p.add_argument("train", help="Training set (.csv, or .tsv/.txt for tab-separated)")
    p.add_argument("test", help="Test set (same formats)")
    p.add_argument("-o", "--out", default="leakage_report.csv", help="Per-molecule report (CSV)")
    p.add_argument("--smiles-col", default="smiles", help="SMILES column in both files")
    p.add_argument("--train-smiles-col", help="Override SMILES column for the training file")
    p.add_argument("--test-smiles-col", help="Override SMILES column for the test file")
    p.add_argument("--radius", type=int, default=2, help="Morgan radius (default: 2)")
    p.add_argument("--n-bits", type=int, default=2048, help="Fingerprint length (default: 2048)")
    p.add_argument("--generic-scaffold", action="store_true",
                   help="Use generic Murcko scaffolds instead of standard ones")
    p.add_argument("--thresholds", type=float, nargs="+", default=list(DEFAULT_THRESHOLDS),
                   help="Similarity cutoffs for the summary (default: 0.4 0.6 0.8 1.0)")
    p.add_argument("--summary-json", help="Also write the summary to this JSON file")
    return p


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)

    train_df, test_df = _read(args.train), _read(args.test)
    train_col = args.train_smiles_col or args.smiles_col
    test_col = args.test_smiles_col or args.smiles_col

    report = assess_leakage(
        _get_col(train_df, train_col, args.train),
        _get_col(test_df, test_col, args.test),
        radius=args.radius,
        n_bits=args.n_bits,
        generic_scaffold=args.generic_scaffold,
    )
    report.to_csv(args.out, index=False)

    summary = summarize(report, args.thresholds)
    for key, val in summary.items():
        print(f"{key:35s} {val:.3f}" if isinstance(val, float) else f"{key:35s} {val}")
    print(f"\nPer-molecule report written to {args.out}")

    if args.summary_json:
        with open(args.summary_json, "w") as fh:
            json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
