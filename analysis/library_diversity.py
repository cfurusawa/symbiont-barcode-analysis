#!/usr/bin/env python3
"""Barcode diversity of the four chromosomal libraries (Extended Data Table 1).

Usage:
    python analysis/library_diversity.py --barcodes DIR
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, ".")
from barcodeflow.measures import expected_barcodes, inverse_simpson  # noqa: E402

# Sequencing run of each chromosomal library, named by its barcode table.
LIBRARIES = {
    "Anc": "BW25113-dmutS_fastp_barcode_barcode.csv",
    "CmL": "CmL0521A_fastp_barcode_barcode.csv",
    "GmL": "GmL1220A_fastp_barcode_barcode.csv",
    "SymD": "11_S11_L001_R1_001.fastp_extracted_barcode_barcode.csv",
}
INOCULUM = 10_000          # cells ingested by a newborn nymph, order of magnitude


def load(path: str) -> tuple[pd.Series, int]:
    table = pd.read_csv(path)
    table.columns = ["sequence", "reads", "cluster"]
    table = table.sort_values("reads", ascending=False)
    representative = table.groupby("cluster", sort=False).agg(
        sequence=("sequence", "first"), reads=("reads", "sum"))
    counts = representative.groupby("sequence")["reads"].sum()
    return counts / counts.sum(), int(counts.sum())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--barcodes", required=True)
    ap.add_argument("--out", default="results")
    args = ap.parse_args()

    files = {os.path.basename(p): p for p in
             glob.glob(os.path.join(args.barcodes, "**", "*.csv"), recursive=True)
             if "old_ver" not in p}

    rows = []
    for name, filename in LIBRARIES.items():
        profile, reads = load(files[filename])
        unique = len(profile)
        effective = inverse_simpson(profile, threshold=0.0)
        expected, sd = expected_barcodes(profile, INOCULUM)
        rows.append({
            "library": name,
            "reads": reads,
            "unique_barcodes": unique,
            "inverse_simpson": round(effective),
            "inverse_simpson_over_unique": round(effective / unique, 3),
            "max_barcode_frequency": float(profile.max()),
            "expected_barcodes_per_inoculum": round(expected),
            "sd_of_expected": round(sd),
        })

    table = pd.DataFrame(rows)
    table.to_csv(f"{args.out}/library_diversity.csv", index=False)
    pd.set_option("display.width", 160)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
