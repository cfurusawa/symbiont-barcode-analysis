#!/usr/bin/env python3
"""Checks on the sample assignment, and the carry-over estimate in the Methods.

Reproduces three statements:

  1. the dominant clones are specific to a single insect;
  2. the level of carry-over between experimental samples;
  3. the choice of egg mass for insects with more than one.

Usage:
    python analysis/validate_assignment.py --barcodes DIR --master FILE
"""
from __future__ import annotations

import argparse
import collections
import itertools
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from barcodeflow.core import DETECTION_THRESHOLD, REGIONS, Dataset  # noqa: E402

DOMINANT = 0.10        # a clone "dominates" an insect above this frequency


def dominant_clones(data: Dataset) -> dict[str, set]:
    """Barcodes above 10% in any midgut region, per insect."""
    out = {}
    for individual in data.individuals:
        peak: dict[str, float] = {}
        for column in REGIONS:
            profile = data.profile(individual, column)
            for barcode, value in profile.items():
                if value > peak.get(barcode, 0.0):
                    peak[barcode] = value
        out[individual] = {b for b, v in peak.items() if v > DOMINANT}
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--barcodes", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--out", default="results")
    args = ap.parse_args()

    data = Dataset(args.barcodes, args.master)
    dominant = dominant_clones(data)
    total = sum(len(v) for v in dominant.values())

    # ---------------------------------------------- 1. specificity
    shared = collections.Counter()
    for a, b in itertools.combinations(data.individuals, 2):
        for barcode in dominant[a] & dominant[b]:
            shared[barcode] += 1
    print(f"Dominant clones (>{DOMINANT:.0%} in a symbiotic organ): {total}")
    print(f"  shared by more than one insect: {len(shared)}")

    # ---------------------------------------------- 2. carry-over
    samples = {}
    for individual in data.individuals:
        for column in REGIONS + ["egg"]:
            if data.run_of(individual, column) is not None:
                samples[(individual, column)] = data.profile(individual, column)

    rows = []
    for (individual, column), profile in samples.items():
        strain = data.series[individual]["strain"]
        for other in data.individuals:
            if other == individual:
                continue
            if data.series[other]["strain"] != strain:
                continue        # different libraries, barcodes cannot be shared
            for barcode in dominant[other]:
                rows.append({"sample": f"{individual}_{column}",
                             "individual": individual,
                             "donor": other,
                             "frequency": float(profile.get(barcode, 0.0))})
    carry = pd.DataFrame(rows)
    carry.to_csv(f"{args.out}/carryover.csv", index=False)

    detected = carry[carry.frequency > 0]
    above = carry[carry.frequency > DETECTION_THRESHOLD]
    print(f"\nCarry-over, over {len(carry):,} combinations of a sample and a "
          f"barcode dominant in another insect of the same strain:")
    print(f"  detected at all        {len(detected):,}  "
          f"(median frequency {detected.frequency.median():.2e})")
    print(f"  above {DETECTION_THRESHOLD:.1%}             {len(above):,}  "
          f"({len(above)/len(carry):.2%})")
    by_insect = collections.Counter(above["individual"])
    print("  insects contributing:", dict(by_insect.most_common()))

    # ---------------------------------------------- 3. egg mass choice
    print("\nEgg mass selection (fraction of reads from clones established in the mother):")
    owner = {}
    for individual in data.individuals:
        for column in REGIONS:
            for barcode, value in data.profile(individual, column).items():
                best = owner.get(barcode)
                if best is None or value > best[1]:
                    owner[barcode] = (individual, value)
    assigned = {b: who for b, (who, v) in owner.items() if v > 1e-4}

    rows = []
    for individual in data.individuals:
        candidates = data.series[individual]["stages"].get("4a", [])
        for run in candidates:
            if data.runs[run]["barcode_file"] not in data._files:
                continue
            profile = data._raw_profile(run)
            own = float(sum(v for b, v in profile.items()
                            if assigned.get(b) == individual))
            rows.append({"individual": individual, "run": run,
                         "maternal_fraction": round(own, 4),
                         "used": run == data.egg[individual]})
    eggs = pd.DataFrame(rows)
    eggs.to_csv(f"{args.out}/egg_selection.csv", index=False)
    multiple = eggs.groupby("individual").filter(lambda g: len(g) > 1)
    print(multiple.to_string(index=False))


if __name__ == "__main__":
    main()
