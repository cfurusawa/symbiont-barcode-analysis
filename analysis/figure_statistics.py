#!/usr/bin/env python3
"""Statistics reported in Fig. 2c, f, i, Fig. 4c and Extended Data Figs. 2 and 6.

Usage:
    python analysis/figure_statistics.py --barcodes DIR --master FILE
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import kruskal, wilcoxon

sys.path.insert(0, ".")
from barcodeflow.core import COLUMNS, REGIONS, Dataset           # noqa: E402
from barcodeflow.corrections import EXCLUDED_FROM_SIMILARITY     # noqa: E402
from barcodeflow.measures import cosine_similarity, richness_ratio  # noqa: E402

STRAINS = ["Anc", "CmL", "GmL", "SymD"]
TRANSITIONS = [("S1", "S2"), ("S2", "S3"), ("S3", "S4"), ("S4", "SE"), ("SE", "egg")]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--barcodes", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--out", default="results")
    args = ap.parse_args()

    data = Dataset(args.barcodes, args.master)
    counts = data.count_table()
    counts.to_csv(f"{args.out}/barcode_counts.csv", index=False)

    # ------------------------------------------- Extended Data Fig. 2
    print("Extended Data Fig. 2  barcodes acquired by newborn nymphs")
    groups = [counts.loc[counts.strain == s, "N"].values for s in STRAINS]
    print(f"  Kruskal-Wallis  P = {kruskal(*groups).pvalue:.3g}")
    for s in STRAINS:
        v = counts.loc[counts.strain == s, "N"]
        print(f"    {s:5s} n={len(v)}  median {v.median():,.0f}  "
              f"range {v.min():,.0f}-{v.max():,.0f}")

    # ------------------------------------ Fig. 2c, f, i and Fig. 4c
    print("\nFig. 2c, f, i and Fig. 4c  decline along the symbiotic organ")
    print("  linear mixed model: count ~ region position, individual as random intercept")
    for s in STRAINS:
        long = counts[counts.strain == s].melt(
            id_vars=["individual"], value_vars=REGIONS,
            var_name="region", value_name="count")
        long["position"] = long["region"].map({c: k + 1 for k, c in enumerate(REGIONS)})
        model = smf.mixedlm("count ~ position", long, groups=long["individual"])
        fit = model.fit(reml=False)
        print(f"    {s:5s} slope {fit.params['position']:+8.2f}  "
              f"P = {fit.pvalues['position']:.3g}")

    # ---------------------------------------------- Extended Data Fig. 6
    keep = [i for i in data.individuals
            if data.series[i]["strain"] not in EXCLUDED_FROM_SIMILARITY]
    print(f"\nExtended Data Fig. 6  (Anc excluded; n = {len(keep)} insects)")

    rows = []
    for individual in keep:
        for upstream, downstream in TRANSITIONS:
            value = richness_ratio(data, individual, upstream, downstream)
            if not np.isnan(value):
                rows.append({"individual": individual,
                             "transition": f"{upstream}->{downstream}",
                             "ratio": value})
    ratios = pd.DataFrame(rows)
    ratios.to_csv(f"{args.out}/richness_ratios.csv", index=False)
    wide = ratios.pivot_table(index="individual", columns="transition",
                              values="ratio").dropna()
    within = wide[["S1->S2", "S2->S3", "S3->S4"]].mean(axis=1)
    print("  (a) ratio of barcode numbers between successive samples")
    for name in [f"{a}->{b}" for a, b in TRANSITIONS]:
        print(f"      {name:10s} median {ratios.loc[ratios.transition == name, 'ratio'].median():.3f}")
    for target in ["S4->SE", "SE->egg"]:
        p = wilcoxon(wide[target], within).pvalue
        print(f"      {target} vs mean of S1-S4 transitions: P = {p:.3g}  "
              f"(n = {len(wide)})")

    # Panel (b) needs at least five barcodes in the end crypts for the
    # similarity to be informative; panel (c) uses every insect with an egg mass.
    for panel, target, columns, label, minimum in [
        ("b", "SE", REGIONS[:4], "similarity to the end crypts", 5),
        ("c", "egg", REGIONS, "similarity to the egg mass", 1),
    ]:
        rows = []
        for individual in keep:
            reference = data.profile(individual, target)
            if reference is None or (reference > 1e-3).sum() < minimum:
                continue
            rows.append({"individual": individual,
                         **{c: cosine_similarity(reference, data.profile(individual, c))
                            for c in columns}})
        table = pd.DataFrame(rows)
        table.to_csv(f"{args.out}/cosine_{panel}.csv", index=False)
        last = columns[-1]
        others = columns[:-1]
        p = wilcoxon(table[last], table[others].mean(axis=1)).pvalue
        print(f"  ({panel}) {label}: {last} vs mean of {others[0]}-{others[-1]}, "
              f"P = {p:.3g}  (n = {len(table)})")


if __name__ == "__main__":
    main()
