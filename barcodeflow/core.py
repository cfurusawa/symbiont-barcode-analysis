"""Sample assignment and barcode profiles.

Reads the per-sample barcode tables produced by Bartender, applies the sample
assignment described in the Methods, and exposes the relative-frequency profile
of every sample used in the analyses.

The assignment is fixed in ``corrections.py`` and is reproduced in Extended Data
Table 4 of the paper.
"""
from __future__ import annotations

import glob
import os
from functools import lru_cache

import numpy as np
import pandas as pd

from .corrections import (ANC3_EXCEPTION, EGG_ASSIGNMENT, REGION_OVERRIDES,
                          SE_SWAPS, STAGE_CODES)

DETECTION_THRESHOLD = 1e-3      # relative frequency; see Methods
REGIONS = ["S1", "S2", "S3", "S4", "SE"]
COLUMNS = ["N"] + REGIONS + ["egg"]


class Dataset:
    """Barcode profiles for every sample, with the published assignment applied."""

    def __init__(self, barcode_dir: str, master_table: str):
        self.barcode_dir = barcode_dir
        self._files = {
            os.path.basename(p): p
            for p in glob.glob(os.path.join(barcode_dir, "**", "*.csv"), recursive=True)
            if "old_ver" not in p
        }
        self._read_master(master_table)
        self._build_assignment()

    # ---------------------------------------------------------------- master
    def _read_master(self, path: str) -> None:
        import openpyxl

        wb = openpyxl.load_workbook(path, data_only=True)

        self.runs = {}
        for row in wb.worksheets[1].iter_rows(min_row=2, values_only=True):
            if isinstance(row[0], (int, float)):
                self.runs[int(row[0])] = {
                    "description": str(row[1]) if row[1] else "",
                    "strain": row[2],
                    "line": row[3],
                    "date": str(row[5])[:10] if row[5] else "",
                    "barcode_file": row[11],
                }

        stages = ["1", "3a", "3b", "3c", "3d", "3e", "4a", "4b", "5"]
        self.series, strain = {}, None
        for row in wb.worksheets[0].iter_rows(min_row=3, values_only=True):
            if row[1]:
                strain = row[1]
            if not row[2]:
                continue
            entry = {}
            for k, stage in enumerate(stages):
                value = row[4 + k]
                if value is None:
                    continue
                nums = [
                    int(x.strip())
                    for x in str(value).replace("，", ",").split(",")
                    if x.strip().isdigit()
                ]
                if nums:
                    entry[stage] = nums
            self.series[str(row[2]).strip()] = {
                "strain": strain,
                "line": str(row[3]).strip(),
                "stages": entry,
            }
        self.individuals = list(self.series)

    # ------------------------------------------------------------ assignment
    def _build_assignment(self) -> None:
        use = {
            i: {col: self.series[i]["stages"][code][0] for col, code in STAGE_CODES.items()}
            for i in self.individuals
        }
        for (individual, column), run in REGION_OVERRIDES.items():
            use[individual][column] = run
        for a, b in SE_SWAPS:
            use[a]["SE"], use[b]["SE"] = use[b]["SE"], use[a]["SE"]
        self.use = use
        self.egg = dict(EGG_ASSIGNMENT)

    def run_of(self, individual: str, column: str) -> int | None:
        if column == "egg":
            return self.egg[individual]
        return self.use[individual][column]

    # -------------------------------------------------------------- profiles
    @lru_cache(maxsize=None)
    def _raw_profile(self, run: int) -> pd.Series:
        """Relative frequency of each barcode in one sequencing run.

        Clusters are collapsed onto their most abundant member sequence, which is
        what makes clusters comparable between samples (cluster identifiers are
        assigned independently in each sample).
        """
        path = self._files[self.runs[run]["barcode_file"]]
        table = pd.read_csv(path)
        table.columns = ["sequence", "reads", "cluster"]
        table = table.sort_values("reads", ascending=False)
        representative = table.groupby("cluster", sort=False).agg(
            sequence=("sequence", "first"), reads=("reads", "sum")
        )
        counts = representative.groupby("sequence")["reads"].sum()
        return counts / counts.sum()

    def profile(self, individual: str, column: str) -> pd.Series | None:
        """Profile of one sample, with the Anc-3 egg exception applied."""
        run = self.run_of(individual, column)
        if run is None:
            return None
        profile = self._raw_profile(run)
        if (individual, column) == ANC3_EXCEPTION["sample"]:
            keep = self._anc3_maternal_barcodes()
            profile = profile.copy()
            profile[~profile.index.isin(keep)] = 0.0
        return profile

    @lru_cache(maxsize=1)
    def _anc3_maternal_barcodes(self) -> frozenset:
        """Barcodes above threshold in any midgut region of the Anc-3 mother."""
        individual = ANC3_EXCEPTION["sample"][0]
        keep: set[str] = set()
        for column in REGIONS:
            p = self._raw_profile(self.run_of(individual, column))
            keep |= set(p.index[p > DETECTION_THRESHOLD])
        return frozenset(keep)

    # ----------------------------------------------------------------- counts
    def barcode_count(self, individual: str, column: str) -> float:
        """Number of barcodes scored as detected (see Methods).

        All detected barcodes are counted for the newborn nymph, because the
        population there is too diverse for any barcode to reach 0.1%.
        """
        profile = self.profile(individual, column)
        if profile is None:
            return np.nan
        if column == "N":
            return int((profile > 0).sum())
        return int((profile > DETECTION_THRESHOLD).sum())

    def count_table(self) -> pd.DataFrame:
        rows = []
        for individual in self.individuals:
            rows.append(
                {
                    "strain": self.series[individual]["strain"],
                    "individual": individual,
                    "line": self.series[individual]["line"],
                    **{c: self.barcode_count(individual, c) for c in COLUMNS},
                }
            )
        return pd.DataFrame(rows)
