"""Sample assignment and the corrections applied to it.

Every decision made about which sequencing run belongs to which sample is
recorded here, so that the assignment used in the paper can be inspected in one
place. It is reproduced in Extended Data Table 4.
"""

# Column of the paper -> stage code in the master table
STAGE_CODES = {
    "N": "1",     # newborn nymph
    "S1": "3a",   # adult midgut M4, region 1
    "S2": "3b",
    "S3": "3c",
    "S4": "3d",
    "SE": "3e",   # end crypts
}

# ---------------------------------------------------------------------------
# Anc-9: the master table lists two runs for S1 and S2. The runs from the later
# experiment are used so that all five midgut regions of this insect come from
# the same sequencing run.
# ---------------------------------------------------------------------------
REGION_OVERRIDES = {
    ("Anc-9", "S1"): 320,
    ("Anc-9", "S2"): 321,
}

# ---------------------------------------------------------------------------
# Samples interchanged during handling, identified from their barcode
# composition (see Methods and validate_assignment.py).
# ---------------------------------------------------------------------------
SE_SWAPS = [("SymD-1", "SymD-5")]
EGG_SWAPS = [("SymD-2", "SymD-5")]

# ---------------------------------------------------------------------------
# Egg mass used for each insect. Where several egg masses were available, the
# one in which the clones established in the mother accounted for the largest
# fraction of reads was used; see select_egg_masses.py, which reproduces this
# choice from the data. ``None`` marks an insect with no usable egg mass.
# ---------------------------------------------------------------------------
EGG_ASSIGNMENT = {
    "Anc-1": 13,   "Anc-2": 293,  "Anc-3": 294,  "Anc-4": 14,   "Anc-5": None,
    "Anc-6": 296,  "Anc-7": 330,  "Anc-8": 327,  "Anc-9": 148,
    "CmL-1": 11,   "CmL-2": 7,    "CmL-3": 316,  "CmL-4": 9,    "CmL-5": 299,
    "CmL-6": 300,  "CmL-7": 301,  "CmL-8": 332,  "CmL-9": 335,
    "GmL-1": 343,  "GmL-2": 317,  "GmL-3": 305,  "GmL-4": 12,   "GmL-5": 337,
    "GmL-6": 307,  "GmL-7": 308,  "GmL-8": 309,  "GmL-9": 310,
    "SymD-1": 311, "SymD-2": 341, "SymD-3": 338, "SymD-4": 314, "SymD-5": 312,
}

# ---------------------------------------------------------------------------
# The single exception applied when counting barcodes. In the egg mass of Anc-3,
# 71 of the 72 barcodes detected were absent from every midgut region of the
# same insect; they are excluded from the count for this sample only.
# ---------------------------------------------------------------------------
ANC3_EXCEPTION = {
    "sample": ("Anc-3", "egg"),
    "rule": "exclude barcodes absent above threshold from all midgut regions "
            "of the same insect",
}

# Strains excluded from the analyses of Extended Data Fig. 6, because only one
# to six barcodes establish per insect, so similarity measures saturate and the
# ratio between adjacent regions turns on a single barcode.
EXCLUDED_FROM_SIMILARITY = ["Anc"]
