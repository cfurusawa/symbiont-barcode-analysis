# Barcode analysis for *Experimental evolution of population bottleneck in symbiosis*

Code for the barcode-sequencing analyses in the paper. It covers the steps
performed by the Furusawa/Maeda group: sample assignment, barcode counting,
diversity and similarity measures, and the statistics reported in the figures.

Analyses performed by other groups (quantitative PCR, colony counts,
fluorescence *in situ* hybridization and the related figures) used the
commercial and open-source software listed in the Methods and are not included
here.

## What this reproduces

| Output | Figure or table in the paper |
|---|---|
| Barcode counts per sample | Fig. 2b, c, e, f, h, i; Fig. 4b, c; Extended Data Figs. 3–5, 11; Extended Data Table 4 |
| Kruskal–Wallis test on barcodes acquired by nymphs | Extended Data Fig. 2 |
| Linear mixed models along the symbiotic organ | Fig. 2c, f, i; Fig. 4c |
| Ratios of barcode numbers and cosine similarities | Extended Data Fig. 6 |
| Library diversity and expected barcodes per inoculum | Extended Data Table 1 |
| Specificity of dominant clones and the carry-over estimate | Methods |
| Choice of egg mass where several were available | Extended Data Table 4 |

## Requirements

Python 3.12 with

```
numpy==2.4.4
scipy==1.17.1
pandas==3.0.2
statsmodels==0.15.0
matplotlib==3.10.8
openpyxl==3.1.5
```

Install with `pip install -r requirements.txt`.

## Input

Two inputs are needed.

**Barcode tables.** One CSV per sequencing run, with the columns
*sequence*, *reads* and *cluster*, as produced by Bartender v1.1. Upstream of
this, reads were filtered with fastp v0.21.0 (`-q 30`), the 15-nt barcode was
located in Read 1 by the pattern `GCCGG[15]TATCT`, and barcodes were clustered
with Bartender (`-d 3 -z -1`). Raw reads are deposited in the Sequence Read
Archive under BioProject PRJNA1551028; see Extended Data Table 4 for the run
corresponding to each sample.

**Sample table.** The laboratory record that maps run numbers to insects and
tissues. Its contents are reproduced in Extended Data Table 4; the assignment
actually used, including every correction, is also written out in
`barcodeflow/corrections.py` so that the analysis can be followed without it.

## Running

```bash
pip install -r requirements.txt
mkdir -p results

python analysis/library_diversity.py   --barcodes path/to/barcodes
python analysis/figure_statistics.py   --barcodes path/to/barcodes --master path/to/sample_table.xlsx
python analysis/validate_assignment.py --barcodes path/to/barcodes --master path/to/sample_table.xlsx
```

Each script writes CSV files to `results/` and prints the values quoted in the
paper.

## Decisions recorded in the code

Three choices affect the counts and are set out in `barcodeflow/corrections.py`.

**Detection threshold.** A barcode counts as detected when its relative
frequency within a sample exceeds 0.1%. The newborn nymphs are the exception:
their symbiont population is too diverse for any barcode to reach that
frequency, so all detected barcodes are counted there.

**Samples interchanged during handling.** The end-crypt samples of SymD-1 and
SymD-5, and the egg masses of SymD-2 and SymD-5, were interchanged; this was
identified from their barcode composition and corrected.

**One excluded set of barcodes.** In the egg mass of Anc-3, 71 of the 72
barcodes detected were absent from every midgut region of the same insect and
are excluded from the count for that sample. No other sample is treated this
way.

`analysis/validate_assignment.py` recomputes the evidence for the second and
third of these from the data.

## Layout

```
barcodeflow/
    core.py          sample assignment, barcode profiles, counts
    corrections.py   every assignment decision, in one place
    measures.py      inverse Simpson, cosine similarity, expected barcodes
analysis/
    library_diversity.py     Extended Data Table 1
    figure_statistics.py     Fig. 2, Fig. 4, Extended Data Figs. 2 and 6
    validate_assignment.py   checks behind the Methods statements
```

## Citing

[Citation to be added on publication.]
