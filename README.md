# How much can five donors tell us? A donor-level reanalysis of PBMC single-cell data after malaria treatment

Exploratory reanalysis, in Python, of the public single-cell RNA-seq dataset from
Dooley et al., *Nature Communications* 2023 (14:7387), GEO accession
[GSE217930](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217930).
The dataset follows PBMCs from donors with uncomplicated *P. falciparum* malaria before and after treatment.
This project asks what can be concluded when the **donor**, not the cell, is treated as the unit of replication.

**Status: exploratory.** There are 5 donors with paired day 0 and day 28 samples (4 for NK cells), so the
smallest p-value an exact test can return is 0.0625 (0.125 with 4 donors). Results are reported as direction and
consistency across donors, not as significance claims.

## Main takeaways

1. On the same cells, a cell-level test calls 20 to 200 times more differentially expressed genes than a donor-level test, and most of the extra genes do not move the same way in every donor.
2. Donor identity explains more expression variance than time point, and there is no evidence that between-donor spread differs between day 0 and day 28.
3. CD8 T cells rose in all 5 donors, and cycling and activated T cells fell in all informative donors. None of this is significant, because with 4 to 5 donors an exact test cannot get below p = 0.0625 to 0.125.
4. Inferred transcription-factor activity shows a signal in NK cells only (lower NF-kB, STAT and AP-1 activity at day 28), and that result rests on 4 donors.
5. Monocytes, cDC and pDC are too sparse per donor to test.

## Questions

1. Which immune cell types change between day 0 (acute infection) and day 28 after treatment?
2. How many differentially expressed genes survive when donors, not cells, are the replicates, compared with a standard cell-level test on the same cells?
3. How much of the expression variance is donor-to-donor versus day, and does between-donor spread differ between day 0 and day 28?
4. Does transcription-factor activity inferred from donor-level statistics point to the NF-kB, IRF, STAT and AP-1 families highlighted in the original paper?

## Data

| | |
|---|---|
| Source | GEO GSE217930, Cell Ranger outputs (10x 3' v3.1), 20 libraries |
| Design | 6 infected donors (ages 6 to 24, natural infection in Sabah, Malaysia) sampled at day 0, 7 and 28 after treatment, plus 2 endemic controls |
| Removed | Child1 day 0 (failed library, also removed in the paper), leaving day 0 n = 5, day 7 n = 6, day 28 n = 6, controls n = 2 |
| Batches | Libraries 1 to 8 and 9 to 20 come from two sequencing batches (different feature lists). Every infected donor has all of its time points inside one batch, so batch is nested in donor and does not confound day. No batch integration was applied. |

Raw data are not included in this repository. Download `GSE217930_RAW.tar` from GEO and extract it into `data/`.

## Pipeline

| Step | What was done |
|---|---|
| Load and merge | 20 libraries read with `scanpy.read_10x_mtx`, genes matched by Ensembl ID, inner join: 122,518 cells x 19,810 genes |
| QC | at least 250 genes, at least 500 UMIs, under 20% mitochondrial reads, complexity (log10 genes / log10 UMIs) above 0.8, haemoglobin genes removed, genes in fewer than 10 cells removed: 106,385 cells x 15,829 genes; 105,952 after dropping Child1 day 0 |
| Clustering | normalise to 10^4, log1p, 3,000 highly variable genes (selected per donor), 30 PCs, Leiden resolution 0.6 (18 clusters). The clustering matrix was normalised and log-transformed twice by accident (see Limitations); all downstream tests use raw counts or a correctly rebuilt log-normalised layer |
| Annotation | own annotation from canonical markers and co-expression checks (not the authors' labels). A cDC2 population was recovered from a monocyte sub-cluster. 2,126 cells co-expressing B+T or monocyte+lymphocyte programmes were flagged as doublet-like and excluded from all downstream tests |
| Cycling cells | scored with a 15-gene proliferation signature; cells with score above 0.5 were split into their own group (364 cells). The threshold was chosen by inspecting the score histogram: 457 / 375 / 314 cells lie above 0.4 / 0.5 / 0.6 |
| Pseudobulk | raw counts summed per donor x day x cell type; samples with fewer than 20 cells dropped |
| Donor-level DE | PyDESeq2, design `~ donor + day`, contrast day 28 vs day 0, only donors with both time points, genes with at least 10 total counts |
| Cell-level DE | Wilcoxon rank-sum on log-normalised cells, same cells and donors as the donor-level test |
| Composition | per-sample cell-type proportions, log-proportion change per donor, exact sign-flip test, donors with at least 5 cells of that type in both samples |
| Variance | per-gene variance partitioning (donor / day / residual) on pseudobulk log-CPM, compared with a label-flip null; donor-to-donor variance ratio day 0 / day 28 with equal cell numbers per sample |
| TF activity | decoupler univariate linear model on the donor-level Wald statistics with the CollecTRI network |

Final cell groups (non-doublet): CD4 T 46,070; CD8 T 27,243; NK 14,529; monocyte 6,871; B 6,365; T activated 794;
plasma 604; pDC 485; cDC 454; cycling 364; platelet 47.

## Results

![Summary panel](results/figures/summary_panel.png)

### 1. Cell-level testing calls far more genes than donor-level testing

Genes with BH-adjusted p < 0.05, day 28 vs day 0, same cells and same donors:

| Cell type | Paired donors | Cells day 0 / day 28 | Cell-level (naive) | Donor-level (pseudobulk) | Donor-level genes also found by cell-level test |
|---|---|---|---|---|---|
| B | 5 | 582 / 1,437 | 545 | 3 | 1 of 3 |
| CD4 T | 5 | 9,601 / 12,991 | 2,699 | 135 | 126 of 135 |
| CD8 T | 5 | 3,734 / 8,620 | 2,115 | 60 | 53 of 60 |
| NK | 4 | 760 / 2,701 | 1,461 | 7 | 6 of 7 |

Before the cycling cells were split off, donor-level counts were B 4, CD4 T 132, CD8 T 56, NK 80. NK lost one donor in
the split version (4 instead of 5 paired donors), so the NK donor-level result is the least stable. PyDESeq2 also warned
about low residual degrees of freedom for NK.

Donor-level testing has little power with 4 to 5 donors, so fewer calls does not mean the extra cell-level genes are
all false. Panel B checks this directly: for each gene, did the log-CPM change go the same way in every donor?

- Genes found **only** by the cell-level test: same direction in all donors for 23.8% (B, 542 genes), 17.3% (CD4 T, 2,573), 16.3% (CD8 T, 2,062) and 34.9% (NK, 1,454). Chance level is 3.1% (6.2% for NK with 4 donors). Above chance, but far from consistent.
- Genes found by **both** tests: 82.5% (CD4 T, 126 genes), 77.4% (CD8 T, 53) and 83.3% (NK, 6); B has only one such gene. Every one of these genes has the same direction in all but at most one donor.
- Nearly every donor-level hit is also a cell-level hit (table above), so the donor-level test mostly acts as a much stricter filter on the same signal.

### 2. Composition changes between day 0 and day 28

![Composition per donor](results/figures/composition_paired.png)

| Cell type | Donors tested | Donors up / down | Mean % day 0 | Mean % day 28 | Exact p (min. possible) |
|---|---|---|---|---|---|
| CD8 T | 5 | 5 / 0 | 21.8 | 31.8 | 0.062 (0.062) |
| Cycling | 4 | 0 / 4 | 0.85 | 0.25 | 0.125 (0.125) |
| T activated | 4 | 0 / 4 | 2.24 | 0.38 | 0.125 (0.125) |
| NK | 5 | 4 / 1 | 4.7 | 11.0 | 0.188 (0.062) |
| CD4 T | 5 | 2 / 3 | 61.5 | 48.3 | 0.250 |
| B | 5 | 2 / 3 | 3.5 | 5.4 | 0.500 |
| Plasma | 4 | 1 / 3 | 0.88 | 0.82 | 0.625 |

No change is significant after BH correction (best adjusted p = 0.29), which was not attainable with this many donors.
CD8 T cells rose in every donor, and cycling and activated T cells fell in every informative donor. Proportions are
compositional, so the CD4 T decrease partly mirrors the CD8 T and NK increase. The B-cell increase comes from two donors
that start near 1.3%.

Monocytes, cDC, pDC and platelets were too sparse to test. At day 0, 707 of the 750 monocytes come from one donor
(Adult2), so a donor-level monocyte comparison was not possible.

### 3. Donor identity explains more variance than time point

Median share of pseudobulk log-CPM variance per gene (well-expressed genes only), with the label-flip null for the day share:

| Cell type | Donors | Donor | Day | Residual | Day share: null mean (max) | p for day (min. possible) |
|---|---|---|---|---|---|---|
| B | 5 | 0.515 | 0.074 | 0.352 | 0.055 (0.080) | 0.188 (0.062) |
| CD4 T | 5 | 0.535 | 0.075 | 0.318 | 0.048 (0.079) | 0.125 (0.062) |
| CD8 T | 5 | 0.537 | 0.067 | 0.324 | 0.050 (0.083) | 0.312 (0.062) |
| NK | 4 | 0.467 | 0.105 | 0.355 | 0.083 (0.095) | 0.125 (0.125) |

Donor explains more variance than day in 84% (NK) to 91% (B) of genes. On the F-ratio scale (median over genes), donor is
1.29 to 1.66 against a null median of 1.0, and day is 0.87 to 0.97 against a null median of 0.55 to 0.59. The day share is
above the label-flip null mean in all four cell types, and in NK it is above the largest value flipping produced, which with
4 donors is the smallest p-value possible (0.125). The other three cell types fall inside the null range. A small day effect
may be present, but it cannot be separated from chance with this many donors.

### 4. Between-donor spread, day 0 versus day 28

Donor-to-donor variance of the 2,000 most-counted genes (log-CPM), day 0 divided by day 28. Every sample was subsampled to the same number of cells (30 resamplings) so that differences in cell number cannot create a difference in noise, and p-values come from flipping day labels within donors.

| Cell type | Donors | Cells per sample | Variance ratio day 0 / day 28 | p (flip test) |
|---|---|---|---|---|
| B | 5 | 26 | 1.10 | 0.710 |
| CD4 T | 5 | 1,151 | 1.50 | 0.442 |
| CD8 T | 5 | 325 | 1.23 | 0.608 |
| NK | 4 | 26 | 1.21 | 0.729 |

All four ratios are slightly above 1 (more spread at day 0), but none differs from what label flipping produces. B and NK were
subsampled to 26 cells per sample because one day 0 sample has only 26 cells of those types, so those two have little power.
This is an absence of evidence, not evidence that variability is unchanged.

### 5. Transcription-factor activity

Univariate linear model (decoupler) on the donor-level Wald statistics, day 28 vs day 0, with the CollecTRI network (1,185 TFs). Negative scores mean lower activity at day 28.

Smallest BH-adjusted p-value over all TFs: B 0.997, CD4 T 0.995, CD8 T 0.175, NK 0.001. Nothing is significant in B, CD4 T or CD8 T.
In NK, TFs from three of the families named in the paper score negative (higher activity at day 0); unadjusted p-values:

| TF | NK score | NK p |
|---|---|---|
| STAT3 | -4.85 | < 0.001 |
| JUN | -3.70 | 0.013 |
| RELA | -3.32 | 0.037 |
| STAT5A | -3.29 | 0.038 |

NFKB1, STAT1, STAT5B and FOSB reach p = 0.06 to 0.08 in the same direction. The IRF family shows no signal in any cell type (all p > 0.5).

Treat this as a hypothesis only. The p-values come from a regression across genes, not across donors, so they do not carry the
uncertainty from having 4 donors; TFs in the network share target genes, so the scores are not independent; and NK is also
the least stable donor-level result above. Expression-based TF activity is an indirect proxy for regulation and says nothing about chromatin.

## Limitations

- Four to five paired donors; all inference is exploratory and exact p-values have a floor of 0.0625 to 0.125.
- Natural infection with an unknown time since infection; day 0 is acute illness, and day 28 is not a healthy baseline (the paper reports differences from endemic controls). The 2 endemic controls are from one batch only and were not used in the day comparisons.
- Annotation and the cycling and doublet cut-offs are my own choices from marker inspection, not the authors' labels; the cycling threshold was set by eye.
- The two sequencing batches were not integrated; within-donor contrasts are unaffected, but between-donor comparisons carry batch.
- TF-activity p-values are gene-level (regression across genes), not donor-level, and TF target sets overlap.
- The Leiden clusters used for annotation were computed on a matrix that had been normalised and log-transformed twice (a cell was re-run after installing `igraph`). Annotation was checked against canonical markers, and pseudobulk counts, DE, composition and variance analyses use raw counts or a rebuilt, correctly log-normalised layer. Notebook 02 reproduces the original clustering (`DOUBLE_TRANSFORM = True`) and ends with an optional check that re-clusters on single normalisation and reports the agreement.
- No genotypes, so no genetic-effect analysis.
- Monocyte, cDC and pDC responses cannot be assessed at donor level in these data.

## Reproduce

```bash
python -m venv .venv
.venv\Scripts\activate        # Windows PowerShell
pip install -r requirements.txt
jupyter lab
```

Run the notebooks in order (`01`, `02`, `03`). Set `data_dir` and `out_dir` in the first cell of each notebook; `out_dir`
needs about 5 GB free, and the `.h5ad` files written there are large and are not tracked in git. Notebook 02 stops with an
error if the clustering differs from the reported run, because the annotation maps cluster numbers to cell types.

Tested with Python 3.13, scanpy 1.12.4, anndata 0.13.4, pydeseq2 0.5.4, decoupler 2.2.0.

## Repository layout

```
notebooks/
  01_load_qc_metadata.ipynb       load 20 libraries, merge, QC, GEO metadata
  02_clustering_annotation.ipynb  clustering, annotation, cDC2, cycling score
  03_donor_level_analysis.ipynb   pseudobulk DE, naive vs donor test, variance, TF activity, composition, figure
results/
  tables/     small CSV outputs (DE results, composition, variance, TF activity)
  figures/    summary panel and composition plot
requirements.txt
```

## Credit

Data and original analysis: Dooley et al., *Nat Commun* 14:7387 (2023); authors' code at
https://github.com/MichelleBoyle/scRNAseq_malaria_2023. This repository is an independent reanalysis.

Author: Mahnoor Asif
