# Stool antigen immunoassays for *Campylobacter* — diagnostic test accuracy meta-analysis

Analysis code and extracted data for a systematic review and bivariate meta-analysis of the
diagnostic accuracy of stool antigen immunoassays for *Campylobacter* spp., stratified by platform
(lateral flow immunoassay vs microplate enzyme immunoassay) and by reference standard (culture vs
PCR).

24 studies, 68 test × reference 2×2 estimates, 13 677 stool specimens.

This repository holds what is needed to **reproduce the calculations**. The figures, the summary
tables and their interpretation are published with the article and are not duplicated here.

> **Status:** manuscript under preparation. This repository will be tagged at submission and the
> article citation added here.

## Contents

```
data/
  DATA_METAANALYSIS.xlsx   extraction workbook (sheet MetaDTA = one row per study × index test × reference)
  patch_datos.csv          documented correction layer applied on top of the workbook (see below)
R/
  01_prep_datos.R          builds the analysis datasets: eligibility, deduplication, comparison labels
  02_analisis_bivariado.R  Reitsma bivariate models, summary points, SROC, prediction regions
  03_extensiones.R         heterogeneity, QUADAS-2, meta-regression, subgroups, sensitivity, leave-one-out, Deeks
  03b_i2_zhou.R            Zhou & Dendukuri I² recomputed from the fitted bivariate model
  04_bayesiano_jags.R      Bayesian bivariate binomial-normal model (Chu & Cole) in JAGS
  05_subgrupos_marca.R     brand-level subgroups, brand-effect test, study characteristics table
  06_bayes_priors.R        prior sensitivity (three specifications)
run_all.sh                 runs the seven scripts in order
sessionInfo.txt            R, JAGS and package versions used for the published run
```

## Running it

Requires R ≥ 4.4 with `readxl`, `dplyr`, `stringr`, `tidyr`, `mada`, `metafor`, `mixmeta`,
`writexl` and `ellipse`, plus [JAGS](https://mcmc-jags.sourceforge.io/) ≥ 4.3 with `rjags` and
`coda` for the Bayesian step. Exact versions are in `sessionInfo.txt`.

```bash
bash run_all.sh
```

Every result table is written to `salidas/` (not tracked): pooled estimates, heterogeneity,
meta-regression, subgroups, sensitivity scenarios, leave-one-out, Deeks regression and the Bayesian
posteriors. `01_prep_datos.R` looks for the workbook in `data/`; override with the `CAMPY_DATA`
environment variable.

The scripts are commented in Spanish and so are the column names of the output tables; the figures
and tables as published are in English.

## Analysis decisions worth knowing before you read the code

**Search and its update.** The protocol search (PROSPERO CRD420251160098) closed in December 2025.
Because manuscript preparation ran long, the search was updated in **PubMed only** to October 2026;
the two eligible reports it identified are included and the update is shown as a separate branch of
the PRISMA flow diagram. The update was not replicated in the other four databases, which is a
limitation of the search process.

**Deduplication.** A study that evaluated several brands on the same platform against the same
reference standard contributes one entry to each pooled comparison, chosen by a fixed market-share
ranking declared in `01_prep_datos.R` — not by its results. The full non-deduplicated set is kept
and used for the brand-level analysis and as a sensitivity scenario.

**Four pooled comparisons.** Platform × reference standard: LFA vs culture (k = 11), microplate EIA
vs culture (k = 11), LFA vs PCR (k = 6), microplate EIA vs PCR (k = 4) — 32 estimates from 17
studies. Composite reference standards, filtration/enrichment culture, rapid membrane EIA and
chemiluminescent assays are described narratively and never pooled.

**Correction layer.** `patch_datos.csv` carries corrections and additions applied on top of the
extraction workbook, each with its reason in the `notes` field, rather than editing the workbook in
place. It covers a recalculated 2×2 table for one study whose original extraction did not match its
declared reference standard, two entries added from the same study, the two reports added by the
search update, and one study reconstructed from its abstract that is flagged
`solo_sensibilidad = 1` and therefore enters only the sensitivity analysis. One of the added
reports publishes a single 2×2 table for two kits of the same platform evaluated on the same
specimens; it is entered once, as a combined entry, not as two brands.

**QUADAS-2 harmonisation.** The extraction workbook rated the same sampling design inconsistently
across studies. `01_prep_datos.R` applies a single rule: non-consecutive, non-random sampling
(`preselected = 1`) is rated **high risk** in the patient-selection domain, which QUADAS-2 is
explicit about. The rule and the entries it affects are documented in the script.

**Continuity correction.** `reitsma()` is used with its default `correction = 0.5` and
`correction.control = "all"`: when any cell of any study in a comparison is zero, 0.5 is added to
all four cells of **every** study in that comparison — not only to the zero. This is active in three
of the four comparisons.

**Heterogeneity.** The `i2` component of `summary.reitsma()` in mada 0.5.12 computes the
between-study term from the fixed-effect standard errors rather than from the between-study
covariance matrix, and uses TP+FP / TN+FN as group sizes rather than TP+FN / FP+TN. `03b_i2_zhou.R`
recomputes the Zhou & Dendukuri statistic from `Psi` and the diseased / non-diseased group sizes;
both versions are reported side by side.

**Multiplicity.** Meta-regression is univariable, one moderator per model, with Benjamini-Hochberg
control applied within each comparison. Moderators that induce the same partition of studies are
grouped and a single representative is reported; the grouping is in the `grupo_colineal` column and
the size of each correction family in `n_familia_BH`.

## Citation

See `CITATION.cff`. Until the article is published, cite this repository by its commit hash.

## License

Code: MIT (`LICENSE`). Data: CC BY 4.0 — reuse freely with attribution to the authors and to the
primary studies, which remain the copyright of their respective publishers.
