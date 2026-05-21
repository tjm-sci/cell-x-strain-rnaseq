# EXP383 FANS RNA-seq: DESeq2-first troubleshooting plan

## Purpose

This document is an ordered troubleshooting runbook for a Codex/code agent working in the repository:

```text
~/exp383_code/cell-x-strain-rnaseq
```

The immediate problem is that the first-pass DESeq2 differential expression analysis, especially in PU1 nuclei, produced alarmingly few hits for prion-infected versus CBH/uninfected contrasts across strains and dpi levels. The current data were generated from FANS-sorted nuclei using nf-core/rnaseq QC with fastp and transcriptome-only Salmon quantification, followed by tximport to a gene-level count object. The unified object was then split into four population-level analyses.

The first troubleshooting phase must focus only on the following six areas:

1. Re-test DESeq2 designs before re-running nf-core.
2. Test whether population-specific covariate adjustment is suppressing signal.
3. Check confounding between retained covariates and biological variables.
4. Distinguish bad-sample effects from covariate-adjustable nuisance effects.
5. Test whether population-split DESeq2 fitting is costing power.
6. Inspect expected positive-control genes directly rather than relying only on total DEG counts.

The nf-core re-run / STAR-Salmon / exon-plus-intron investigation comes after the quick DESeq2 troubleshooting pass.

---

## Global instructions for the Codex agent

### Do not overwrite current primary results

The existing outputs under:

```text
results/dea/03_deseq2_dea/
```

are the current first-pass DEA results. Preserve them.

All troubleshooting outputs must be written to new directories with names that clearly mark them as exploratory/troubleshooting analyses.

Recommended root:

```text
results/dea/troubleshooting/
```

Suggested substructure:

```text
results/dea/troubleshooting/
├── 01_deseq2_design_ablation/
├── 02_covariate_confounding/
├── 03_sample_qc_sensitivity/
├── 04_combined_population_model/
├── 05_positive_control_gene_audit/
├── 06_summary_reports/
└── 99_nfcore_star_salmon_planning/
```

### Reuse existing code aggressively

Do not build a parallel pipeline from scratch.

The existing DESeq2 scripts already contain the required mechanics for loading input objects, fitting models, exporting contrast tables, generating MA/volcano plots, saving `dds.rds`, saving `vst.rds`, exporting `resultsNames(dds)`, and writing model QC outputs.

Copy and edit the existing code where appropriate:

```text
scripts/03a_deseq2_dea_utilities.R
scripts/03b_deseq2_dea_fit_models.R
scripts/03c_deseq2_dea_export_results.R
```

Recommended approach:

```text
scripts/troubleshooting/
├── 04a_deseq2_troubleshooting_utilities.R
├── 04b_deseq2_design_ablation.R
├── 04c_covariate_confounding_audit.R
├── 04d_sample_qc_sensitivity.R
├── 04e_combined_population_model.R
├── 04f_positive_control_gene_audit.R
└── 04g_write_troubleshooting_summary.R
```

These scripts may source the existing utility scripts rather than duplicating functions unnecessarily.

### Retain comparable output formats

Where possible, preserve the same output file formats and naming conventions used in the current DEA outputs:

```text
contrast_summary.tsv
model_fit_summary.tsv
model_formula.txt
design_matrix.tsv
results_names.txt
sample_model_qc.tsv
vst_pca_scores.tsv
tables/all_genes/*.tsv
tables/significant_genes/*.tsv
tables/significant_lfc_threshold/*.tsv
plots/ma/*.png
plots/volcano_raw_p/*.png
plots/volcano_adjusted_p/*.png
```

This makes it easy to compare troubleshooting outputs to:

```text
results/dea/03_deseq2_dea/<POPULATION>/
```

### Treat covariates as population-specific

The covariate screen was run after splitting the data by population. Therefore, each population has its own retained covariates and recommended design formula:

```text
results/dea/02_covariate_screening/NeuN/covariate_screening/
results/dea/02_covariate_screening/PU1/covariate_screening/
results/dea/02_covariate_screening/SOX10/covariate_screening/
results/dea/02_covariate_screening/SOX2/covariate_screening/
```

In particular, read these per population:

```text
recommended_design_formula.txt
selected_covariates.tsv
selected_covariate_biological_associations.tsv
selected_numeric_covariates_by_biological_variables.tsv
selected_categorical_covariates_by_biological_variables.tsv
```

Do not assume the same formula across populations.

### Use the label-corrected input object

The repo contains both original and label-corrected input stages. Prefer the label-corrected object unless the existing DESeq2 scripts clearly use something else:

```text
results/dea/01e_sample_label_correction/exp383_salmon_gene_input_label_corrected.rds
results/dea/01e_sample_label_correction/sample_metadata_label_corrected.tsv
```

Also refer back to:

```text
results/dea/01_build_salmon_gene_inputs/
```

for tximport-derived gene input, gene annotation, tx2gene, and sample QC metrics.

---

# Phase 1 — DESeq2-first troubleshooting

## 1. Design ablation: test whether the retained covariates are suppressing biological signal

### Rationale

The current first-pass DESeq2 model uses empirically retained covariates from population-level screening. This is sensible as a starting point, but several retained covariates may absorb biological signal if they are correlated with inoculum, dpi, group assignment, or sample quality changes caused by disease.

Examples of potentially risky covariates include:

```text
date_nuc_prep_days
filtered_total_sequences_millions
salmon_num_mapped_millions
salmon_percent_mapped
rrna_alignment_percent
fastp_percent_adapter
filtered_percent_gc
inoculation_batch
```

The agent should run a controlled ladder of DESeq2 designs for each population, with emphasis on PU1 first.

### Target output directory

```text
results/dea/troubleshooting/01_deseq2_design_ablation/
```

Suggested per-population structure:

```text
results/dea/troubleshooting/01_deseq2_design_ablation/
├── PU1/
│   ├── design_01_minimal_group/
│   ├── design_02_known_batch_group/
│   ├── design_03_date_group/
│   ├── design_04_selected_non_qc_covariates_group/
│   ├── design_05_selected_qc_covariates_group/
│   └── design_06_original_screened_design/
├── NeuN/
├── SOX10/
└── SOX2/
```

### Required models

For each population, run at least the following DESeq2 models.

#### Model 1: minimal biological model

Use this as the baseline.

```r
design(dds) <- ~ group_assignment
```

or, if the current scripts use `group_id`, use the existing biological grouping variable consistently:

```r
design(dds) <- ~ group_id
```

The agent must inspect the existing metadata and current model scripts to determine whether the canonical variable is `group_assignment`, `group_id`, or another equivalent variable.

#### Model 2: known batch plus biology

If present and not perfectly confounded:

```r
design(dds) <- ~ inoculation_batch + group_id
```

This model tests whether a clearly known experimental batch explains variance without absorbing as much signal as several continuous QC covariates.

#### Model 3: lab-date plus biology

For populations where `date_nuc_prep_days` or equivalent was retained:

```r
design(dds) <- ~ date_nuc_prep_days_z + group_id
```

Use the actual z-scored variable name present in the sample metadata. Do not invent a new column if the pipeline already created one.

#### Model 4: selected non-QC covariates plus biology

Use selected covariates that are plausibly experimental/technical but not direct sequencing depth or mapping metrics.

Examples:

```text
date_nuc_prep_days_z
inoculation_batch
```

Avoid including mapping depth, percent mapped, total reads, rRNA percentage, or GC percentage in this model.

#### Model 5: selected QC covariates plus biology

Include sequencing/QC covariates from the selected population-specific covariate list.

Examples:

```text
salmon_percent_mapped
salmon_num_mapped_millions
filtered_total_sequences_millions
rrna_alignment_percent
fastp_percent_adapter
filtered_percent_gc
```

This is a stress test for whether QC adjustment collapses expected disease signal.

#### Model 6: original screened design

Reproduce the existing screened design from:

```text
results/dea/02_covariate_screening/<POPULATION>/covariate_screening/recommended_design_formula.txt
```

This should match the current first-pass model as closely as possible.

### Required comparisons

For every population and every design model:

1. Fit DESeq2.
2. Export all current infected-vs-CBH contrasts:
   - `RML_60_vs_CBH_60`
   - `RML_90_vs_CBH_90`
   - `RML_120_vs_CBH_120`
   - `ME7_60_vs_CBH_60`
   - `ME7_90_vs_CBH_90`
   - `ME7_120_vs_CBH_120`
   - `22L_60_vs_CBH_60`
   - `22L_90_vs_CBH_90`
   - `22L_120_vs_CBH_120`
3. Export the same tables and plots as the current DEA scripts.
4. Create a combined comparison table:

```text
results/dea/troubleshooting/01_deseq2_design_ablation/design_ablation_contrast_summary.tsv
```

Minimum columns:

```text
population
design_id
design_label
model_formula
contrast
n_genes_tested
n_nominal_p_lt_0_05
n_padj_lt_0_1
n_padj_lt_0_05
n_padj_lt_0_01
median_abs_lfc
median_lfcSE
n_genes_with_padj_NA
n_genes_with_pvalue_NA
```

### Key decision logic

After running the design ladder, classify each contrast into one of these patterns:

```text
A. Signal present in minimal model but lost in screened model
B. Signal absent in all models
C. Signal only appears after batch/date adjustment
D. Signal present at raw p-value level but weak after FDR
E. Signal present only in later dpi or only for specific strains
```

Pattern A strongly suggests over-adjustment or covariate/biology confounding.

Pattern B points toward quantification, power, biology, or sample quality.

Pattern C supports retaining some adjustment.

Pattern D supports moving to pathway/rank-based analysis later, but do not implement that yet in this phase unless requested.

---

## 2. Covariate confounding audit

### Rationale

Before deciding that a selected covariate belongs in the final DESeq2 model, quantify whether it is associated with the biological design. A covariate that is strongly associated with `group_id`, inoculum, dpi, or infected-vs-CBH status can remove true biological signal.

This audit is especially important because the retained covariates differ across the four population-specific covariate screens.

### Target output directory

```text
results/dea/troubleshooting/02_covariate_confounding/
```

Suggested structure:

```text
results/dea/troubleshooting/02_covariate_confounding/
├── PU1/
│   ├── tables/
│   └── plots/
├── NeuN/
├── SOX10/
└── SOX2/
```

### Inputs

For each population:

```text
results/dea/02_covariate_screening/<POPULATION>/covariate_screening/selected_covariates.tsv
results/dea/02_covariate_screening/<POPULATION>/covariate_screening/recommended_design_formula.txt
results/dea/02_covariate_screening/<POPULATION>/covariate_screening/selected_covariates_by_biological_variables_data/
```

Also use the sample metadata used for DESeq2 model fitting.

### Required plots

For each retained numeric covariate, generate boxplots/jitter plots against:

```text
group_id or group_assignment
inoculum
dpi
infected_status
```

Use the actual column names in the metadata.

Example plot names:

```text
PU1/plots/date_nuc_prep_days_by_group_id.png
PU1/plots/salmon_percent_mapped_by_group_id.png
PU1/plots/rrna_alignment_percent_by_inoculum.png
```

For retained categorical covariates, generate stacked bar plots and contingency heatmaps against:

```text
group_id or group_assignment
inoculum
dpi
infected_status
```

### Required statistical summaries

For numeric covariates, fit simple models:

```r
lm(covariate ~ group_id, data = coldata_population)
lm(covariate ~ inoculum, data = coldata_population)
lm(covariate ~ dpi, data = coldata_population)
lm(covariate ~ infected_status, data = coldata_population)
```

For categorical covariates, use contingency tables and Fisher's exact tests or chi-squared tests as appropriate.

Export:

```text
results/dea/troubleshooting/02_covariate_confounding/covariate_biology_association_summary.tsv
```

Minimum columns:

```text
population
covariate
covariate_type
biological_variable
test_type
p_value
effect_size_or_r2
notes
```

### Required classification

Create:

```text
results/dea/troubleshooting/02_covariate_confounding/covariate_risk_classification.tsv
```

Minimum columns:

```text
population
covariate
retained_in_original_screen
included_in_original_formula
risk_class
reason
recommended_handling
```

Use these risk classes:

```text
low_risk
moderate_risk
high_risk
do_not_include_without_strong_justification
```

Suggested rules:

- Mark direct library size / mapped-read count covariates as high-risk unless there is a very clear reason to include them.
- Mark covariates strongly associated with `group_id`, inoculum, dpi, or infected status as high-risk.
- Mark known experimental batch variables as moderate risk if partially confounded and low risk if well balanced.
- Mark covariates with no obvious biological association and clear technical interpretation as lower risk.

---

## 3. Sample QC sensitivity: decide whether poor samples should be excluded rather than covariate-adjusted

### Rationale

A poor-quality sample should not automatically be rescued by adding QC covariates. Extreme samples can inflate dispersion, increase LFC standard errors, degrade PCA structure, and reduce DEG discovery.

The goal here is not to invent an exclusion scheme after seeing desired results. The goal is to determine whether a small number of technical outliers are dominating the weak DESeq2 signal.

### Target output directory

```text
results/dea/troubleshooting/03_sample_qc_sensitivity/
```

Suggested structure:

```text
results/dea/troubleshooting/03_sample_qc_sensitivity/
├── sample_qc_audit_all_populations.tsv
├── exclusion_sets/
├── PU1/
├── NeuN/
├── SOX10/
└── SOX2/
```

### Inputs

Use available QC and model files:

```text
results/dea/01_build_salmon_gene_inputs/sample_qc_metrics.tsv
results/dea/03_deseq2_dea/<POPULATION>/model_fit/sample_model_qc.tsv
results/dea/03_deseq2_dea/<POPULATION>/model_fit/vst_pca_scores.tsv
results/dea/03_deseq2_dea/<POPULATION>/model_fit/filter_summary.tsv
results/dea/03_deseq2_dea/<POPULATION>/model_fit/size_factors.png
results/dea/03_deseq2_dea/<POPULATION>/model_fit/sample_distance_heatmap.png
```

Also use MultiQC-derived sample-level metrics if already incorporated into `sample_qc_metrics.tsv`.

### Required QC audit table

Create:

```text
results/dea/troubleshooting/03_sample_qc_sensitivity/sample_qc_audit_all_populations.tsv
```

Minimum columns:

```text
sample_id
animal_id
population
group_id
inoculum
dpi
infected_status
salmon_percent_mapped
salmon_num_mapped_millions
rrna_alignment_percent
filtered_total_sequences_millions
n_expressed_genes
size_factor
pca1
pca2
pca_outlier_flag
sample_distance_outlier_flag
qc_outlier_flag
outlier_reason
```

Only include columns that exist; do not fail if some are absent. Log absent columns in a plain text note:

```text
results/dea/troubleshooting/03_sample_qc_sensitivity/missing_qc_columns.txt
```

### Define objective exclusion sets

Create at least these sample sets:

```text
set_00_no_exclusion
set_01_extreme_qc_outliers_only
set_02_pca_distance_outliers_only
set_03_combined_extreme_qc_and_pca_outliers
```

The thresholds must be written to:

```text
results/dea/troubleshooting/03_sample_qc_sensitivity/exclusion_sets/exclusion_thresholds.tsv
```

Suggested thresholds:

- extreme low Salmon mapping: below Q1 - 3 × IQR within population
- extreme high rRNA percentage: above Q3 + 3 × IQR within population
- extreme size factor: outside Q1 - 3 × IQR / Q3 + 3 × IQR within population
- PCA outlier: robust distance in first 2-5 PCs, threshold stated explicitly
- sample-distance outlier: unusually large mean distance to all other samples within population

If any threshold removes too many samples or creates empty/imbalanced groups, do not use that exclusion set for DESeq2. Record this in the summary.

### Required DESeq2 sensitivity fits

For each population, fit at least:

```text
set_00_no_exclusion + minimal design
set_00_no_exclusion + original screened design
set_03_combined_extreme_qc_and_pca_outliers + minimal design
set_03_combined_extreme_qc_and_pca_outliers + safer selected design
```

The safer selected design should come from the covariate-risk audit: avoid high-risk QC covariates and include only defensible low/moderate-risk covariates.

Export the same contrast results as in Phase 1.

### Required summary table

Create:

```text
results/dea/troubleshooting/03_sample_qc_sensitivity/sample_qc_sensitivity_contrast_summary.tsv
```

Minimum columns:

```text
population
exclusion_set
n_samples_removed
samples_removed
design_label
contrast
n_genes_tested
n_nominal_p_lt_0_05
n_padj_lt_0_1
n_padj_lt_0_05
median_lfcSE
notes
```

---

## 4. Combined-population DESeq2 sensitivity model

### Rationale

The current primary strategy splits the unified tximport gene object into four population-specific analyses. That is biologically defensible because the project is not asking for direct between-population differential expression.

However, splitting may reduce power for dispersion estimation and increase instability when each population has its own selected covariates. A combined-population sensitivity model can test whether the weak signal is partly a modelling-power issue.

This is a sensitivity analysis only. It does not automatically replace the population-specific primary analysis.

### Target output directory

```text
results/dea/troubleshooting/04_combined_population_model/
```

### Required model designs

Fit at least one simple combined model that avoids the population-specific covariate problem:

```r
design(dds_all) <- ~ population + group_id + population:group_id
```

If `population` and `group_id` have different names in the metadata, use the actual names. The model should allow population-specific biological effects through the interaction.

Alternatively, create a single combined factor:

```r
coldata$population_group <- interaction(coldata$population, coldata$group_id, sep = "_")
design(dds_all) <- ~ population_group
```

The `population_group` formulation is often easier for explicit contrasts.

### Required contrasts

Extract the same infected-vs-CBH contrasts within each population:

```text
PU1_RML_60_vs_PU1_CBH_60
PU1_RML_90_vs_PU1_CBH_90
PU1_RML_120_vs_PU1_CBH_120
...
NeuN_RML_60_vs_NeuN_CBH_60
...
SOX10_...
SOX2_...
```

Focus initial interpretation on PU1, but export all populations if the model fits successfully.

### Required output

Create:

```text
results/dea/troubleshooting/04_combined_population_model/combined_population_contrast_summary.tsv
```

Minimum columns:

```text
population
contrast
combined_model_n_padj_lt_0_05
split_original_n_padj_lt_0_05
split_minimal_n_padj_lt_0_05
combined_model_median_lfcSE
split_original_median_lfcSE
split_minimal_median_lfcSE
notes
```

Also export full contrast tables under:

```text
results/dea/troubleshooting/04_combined_population_model/contrast_results/tables/all_genes/
```

### Decision logic

If the combined model recovers expected PU1 signal while the split model does not, treat this as evidence that population-level splitting and/or population-specific covariate adjustment is reducing power.

If both combined and split models are weak, move attention to quantification strategy, sample quality, and nuclear RNA biology.

---

## 5. Positive-control gene audit

### Rationale

Do not evaluate the pipeline only by total DEG count. A low DEG count could reflect weak power, heavy multiple-testing burden, over-adjustment, quantification failure, or true absence of signal.

The agent must inspect expected genes directly and classify why they are or are not significant.

### Target output directory

```text
results/dea/troubleshooting/05_positive_control_gene_audit/
```

### Positive-control gene panels

Create a version-controlled input file:

```text
config/troubleshooting_positive_control_genes.tsv
```

Minimum columns:

```text
gene_symbol
panel
expected_population
expected_direction
notes
```

Start with this PU1/microglial/prion-response panel:

```text
Aif1
C1qa
C1qb
C1qc
Tyrobp
Trem2
Apoe
Lpl
Cst7
Itgax
Clec7a
Cx3cr1
P2ry12
Tmem119
Hexb
Spp1
Lgals3
B2m
H2-Ab1
Stat1
Irf7
Isg15
Ifit1
Ifit3
Oasl2
Gfap
Serpina3n
```

Notes:

- `Gfap` and `Serpina3n` are not expected PU1 markers but are useful contamination / cross-population sanity checks.
- Add literature-specific genes later if available, but do not block this audit on perfect curation.
- Map gene symbols to Ensembl IDs using the existing gene annotation table:

```text
results/dea/01_build_salmon_gene_inputs/gene_annotation.tsv
```

### Required audit outputs

For every design tested in Phase 1 and every PU1 contrast, export:

```text
results/dea/troubleshooting/05_positive_control_gene_audit/PU1_positive_control_gene_audit.tsv
```

Minimum columns:

```text
population
design_id
design_label
contrast
gene_symbol
gene_id
baseMean
log2FoldChange
lfcSE
stat
pvalue
padj
padj_is_na
pvalue_is_na
detected_in_count_matrix
passed_prefilter
classification
```

Use these classifications:

```text
not_detected
filtered_out
detected_low_baseMean
large_lfc_high_se
nominal_only
fdr_significant
small_lfc
opposite_direction
missing_gene_mapping
```

Also generate normalized-count plots for positive-control genes:

```text
results/dea/troubleshooting/05_positive_control_gene_audit/plots/PU1/<gene_symbol>_normalised_counts_by_group.png
```

Use DESeq2 normalized counts or VST values consistently and state which was used in the plot subtitle/caption.

### Required summary outputs

Create:

```text
results/dea/troubleshooting/05_positive_control_gene_audit/positive_control_summary_by_design.tsv
```

Minimum columns:

```text
population
design_id
contrast
n_positive_control_genes_in_matrix
n_not_detected
n_filtered_out
n_nominal_only
n_fdr_significant
n_large_lfc_high_se
n_opposite_direction
median_positive_control_lfc
median_positive_control_lfcSE
```

### Decision logic

Use this audit to answer:

1. Are expected PU1 genes absent from the Salmon/tximport count matrix?
2. Are they present but lowly counted?
3. Are they changing in the expected direction but missing FDR?
4. Are they losing LFC only after covariate adjustment?
5. Are they high-LFC but high-SE, suggesting dispersion/sample-quality problems?
6. Are astrocyte or other contamination markers unexpectedly driving signal in PU1?

---

## 6. Filtering and DESeq2 result-state audit

### Rationale

Some expected genes may have `padj = NA` or `pvalue = NA`. This can happen through independent filtering, outlier handling, all-zero rows, or model-fitting issues.

The agent should make the result state explicit for every contrast and especially for the positive-control genes.

### Target output directory

```text
results/dea/troubleshooting/06_summary_reports/
```

### Required checks

For each population, design, and contrast:

```r
summary(res)
table(is.na(res$pvalue), is.na(res$padj))
```

Run diagnostic-only versions of selected contrasts with:

```r
results(
  dds,
  contrast = c("group_id", "RML_120", "CBH_120"),
  independentFiltering = FALSE,
  cooksCutoff = FALSE
)
```

Use the actual contrast variable names from the fitted model.

Do not treat this diagnostic run as the final result. Its purpose is to determine whether genes disappear because of filtering/outlier rules.

### Required output

Create:

```text
results/dea/troubleshooting/06_summary_reports/deseq2_result_state_audit.tsv
```

Minimum columns:

```text
population
design_id
contrast
n_total_rows
n_all_zero_or_not_tested
n_pvalue_na
n_padj_na
n_padj_na_but_pvalue_present
n_nominal_p_lt_0_05
n_padj_lt_0_05
n_nominal_p_lt_0_05_diagnostic_no_filter_no_cooks
n_padj_lt_0_05_diagnostic_no_filter_no_cooks
notes
```

Also create a human-readable markdown summary:

```text
results/dea/troubleshooting/06_summary_reports/deseq2_troubleshooting_interim_summary.md
```

This summary must include:

```text
- Whether the minimal model recovers more signal than the screened model.
- Whether retained covariates are confounded with group, dpi, inoculum, or infection.
- Whether sample outliers materially change results.
- Whether combined-population modelling improves power.
- What happened to the positive-control gene panel.
- Whether many expected genes are lost through filtering / Cook's cutoff / NA padj.
- Whether the evidence now points mainly to DESeq2 modelling or upstream quantification.
```

---

# Phase 2 — nf-core / quantification follow-up after DESeq2 troubleshooting

## 7. STAR-Salmon / genome-aligned re-run planning

### Rationale

The current data were processed using transcriptome-only Salmon quantification and showed low mapping rates, which is a major concern for FANS-sorted nuclei. Nuclear RNA contains substantial unspliced/intronic signal. Transcriptome-only quantification against mature transcript annotations may undercount nuclear signal.

Therefore, after the DESeq2-first troubleshooting pass, plan a genome-aligned re-run using nf-core/rnaseq with STAR-Salmon or an equivalent genome-aligned strategy. The goal is to recover exonic and intronic signal, or at minimum to quantify how much nuclear signal is being missed by transcriptome-only Salmon.

### Storage-aware execution plan

Available storage for the STAR-Salmon troubleshooting phase:

```text
1TB boot drive
2TB SSD
4TB HDD containing the current raw/trimmed FASTQ files
```

Use the disks as follows:

```text
4TB HDD  = read-only FASTQ input store
2TB SSD  = Nextflow STAR-Salmon work/scratch directory
1TB boot = repo, scripts, logs, small final results only
```

Do not run all 288 samples through STAR-Salmon initially. Run in this order:

```text
1. Tiny PU1 smoke test: 4-8 samples.
2. PU1-only pilot: all PU1 samples.
3. Review peak work-dir size, final outdir size, mapping rescue, and positive-control gene recovery.
4. Only then consider NeuN, SOX10, and SOX2 one population at a time.
```

Do not save BAM/intermediate alignment files for the pilot unless explicitly required:

```yaml
save_align_intermeds: false
save_unaligned: false
save_trimmed: false
save_non_ribo_reads: false
save_reference: false
```

The STAR-Salmon pilot command should use an external work directory on the 2TB SSD, for example:

```bash
nextflow run nf-core/rnaseq \
  -profile singularity \
  -params-file params/nfcore_rnaseq_star_salmon_troubleshooting.yaml \
  --input metadata/batches/exp383_nfcore_rnaseq_PU1.csv \
  --outdir results/nfcore_rnaseq_star_salmon_troubleshooting/PU1 \
  -work-dir /mnt/2tb_ssd/exp383_nextflow_work/star_salmon_PU1
```

Replace `/mnt/2tb_ssd/` with the actual SSD mount point.

Add a storage audit file for the pilot:

```text
results/dea/troubleshooting/99_nfcore_star_salmon_planning/star_salmon_storage_audit.tsv
```

Minimum columns:

```text
population
n_samples
work_dir_peak_gb
final_outdir_gb
save_align_intermeds
runtime
mean_star_uniquely_mapped_percent
mean_salmon_mapping_percent
notes
```

If PU1-only STAR-Salmon already approaches the available 2TB SSD limit, do not proceed to the remaining populations until work-dir cleanup, scratch relocation, or a reduced output strategy has been agreed.

### Target output directory for planning notes

```text
results/dea/troubleshooting/99_nfcore_star_salmon_planning/
```

### Required planning outputs

Create:

```text
results/dea/troubleshooting/99_nfcore_star_salmon_planning/nfcore_star_salmon_rerun_plan.md
```

Include:

1. Existing nf-core inputs:
   ```text
   scripts/00_run_nfcore_rnaseq.sh
   params/nfcore_rnaseq_fastp_salmon_rrna.yaml
   metadata/batches/exp383_nfcore_rnaseq_PU1.csv
   metadata/batches/exp383_nfcore_rnaseq_NeuN.csv
   metadata/batches/exp383_nfcore_rnaseq_SOX10.csv
   metadata/batches/exp383_nfcore_rnaseq_SOX2.csv
   resources/reference/Mus_musculus.GRCm39.115.gtf.gz
   resources/reference/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
   ```

2. Proposed new parameters file:
   ```text
   params/nfcore_rnaseq_star_salmon_troubleshooting.yaml
   ```

3. Proposed new logs/output directories:
   ```text
   logs/nextflow/troubleshooting_star_salmon/
   results/nfcore_rnaseq_star_salmon_troubleshooting/
   ```

4. A note that this re-run should not overwrite existing nf-core results.

### STAR-Salmon and tximport-specific requirements

The agent must include the following in the planning document.

#### STAR-Salmon mode

The re-run should use STAR-based genome alignment with Salmon quantification if supported by the current nf-core/rnaseq version and local reference/index setup.

The purpose is not merely to get Salmon output again, but to obtain genome-alignment-level evidence about:

```text
exonic assignment
intronic assignment
intergenic/unassigned reads
gene-body coverage
splice-aware mapping rates
```

#### Exon-plus-intron counting

The agent should evaluate options for producing gene-level counts that include intronic reads. Potential strategies include:

```text
featureCounts against exon annotations
featureCounts against gene-body annotations
STARsolo-style exon+intron counting if appropriate
custom GTF transformation from exon features to gene-body intervals
separate exon-only versus exon+intron count matrices
```

The final choice should be documented before implementation.

#### tximport setting must be correct

If Salmon transcript-level output is used again and imported through tximport, the tximport settings must be explicit and appropriate.

The agent must check and record:

```r
txi$countsFromAbundance
```

The agent must also check whether the DESeq2 object was created via:

```r
DESeqDataSetFromTximport(txi, colData = sample_table, design = ...)
```

or by using a matrix constructor such as:

```r
DESeqDataSetFromMatrix(...)
```

If using original Salmon estimated counts, prefer `DESeqDataSetFromTximport()` so that length offsets are handled correctly.

If using `countsFromAbundance = "lengthScaledTPM"` or another scaled option, document why and ensure downstream DESeq2 construction is consistent with that choice.

Do not mix tximport modes across primary and troubleshooting analyses without clearly naming the outputs.

### Required post-rerun comparison plan

Once genome-aligned / exon-plus-intron counts exist, compare at least:

```text
A. Existing transcriptome-only Salmon/tximport gene counts
B. STAR/featureCounts exon-only gene counts
C. STAR/featureCounts exon+intron or gene-body counts
D. STAR-Salmon transcript-derived gene counts, if retained
```

For each count source, rerun:

```text
minimal DESeq2 design
safer selected design from Phase 1/2
positive-control gene audit
```

Write these outputs under a new root, for example:

```text
results/dea/troubleshooting/star_salmon_count_source_comparison/
```

---

# Final deliverables expected from the agent

At the end of this troubleshooting pass, produce the following files:

```text
results/dea/troubleshooting/01_deseq2_design_ablation/design_ablation_contrast_summary.tsv
results/dea/troubleshooting/02_covariate_confounding/covariate_biology_association_summary.tsv
results/dea/troubleshooting/02_covariate_confounding/covariate_risk_classification.tsv
results/dea/troubleshooting/03_sample_qc_sensitivity/sample_qc_audit_all_populations.tsv
results/dea/troubleshooting/03_sample_qc_sensitivity/sample_qc_sensitivity_contrast_summary.tsv
results/dea/troubleshooting/04_combined_population_model/combined_population_contrast_summary.tsv
results/dea/troubleshooting/05_positive_control_gene_audit/PU1_positive_control_gene_audit.tsv
results/dea/troubleshooting/05_positive_control_gene_audit/positive_control_summary_by_design.tsv
results/dea/troubleshooting/06_summary_reports/deseq2_result_state_audit.tsv
results/dea/troubleshooting/06_summary_reports/deseq2_troubleshooting_interim_summary.md
results/dea/troubleshooting/99_nfcore_star_salmon_planning/nfcore_star_salmon_rerun_plan.md
```

Also produce or update these troubleshooting scripts:

```text
scripts/troubleshooting/04a_deseq2_troubleshooting_utilities.R
scripts/troubleshooting/04b_deseq2_design_ablation.R
scripts/troubleshooting/04c_covariate_confounding_audit.R
scripts/troubleshooting/04d_sample_qc_sensitivity.R
scripts/troubleshooting/04e_combined_population_model.R
scripts/troubleshooting/04f_positive_control_gene_audit.R
scripts/troubleshooting/04g_write_troubleshooting_summary.R
```

---

# Interpretation rules

Use these rules when writing the interim summary.

## If minimal models recover signal but screened models do not

Conclude that the current covariate strategy is likely over-adjusting, especially if retained covariates are associated with group, inoculum, dpi, or infected status.

Recommended next action:

```text
Move toward a reduced, biologically defensible design rather than the maximally screened design.
```

## If expected genes show LFCs in the correct direction but weak FDR

Conclude that the signal exists but may be underpowered or dispersed.

Recommended next action:

```text
Proceed to pathway/rank-based analysis later, and evaluate sample QC / dispersion.
```

## If expected genes are absent or extremely low in the count matrix

Conclude that transcriptome-only Salmon may be failing for nuclear RNA.

Recommended next action:

```text
Prioritise STAR/genome-aligned exon+intron counting.
```

## If expected genes are present in PU1 but only after outlier exclusion

Conclude that sample quality is materially affecting inference.

Recommended next action:

```text
Define objective exclusion criteria, document them, and avoid ad hoc post hoc removal.
```

## If combined-population modelling improves signal

Conclude that population-level splitting may be costing power.

Recommended next action:

```text
Consider using the combined model as a sensitivity analysis or possibly as the primary model if it remains interpretable.
```

## If nothing improves in DESeq2-first troubleshooting

Conclude that the most likely next issue is upstream quantification, especially given FANS nuclei, transcriptome-only Salmon, and low mapping rates.

Recommended next action:

```text
Move to nf-core STAR-Salmon / genome-aligned exon+intron count generation.
```
