# How the pipeline works

What each step does and why. For how to run it see [usage.md](usage.md); for what it
writes see [output.md](output.md).

![nf-rare-var-assoc pipeline](draw_io_diagrams/nf_rare_var_assoc_metromap.drawio.png "Pipeline diagram")

## The problem this pipeline addresses

In a rare-variant association study, no single variant occurs often enough to be tested
on its own. Variants are instead grouped -- usually by gene -- and each group is tested as
a whole. That makes the result depend heavily on decisions taken long before the
statistical test: which variants survive quality control, which are assigned to which
gene, how a genotype's uncertainty is represented, and how population structure is
handled. A mistake in any of these produces a result that looks perfectly plausible.

This pipeline performs all of those steps in one workflow, with each recording what it
did, so the result can be traced back through every decision that produced it.

## 1. Preparing the VCF

Delegated to [nf-prepare-vcf](https://github.com/psuszyns/nf-prepare-vcf), and skipped
when `--skip_preparation true` is set.

- **Sorting, splitting multi-allelic sites and removing exact duplicates.** A site with
  several alternative alleles becomes one record per allele, so that each can be tested
  and annotated separately.
- **Chromosome renaming and variant identifiers.** Every variant is given the identifier
  `CHROM_POS_REF_ALT`. This happens *after* the multi-allelic split, so that identifiers
  stay unique. All later steps match variants to genes by this identifier, so it has to
  be assigned consistently.
- **Left-alignment is deliberately switched off** (`--do-not-normalize`). Left-aligning
  indels would shift their positions and therefore change their identifiers, breaking
  the correspondence with the annotation files. This is a deliberate departure from the
  usual convention, made so that variant identifiers stay stable through the pipeline.
- **VEP annotation.** Predicted consequences are written into the VCF's `CSQ` field, and
  are what the gene grouping step later reads.
- **Dosage computation from genotype likelihoods.** A hard genotype call throws away how
  confident the caller was. Instead, a dosage -- the expected number of alternative
  alleles, between 0 and 2 -- is computed from the `PL` genotype likelihoods and stored
  in a `DS` field. The computation also corrects calls where every likelihood is zero,
  which some callers emit for homozygous-reference sites, using the genotype quality
  instead. Only genotypes above a minimum quality contribute.
- **Gene ranges and conversion to PLINK format.**

## 2. Variant and genotype quality filtering

Two filters, applied to the VCF:

- **Per-variant**, on site quality and on the average genotype quality and read depth
  across samples. A variant failing these is removed entirely.
- **Per-genotype**, on that genotype's own quality and depth. A genotype failing these is
  **set to missing rather than removed**, so one poorly covered sample does not cost every
  other sample the variant.

The thresholds are the `--filter_vcf_*` parameters in [usage.md](usage.md).

## 3. Sample quality control

- **Sex imputation** from chromosome X heterozygosity. This runs as a separate PLINK2
  invocation, because PLINK2 requires the pseudo-autosomal region to be split first, and
  it cannot do both in one call. The imputed sex is used from then on, rather than any
  reported sex, so the pipeline does not depend on sample metadata being correct.
- **Per-phenotype missingness filtering.** Missingness is measured **within each phenotype
  group separately** and the results intersected, rather than across all samples at once.
  A variant that is well covered in controls but poorly covered in cases has an
  acceptable overall missingness while being a serious source of false positives; measuring
  each group separately catches it.
- **Removal of inbreeding-coefficient outliers.** The inbreeding coefficient measures
  observed against expected heterozygosity. Samples more than
  `--inbreeding_outliers_range_stds` standard deviations from the mean are removed --
  usually indicating contamination, sample mix-up, or unusually poor quality.

## 4. Population structure

Ancestry differences between cases and controls create associations that have nothing to
do with the phenotype, and are the most common cause of false results in association
studies. The correction has four stages:

1. **Linkage-disequilibrium pruning** to an approximately independent set of common
   variants, excluding regions of unusually strong linkage disequilibrium listed in
   `--hild_path` -- the major histocompatibility complex and a small number of large
   inversions, whose structure would otherwise dominate the components.
2. **Relatedness estimation**, to choose a set of mutually unrelated samples.
3. **Principal-component analysis** on that unrelated set.
4. **Projection** of *every* sample, including the related ones, onto those components.

The resulting scores become covariates in both association steps.

**Related samples are not removed from the analysis.** The relatedness cutoff
(`--plink2_king_cutoff_threshold_pca`) only decides which samples the components are
*computed from*, because principal components are distorted if close relatives dominate
the sample. All samples are then projected onto those components and all of them enter
the association test, where REGENIE's whole-genome model accounts for the remaining
relatedness. Describing this step as relatedness filtering would be wrong.

## 5. The four quality-control paths, and why the thresholds differ

After the shared steps above, the data splits into four paths, each with its own PLINK2
thresholds. This is deliberate: the four uses have genuinely different requirements, and
a single set of thresholds would be wrong for at least three of them.

| Path | What it feeds | What it needs |
|---|---|---|
| Inbreeding-coefficient filtering | sample quality control | common, independent variants -- the coefficient is meaningless on rare ones |
| Principal-component analysis | covariates | common, independent variants, pruned harder and over a wider window |
| REGENIE step 1 | the whole-genome model | common, high-quality variants; **rare variants contribute nothing here** |
| REGENIE step 2 | the association tests | **rare variants, kept** -- this is the only path where they matter, so frequency filters must not be applied |

The essential point is the contrast between the last two. REGENIE step 1 fits a
whole-genome model to capture relatedness and background genetic effects, and needs
common variants for that; the default `--plink2_makepgen_3_options` accordingly applies a
minor-allele-frequency floor. Step 2 performs the actual rare-variant tests, so applying
that same floor would remove exactly the variants the study is about. Its settings
(`--plink2_makepgen_4_options`, `--plink2_makepgen_5_options`) filter only on missingness.

Missingness thresholds also differ between the two: step 1 can afford to be strict,
because it only needs enough common variants to fit a model, while step 2 must be looser
or rare variants -- which are inherently harder to genotype -- would be lost.

## 6. Grouping variants into genes

Association testing needs to know which variants belong to which gene, and which of them
are worth testing together. The VEP consequences written into the VCF during preparation
are parsed, and three files are produced per dataset:

- a **variant set list**, assigning variants to genes,
- an **annotation file**, giving each variant its impact group,
- an **allele-frequency file**, used to split each gene's variants into frequency bands.

Impact groups come from `--input_masks`. The default defines `Mask_High` (variants
predicted to disrupt the protein: stop gained or lost, frameshift, splice site),
`Mask_Mod` (missense and in-frame insertions or deletions) and `Mask_HighMod` (both
together). Each group is tested separately, so a signal driven only by clearly disruptive
variants can be distinguished from one that also needs the milder ones.

Groups are **binary tiers**, not continuous weights: a variant either belongs to a group
or it does not. Some other tools weight each variant by a continuous functional score
instead. The advantage of tiers is that they are easy to interpret and apply equally to
indels; the disadvantage is that they cannot express fine gradations of predicted effect.

## 7. Association testing

REGENIE runs in two steps.

**Step 1** fits a whole-genome model on the common-variant set, producing a per-sample
prediction that captures relatedness and the polygenic background. This is what allows
related samples to be kept.

**Step 2** tests each gene, in each impact group, at each allele-frequency threshold in
`--aaf-bins`, using:

- a **burden test**, which sums the variants in a group and asks whether the total differs
  between cases and controls. Powerful when the variants act in the same direction.
- **SKAT-O**, which combines a burden test with a variance-component test. Retains power
  when some variants raise and others lower risk, which a burden test cancels out.
- **Firth correction**, which keeps p-values valid when cases are few or a variant is seen
  almost exclusively in one group -- the situation ordinary logistic regression handles
  badly, and the normal situation in rare-variant analysis.

With `--use_dosage true`, the dosages computed during preparation are used instead of hard
genotype calls, so genotype uncertainty carries through into the test.

## 8. Reporting and provenance

- **The main report** collects the Manhattan plot, QQ plot with genomic inflation factor,
  the strongest associations, phenotype summaries, principal-component plots and the
  data-quality figures into one HTML file.
- **Exploratory figures** cover genotype quality, read depth, dosage, missingness,
  heterozygosity and allele frequency, plotted separately for cases and controls. The
  case/control comparisons matter most: a systematic difference in coverage between the
  two groups produces association signals that are entirely technical.
- **Per-variant carrier tables** record, for every variant in a gene group, its allele
  counts and frequencies in cases and controls and the identifiers of the samples carrying
  it. This is what turns a significant gene into something that can be examined.
- **A data-flow record.** Every step writes how many samples and variants it received and
  returned, plus the parameters it used. These are assembled into a diagram and a text
  report, so any loss of data can be traced to the step that caused it.

Reporting is by far the slowest part of the pipeline. `--skip_reporting true` disables all
of it when only the association results are needed.

## Design decisions worth knowing about

**Left-alignment is off.** Standard practice, but it would change indel positions and
therefore variant identifiers, breaking the link between the VCF and the gene groupings.
Identifier stability was chosen over convention.

**Dosages are computed, imported late, and not used everywhere.** Dosages are exported and
then re-imported immediately before the association test rather than being carried through
the whole pipeline. Importing them at the start breaks the intermediate steps:
heterozygosity is computed incorrectly and sex imputation fails, because those calculations
assume discrete genotypes.

**Sex is imputed rather than read.** The pipeline derives sex from chromosome X
heterozygosity instead of trusting sample metadata, and acts on the result. On exome data,
where chromosome X is sparsely covered, the estimate is noisier and more samples fall
between the two thresholds -- widen them with `--plink2_makepgen_2_options` if this
happens.

**The focus is coding variation.** The impact groups are built from protein-coding
consequences. Non-coding variants are annotated but are not grouped into testable sets.

**Relatedness is modelled, not filtered.** See section 4.
