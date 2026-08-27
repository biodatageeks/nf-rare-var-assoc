# Comparison against RICOPILI combined with STAARpipeline

No single existing tool covers the same range of steps as nf-rare-var-assoc. This
comparison therefore builds a complete alternative out of two tools: RICOPILI performs
quality control and computes principal components, and STAARpipeline performs the
association test using its own functional annotations. Together they are treated as one
method and compared against this pipeline on the same simulated datasets.

The combined method receives only three things: the raw VCF, the phenotype files, and
the list of causal genes used for scoring. Everything else -- quality control, population
structure, annotation and the statistical test -- is done by RICOPILI and STAARpipeline
themselves.

The scripts that are not specific to these two tools (running this pipeline, scoring
results, comparing two methods statistically) are in
[../benchmark-common/](../benchmark-common/).

## Contents of this directory

| File | What it is |
|---|---|
| `Dockerfile.ricopili` | Builds a RICOPILI container. Upstream publishes no image. Three stages: unpack and prune the dependency archive, fetch `rp_bin` from git at a fixed commit, and assemble a runtime on `condaforge/miniforge3` using upstream's own pinned conda environment. |
| `Dockerfile.favorannotator` | Builds a FAVORannotator container. Upstream publishes no working image either. Based on `bioconductor/bioconductor_docker:RELEASE_3_20`, with gdsfmt, SeqArray, SeqVarTools, readr, dplyr, stringr, stringi, the `xsv` command-line tool, and a copy of `FAVORdatabase_chrsplit.csv`. |
| `convert_vcf_to_gds.R` | Converts a VCF into a SeqArray GDS file. Uses upstream's own conversion settings; paths are command-line arguments rather than constants in a configuration file. |
| `favorannotator_csv_essential.R` | Adds FAVOR functional annotations to a GDS file, producing an annotated GDS ("aGDS"). Upstream's annotation logic is unchanged; paths are parameters, the per-run database re-download is removed, and preflight checks are added. |
| `staar_gene_centric_coding.R` | Runs the association test: fits the STAAR null model and calls `Gene_Centric_Coding` over the annotated GDS, writing one CSV per functional category. Contains the mapping from FAVOR annotation channels to STAAR's `Annotation_name_catalog` (20 channels). `--use-spa` selects between the STAAR-O and STAAR-B tests; `--grm` accepts a relatedness matrix. |
| `genesis_pcair_pcrelate.R` | Computes relatedness-robust principal components (PC-AiR) and a relatedness matrix (PC-Relate) from a PLINK bed file, using GENESIS. Needed for the "Full" version described below. |
| `run_chain_ricopili_staar.sh` | Runs the whole combined method for one or more datasets, in both versions. This is the main entry point. |
| `run_staar_eval.sh` | Scores the association results of both versions. |
| `staar_to_eval.py` | Converts STAAR result CSVs into the REGENIE-shaped table that the scoring step expects. Maps `ID=<SYMBOL>.<category>` and `LOG10P=-log10(STAAR-O p-value)`, takes gene positions from `genes_info_hgnc.tsv`, and leaves BETA and A1FREQ empty. |
| `genes_info_hgnc.tsv` | Gene coordinates used by `staar_to_eval.py`: 18,445 rows (1000 on chromosome 12, 427 on chromosome 22, none on chromosome X). Exported from STAARpipeline's own internal data, see below. |

## Requirements

- **podman** or docker. All commands below use podman. On SELinux systems every mount
  needs `:z`, and rootless podman needs `--userns=keep-id` to write into directories
  owned by your user.
- **Disk space: roughly 250 GB.** The FAVOR annotation database is 23 GB compressed and
  expands about tenfold: 53 GB for chromosome 22 and about 196 GB for chromosome 12.
- **A reference genome with `chr`-prefixed chromosome names**, plus its `.fai` and `.gzi`
  index files. This work used the Genome in a Bottle GRCh38 file. An Ensembl-named
  reference (chromosome `22` rather than `chr22`) will not match the VCF.
- The raw VCF, the simulated datasets, and the 1000 Genomes pedigree file
  `assets/integrated_call_samples_v3.20250704.ALL.ped` from this repository.

Working data is kept outside the repository. The scripts default to this machine's
layout but every location is an environment variable, so the same scripts run elsewhere
unchanged -- set `DATA` and `RVA_REPO` and the rest follow (`REFDIR`, `DATASETS_DIR`,
`ARM`, `COMMON`, `INPUT_VCF`, `FAVOR_DB` and `RUN_DIR` can also be overridden
individually).

## One-time setup

```bash
DATA=/data/doktorat/biodatageeks/article_on_nf_rare_var_assoc
F=$DATA/tools_comparison/favor                  # annotation database and working files
D=$DATA/tools_comparison/ricopili               # RICOPILI dependency archive
ARM=/data/git/doktorat_pw/wum_pims/nf-rare-var-assoc/docs/tool-comparison/chained-benchmark
REFDIR=/data/doktorat/biodatageeks/genome_in_a_bottle/reference
REF=GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
mkdir -p $F/db $F/work $D
```

### 1. Build the two containers

```bash
cd $ARM
podman build -f Dockerfile.favorannotator -t favorannotator:1.0.0 . > /tmp/favor_build.log 2>&1
echo "BUILD_EXIT=$?"          # check this directly; see the notes about pipes below

podman run --rm localhost/favorannotator:1.0.0 bash -lc \
  'xsv --version; Rscript -e "library(gdsfmt);library(SeqArray);library(readr)"; ls /opt/favorannotator/'
```

RICOPILI needs its dependency archive first. Upstream publishes a checksum next to it;
use it.

```bash
cd $D
curl -sSLO https://personal.broadinstitute.org/braun/sharing/ricopili_dependencies_0225b.md5.cksum
curl -SL -C - --speed-limit 102400 --speed-time 60 \
     -O https://personal.broadinstitute.org/braun/sharing/ricopili_dependencies_0225b.tar.gz
md5sum -c ricopili_dependencies_0225b.md5.cksum      # must report OK

# the build context is the archive directory, not the repository
podman build -f $ARM/Dockerfile.ricopili -t ricopili:2025_Feb_20.001 "$D" > /tmp/rico_build.log 2>&1
echo "BUILD_EXIT=$?"

podman run --rm localhost/ricopili:2025_Feb_20.001 bash -lc '
  smartpca 2>&1 | head -3; /opt/rp_dep/plink/plink --version; R --version | head -1
  preimp_dir --help | head -5'
```

The RICOPILI image may be copied directly between machines with `podman save` and
`podman load`, but **must not be pushed to any container registry**: it contains
`rp_bin`, which makes the image a CC BY-NC-SA derivative.

Two more containers are pulled ready-made:

```bash
podman pull docker.io/zilinli/staarpipeline:0.9.7      # the association engine
podman pull docker.io/uwgac/topmed-roybranch:latest    # GENESIS, for the Full version
podman pull docker.io/psuszynski/bioinf_combo:1.5.1    # bcftools and friends
```

### 2. Download and extract the FAVOR annotation database

Downloads from Harvard Dataverse stall silently, so a retry loop with a minimum speed is
required. **Do not add `curl --retry`**: it silently truncates a `-C -` resume and can
destroy hours of progress.

```bash
cd $F
echo "0da5720d1c0d7a6e028a577ca697cfdf  chr22.tar.gz" > chr22.tar.gz.md5
echo "709dd56d9b5654009f161d010fadd597  chr12.tar.gz" > chr12.tar.gz.md5

for N in 22:6170504 12:6170520; do          # chromosome:Dataverse file id
  CHR=${N%%:*}; ID=${N##*:}
  for i in $(seq 1 60); do
    before=$(stat -c %s chr$CHR.tar.gz 2>/dev/null || echo 0)
    curl -SL -C - --speed-limit 102400 --speed-time 60 \
      -o chr$CHR.tar.gz "https://dataverse.harvard.edu/api/access/datafile/$ID" && break
    after=$(stat -c %s chr$CHR.tar.gz 2>/dev/null || echo 0)
    echo "chr$CHR attempt $i: $before -> $after bytes"
    [ "$after" -lt "$before" ] && echo "WARNING: file shrank, the resume is broken"
    sleep 15
  done
done
md5sum -c chr22.tar.gz.md5 chr12.tar.gz.md5     # both must report OK
```

Expected sizes: chromosome 22 is 5,574,054,308 bytes, chromosome 12 is 19,323,745,439
bytes. At 2-4 MB/s chromosome 12 takes roughly two hours.

The archives contain the authors' own cluster paths, so seven leading path components
have to be stripped:

```bash
df -h $F                                                     # check space first
tar -xzf $F/chr22.tar.gz -C $F/db --strip-components=7        # ~3 min  -> 53 GB
tar -xzf $F/chr12.tar.gz -C $F/db --strip-components=7        # ~9 min  -> ~196 GB
ls -1 $F/db/chr{12,22}_*.csv     | wc -l    # 11
ls -1 $F/db/chr{12,22}_*.csv.idx | wc -l    # 11
```

**Keep the `.idx` files.** They are `xsv` indexes, and they are the reason the annotation
step needs only about 53 MB of memory against a 23 GB table.

### 3. Build the annotated GDS file

This is done once per chromosome and is independent of the datasets. The annotated file
is small (a few MB) and is the only thing the association step needs afterwards -- the
large database is required at annotation time only, which makes the association step
portable to any machine.

```bash
# a) cut the chromosome out of the raw VCF: a plain region subset, no filtering
FIXTURE=20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr_12_22_X.recalibrated_variants.exome.vcf.gz
podman run --rm -v "$DATA/tools_comparison":/d:z docker.io/psuszynski/bioinf_combo:1.5.1 bash -lc "
  bcftools view -r chr22 -Oz -o /d/favor/work/fixture_chr22.vcf.gz /d/$FIXTURE &&
  bcftools index -t /d/favor/work/fixture_chr22.vcf.gz"

# b) preprocess exactly as FAVORannotator's own Scripts/UTL/preProcessingVCF.sh does.
#    This is required: the annotator fails on multi-allelic sites.
podman run --rm -v "$F/work":/w:z -v "$REFDIR":/ref:z,ro \
    docker.io/psuszynski/bioinf_combo:1.5.1 bash -lc "
  set -e
  bcftools annotate -x ^FORMAT/GT --threads 4 -Oz -o /w/c22.gt.vcf.gz       /w/fixture_chr22.vcf.gz
  bcftools norm -m -any          --threads 4 -Oz -o /w/c22.gt.bk.vcf.gz     /w/c22.gt.vcf.gz
  bcftools norm -f /ref/$REF     --threads 4 -Oz -o /w/c22.gt.bk.nm.vcf.gz  /w/c22.gt.bk.vcf.gz
  bcftools index -t /w/c22.gt.bk.nm.vcf.gz"

# c) VCF -> GDS
podman run --rm -v "$F/work":/work:z localhost/favorannotator:1.0.0 \
  Rscript /opt/favorannotator/convert_vcf_to_gds.R /work/c22.gt.bk.nm.vcf.gz /work/c22.agds 4

# d) add the annotations, in place
podman run --rm -v "$F/work":/work:z -v "$F/db":/favordb:z,ro localhost/favorannotator:1.0.0 \
  Rscript /opt/favorannotator/favorannotator_csv_essential.R \
          /work/c22.agds 22 /favordb /work/anno_out/
```

Notes on this sequence:

- Removing everything except the GT field in step (b) matters a great deal. Without it
  the GDS is 250 MB and takes 3m25s (the PL field alone is 133 MB); with it the same
  chromosome is 3 MB and takes 11 seconds.
- Step (b) also splits multi-allelic sites. This is not optional and not a substitution
  of our own: it is upstream's prescribed preprocessing. About 10% of raw chromosome-22
  records are multi-allelic, and the annotator writes its intermediate file without
  quoting, so a multi-allelic ALT such as `22-15528754-G-T,A` produces a malformed row
  and the tool fails with "found record with 2 fields, but the previous record has 1
  fields".
- `bcftools norm -f` reports `mismatch_removed: 0` when the reference is correct. A
  non-zero count means the wrong reference file.
- Step (d) **modifies the GDS in place**, which is why step (c) writes straight to the
  final `.agds` name. To repeat the annotation, re-run step (c) first; otherwise the
  annotation node is added to an already-annotated file.

## Running the comparison

`run_chain_ricopili_staar.sh` performs, for each dataset in turn: building a PLINK
dataset from the raw VCF, RICOPILI quality control (`preimp_dir`), RICOPILI principal
components (`pcaer`), GENESIS structure estimation, export back to VCF, and the STAAR
association test for both versions described in the next section. It writes one
REGENIE-shaped result table per dataset per version, plus a table recording how many
variants and samples survived each quality-control step.

```bash
ARM=$RVA_REPO/docs/tool-comparison/chained-benchmark

# quick check that everything is wired up: one dataset, one chromosome, 20 genes, ~7 minutes
DATASET_IDXS="18" CHRS="22" VERSIONS="filtered full" MAX_GENES=20 CLEANUP=false \
  bash "$ARM/run_chain_ricopili_staar.sh" 2>&1 | tee ~/smoke_run18.log
echo "EXIT=${PIPESTATUS[0]}"

# a real run
nohup env DATASET_IDXS="4 18 27" CHRS="12 22" VERSIONS="filtered full" \
      CLEANUP=false SCORE=false THREADS=4 \
      bash "$ARM/run_chain_ricopili_staar.sh" > ~/run.log 2>&1 &

tail -f ~/run.log | grep -E '^\[[A-G]\]|^\[run_|VERDICT|ERROR|eval table'
```

`MAX_GENES` **must be left unset for any run whose results will be scored** -- it exists
only for quick checks. `CHRS` defaults to `"12 22"`.

**Expected cost:** roughly 2 seconds per gene per version at this sample size, so one
fully scored dataset is about (427 genes on chromosome 22 + 1000 on chromosome 12) x 2
seconds x 2 versions, or about 95 minutes of association testing, plus about 7 minutes
for quality control, structure estimation and export.

## The two versions: removing relatives against modelling them

RICOPILI's `pcaer` removes one member of every pair of samples whose estimated
relatedness exceeds its threshold. On exome data this removes far more than intended:
**902 of 1,356 samples, leaving 454.**

The cause is a shortage of usable markers, not genuine relatedness. `pcaer` prunes
markers with a hardcoded minor-allele-frequency floor of 0.05, and only 4,073 of the
76,634 variants that pass quality control reach that frequency. Pruning those to
independence leaves about 1,220 markers, and relatedness estimated from so few markers
is very noisy, so tens of thousands of sample pairs cross the threshold. The 0.05 floor
is not a command-line option, so changing it would mean modifying RICOPILI's source.

Rather than modify RICOPILI, the comparison runs **two versions of the combined method
for every dataset** and reports both:

| Version | Samples entering the association test | Principal components | How relatedness is handled |
|---|---|---|---|
| **Filtered** | 454 unrelated | RICOPILI's smartpca (`menv.mds`) | relatives are removed -- RICOPILI as shipped |
| **Full** | all ~1,356 | GENESIS PC-AiR, covering all samples | relatives are kept, and a GENESIS PC-Relate relatedness matrix is included in the null model |

The Full version is not an unfair advantage. STAARpipeline is designed to be given a
relatedness matrix, and providing one is the single substitution this comparison makes,
because RICOPILI has no equivalent concept. PC-AiR is used instead of RICOPILI's own
projected components because it covers every sample (RICOPILI's projection leaves 97
without components) and because GENESIS is being used for the relatedness matrix anyway.
This also matches how nf-rare-var-assoc works: relatedness affects only which samples
principal components are computed from, while the association test keeps everyone.

Comparing Filtered against Full isolates exactly one thing: whether **removing**
relatives or **modelling** them works better.

The principal components and the relatedness matrix are computed on all available
markers (chromosomes 12, 22 and X), independently of which chromosomes the association
test covers. They are genome-wide covariates.

### The relatedness matrix must be sparse

STAAR's `fit_nullmodel` leaves `use_sparse` unset, which sends a dense matrix down
GMMAT's dense code path. STAAR's C++ code then rejects that path's projection matrix, and
because the call is wrapped in `try()`, **every gene silently returns nothing while the
run reports zero errors**. A sparse matrix is therefore the only working option, and it
is also what STAARpipeline's own tutorial uses.

Making it sparse is difficult on exome data, for the same reason as above. Measured on
one dataset (3,149 pruned markers, 918,690 sample pairs), thresholding the relatedness
graph at the conventional fourth-degree cutoff leaves it fully connected and still
dense. Only a second-degree cutoff breaks it into usable blocks:

| Relatedness cutoff | Pairs kept | Connected groups | Largest group | Usable |
|---|---|---|---|---|
| 0.0221 (4th degree, the conventional value) | 79,075 | 1 | 1,356 | no |
| 0.0442 (3rd degree) | 4,304 | 11 | 1,346 | no |
| **0.0884 (2nd degree)** | **571** | **275** | **12** | **yes** |
| 0.1768 (1st degree) | 544 | 817 | 12 | yes |

The second-degree result matches the data's real structure: 544 of the 571 retained
pairs are first-degree, consistent with the roughly 600 trios in this 1000 Genomes
release. Two scaling traps: `pcrelateToMatrix` applies its threshold **after** doubling
the kinship values, so a second-degree cutoff is passed as `2^-2.5` (0.1768), not
`2^-3.5`; and STAAR's own `kins_cutoff` default of 0.022 is on that same doubled scale,
which is why it keeps everything.

## Quality-control settings

RICOPILI's defaults were designed for genotyping arrays and are five to ten times
stricter than what nf-rare-var-assoc applies. On exome data with few cases they remove
so many variants that some datasets produce no association result at all. Four
thresholds are therefore set explicitly, using **RICOPILI's own documented command-line
options** -- nothing is patched:

| Threshold | RICOPILI default | Used here | Matching setting in nf-rare-var-assoc |
|---|---|---|---|
| `--geno` (variant missingness) | 0.02 | **0.10** | `plink2_makepgen_4_options --geno 0.100` |
| `--mind` (sample missingness) | 0.02 | **0.20** | `plink2_makepgen_5_options --mind 0.200` |
| `--pre_geno` (pre-filter) | 0.05 | **0.35** | `plink2_missing_per_pheno_options --geno 0.350` |
| `--midi` (differential missingness) | 0.02 | **0.35** | the same step, see below |

Everything else is left as shipped: Hardy-Weinberg thresholds, the inbreeding-coefficient
threshold, no minor-allele-frequency floor, and all of `pcaer`. Set
`PREIMP_GENO`, `PREIMP_MIND`, `PREIMP_PRE_GENO` and `PREIMP_MIDI` back to
`0.02/0.02/0.05/0.02` to reproduce a run with RICOPILI's defaults.

**`--pre_geno` is easy to get wrong, and raising `--geno` alone would have had no
effect.** RICOPILI computes the pre-filter list before any sample filtering and then
excludes those variants from every later step, so they never return. The effective
variant cut is therefore `min(pre_geno, geno)`. Left at 0.05 against `--geno 0.10`, the
entire 5-10% missingness band would still have been removed while the log reported 0.10.
RICOPILI's own defaults keep `pre_geno` looser than `geno` for exactly this reason, and
0.35 restores that ordering. Both run scripts abort before doing any work if
`PREIMP_PRE_GENO` is smaller than `PREIMP_GENO`.

**`--midi` is not an exact equivalent of anything in this pipeline.** RICOPILI limits the
*difference* in missingness between cases and controls; nf-rare-var-assoc splits samples
by phenotype, applies `--geno 0.350` within each group, and intersects the results, which
limits each group's *absolute* rate. A variant missing in 0% of controls and 34% of cases
passes ours and fails `--midi 0.02` by a factor of seventeen. 0.35 is the largest
difference our own filter can allow through, so it is the loosest setting still implied
by this pipeline rather than an invented one. The default is severe here for a reason
worth knowing: in the most imbalanced datasets there are about eleven cases, so **a
single missing genotype in one case is a 9% difference**, four and a half times over
RICOPILI's default threshold.

One difference cannot be removed by any setting: RICOPILI computes `pre_geno` across all
samples before any sample removal and `--geno` on the remaining set, in one pass over a
single dataset, whereas nf-rare-var-assoc computes missingness within each phenotype
group. This follows from RICOPILI producing one quality-controlled dataset where this
pipeline produces separate sets for the two REGENIE steps.

The settings actually used are echoed in each dataset's banner and recorded as
`preimp_qc_args` on the first line of `retention/run_<N>.tsv`, so results from different
settings identify themselves.

## Scoring the results and producing the figures

The association results are scored against the complete list of causal genes, including
genes on chromosome X and non-coding genes that STAAR cannot test -- these count as
misses. This is the honest whole-method number and it shares an axis with the other
comparisons.

Everything downstream of the association tables is cheap to regenerate (seconds to a
couple of minutes) and rescans whatever tables are on disk, so newly finished datasets
are picked up automatically and no list needs editing.

```bash
ARM=/data/git/doktorat_pw/wum_pims/nf-rare-var-assoc/docs/tool-comparison/chained-benchmark
COMMON=/data/git/doktorat_pw/wum_pims/nf-rare-var-assoc/docs/tool-comparison/benchmark-common
T=/data/doktorat/biodatageeks/article_on_nf_rare_var_assoc/tools_comparison
PY=/data/git/playground_all/python/.venv/bin/python   # needs pandas, scipy, matplotlib

# 1. score both versions
bash "$ARM/run_staar_eval.sh"
# -> runs/ricopili_staar_{filtered,full}_eval/results/compute_score/*_auc_summary.csv
#    (VERSIONS=full bash "$ARM/run_staar_eval.sh" for just one)

# 2. compare each version against this pipeline
for V in filtered full; do
  $PY "$COMMON/pairwise_compare.py" --runs "$T/runs" --missing drop \
    --arm-a nf_rare_var_assoc nf_rare_var_assoc_eval \
    --arm-b "ricopili_staar_$V" "ricopili_staar_${V}_eval" \
    --out "$T/runs/pairwise_ricopili_staar/vs_reference_$V"
done

# 3. draw the three-method bar chart
cd "$COMMON" && $PY three_method_ap_bar.py --runs "$T/runs"
```

`--missing drop` compares only the datasets a version actually produced; `--missing zero`
instead counts a missing result as a score of zero, so a failure counts against the
method that failed. Which is appropriate depends on what you are measuring.

The bar chart uses the Full version. To switch it to Filtered, or to change which
datasets appear, edit the `METHODS` list and `all_idx` in
`benchmark-common/three_method_ap_bar.py`.

After a refresh, check that the scoring log prints one `_auc_summary.csv` per dataset per
version, that `pairwise_compare.py` prints `coverage: ... produced N/30`, and that the
bar chart log prints `RICOPILI + STAAR ... datasets=N`. If a dataset is missing, check
that its `*_dataset_idx_<N>_step2_Y1.regenie` table exists under
`regenie_per_dataset/<version>/`.

### Regenerating `genes_info_hgnc.tsv`

The gene coordinates come from STAARpipeline's own internal data file, which is a plain
`.rda` and can be read without installing the package:

```r
e <- new.env(); load("<STAARpipeline>/R/sysdata.rda", envir = e)
g <- get("genes_info", envir = e)
g <- data.frame(hgnc_symbol=as.character(g[[1]]), chromosome_name=as.character(g[[2]]),
                start_position=as.integer(g[[3]]), end_position=as.integer(g[[4]]))
write.table(g, "genes_info_hgnc.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```

## What this comparison does and does not cover

**Indels are not tested.** FAVOR annotates single-nucleotide variants richly and indels
barely: across the 42,476 chromosome-22 variants, all of the continuous annotation
channels are populated for about 97% of single-nucleotide variants and for **0%** of
indels. The indel rows are present in the database and do join correctly -- they simply
carry no annotation scores. More decisively, STAAR reads its annotation weights only when
`variant_type` is `"SNV"`; setting it to `"variant"` to include indels makes STAAR apply
no annotation weights to anything, reducing it to plain unweighted tests. There is no
mode that weights single-nucleotide variants and also tests indels. Since annotation
weighting is STAAR's entire purpose, `"SNV"` is the only setting that represents the tool
fairly, and the cost is that indels are not tested. nf-rare-var-assoc tests indels
natively through its VEP consequence tiers.

**Only coding regions are tested.** The comparison runs gene-centric coding tests only,
for runtime reasons. In one examined dataset, 2 of the 10 causal variants on chromosome
22 are in 3' untranslated regions and cannot be reached by any coding category. Scoring
against the complete causal list already accounts for this honestly, but it is a real
limitation of the coding-only setting.

**Correcting for population structure is mandatory, not optional.** Both tools
independently measured a genomic inflation factor above 8 on this data with no structure
correction (STAAR 8.33 on genes, RICOPILI 8.479 on single variants). The inflation factor
compares observed test statistics against what would be expected by chance; a value near
1.0 indicates a well-calibrated test. Scoring an uncorrected run against a corrected one
would be grossly unfair. The inflation factor is therefore measured for every dataset in
every method, every time -- it is free and it is the cheapest evidence that the
comparison is sound.

**Sex information comes from the public pedigree, not from the genotypes.** PLINK2
refuses to import chromosome X without sex information, so the conversion cannot happen
at all without it. It is taken from `assets/integrated_call_samples_v3.20250704.ALL.ped`,
which ships with the same data release and is available to any user of this VCF. This
gives the combined method slightly more than nf-rare-var-assoc uses, since this pipeline
imputes sex from genotypes rather than reading it. The alternative -- dropping chromosome
X -- would have disabled RICOPILI's sex check, one of the capabilities being measured.

**RICOPILI's own quality report flags exome data as problematic** on several checks:
too few variants (its threshold assumes genome-wide array coverage), an extreme
case/control ratio (deliberate in these simulated datasets), and a 16.8% sex-check
mismatch rate. The last is worth understanding: 241 of 1,416 samples are flagged, almost
all of them women whose X-chromosome inbreeding estimate falls in the ambiguous range --
which is what an exome's sparse, rare-variant-heavy X coverage does to a statistic
designed for arrays. RICOPILI only warns: **2 samples were actually excluded**, the other
239 were merely reported.

**Most of what RICOPILI's variant filtering removes is not a quality judgement.** Of the
65,072 variants removed from one dataset, **60,290 (93%) were monomorphic** -- they have
no minor allele in that dataset's sample subset, which is an inevitable consequence of
selecting 1,416 samples out of 3,202. The genuine quality filters together remove about
2,750 variants (4%). RICOPILI applies no minor-allele-frequency floor to the association
data, so rare variants survive intact, as they must for a rare-variant test; the 0.05
floor is confined to `pcaer`'s pruning and never touches the data used for association.

## Practical notes

Things that are easy to get wrong and cost real time.

**Containers and builds**

- **Never read a build's exit code through a pipe.** `podman build ... | tail` reports
  `tail`'s status. Redirect to a log file, check `$?`, then read the log.
- **`R -e "install.packages(...)"` exits successfully even when the install fails.** The
  first FAVORannotator image built "successfully" with no gdsfmt in it. Every install
  step in the Dockerfile is now followed by a check that exits non-zero.
- **`podman build` in a detached background shell had no DNS on this machine**, while
  `podman run` and a foreground build both did. If a build fails on "Could not resolve
  host", try it in the foreground before suspecting the Dockerfile.
- The GENESIS image runs as its own user (uid 2049). Under `--userns=keep-id` that maps
  to an ID with no write access to your working directory, so that one container must be
  run as your own uid.

**LaTeX inside RICOPILI**

RICOPILI renders its quality report with pdflatex during the quality-control step, so a
LaTeX failure aborts the whole step after all the work is done.

- conda's texlive cannot build its own formats and must be shadowed by Debian's texlive.
- Shadowing `pdflatex` alone is not enough: `pcaer` assembles plots with `pdfjam`, which
  finds packages through `kpsewhich`. The whole Debian toolchain has to come first on
  PATH.
- `cm-super` is required, not optional -- without it pdflatex fails with "Font tcrm1000
  at 600 not found".

**PLINK and file conversion**

- PLINK2 will not import chromosome X without sex information.
- `--split-par` produces a file PLINK 1.9 refuses to read ("the .bim file has a split
  chromosome"). Running `plink --make-bed` once re-sorts it.
- Use `--new-id-max-allele-len 1000`, not `100 truncate`: truncation produces duplicate
  identifiers for long indels at the same position, and `pcaer`'s merge step then fails
  with "Duplicate ID". Long identifiers are safe -- `pcaer` shortens anything over 30
  characters itself.
- A PLINK bed file does not record which allele was the reference, and PLINK 1.9 makes
  the minor allele A1. Exporting back to VCF therefore needs `--keep-allele-order` on the
  re-sort and `--ref-allele` at export, taking the reference from the variant identifier.
  Without this, `bcftools norm -f` fails on a reference mismatch.
- Export with `id-paste=iid`. PLINK2 otherwise names samples `FID_IID`, and `preimp_dir`
  rewrites the family ID.
- `bcftools annotate -x ^FORMAT/GT` fails on a PLINK-exported VCF, because it is already
  genotype-only and bcftools errors when the operation would remove nothing.
- **STAAR keeps only variants whose FILTER column is exactly `PASS`**, and PLINK writes
  `.`. Without marking records as PASS after quality control, every gene tests zero
  variants silently.

**RICOPILI's own behaviour**

- **`preimp_dir` always runs twice.** The first call writes a template name file and
  exits asking you to edit it; the study name column must be filled in with five
  alphanumeric characters before the second call does any work.
- **A failed run leaves state that blocks the retry.** RICOPILI records progress and
  refuses to repeat a step that "has been done repeatedly without any progress". After a
  failure, start from a clean directory rather than re-running in place.
- The platform guesser finds nothing and prints "Error, something went wrong with plague"
  followed by "Warning: Unknown platform". This is harmless: our variant identifiers are
  `chr:pos:ref:alt`, so no genotyping array matches. The run continues normally with an
  empty platform suffix in the dataset name.
- EIGENSOFT is **not** in the 2025 dependency archive, although `pcaer` requires
  `smartpca`. It comes from conda instead, via `rp_env_0225b.yaml`, published in the same
  Broad directory as the dependency archive with its own checksum. Note that upstream's
  README never names this file; it was found by listing that directory.

**R and logging**

- **STAAR returns list-matrices, which `write.csv` cannot encode.** Its result tables are
  built with `rbind()` over mixed-type vectors, so `as.data.frame()` produces list
  columns and writing fails with "unimplemented type 'list' in 'EncodeElement'".
  `staar_gene_centric_coding.R` unpacks each column separately, preserving its type. Do
  not convert with `as.character` instead: that would round-trip p-values through text.
- `info.import=NULL` and `fmt.import=NULL` in SeqArray import **all** fields, not none.
  Use `character(0)` for none.
- **`tee` in `podman run ... | tee /work/x.log` runs on the host**, where `/work` usually
  does not exist. Put the whole pipeline inside the container instead.
- `grep` buffers by block when piped, so a long run's log looks empty for minutes. Use
  `grep --line-buffered`. Relatedly, `staar_gene_centric_coding.R` writes its CSVs only
  after the full gene loop finishes, so an empty output directory during a run is
  expected.
- **A silent empty result is the failure mode to watch for.** STAAR wraps its test call
  in `try()`, so a failure inside it produces no error and no results. The run script's
  check that p-values are non-degenerate is what catches this. Keep it.
