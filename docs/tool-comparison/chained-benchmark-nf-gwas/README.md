# Comparison against RICOPILI combined with nf-gwas

The second combined comparison, and the one where both halves are complete pipelines that
a user configures and runs. RICOPILI performs quality control and computes principal
components; [nf-gwas](https://github.com/genepi/nf-gwas) runs the gene-based association
test with REGENIE. Together they are treated as one method and compared against
nf-rare-var-assoc on the same simulated datasets.

Unlike the [comparison against RICOPILI combined with STAARpipeline](../chained-benchmark/),
no analysis script of ours sits between the two tools. One thing is still borrowed: the
gene groupings, which nf-gwas cannot build itself. That is the only input this method
receives from nf-rare-var-assoc besides the raw VCF, the phenotypes and the causal genes
used for scoring.

The scripts that are not specific to these tools are in
[../benchmark-common/](../benchmark-common/). **The setup instructions and practical notes
in the [RICOPILI + STAARpipeline README](../chained-benchmark/README.md) apply here
unchanged and are not repeated** -- in particular building the RICOPILI container, how
`preimp_dir` and `pcaer` behave, and the quality-control settings.

## Contents of this directory

`run_chain_ricopili_nf_gwas.sh` runs the whole method for one or more datasets. Its
first three stages are copied verbatim from the STAARpipeline comparison's script, so
that each comparison keeps the exact code that produced its own numbers.

| Stage | What happens | Which tool |
|---|---|---|
| A | raw VCF -> sample subset -> genotypes only -> split multi-allelic sites -> PLINK bed | RICOPILI's entry point |
| B | `preimp_dir`: sample and variant quality control | RICOPILI |
| C | `pcaer`: principal components, with relatives projected back onto them | RICOPILI |
| D | quality-controlled bed -> association VCF, prediction bed, principal-component covariates; then the identifier check | conversion |
| E | nf-gwas: REGENIE steps 1 and 2, gene-based, on the borrowed gene groupings | nf-gwas |
| F | rename the result table; check inflation factor and that p-values are not degenerate | this comparison |

Settings: `DATASET_IDXS`, `EXPORT_ONLY=true` (stop after stage D -- the cheap way to check
the identifier convention), `SCORE=true` (score at the end), `NPCS_COVAR` (default 4),
`MIN_ID_MATCH` (default 0.5). **The script deletes nothing.** It refuses to start in a
working directory that is not clean, and prints the sizes of the disposable directories
when a dataset finishes.

Relatives are kept and REGENIE models them, so unlike the STAARpipeline comparison there
is only one version, not a Filtered/Full pair.

## Converting RICOPILI's output for nf-gwas

Stage D is the only genuinely new part, and four details in it fail quietly if they are
got wrong.

### Normalising against the reference genome

The borrowed gene groupings key variants by `CHROM_POS_REF_ALT`, and nf-prepare-vcf builds
those keys after left-aligning and trimming indels against a reference FASTA. The
association VCF therefore has to be normalised against the same FASTA
(`GRCh38_full_analysis_set_plus_decoy_hla.fa`, `REFDIR`/`REF` in the script) before its
identifiers are assigned, or every indel key lands a few bases away from the mask entry
that names it. The FASTA is chr-prefixed, which is why the export writes `chr12` and the
prefix is stripped in the step after the normalisation rather than by plink2.

`bcftools norm` exits on a REF/ALT mismatch, so a wrong reference stops the run instead of
quietly dropping variants — provided the reference allele has been restored first.

### Restoring the reference allele

A PLINK bed file does not record which allele was the reference, and `preimp_dir`'s PLINK
1.9 steps make the minor allele A1. The variant identifiers assigned in stage A
(`chr:pos:REF:ALT`) carry the answer, so the reference allele is restored from them.
On one dataset this rotated 1,099 sets of allele codes.

### Bringing the pseudo-autosomal region back onto chromosome X, with `--merge-x`

The gene groupings name chromosome X variants plain `X`, so the pseudo-autosomal region
has to be merged back. `--merge-par` is the obvious flag and **does nothing here**:

```
Warning: --merge-par had no effect (no PAR1/PAR2 chromosome codes present).
```

Stage A splits the region with plink2 (codes `PAR1`/`PAR2`) but then re-sorts with PLINK
1.9, which knows only the older `XY` code. By the time the quality-controlled bed exists
the region sits under `XY`, and `--merge-par` finds nothing. `--merge-x` handles the `XY`
code, but it needs `--sort-vars`, which in turn works only with `--make-pgen`. That is why
the export is two plink2 passes rather than one.

Left unfixed this silently loses about 934 variants per dataset to a chromosome name that
matches nothing -- a warning, not an error, and invisible afterwards.

### Keeping chromosome X diploid in males

plink2 refuses to merge and export in one pass, and says why:

```
Warning: --merge-x should not be used in the same run as VCF export; this causes some
ploidies to be wrong.  Instead, use --merge-x + --sort-vars + --make-[b]pgen in one run,
and follow up with --split-par + --export vcf.
```

Following that advice literally produces a correct file, with male non-pseudo-autosomal X
haploid. **But that is not what this comparison should produce.** The input VCF codes male
chromosome X as diploid, as the 1000 Genomes data does, so nf-rare-var-assoc and the
earlier nf-gwas comparison both test chromosome X on diploid male dosages. Exporting
haploid would halve the male contribution to every chromosome-X gene test and mix a pure
encoding difference into a comparison that is meant to be about pipelines.

The second pass therefore exports with `--update-sex` marking every sample female, which
reproduces the input's encoding exactly. Nothing else uses that file: RICOPILI's quality
control already ran on the true sex, and the true sex is kept in the prediction genotypes
that REGENIE step 1 reads.

```bash
# restore REF from our own variant identifiers
awk -F'\t' '{n=split($2,a,":"); if (n==4 && a[3]!="") print $2"\t"a[3]}' qc.bim > ref_alleles.txt

# pass 1: pseudo-autosomal region back onto X, reference restored, sorted
plink2 --bfile <qc> --merge-x --sort-vars --ref-allele force ref_alleles.txt 2 1 \
       --make-pgen --out xmerged

# pass 2: diploid X, chromosome names matching the reference FASTA, sample name = IID
plink2 --pfile xmerged --update-sex sex_diploid.txt --output-chr chrM \
       --export vcf bgz id-paste=iid --out assoc.export

# normalise against the reference, strip the chr prefix, then assign identifiers
bcftools norm --fasta-ref <ref.fa> -m -any --rm-dup exact -Ou assoc.export.vcf.gz \
  | bcftools annotate --rename-chrs chr_map.txt -Ou \
  | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o assoc.vcf.gz

# prediction genotypes: same bed, family ID rewritten to the sample id, TRUE sex kept
plink2 --bfile <qc> --update-ids update_ids.txt --output-chr MT --make-bed --out prediction
```

### nf-gwas's own plink2 is five years old, and this comparison depends on it

A current plink2 (2.00a5.10) refuses the association VCF outright:

```
Error: chrX is present in the input file, but no sex information was provided;
rerun this import with --psam or --update-sex.
```

nf-gwas's import step runs `plink2 --vcf <file> dosage=DS --double-id` with no sex
information anywhere. It works only because the `quay.io/genepi/nf-gwas:v1.0.9` image
ships PLINK v2.00a2.3 (January 2020), which predates that check. Worth knowing before
anyone upgrades the image.

### The identifier check

Before running nf-gwas, stage D verifies that the exported variants actually match the
borrowed gene groupings. On one dataset:

| Check | Value |
|---|---|
| variants exported | 76,634 |
| named in the borrowed annotations | 51,135 (**0.667**); the threshold is 0.5 |
| by chromosome | 12: 0.905, 22: 0.748, **X: 0.052** |
| grouping variants surviving RICOPILI's quality control | 69,225 / 127,963 (**0.541**) |
| genes with at least one surviving variant | 2,473 / 2,706 |
| genes with at least two surviving variants | 2,179 / 2,706 |

The chromosome-X fraction is low for a reason unrelated to identifiers: the borrowed
groupings cover coding regions only and name just 2,064 chromosome-X variants against
81,881 on chromosome 12, while RICOPILI's output carries 17,136 chromosome-X variants. The
identifiers themselves are correct -- 896 of them match. Because the overall fraction is
dominated by the autosomes and would pass the threshold even if chromosome X were named
wrongly, the check **also requires every exported chromosome to match something**.

The last three rows are worth recording per dataset, and the script writes them to
`retention/run_<N>.tsv`: RICOPILI's quality control costs the average gene 46% of its
named variants, and 233 genes (8.6%) lose every variant and become untestable.

## REGENIE step 1 and the common-variant supply

nf-gwas always runs its own quality control on the prediction genotypes -- there is no way
to switch it off -- so this method applies RICOPILI's filtering and nf-gwas's filtering in
sequence. That is exactly what combining the two tools gives a user, so it is left as it
is.

The thresholds are a direct copy of nf-rare-var-assoc's own tuned settings for this data,
so the step-1 quality-control layer is identical across all the compared methods by
construction:

| nf-rare-var-assoc setting | nf-gwas setting |
|---|---|
| `--geno 0.250` | `qc_geno 0.25` |
| `--hwe 1e-9 0.01` | `qc_hwe 1e-9` |
| `--mac 16` | `qc_mac 16` |
| `--maf 0.045` | `qc_maf 0.045` |
| `--mind 0.150` | `qc_mind 0.15` |

These come from the benchmark's own configuration, **not** from `nextflow.config`'s
defaults (`--mac 70 --maf 0.01`), which this benchmark never used.

The minor-allele-frequency floor removes about 94% of variants, which looks alarming and
is not. It applies only to the prediction genotypes that REGENIE step 1 fits on; the
association test still sees all 76,634 variants through the gene groupings. A
three-chromosome exome simply contains very little common variation:

| Minor-allele-frequency cutoff | Variants |
|---|---|
| >= 0.045 (what is applied) | 4,263 |
| >= 0.02 | 6,273 |
| >= 0.01 | 8,933 |
| >= 0.005 | 12,533 |

RICOPILI removes 46% of variants overall but barely touches this pool, because most of
what it removes is monomorphic. All the compared methods fit step 1 on the same order of
magnitude of markers (about 4,000), so none is starved relative to the others.

## Quality-control settings

`run_chain_ricopili_nf_gwas.sh` runs `preimp_dir` with `--geno 0.10 --mind 0.20
--pre_geno 0.35 --midi 0.35` rather than RICOPILI's defaults, matching
nf-rare-var-assoc's thresholds. These are `preimp_dir`'s own command-line options; nothing
is patched, and everything else is left as shipped.

**The full reasoning -- including why raising `--geno` alone has no effect, and why
`--midi` is not an exact equivalent of anything in this pipeline -- is in one place: the
quality-control section of [`../chained-benchmark/README.md`](../chained-benchmark/README.md).**
Both combined comparisons must run RICOPILI at identical thresholds or their
quality-control halves stop being the same pipeline, which is why the settings block is
duplicated verbatim in both scripts rather than shared.

## Running it

```bash
export DATA=/data/runs/tools_comparison
export RVA_REPO=/data/git/nf-rare-var-assoc
export NFGWAS_REPO=/data/git/nf-gwas
ARM=$RVA_REPO/docs/tool-comparison/chained-benchmark-nf-gwas
```

### Prerequisites

Beyond the RICOPILI container and bcftools image from the other comparison:

```bash
git clone https://github.com/genepi/nf-gwas.git "$NFGWAS_REPO"
git -C "$NFGWAS_REPO" checkout 6a5de44          # v1.0.9-7-g6a5de44
podman pull quay.io/genepi/nf-gwas:v1.0.9
nextflow -v || curl -s https://get.nextflow.io | bash   # needs Java 17+
```

The borrowed gene groupings (about 227 MB) must be present. They are the
`*.annotations` and `*.setlist` files produced by a run of nf-rare-var-assoc with
`--publish_intermediate true`, under
`runs/nf_rare_var_assoc/results/bcftools/`.

Disk: about 60 GB for 30 datasets (each dataset's working directory is roughly 1.8 GB,
plus about 0.6 GB of results). Check `nproc` and `free -g` too: nf-gwas requests 8 CPUs
and 16 GB per REGENIE process, so on a smaller machine add `NFGWAS_MAX_CPUS=4
NFGWAS_MAX_MEM=12.GB` or Nextflow will abort.

### A single dataset first

```bash
RUN_DIR=$DATA/runs/ricopili_nf_gwas_smoke DATASET_IDXS=18 THREADS=4 \
  bash "$ARM/run_chain_ricopili_nf_gwas.sh" 2>&1 | tee ~/smoke_run18.log
echo "EXIT=${PIPESTATUS[0]}"
```

Expect 141,706 -> 76,634 variants, 1,416 -> 1,356 samples, 76,634 exported with 51,135
named in the groupings (0.667), 1,320 genes tested, and an inflation factor near 1.19.
Overriding `RUN_DIR` keeps the real run's directory clean, since the script refuses to
start in a directory that is not empty.

### The full run

**Launch with `setsid`, not a bare `nohup ... &`.** See the note about SIGTTIN below; the
alternative stops the whole run within seconds and takes an hour to notice.

```bash
setsid nohup env DATA=$DATA RVA_REPO=$RVA_REPO NFGWAS_REPO=$NFGWAS_REPO \
          THREADS=4 SCORE=false \
     bash "$ARM/run_chain_ricopili_nf_gwas.sh" < /dev/null > ~/nfgwas_all.log 2>&1 &

tail -f ~/nfgwas_all.log | grep -E '^ chain run|^\[[A-F]\]|named in masks|lambda|eval table|ERROR'
```

Budget about 14 minutes per dataset on a modern workstation.

**Failures are tolerated and the run resumes.** Each dataset runs in its own subshell, so
one dataset's failure is recorded and skipped rather than fatal. A dataset whose result
table already exists is skipped untouched, so **re-running the identical command
resumes**: finished datasets are not redone and nothing is deleted. At the end the script
prints which datasets produced no table; scoring simply omits them.

Finished datasets leave about 1.8 GB each that nothing downstream needs. The script never
removes anything, so delete them yourself when space is short:

```bash
rm -rf $DATA/runs/ricopili_nf_gwas/{work,nfwork}/run_<N>
```

## Scoring and figures

Scoring is wired into the script's `SCORE=true` path, which scores every dataset and then
compares this method against nf-rare-var-assoc and against the earlier nf-gwas
comparison. Both steps can also be run by hand on tables that already exist:

```bash
COMMON=/data/git/doktorat_pw/wum_pims/nf-rare-var-assoc/docs/tool-comparison/benchmark-common
T=/data/doktorat/biodatageeks/article_on_nf_rare_var_assoc/tools_comparison

EVAL_REPO=/data/git/doktorat_pw/wum_pims/nf-eval-gene-assoc \
EVAL_RUN_DIR="$T/runs/ricopili_nf_gwas_eval" \
EVAL_PROJECT=ricopili_nf_gwas \
INPUT_VCF="$T/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr_12_22_X.recalibrated_variants.exome.vcf.gz" \
SKIP_PREP=true \
REGENIE_GLOB="$T/runs/ricopili_nf_gwas/regenie_per_dataset/ricopili_nf_gwas_dataset_idx_*_step2_Y1.regenie" \
CAUSAL_SNPLIST_GLOB="$T/datasets/run_*/select_genes/*_dataset_idx_*_in_*.snplist" \
CAUSAL_GENES_GLOB="$T/datasets/run_*/select_genes/*_genes_dataset_idx_*.txt" \
  bash "$COMMON/run_eval.sh"

PY=/data/git/playground_all/python/.venv/bin/python   # needs pandas, scipy, matplotlib
$PY "$COMMON/pairwise_compare.py" --runs "$T/runs" --missing zero \
  --arm-a nf_rare_var_assoc nf_rare_var_assoc_eval \
  --arm-b ricopili_nf_gwas  ricopili_nf_gwas_eval \
  --out "$T/runs/pairwise_ricopili_nf_gwas/vs_reference"
```

Each output directory gets `pairwise_per_dataset.csv`, `pairwise_summary.csv` and the
figures described in [../benchmark-common/](../benchmark-common/).

**Choosing `--missing`.** A dataset that one method could not analyse at all needs a
decision. `--missing zero` gives the method that produced nothing a score of zero -- which
is exactly what the recall-scaling formula gives for an empty result -- and keeps the
dataset in the comparison. `--missing drop` excludes it entirely. Dropping such a dataset
*rewards* a method for failing to run, so `--missing zero` is the honest default and
`--missing drop` is the lenient sensitivity check. Either way the tool prints how many
datasets each method actually produced a result for, which is itself worth reporting.

If you are deciding how many datasets to generate, run `project_power.py` on an existing
`pairwise_per_dataset.csv` first. When the difference between two methods is small,
statistical power grows very slowly: at the effect size measured here (Cohen's d_z around
0.25), reaching 80% power would need roughly 130 datasets, and going from 30 to 60 raises
the chance of reaching p<0.05 only to about 46%.

## How the inflation factor is computed

`lambda_<test>` in each `retention/run_<N>.tsv` is `median(CHISQ) / 0.4549364`, where
CHISQ is REGENIE's one-degree-of-freedom statistic for that test (`ADD`, `ADD-SKAT`,
`ADD-SKATO`) and 0.4549364 is the median of the null chi-square distribution with one
degree of freedom. It is the standard measure, but computed here over roughly 1,300 gene
tests rather than millions of variants, so it is noisy and, for rare-variant burden tests,
structurally conservative. **Read it as a sanity check, not as a calibrated statement
about false discovery.**

Values track the number of cases almost monotonically, from about 1.2-1.3 in datasets with
several hundred cases down to 0.5-0.95 in datasets with a dozen. The low values are a
power artefact, not a problem: with few cases the burden statistic is discrete and
conservative, which pulls the value below 1. It is not corrected, because it is a true
property of those datasets.

What matters is the upper end. Without any structure correction both tools measured
inflation above 8 on this data; with RICOPILI's principal components the maximum here is
about 1.29 for burden tests, so the structure correction is working in every dataset.

Note that REGENIE writes chromosome X as `23` in its output. The scoring step handles
this.

## Datasets that legitimately cannot be analysed

At least one simulated dataset cannot be analysed by any burden method. The datasets share
a fixed control pool of about 1,333 samples and draw a variable number of cases; one draw
has 11 cases, of which quality control removes enough to fall below REGENIE's minimum of
10, and REGENIE stops with *"all phenotypes have less than 10 cases"*. The thinnest
dataset that still runs has 12 cases. This is a property of the data, not a failure of
either tool, and the failure-tolerant loop is what lets a full run carry past it
unattended.

## Practical note: `nohup ... &` alone stops the run within seconds

Symptom: the log stops right after nf-gwas prints its parameter banner, no container ever
starts, the CPU is idle, and nothing reports an error. `ps` gives it away -- **every
process in the job is in state `T` (stopped)** -- and two extra children are present:

```
16590  sh -c stty -icanon min 1 -icrnl -inlcr < /dev/tty
16591  stty -icanon min 1 -icrnl -inlcr
```

Nextflow's terminal handling runs `stty` against `/dev/tty` explicitly, which bypasses the
`/dev/null` standard input that `nohup` provides. A background process group that reads
the controlling terminal is sent SIGTTIN, and the kernel stops the entire group -- shell
script, Java virtual machine and all. `kill -CONT` does not fix it: the read is retried
and the signal fires again.

`setsid` fixes it permanently: with no controlling terminal, opening `/dev/tty` simply
fails and Nextflow carries on. This matters more here than usual, because this script
launches Nextflow once per dataset and so has one chance to freeze per dataset.

To recover, kill the process group, remove the half-finished dataset's directories, and
relaunch under `setsid`:

```bash
kill -CONT -<PGID>; kill -TERM -<PGID>          # PGID from `ps -o pid,pgid,stat,args`
rm -rf $RUN_DIR/work/run_<N> $RUN_DIR/nfwork/run_<N>
```

## What is written where

None of this is committed to the repository.

| Path | What | Keep? |
|---|---|---|
| `<DATA>/runs/ricopili_nf_gwas/work/run_<N>/` | quality control, principal components, exported genotypes, conversion intermediates | disposable once the result table exists |
| `<DATA>/runs/ricopili_nf_gwas/nfwork/run_<N>/` | Nextflow working directory | disposable |
| `<DATA>/runs/ricopili_nf_gwas/results/run_<N>/` | nf-gwas output, including logs, trace and report | keep -- the logs carry the quality-control and step-1 numbers |
| `<DATA>/runs/ricopili_nf_gwas/regenie_per_dataset/` | the result tables that get scored | keep |
| `<DATA>/runs/ricopili_nf_gwas/retention/run_<N>.tsv` | per-step retention, identifier check, gene attrition, inflation factor | keep |
| `<DATA>/runs/ricopili_nf_gwas_eval/` | scores | keep |
| `<DATA>/runs/pairwise_ricopili_nf_gwas/` | comparison tables and figures | keep |

Containers used, none built here: `localhost/ricopili:2025_Feb_20.001`,
`docker.io/psuszynski/bioinf_combo:1.5.1`, `quay.io/genepi/nf-gwas:v1.0.9`.
