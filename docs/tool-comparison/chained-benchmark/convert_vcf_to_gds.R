#############################################################################
# convert_vcf_to_gds.R -- VCF -> SeqArray GDS, the first step of STAAR's own recipe
#
# Derived from zhouhufeng/FAVORannotator @ d65f565 Scripts/CSV/convertVCFtoGDS.r.
# Upstream sources a config.R holding one hardcoded vcf/gds path pair per
# chromosome; this takes them as arguments instead. seqVCF2GDS options are kept
# byte-for-byte identical to upstream (genotype.var.name="GT",
# info.import=NULL, fmt.import=NULL, ignore.chr.prefix="chr") so the aGDS we
# build matches what FAVORannotator and STAARpipeline expect.
#
# NOTE on the NULL arguments: in SeqArray, info.import=NULL / fmt.import=NULL
# import ALL INFO/FORMAT fields -- character(0) is what imports none. So the
# GDS carries the full FORMAT payload (AB/AD/DP/GQ/PL/PGT/PID); on the chr22
# fixture, PL alone is 133 MB of the 250 MB file. STAAR itself reads only hard
# calls, so this is wasted space rather than a correctness problem, and we keep
# upstream's behaviour rather than diverge from it.
#
# Usage: Rscript convert_vcf_to_gds.R <in.vcf[.gz]> <out.gds> [threads]
#############################################################################
suppressPackageStartupMessages({
    library(gdsfmt)
    library(SeqArray)
})

args <- commandArgs(TRUE)
if (length(args) < 2) {
    stop("usage: convert_vcf_to_gds.R <in.vcf[.gz]> <out.gds> [threads]")
}
vcf.fn <- args[1]
gds.fn <- args[2]
threads <- if (length(args) >= 3) as.integer(args[3]) else 1L

cat("vcf    :", vcf.fn, "\n")
cat("gds    :", gds.fn, "\n")
cat("threads:", threads, "\n")

start_time <- Sys.time()

seqVCF2GDS(vcf.fn, gds.fn,
           header = NULL,
           genotype.var.name = "GT",
           info.import = NULL,
           fmt.import = NULL,
           ignore.chr.prefix = "chr",
           raise.error = TRUE,
           parallel = threads,
           verbose = TRUE)

genofile <- seqOpen(gds.fn, readonly = TRUE)
cat("GDS built\n")
cat("variants:", length(seqGetData(genofile, "variant.id")), "\n")
cat("samples :", length(seqGetData(genofile, "sample.id")), "\n")
print(genofile)
seqClose(genofile)

cat("elapsed:", format(Sys.time() - start_time), "\n")
