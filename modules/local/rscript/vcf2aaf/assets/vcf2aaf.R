# podman run --rm -v $PWD/:/work_dir/:Z -w /work_dir rocker/r-ver:4.4.2 Rscript /work_dir/vcf_to_aaf.R /work_dir/gnomad.exomes.v4.1.sites.chr1_head1k.vcf /work_dir AF_nfe AF

# required.packages <- c("data.table", "stringr", "dplyr", "R.utils")
# new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)


library("data.table")
library("stringr")
library("dplyr")
library("R.utils")


args <- commandArgs(trailingOnly = TRUE)

if (length(args)<=0) {
  stop("Provide the path to the vcf file", call.=FALSE)
}
if (length(args)<=1) {
  stop("Provide the path to the output aaf file", call.=FALSE)
}
if (length(args)<=2) {
  stop("Provide the name of the INFO tag to extract", call.=FALSE)
}
if (length(args)<=3) {
  stop("Provide the name of the INFO tag to extract when the normal tag name is absent", call.=FALSE)
}

vcf_path <- args[1]
out_aaf_path <- args[2]
tag_name <- args[3]  # tag_name = "AF_nfe"
default_tag_name <- args[4]  # default_tag_name = "AF"

vars <- fread(vcf_path, skip="#CHROM", sep="\t")

chrom_num_col <- stringr::str_extract(vars[["#CHROM"]], regex("(chr[0-9XY]+)", ignore_case = T), group = 1)
pos_col <- paste(chrom_num_col, vars$POS, vars$REF, vars$ALT, sep="_")

af_col <- stringr::str_extract(vars$INFO, regex(paste0(default_tag_name, "=([\\.0-9]+)"), ignore_case = T), group = 1)
af_nfe_col <- stringr::str_extract(vars$INFO, regex(paste0(tag_name, "=([\\.0-9]+)"), ignore_case = T), group = 1)

# should we do: coalesce(na_if(af_nfe_col, "0"), af_col, "0") ?
aafs <- data.table(pos=pos_col, af=coalesce(af_nfe_col, af_col, "0"))
fwrite(aafs, out_aaf_path, sep="\t", col.names=F, quote=F)
