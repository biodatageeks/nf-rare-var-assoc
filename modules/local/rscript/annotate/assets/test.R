required.packages <- c("data.table")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


library("data.table")

args <- commandArgs(trailingOnly = TRUE)

if (length(args)<=0) {
  stop("Provide the path to the input fam file", call.=FALSE)
}
if (length(args)<=1) {
  stop("Provide the path to the input controls file", call.=FALSE)
}
if (length(args)<=2) {
  stop("Provide the path to the input cases file", call.=FALSE)
}
if (length(args)<=3) {
  stop("Provide the path to the input vcf file", call.=FALSE)
}
if (length(args)<=4) {
  stop("Provide the path to the input sample file", call.=FALSE)
}
if (length(args)<=5) {
  stop("Provide the path to the output fam file", call.=FALSE)
}
if (length(args)<=6) {
  stop("Provide the path to the output sample file", call.=FALSE)
}
if (length(args)<=7) {
  stop("Provide the path to the output phenotype file", call.=FALSE)
}
if (length(args)<=8) {
  stop("Provide the path to the output annotations file", call.=FALSE)
}
if (length(args)<=9) {
  stop("Provide the path to the output masks file", call.=FALSE)
}
if (length(args)<=10) {
  stop("Provide the path to the output setlist file", call.=FALSE)
}

r_in_fam_path <- args[1]
r_in_controls_path = args[2]
r_in_cases_path = args[3]
r_in_vcf_path <- args[4]
r_in_sample_path <- args[5]

r_out_fam_path <- args[6]
r_out_sample_path <- args[7]
r_out_phenotype_path <- args[8]
r_out_annotations_path <- args[9]
r_out_masks_path <- args[10]
r_out_setlist_path <- args[11]


print(r_out_phenotype_path)


df1 <- data.table(mask="Mask1", csq="missense_variant")
fwrite(df1, r_out_fam_path, sep="\t", quote=F, col.names=T)
fwrite(df1, r_out_sample_path, sep="\t", quote=F, col.names=T)
fwrite(df1, r_out_phenotype_path, sep="\t", quote=F, col.names=T)
fwrite(df1, r_out_annotations_path, sep="\t", quote=F, col.names=T)
fwrite(df1, r_out_masks_path, sep="\t", quote=F, col.names=T)
fwrite(df1, r_out_setlist_path, sep="\t", quote=F, col.names=T)
