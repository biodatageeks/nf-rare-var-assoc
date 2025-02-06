# required.packages <- c("data.table", "R.utils")
# new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)


library("data.table")
library("R.utils")

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

filter_annotations_threshold <- 50
if (length(args) > 11) {
  filter_annotations_threshold <- as.numeric(args[12])
}
cat(paste0("filter_annotations_threshold = ", filter_annotations_threshold, "\n"))


tmp_dir = "./"

study_name <- "pims"

fam <- fread(r_in_fam_path, header=F)
controls <- fread(r_in_controls_path, header=F)
cases <- fread(r_in_cases_path, header=F)
fam$V6[which(gsub("\\.hg38", "",fam$V1) %in% controls$V1)] <- 0
fam$V6[which(gsub("\\.hg38", "",fam$V1) %in% cases$V1)] <- 1
fwrite(fam, r_out_fam_path, sep="\t", quote=F, col.names=F)
pheno <- fam[,c(1:2, 6)]
colnames(pheno) <- c("FID", "IID", "Y1")
fwrite(pheno, r_out_phenotype_path, sep="\t", quote=F, col.names=T)


dd <- fread(r_in_vcf_path, skip="#CHROM")
dd[, c("A", "CSQ") := tstrsplit(INFO, "CSQ=", fixed=TRUE)]
dd2 <- dd[,"CSQ",with=F]
csq_file <- paste0(tmp_dir, "csq.tsv")
out_csq_file <- paste0(tmp_dir, "csq_split.tsv")
fwrite(dd2, csq_file, sep="\t", col.names=F, quote=F)
cmd <- paste0("cat ", csq_file, " | awk '{split($0,a,\"|\"); print a[2],a[3],a[4]}' > ", out_csq_file)
system(cmd)
dd3 <- fread(out_csq_file, fill=T, header=F)
setnames(dd3, c("Consequence", "Impact", "Symbol"))
dd[,Consequence:=dd3$Consequence]
dd[,Impact:=dd3$Impact]
dd[,Symbol:=dd3$Symbol]

cn <- colnames(dd)
cn[1] <- "CHROM"
setnames(dd, cn)
dd[,key := paste0(CHROM, '_', POS, '_', REF,  '_',ALT )]

dd_final <- dd[,c("key","Symbol", "Consequence"),with=F]
anno_file <- r_out_annotations_path
fwrite(dd_final, anno_file, sep="\t", col.names=F, quote=F )

masks_file <- r_out_masks_path
# masks <- data.table(mask="Mask1", csq="missense_variant")
masks <- data.table(mask=c("Mask.High", "Mask.Mod", "Mask.HighMod"), csq=c(
    "stop_gained,stop_lost,start_lost,stop_gained&splice_region_variant,stop_gained&frameshift_variant,stop_gained&NMD_transcript_variant,frameshift_variant,frameshift_variant&NMD_transcript_variant,frameshift_variant&splice_region_variant,splice_donor_variant,splice_acceptor_variant,splice_donor_variant&NMD_transcript_variant,splice_donor_variant&splice_donor_region_variant&intron_variant,splice_acceptor_variant&non_coding_transcript_variant,splice_donor_variant&non_coding_transcript_variant",
    "missense_variant,missense_variant&splice_region_variant,missense_variant&splice_region_variant&NMD_transcript_variant,missense_variant&NMD_transcript_variant,inframe_insertion,inframe_deletion,inframe_deletion&splice_region_variant",
    "stop_gained,stop_lost,start_lost,stop_gained&splice_region_variant,stop_gained&frameshift_variant,stop_gained&NMD_transcript_variant,frameshift_variant,frameshift_variant&NMD_transcript_variant,frameshift_variant&splice_region_variant,splice_donor_variant,splice_acceptor_variant,splice_donor_variant&NMD_transcript_variant,splice_donor_variant&splice_donor_region_variant&intron_variant,splice_acceptor_variant&non_coding_transcript_variant,splice_donor_variant&non_coding_transcript_variant,missense_variant,missense_variant&splice_region_variant,missense_variant&splice_region_variant&NMD_transcript_variant,missense_variant&NMD_transcript_variant,inframe_insertion,inframe_deletion,inframe_deletion&splice_region_variant"
))
fwrite(masks , masks_file, sep="\t", col.names=F, quote=F )


unique_genes <- sort(unique(dd_final$Symbol))

library(parallel)
system.time(res <- lapply(unique_genes, function(gene){
    print(gene)
    data.table(symbol=gene,
               chrom=dd$CHROM[dd$Symbol==gene][1],
               pos=dd$POS[dd$Symbol==gene][1],
               variants= paste(dd$key[dd$Symbol==gene],collapse=","))

}))

setlist <- rbindlist(res)
set_list_file <- r_out_setlist_path
fwrite(setlist , set_list_file, sep="\t", col.names=F, quote=F )



# cmd <- paste0("cp ", plink_out_dir, study_name, ".bim", " ", tmp_dir, study_name, "-with-chrM.bim")
# system(cmd)
cmd <- paste0("cp ", r_out_setlist_path, " ", tmp_dir, study_name, ".setlist-with-chrM")
system(cmd)
cmd <- paste0("cp ", r_out_annotations_path, " ", tmp_dir, study_name, ".annotations_with_semicolon")
system(cmd)



# cmd <- paste0("grep -v chrM ", tmp_dir, study_name, "-with-chrM.bim", " > ", tmp_dir, study_name, ".bim")
# system(cmd)
# cmd <- paste0("grep -v chrM ", tmp_dir, study_name, "-with-chrM.setlist", " | grep -v chrY | sed -r 's/:/_/g' | sed 1d > ", r_out_setlist_path)
cmd <- paste0("grep -v chrM ", tmp_dir, study_name, ".setlist-with-chrM", " | grep -v chrY | sed -r 's/:/_/g' > ", r_out_setlist_path)
system(cmd)
cmd <- paste0("grep -v chrM ", tmp_dir, study_name, ".annotations_with_semicolon", " | sed -r 's/:/_/g' > ", r_out_annotations_path)
system(cmd)


anno <- fread(r_out_annotations_path, header=F)
cat(paste0("read ", nrow(anno), " rows from ", r_out_annotations_path, "\n"))
cat("sort(table(anno$V3)):\n")
print(sort(table(anno$V3)))

anno$V3[anno$V3 %in% names(head(sort(table(anno$V3)), filter_annotations_threshold))] <- "NULL"
fwrite(anno, r_out_annotations_path, sep="\t", col.names=F, quote=F )




# need to set missing to 0 in *.sample file
dd <- fread(r_in_sample_path)
dd$missing <- 0
fwrite(dd, file=r_out_sample_path, sep="\t", quote=F)

# remove multiallelic
dd <- fread(r_out_annotations_path, header=F)
cat(paste0("read ", nrow(dd), " rows from ", r_out_annotations_path, "\n"))

dd2 <- dd[!grepl(",", dd$V1),]
cat(paste0("writing ", nrow(dd2), " rows to ", r_out_annotations_path, "\n"))
fwrite(dd2, r_out_annotations_path, sep="\t", col.names=F, quote=F)
