
library("data.table")

args <- commandArgs(trailingOnly = TRUE)

if (length(args)<=0) {
  stop("Provide the path to the regenie step2_masks.snplist input file", call.=FALSE)
}
if (length(args)<=1) {
  stop("Provide the path to the regenie step2_Y1.regenie input file", call.=FALSE)
}
if (length(args)<=2) {
  stop("Provide the path to the vcf input file", call.=FALSE)
}
if (length(args)<=3) {
  stop("Provide the path to the phenotypes input file (produced by annotate.R)", call.=FALSE)
}
if (length(args)<=4) {
  stop("Provide the path to the annotations input file (produced by annotate.R)", call.=FALSE)
}
if (length(args)<=5) {
  stop("Provide the path to the annotated_snps.tsv output file", call.=FALSE)
}
if (length(args)<=6) {
  stop("Provide the path to the res_log10p_1_annotated.tsv output file", call.=FALSE)
}
if (length(args)<=7) {
  stop("Provide the path to the annotated_snps_with_sample_ids.tsv output file", call.=FALSE)
}

r_in_regenie_step2_masks_snplist_path <- args[1]
r_in_regenie_step2_Y1_regenie_path <- args[2]
r_in_vcf_path <- args[3]
r_in_phenotype_path <- args[4]
r_in_annotations_path <- args[5]

r_out_annotated_snps_tsv_path <- args[6]
r_out_res_log10p_1_annotated_tsv_path <- args[7]
r_out_annotated_snps_with_sample_ids_tsv_path <- args[8]


vars <- fread(vcf_path, skip="#CHROM", sep="\t")

masks <- fread(r_in_regenie_step2_masks_snplist_path, header=F)
umaskvars <- unique(unlist(strsplit(masks$V4,",")))

pheno <- fread(r_in_phenotype_path)
anno <- fread(r_in_annotations_path, header=F)
sel <- anno[which(anno$V1 %in% umaskvars),]
print(paste0("nrow(sel) = ", nrow(sel)))

cat("Creating annotated_snps.tsv ")
resList <- lapply( (1:nrow(sel)), function(i){
        if (i %% (nrow(sel) %/% 100) == 0) {
            cat(".")
        }
        selV <- sel$V1[i]
        vars1 <- t(vars[which(vars$ID == selV),10:ncol(vars)])
        samples <- rownames(vars1)
        ref <- samples[grep("0/0:", vars1)]
        htz <- samples[grep("0/1:", vars1)]
        hmz <- samples[grep("1/1:", vars1)]
        nas <- setdiff(samples, union(union(ref, htz), hmz))
        pheno[match(pheno$FID, samples),]
        cases <- pheno$FID[which(pheno$Y1 == 1)]
        controls <- pheno$FID[which(pheno$Y1 == 0)]
        cases_htz <- intersect(cases, htz)
        cases_hmz <- intersect(cases, hmz)
        cases_nas <- intersect(cases, nas)
        controls_nas <- intersect(controls, nas)
        controls_htz <- intersect(controls, htz)
        controls_hmz <- intersect(controls, hmz)
        cases_ac <-  (length(cases_htz) + 2*(length(cases_hmz)))
        cases_af <-  cases_ac/ (2*(length(cases)-length(cases_nas)))
        controls_ac <-  (length(controls_htz) + 2*(length(controls_hmz)))
        controls_af <-  controls_ac/  (2*(length(controls)-length(controls_nas)))
        data.table(cbind(selV, cases_af=cases_af, controls_af=controls_af,cases_ac=cases_ac, controls_ac=controls_ac, cases_na=length(cases_nas), controls_na=length(controls_nas)))
})

anno_snps <- rbindlist(resList)
fwrite(anno_snps, r_out_annotated_snps_tsv_path)
cat(" - done!\n")

# anno_snps <- fread(r_out_annotated_snps_tsv_path)
dd <- fread(r_in_regenie_step2_Y1_regenie_path)
dd <- dd[order(dd$LOG10P,decreasing=T),]
sel_dd <- dd[which(dd$LOG10P > 1.5),]
# cbind(sel_dd[i,],anno_snps[which(anno_snps$selV %in% strsplit(masks$V4[which(masks$V1 == sel_dd$ID[i])],",")[[1]]),])
anno <- fread(r_in_annotations_path, header=F)
# cbind(anno_snps, anno[anno$V1])
anno_snps2 <- cbind(anno[match( anno_snps$selV, anno$V1),], anno_snps)


cat("Creating res_log10p_1_annotated.tsv ..")
sel_dd <- dd[which(dd$LOG10P > 1),]

ll <- lapply(1:nrow(sel_dd), function(i){
        anno_snps2$selV %in% strsplit(masks$V4[which(masks$V1 == sel_dd$ID[i])],",")
        cbind(sel_dd[i,],anno_snps2[which(anno_snps2$selV %in% strsplit(masks$V4[which(masks$V1 == sel_dd$ID[i])],",")[[1]]),])
})
final <- rbindlist(ll)

fwrite(final, r_out_res_log10p_1_annotated_tsv_path)
cat(" - done!\n")

# sel <- anno[which(anno$V1 %in% umaskvars),]
# sel <- sel[which(sel$V2 == "KIAA1549"),]
# print(paste0("for KIAA1549  nrow(sel) = ", nrow(sel)))
#
# cat("Creating annotated_snps_with_sample_ids.tsv (1) ")
# resList <- lapply( (1:nrow(sel)), function(i){
#         selV <- sel$V1[i]
#         vars1 <- t(vars[which(vars$ID == selV),10:ncol(vars)])
#         samples <- rownames(vars1)
#         ref <- samples[grep("0/0:", vars1)]
#         htz <- samples[grep("0/1:", vars1)]
#         hmz <- samples[grep("1/1:", vars1)]
#         nas <- setdiff(samples, union(union(ref, htz), hmz))
#         pheno[match(pheno$FID, samples),]
#         cases <- pheno$FID[which(pheno$Y1 == 1)]
#         controls <- pheno$FID[which(pheno$Y1 == 0)]
#         cases_htz <- intersect(cases, htz)
#         cases_hmz <- intersect(cases, hmz)
#         cases_nas <- intersect(cases, nas)
#         controls_nas <- intersect(controls, nas)
#         controls_htz <- intersect(controls, htz)
#         controls_hmz <- intersect(controls, hmz)
#         cases_ac <-  (length(cases_htz) + 2*(length(cases_hmz)))
#         cases_af <-  cases_ac/ (2*(length(cases)-length(cases_nas)))
#         controls_ac <-  (length(controls_htz) + 2*(length(controls_hmz)))
#         controls_af <-  controls_ac/  (2*(length(controls)-length(controls_nas)))
#
#         data.table(cbind(selV, cases_af=cases_af, controls_af=controls_af,cases_ac=cases_ac, controls_ac=controls_ac, cases_na=length(cases_nas),cases_htz=paste0(cases_htz,collapse=";"),cases_hmz=paste0(cases_hmz,collapse=";"), controls_na=length(controls_nas)))
# })
#
# res <- rbindlist(resList)
# cat(" - done!\n")
#
# sel <- sel[which(sel$V2 == "LRIG1"),]
# sel <- anno[which(anno$V1 %in% umaskvars),]
# sel <- sel[which(sel$V2 == "LRIG1"),]
# print(paste0("for LRIG1  nrow(sel) = ", nrow(sel)))
#
# cat("Creating annotated_snps_with_sample_ids.tsv (2) ")
# resList <- lapply( (1:nrow(sel)), function(i){
#         selV <- sel$V1[i]
#         vars1 <- t(vars[which(vars$ID == selV),10:ncol(vars)])
#         samples <- rownames(vars1)
#         ref <- samples[grep("0/0:", vars1)]
#         htz <- samples[grep("0/1:", vars1)]
#         hmz <- samples[grep("1/1:", vars1)]
#         nas <- setdiff(samples, union(union(ref, htz), hmz))
#         pheno[match(pheno$FID, samples),]
#         cases <- pheno$FID[which(pheno$Y1 == 1)]
#         controls <- pheno$FID[which(pheno$Y1 == 0)]
#         cases_htz <- intersect(cases, htz)
#         cases_hmz <- intersect(cases, hmz)
#         cases_nas <- intersect(cases, nas)
#         controls_nas <- intersect(controls, nas)
#         controls_htz <- intersect(controls, htz)
#         controls_hmz <- intersect(controls, hmz)
#         cases_ac <-  (length(cases_htz) + 2*(length(cases_hmz)))
#         cases_af <-  cases_ac/ (2*(length(cases)-length(cases_nas)))
#         controls_ac <-  (length(controls_htz) + 2*(length(controls_hmz)))
#         controls_af <-  controls_ac/  (2*(length(controls)-length(controls_nas)))
#
#         data.table(cbind(selV, cases_af=cases_af, controls_af=controls_af,cases_ac=cases_ac, controls_ac=controls_ac, cases_na=length(cases_nas),cases_htz=paste0(cases_htz,collapse=";"),cases_hmz=paste0(cases_hmz,collapse=";"), controls_na=length(controls_nas)))
# })
#
# res <- rbindlist(resList)
# cat(" - done!\n")
#
# selV <- sel$V1[12]
# vars1 <- t(vars[which(vars$ID == selV),10:ncol(vars)])

sel <- anno[which(anno$V1 %in% umaskvars),]
print(paste0("for all umaskvars  nrow(sel) = ", nrow(sel)))


cat("Creating annotated_snps_with_sample_ids.tsv (3) ")
resList <- lapply( (1:nrow(sel)), function(i){
        if (i %% (nrow(sel) %/% 100) == 0) {
            cat(".")
        }
        selV <- sel$V1[i]
        vars1 <- t(vars[which(vars$ID == selV),10:ncol(vars)])
        samples <- rownames(vars1)
        ref <- samples[grep("0/0:", vars1)]
        htz <- samples[grep("0/1:", vars1)]
        hmz <- samples[grep("1/1:", vars1)]
        nas <- setdiff(samples, union(union(ref, htz), hmz))
        pheno[match(pheno$FID, samples),]
        cases <- pheno$FID[which(pheno$Y1 == 1)]
        controls <- pheno$FID[which(pheno$Y1 == 0)]
        cases_htz <- intersect(cases, htz)
        cases_hmz <- intersect(cases, hmz)
        cases_nas <- intersect(cases, nas)
        controls_nas <- intersect(controls, nas)
        controls_htz <- intersect(controls, htz)
        controls_hmz <- intersect(controls, hmz)
        cases_ac <-  (length(cases_htz) + 2*(length(cases_hmz)))
        cases_af <-  cases_ac/ (2*(length(cases)-length(cases_nas)))
        controls_ac <-  (length(controls_htz) + 2*(length(controls_hmz)))
        controls_af <-  controls_ac/  (2*(length(controls)-length(controls_nas)))

        data.table(cbind(selV, cases_af=cases_af, controls_af=controls_af,cases_ac=cases_ac, controls_ac=controls_ac, cases_na=length(cases_nas),cases_htz=paste0(cases_htz,collapse=";"),cases_hmz=paste0(cases_hmz,collapse=";"), controls_na=length(controls_nas)))
})

res <- rbindlist(resList)


fwrite(res, r_out_annotated_snps_with_sample_ids_tsv_path)
cat(" - done!\n")
