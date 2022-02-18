rm(list=ls())

library(UpSetR)
library(RColorBrewer)
library(easyPubMed)


#####  Paths  #####

path_base <- "/home/ros259/R_projects/450K_Batch_effects"
path_data <- file.path(path_base, "data")
path_save_rdata <- file.path(path_base, "data", "rdata")
path_results <- file.path(path_base, "results")
path_idats <- file.path(path_base, "data", "idats")
path_code <- file.path(path_base, "scripts", "R")
path_snps <- "/home/ros259/R_projects/utility_scripts"

#####  Functions  #####

myplot <- function(b, title="", col="black") {
    plot(b, cex=0.5, cex.main=0.5,
         ylim=c(0,1), xaxt="n",
         las=2, cex.axis=0.6,
         xlab="", ylab="",
         col=col,
         main=title)
}


myplot2 <- function(p, this, correction, title="", cex=0.5, ...) {
    
    myorder <- order(betas[[this]]$pd$array_num, paste(betas[[this]]$pd$row, betas[[this]]$pd$col))
    
    plot(betas[[this]][[correction]][p, myorder], cex=cex, cex.main=0.7,
         ylim=c(0,1), xaxt="n",
         las=2, cex.axis=0.6,
         xlab="", ylab="",
         col=betas[[this]]$pd$batch_col[myorder], main=title, ...) 
}


publication_plot <- function(myprobes, these = prefix_alts) {
    
    toplot <- function(p, this, correction, cex=0.5, title="") {
        
        myorder <- order(betas[[this]]$pd$array_num, paste(betas[[this]]$pd$row, betas[[this]]$pd$col))
        
        plot(betas[[this]][[correction]][p, myorder], cex=cex, cex.main=0.7,
             ylim=c(0,1), xaxt="n",
             las=2, cex.axis=0.6,
             xlab="", ylab="",
             col=betas[[this]]$pd$batch_col[myorder], main=title) 
    }
    
    myiqrs <- list()
    mymodal <- list()
    
    for(prefix in these) {
        myiqrs[[prefix]] <- betas[[prefix]]$iqrs[betas[[prefix]]$iqrs$is_biggest == TRUE, ]
        mymodal[[prefix]] <- betas[[prefix]]$modal
    }

    for(p in myprobes) {
        #if(!p %in% rownames(betas$`450K`)) next
        
        for(prefix in these) {
            cex=0.5
            
            if(prefix %in% c("EPIC-Italy", "NOVI")) cex=0.3
            idx <- match(p, myiqrs[[prefix]]$probe)
            toplot(p, prefix, correction="original", cex=cex,
                   title=paste(p, ": ", prefix, ",",
                               " SD=", round(var(betas[[prefix]]$original[idx, ]), 3), sep=""))
            toplot(p, prefix, correction="harman", cex=cex,
                   title=paste("Harman,",
                               " SD=", round(var(betas[[prefix]]$harman[idx, ]), 3),
                               ", LVR=", round(log2(myiqrs[[prefix]]$var_ratio_harman[idx]), 2),
                               ", Shift=", round(betas[[prefix]]$meandiffs_harman[idx], 4), sep=""))  # \U0394x\U0304
            toplot(p, prefix, correction="combat", cex=cex,
                   title=paste("ComBat,",
                               " SD=", round(var(betas[[prefix]]$harman[idx, ]), 3),
                               ", LVR=", round(log2(myiqrs[[prefix]]$var_ratio_combat[idx]), 2),
                               ", Shift=", round(betas[[prefix]]$meandiffs_combat[idx], 4), sep=""))
        }
    }
}


probe_cutoffs <- function(prefixes, ps) {
  
  vr_common <- list()
  is_likely_bad <- list()
  is_likely_bio <- list()
  is_surely_bad <- list()
  
  for(prefix in prefixes) {
    stopifnot(sum(ps %in% rownames(var_ratios[[prefix]])) == length(ps))
    vr_common[[prefix]] <- var_ratios[[prefix]][ps, ]
    is_likely_bad[[prefix]] <- vr_common[[prefix]]$var_ratio_combat < log2(1 / 1.5) & vr_common[[prefix]]$meandiffs_combat > 0.01
    is_surely_bad[[prefix]] <- vr_common[[prefix]]$var_ratio_combat < log2(1 / 1.5) & vr_common[[prefix]]$meandiffs_combat > 0.02
    is_likely_bio[[prefix]] <- vr_common[[prefix]]$var_ratio_combat > log2(1.5) & vr_common[[prefix]]$meandiffs_combat > 0.01
  }
  
  to_return <- list("likely_bad"=do.call("cbind", is_likely_bad),
                    "surely_bad"=do.call("cbind", is_surely_bad),
                    "likely_bio"=do.call("cbind", is_likely_bio))
  
  for(i in 1:length(to_return)) {
    rownames(to_return[[i]]) <- ps
  }
  to_return
}


common_probes_figure <- function(x) {
  
  m <- matrix(as.integer(x),
              nrow=nrow(x),
              ncol=ncol(x),
              dimnames = dimnames(x))
  upset(as.data.frame(m), text.scale = 0.8)
}


#####  Load files  #####

prefix_alts <- c("EpiSCOPE", "EPIC-Italy", "BodyFatness", "NOVI", "URECA")
names(prefix_alts) <- c("GSE89278", "GSE51032", "EPIC", "GSE128821", "GSE132181")

load(file=file.path(path_save_rdata, "EpiSCOPE_and_BodyFatness_noob_beta_diffs.RData"))
all_betas <- betas
rm(betas)

for(prefix in c("EPIC-Italy", "NOVI", "URECA")) {
    load(file=file.path(path_save_rdata, paste(prefix, "noob_beta_diffs.RData", sep="_")))
    all_betas[[prefix]] <- betas[[prefix]]
    rm(betas)
}

betas <- all_betas
rm(all_betas)
gc()

var_ratios = list()
for(prefix in prefix_alts) {
    
    is_biggest = betas[[prefix]]$iqrs$is_biggest
    var_ratios[[prefix]] = betas[[prefix]]$iqrs[is_biggest, c("var_ratio_combat", "var_ratio_harman")]
    rownames(var_ratios[[prefix]]) = betas[[prefix]]$iqrs$probe[is_biggest]
    rm(is_biggest)
    var_ratios[[prefix]] = log2(var_ratios[[prefix]])
    var_ratios[[prefix]]$meandiffs_combat = betas[[prefix]]$meandiffs_combat[rownames(var_ratios[[prefix]])]
    var_ratios[[prefix]]$meandiffs_harman = betas[[prefix]]$meandiffs_harman[rownames(var_ratios[[prefix]])]
    var_ratios[[prefix]]$maxdiffs_combat = betas[[prefix]]$maxdiffs_combat[rownames(var_ratios[[prefix]])]
    var_ratios[[prefix]]$maxdiffs_harman = betas[[prefix]]$maxdiffs_harman[rownames(var_ratios[[prefix]])]
}


# Italy colour fix
batch_cols <- rep(brewer.pal(10, 'Paired'), 15)[1:length(levels(betas$`EPIC-Italy`$pd$harman_batch))]
betas$`EPIC-Italy`$pd$batch_col <- batch_cols[betas$`EPIC-Italy`$pd$harman_batch]

# EPIC row, col fix
betas$BodyFatness$pd$row <- sub('\\d+_R(\\d+)C\\d+', '\\1', as.character(betas$BodyFatness$pd$Basename))
betas$BodyFatness$pd$col <- sub('\\d+_R\\d+C(\\d+)', '\\1', as.character(betas$BodyFatness$pd$Basename))


#####  Save lists for publication  #####

modal_pub_cols <- c("probe_id", "chr", "start", "end", "strand", "num_opt_clusters", "distanceSnp", "distanceRsid",
  "xhyb_47", "xhyb_48", "xhyb_49", "xhyb_50", "hwe_p", "sex_p", "batch_p", "superbatch_p", "cells_p")
  

for(prefix in prefix_alts) {
    to_save <- betas[[prefix]]$modal[, modal_pub_cols]
    to_save$chr <- factor(to_save$chr,
                          levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                                   "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"), ordered = TRUE)
    to_save <- to_save[order(to_save$chr, to_save$start, to_save$strand), ]
    p_cols <- grep("_p$", colnames(to_save))
    to_save[, p_cols] <- signif(to_save[, p_cols], 4)
    to_save <- cbind(to_save, round(var_ratios[[prefix]][rownames(to_save), ], 4))
    write.csv(to_save, file = file.path(path_results, paste(prefix, "modal_probes.csv", sep="_")), row.names = FALSE)
}
rm(to_save, modal_pub_cols)


# Modal stats

sink(file = file.path(path_results, "Modal_probe_stats.txt"))

for(prefix in prefix_alts) {
  cat(rep("*", 30), "\n    ", prefix, "\n", rep("*", 30), "\n", sep="")
  
  these_cols <- colnames(betas[[prefix]]$modal[grep("i[sn]_", colnames(betas[[prefix]]$modal))])
  for(mycol in these_cols) {
    cat(rep("-", 20), "\n", sep="")
    cat(prefix, mycol, "\n")
    mytable <- table(betas[[prefix]]$modal[, mycol])
    print(mytable)
    print(prop.table(mytable))
  }
  
  cat(rep("-", 20), "\n", sep="")
  x=table(sex=betas[[prefix]]$modal$is_imprinted,
          "chrX"=betas[[prefix]]$modal$is_chrX,
          "chrY"=betas[[prefix]]$modal$chr=="chrY")
  print(x)
  print(prop.table(x))
  
  cat(rep("-", 20), "\n", sep="")
  x=table(cell=betas[[prefix]]$modal$is_celltype,
          "nclust"=betas[[prefix]]$modal$num_opt_clusters)
  print(x)
  print(prop.table(x))
  
  cat(rep("-", 20), "\n", sep="")
  x=table(sex=betas[[prefix]]$modal$is_imprinted,
          in_hwe=!is.na(betas[[prefix]]$modal$hwe_p) & betas[[prefix]]$modal$hwe_p > 0.001)
  print(x)
  print(prop.table(x))
  
  cat(rep("-", 20), "\n", sep="")
  x=table(batch=betas[[prefix]]$modal$is_batchy,
          superbatch=betas[[prefix]]$modal$is_superbatchy)
  print(x)
  print(prop.table(x))
}
sink()
rm(these_cols)

# LVR/mean diff lists

cols_to_capture <- c("var_ratio_combat", "var_ratio_harman", "meandiffs_combat", "meandiffs_harman")


all_probes <- unique(unlist(lapply(var_ratios, rownames)))
ref_colnames <- unlist(lapply(names(var_ratios), function(x) paste(x, cols_to_capture, sep="_")))

ref <- matrix(rep(NA_real_, length(ref_colnames) * length(all_probes)), ncol=length(ref_colnames), nrow=length(all_probes),
              dimnames = list("probe_id"=all_probes,
                              "dataset"=ref_colnames)
              )

for(prefix in prefix_alts) {
  
  this_m <- var_ratios[[prefix]][, cols_to_capture]
  my_matches <- match(rownames(this_m), rownames(ref))
  stopifnot(sum(is.na(my_matches)) == 0)
  for(this in cols_to_capture) {
    ref[my_matches, paste(prefix, this, sep="_")] <- this_m[, this]
  }
}
rm(this_m, my_matches)

ref <- round(ref, 4)

lvr.combat <- ref[, grepl("var_ratio_combat", colnames(ref))]
md.combat <- ref[, grepl("meandiffs_combat", colnames(ref))]
lvr.harman <- ref[, grepl("var_ratio_harman", colnames(ref))]
md.harman <- ref[, grepl("meandiffs_harman", colnames(ref))]

write.csv(ref, file = file.path(path_results, "Infinium_reference.csv"))
save(lvr.combat, md.combat, lvr.harman, md.harman,
     file = file.path(path_results, "Infinium5.rdata"), compress = TRUE, compression_level = 9)
rm(ref, all_probes, ref_colnames, cols_to_capture, lvr.combat, md.combat, lvr.harman, md.harman)


######  Upset plots  ######

for(prefix in prefix_alts) {
  print(prefix)
  pdf(file=file.path(path_results, paste(prefix, "Figure_upset.pdf", sep="_")), width=7, height = 7)
  modal_res = data.frame("Crosshyb"=as.integer(betas[[prefix]]$modal$is_xhyb),
                         "Batch"=as.integer(betas[[prefix]]$modal$is_batchy),
                         "Super"=as.integer(betas[[prefix]]$modal$is_superbatchy),
                         "Gender"=as.integer(betas[[prefix]]$modal$is_imprinted),
                         "SNP"=as.integer(betas[[prefix]]$modal$is_snp),
                         "SNP_10bp"=as.integer(betas[[prefix]]$modal$is_snp_10bp),
                         "ChrX"=as.integer(betas[[prefix]]$modal$is_chrX),
                         "ChrY"=as.integer(betas[[prefix]]$modal$is_chrY),
                         "Cell"=as.integer(betas[[prefix]]$modal$is_celltype),
                         "HWE"=as.integer(betas[[prefix]]$modal$in_hwe))
  modal_res$Unknown = as.integer(rowSums(modal_res, na.rm=TRUE) == 0)
  print(upset(modal_res, nsets=11, nintersects = 24, text.scale = 0.9, order.by = "freq"))
  dev.off()
}


#####  Batchy probes plot  ####

bio_probes <- c("cg25465065", "cg15544633", "cg00455876", "cg15410402")
batchy_probes <- c("cg01381374", "cg22256960", "cg27298252", "cg04294190")  # "cg17758324", "cg18346402", "cg13701509", "cg26188290", "cg17765025")
dodgy_ewas <- c("cg11963436", "cg18368637", "cg22385669")


pdf(file.path(path_results, paste("Figure06_bio_correction_probes_figure.pdf", sep="_")))
par(mfrow=c(8, 3),  mar=c(0.1, 1.8, 1.0, 0.1))
publication_plot(bio_probes, these = c("EpiSCOPE", "BodyFatness"))
dev.off()


pdf(file.path(path_results, paste("Figure07_technical_correction_probes_figure.pdf", sep="_")))
par(mfrow=c(8, 3),  mar=c(0.1, 1.8, 1.0, 0.1))
publication_plot(batchy_probes, these = c("EpiSCOPE", "BodyFatness"))
dev.off()


##### TIFF plot  ####

for(prefix in prefix_alts) {
  
  tiff(file.path(path_results, paste(prefix, "Probe_correction_Variance.tiff", sep="_")), compression = "lzw", width = 960, height = 960)  # , res=300
  par(mfrow=c(3,2), mar=c(2, 3, 2, 2) + 0.1)
  
  myvr = var_ratios[[prefix]]
  # sanity check  
  stopifnot(identical(names(betas[[prefix]]$meandiffs_combat), rownames(myvr)))
  stopifnot(identical(names(betas[[prefix]]$cluster_props), rownames(myvr)))
  
  num_clust = sapply(betas[[prefix]]$cluster_props, length)
  
  these_probes <- list("All" = names(num_clust),
                       "Modal" = names(num_clust)[num_clust > 1],
                       "Imprinted" = betas[[prefix]]$modal$probe_id[betas[[prefix]]$modal$is_imprinted],
                       "SNP" = betas[[prefix]]$modal$probe_id[betas[[prefix]]$modal$is_snp],
                       "Cell" = betas[[prefix]]$modal$probe_id[betas[[prefix]]$modal$is_celltype],
                       "Batch" = betas[[prefix]]$modal$probe_id[betas[[prefix]]$modal$is_batchy | betas[[prefix]]$modal$is_superbatchy])
  
  probe_cols <- c("All"="#80808020", # grey
                  "Modal"="#7570b330", # purple
                  "Imprinted" = "#66a61e30", # light green
                  "SNP" = "#1b9e7730", # green
                  "Cell" = "#d95f0230", # orange
                  "Batch"="#e7298a30") # crimson 
  
  myxlim = range(myvr$var_ratio_combat)
  myylim = range(betas[[prefix]]$meandiffs_combat)

  for(i in 1:length(these_probes)) {
    
    is_p <- rownames(myvr) %in% these_probes[[i]]
    set_name <- names(these_probes)[i]
    
    my_xaxt="s"
    my_yaxt="s"
    if(i %in% 1:4) {
      my_xaxt="n"
    }
    if(i %in% c(2, 4, 6)) {
      my_yaxt="n"
    }
        
    plot(myvr$var_ratio_combat[!is_p],
         betas[[prefix]]$meandiffs_combat[!is_p],
         xlim=myxlim,
         ylim=myylim,
         xlab="",
         ylab="",
         xaxt= my_xaxt,
         yaxt= my_yaxt,
         cex=0.1,
         cex.axis=2.5,
         col="#99999920",
         main="")
    abline(v=c(log2(1/1.5), log2(1.5)), h=0.01, col="black")
    points(myvr$var_ratio_combat[is_p],
           betas[[prefix]]$meandiffs_combat[is_p],
           cex=0.1, col=probe_cols[set_name])
  }
  dev.off()
}

rm(myxlim, myylim, myvr, num_clust, these_probes, is_p, set_name)  # is_over_cutoff 


#####  Commonality  #####

for(prefix in prefix_alts) {
  cat(prefix, ncol(betas[[prefix]]$original), "\n")
}


# Collect the full set of modal probes associated with cell composition
cell_ps <- character()
for(prefix in prefix_alts) {
  these_cell_ps <- betas[[prefix]]$modal$probe_id[betas[[prefix]]$modal$is_celltype]
  cell_ps <- union(cell_ps, these_cell_ps)
}


#####  bad probes  #####

ps_common_csiro <- intersect(rownames(var_ratios$BodyFatness), rownames(var_ratios$EpiSCOPE))
ps_to_remove <- intersect(ps_common_csiro, cell_ps)
ps_common_csiro_notcell <- setdiff(ps_common_csiro, ps_to_remove)
rm(ps_to_remove)

ps_common_epic <- intersect(rownames(var_ratios$BodyFatness), intersect(rownames(var_ratios$NOVI), rownames(var_ratios$URECA)))
ps_to_remove <- intersect(ps_common_epic, cell_ps)
ps_common_epic_notcell <- setdiff(ps_common_epic, ps_to_remove)
rm(ps_to_remove)

stopifnot(identical(rownames(var_ratios$EpiSCOPE), rownames(var_ratios$`EPIC-Italy`)))
ps_450K <- rownames(var_ratios$EpiSCOPE)
ps_to_remove <- intersect(ps_450K, cell_ps)
ps_450K_notcell <- setdiff(ps_450K, ps_to_remove)
rm(ps_to_remove)

ps_common <- intersect(rownames(var_ratios$EpiSCOPE), intersect(rownames(var_ratios$BodyFatness), intersect(rownames(var_ratios$NOVI), rownames(var_ratios$URECA))))
ps_to_remove <- intersect(ps_common, cell_ps)
ps_common_notcell <- setdiff(ps_common, ps_to_remove)
rm(ps_to_remove)

is_cutoff_csiro <- probe_cutoffs(prefixes=c("EpiSCOPE", "BodyFatness"), ps_common_csiro_notcell)
is_cutoff_epic <- probe_cutoffs(prefixes=c("BodyFatness", "NOVI", "URECA"), ps_common_epic_notcell)
is_cutoff_450K <- probe_cutoffs(prefixes=c("EpiSCOPE", "EPIC-Italy"), ps_450K_notcell)
is_cutoff <- probe_cutoffs(prefixes=prefix_alts, ps_common_notcell)

snp_probes = betas[[prefix]]$modal$probe_id[betas[[prefix]]$modal$is_snp]
snp10_probes = betas[[prefix]]$modal$probe_id[betas[[prefix]]$modal$is_snp_10bp]
over_cutoff_probes = rownames(var_ratios[[prefix]])[var_ratios[[prefix]]$var_ratio_combat > log2(1.5) & var_ratios[[prefix]]$meandiffs_combat > 0.01]

# Upset plots

for(this in names(is_cutoff)) {
  pdf(file=file.path(path_results, paste(this, "overlap_upset.pdf", sep="_")), width=7, height = 5)
  print(common_probes_figure(is_cutoff[[this]]))
  dev.off()
}

#####  Probe overlaps with other sets  #####

# https://rdrr.io/github/markgene/maxprobes/src/R/cross-reactive.R
# Benton, https://github.com/sirselim

# BOWTIE Multi-mappers
mm_ps <- scan(file=file.path(path_data, "HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt"),
              what="character")

# Cross-reactive probes
xrct_probes <- read.csv(file=file.path(path_data, "48639-non-specific-probes-Illumina450k.csv"),
                        stringsAsFactors=FALSE)

tms <- list()
load(file=file.path(path_save_rdata, "450K_probe_tm_info.RData"))
tms$'450K' <- probes
rm(probes)
load(file=file.path(path_save_rdata, "EPIC_probe_tm_info.RData"))
tms$EPIC <- probes
rm(probes)

# 29,233 probes from Chen
# 33,457 probes from Benton

# Bose ICC set
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-312

#icc <- read.csv(file=file.path(path_data, "ICC_values.csv.gz"),
#                stringsAsFactors=FALSE)
#table(icc$ilmnid.CpG.site. %in% common_probes)
#table(common_probes %in% icc$ilmnid.CpG.site.)
#icc_probes <- icc$ilmnid.CpG.site.[icc$ICC.value < 0.1]

# TypeI worse with nCpG <=2
# TypeII worse with nCpG == 0

common_probes <- rownames(is_cutoff$likely_bad)

tms_trim <- lapply(tms, function(x) x[rownames(x) %in% common_probes, ])

sets <- list("All"=common_probes,
             "Chen"=xrct_probes$TargetID[xrct_probes$TargetID %in% common_probes],
             "Benton"=mm_ps[mm_ps %in% common_probes],
             "TypeI"=rownames(tms_trim$'450K')[tms_trim$'450K'$Type == "typeI"],
             "TypeI_nCpG_2orless"=rownames(tms_trim$'450K')[tms_trim$'450K'$Type == "typeI" & tms_trim$'450K'$nCpG <= 2],
             "TypeII"=rownames(tms_trim$'450K')[tms_trim$'450K'$Type == "typeII"],
             "TypeII_nCpG_is0"=rownames(tms_trim$'450K')[tms_trim$'450K'$Type == "typeII" & tms_trim$'450K'$nCpG == 0],
             "Tm_less_70"=rownames(tms_trim$'450K')[tms_trim$'450K'$tm_avg < 70])
             #"nCpG"=rownames(tms_trim$'450K')[tms_trim$'450K'$nCpG])

sources <- data.frame("Set"=names(sets),
                      "n"=sapply(sets, length),
                      "bad_0"=NA_real_,
                      "bad_1"=NA_real_,
                      "bad_2"=NA_real_,
                      "bad_3"=NA_real_,
                      "bad_4"=NA_real_,
                      "bad_5"=NA_real_,
                      "pval"=NA_real_,
                      stringsAsFactors = FALSE)
rownames(sources) <- names(sets)
bad_cols <- c("bad_0", "bad_1", "bad_2","bad_3", "bad_4", "bad_5")

for(s in names(sets)) {
  
  if(s == "All") {
    sources_prop <- sources
    sources[s, bad_cols] <- as.vector(table(rowSums(is_cutoff$likely_bad)))
    sources_prop[s, bad_cols] <- round(as.vector(prop.table(table(rowSums(is_cutoff$likely_bad)))) * 100, 1)   # percentage rounded to 3 dp
  } else {
    x <- table(rowSums(is_cutoff$likely_bad), sets$All %in% sets[[s]])
    res = chisq.test(x=rowSums(is_cutoff$likely_bad), y=sets$All %in% sets[[s]])
    sources[s, bad_cols] <- as.vector(res$observed[, 2])
    sources_prop[s, bad_cols] <- round(as.vector(prop.table(res$observed, margin=2)[, 2]) * 100, 1)  # percentage rounded to 3 dp
    sources[s, c("pval")] <- res$p.value
    sources_prop[s, c("pval")] <- res$p.value
  }
}

sources_table <- sources
sources_table[, c("n", bad_cols)] <- format(sources_table[, c("n", bad_cols)], big.mark=",")

for(bad_col in bad_cols) {
  sources_table[, bad_col] <- paste(sources_table[, bad_col], " (", sources_prop[, bad_col], ")", sep="")
}
sources_table$pval[!is.na(sources_table$pval) & sources_table$pval < 2.2e-16] <- "<2.2e-16"
write.csv(sources_table, file = file.path(path_results, "Probewise_factors_Table.csv"))


#####  Melting Temp  #####

tiff(file=file.path(path_results, "Figure_Tms.tiff"), height=170, width=170, units="mm", res=600, compression="lzw")
par(mfrow=c(5,2), mar=c(3, 4, 1, 2) + 0.1)

for(prefix in prefix_alts) {
  
  epic_set <- c("BodyFatness", "NOVI", "URECA")
  if(prefix %in% epic_set) {
    this_tm <- "EPIC"
  } else {
    this_tm <- "450K"
  }
  
  these_probes <- rownames(is_cutoff$likely_bad)[is_cutoff$likely_bad[, prefix]]
  
  my_tm <- tms[[this_tm]][these_probes, ]
  my_md <- var_ratios[[prefix]][these_probes, "meandiffs_combat"]
  is_typeI <- my_tm$Type == "typeI"
  is_typeII <- my_tm$Type == "typeII"
  
  myxlab = ""
  if(prefix == "URECA") myxlab = "CpG probe melting temperature (\u00B0C)"
  if(prefix == "BodyFatness") prefix <- "BFiN"
  
  smoothScatter(my_tm$tm_avg[is_typeI], my_md[is_typeI],
                nbin=512,
                xlab=myxlab,
                ylab="",
                #ylab="Correction (Mean delta Beta)",
                colramp = colorRampPalette(brewer.pal(9, "Blues")),
                cex.lab=0.7,
                cex.axis=0.8)
  mtext(paste(prefix, "- I"), side=3, cex=0.8)
  
  smoothScatter(my_tm$tm_avg[is_typeII], my_md[is_typeII],
                nbin=512,
                xlab=myxlab,
                ylab="",
                colramp = colorRampPalette(brewer.pal(9, "Reds")),
                cex.lab=0.7,
                cex.axis=0.8)
  mtext(paste(prefix, "- II"), side=3, cex=0.8)
}
dev.off()

