rm(list=ls())

library(minfi)


#####  Paths  #####

path_base <- "/home/ros259/R_projects/450K_Batch_effects"
path_data <- file.path(path_base, "data")
path_save_rdata <- file.path(path_base, "data", "rdata")
path_results <- file.path(path_base, "results")
path_peapod <- "/home/ros259/R_projects/SGA Infants"
path_rgset <- file.path(path_peapod, "data", "rdata")


#####  Load RData  #####

load(file=file.path(path_save_rdata, "GSE89278_RGSet.RData"))
rg <- list(r_450K=getRed(RGSet),
           g_450K=getGreen(RGSet),
           detp_450K=detP,
           pd_450K=pData(RGSet))
rm(RGSet, detP)

load(file=file.path(path_rgset, "SGA_Infants_RGSet.RData"))

# Remove 4 arrays
is_okay <- !is.na(pData(rgSet)$gender)
rgSet <- rgSet[, is_okay]
detP <- detP[, colnames(rgSet)]

rg$r_EPIC <- getRed(rgSet)
rg$g_EPIC <-getGreen(rgSet)
rg$detp_EPIC <- detP
rg$pd_EPIC <- pData(rgSet)
rm(rgSet, detP)


#####  Make summary stats  #####

rg$pd_450K$detected_probes_05 <- colSums(rg$detp_450K <= 0.05)
rg$pd_450K$detected_probes_05_pct <- rg$pd_450K$detected_probes_05 / nrow(rg$detp_450K) * 100
rg$pd_450K$signal_avg_grn <- colMeans(rg$g_450K)
rg$pd_450K$signal_avg_red <- colMeans(rg$r_450K)

rg$pd_EPIC$detected_probes_05 <- colSums(rg$detp_EPIC <= 0.05)
rg$pd_EPIC$detected_probes_05_pct <- rg$pd_EPIC$detected_probes_05 / nrow(rg$detp_EPIC) * 100
rg$pd_EPIC$signal_avg_grn <- colMeans(rg$g_EPIC)
rg$pd_EPIC$signal_avg_red <- colMeans(rg$r_EPIC)


#####  Cex settings  #####

# Make consistent cex pattern

mycex <- c("axis"=0.9, "lab"=1.0, "main"=1.0, "cex"=0.2)


#####  Red-Green by Slide  #####

pdf(file=file.path(path_results, "Figure02_Red_green_signals.pdf"), paper = "a4")
par(mfcol=c(3, 2), mar=c(2.5, 4, 3, 2) + 0.1)

ylims <- list("green"=range(rg$pd_450K$signal_avg_grn, rg$pd_EPIC$signal_avg_grn),
              "red"=range(rg$pd_450K$signal_avg_red, rg$pd_EPIC$signal_avg_red),
              "pos"=range(rg$pd_450K$detected_probes_05_pct, rg$pd_EPIC$detected_probes_05_pct))


# 450K

batches_450K <- c(1, 5, 9, 13, 17, 25)
bad_cols_450K <- c("white", "darkgoldenrod2")[factor(unique(rg$pd_450K$array_num) %in% batches_450K)]

boxplot(signal_avg_grn ~ array_num, data=rg$pd_450K,
        cex.axis=mycex["axis"], cex.lab=mycex["lab"], cex.main=mycex["main"],
        col=bad_cols_450K,
        ylim=ylims$green,
        las=2, xlab="", ylab="Green Signal (Average)", main="Green Signal (Cy3)")
abline(v=batches_450K - 0.5, col="grey80")
mtext(LETTERS[1:length(batches_450K)], at = batches_450K, cex = 0.6, col = "grey30")

boxplot(signal_avg_red ~ array_num, data=rg$pd_450K,
        cex.axis=mycex["axis"], cex.lab=mycex["lab"], cex.main=mycex["main"],
        col=bad_cols_450K,
        ylim=ylims$red,
        las=2, xlab="", ylab="Red Signal (Average)", main="Red Signal (Cy5)")
abline(v=batches_450K - 0.5, col="grey80")
mtext(LETTERS[1:length(batches_450K)], at = batches_450K, cex = 0.6, col = "grey30")

boxplot(detected_probes_05_pct ~ array_num, data=rg$pd_450K,
        las=2, cex.axis=mycex["axis"], cex.lab=mycex["lab"], cex.main=mycex["main"],
        col=bad_cols_450K,
        ylim=ylims$pos,
        xlab="Batch (slide)", ylab="Detected probes (p <= 0.05)",
        main="Positively detected probes")
abline(v=batches_450K - 0.5, col="grey80")
mtext(LETTERS[1:length(batches_450K)], at = batches_450K, cex = 0.6, col = "grey30")

# EPIC

batches_EPIC <- c(1, 13)
bad_cols_EPIC <- c(rep("white", 12), rep("darkgoldenrod2", 10))

boxplot(signal_avg_grn ~ array_num, data=rg$pd_EPIC,
        cex.axis=mycex["axis"], cex.lab=mycex["lab"], cex.main=mycex["main"],
        col=bad_cols_EPIC,
        ylim=ylims$green,
        las=2, xlab="", ylab="Green Signal (Average)", main="Green Signal (Cy3)")
abline(v=batches_EPIC - 0.5, col="grey80")
mtext(LETTERS[1:length(batches_EPIC)], at = batches_EPIC, cex = 0.6, col = "grey30")

boxplot(signal_avg_red ~ array_num, data=rg$pd_EPIC,
        cex.axis=mycex["axis"], cex.lab=mycex["lab"], cex.main=mycex["main"],
        col=bad_cols_EPIC,
        ylim=ylims$red,
        las=2, xlab="", ylab="Red Signal (Average)", main="Red Signal (Cy5)")
abline(v=batches_EPIC - 0.5, col="grey80")
mtext(LETTERS[1:length(batches_EPIC)], at = batches_EPIC, cex = 0.6, col = "grey30")

boxplot(detected_probes_05_pct ~ array_num, data=rg$pd_EPIC,
        las=2, cex.axis=mycex["axis"], cex.lab=mycex["lab"], cex.main=mycex["main"],
        col=bad_cols_EPIC,
        ylim=ylims$pos,
        xlab="Batch (slide)", ylab="Detected probes (p <= 0.05)",
        main="Positively detected probes")
abline(v=batches_EPIC - 0.5, col="grey80")
mtext(LETTERS[1:length(batches_EPIC)], at = batches_EPIC, cex = 0.6, col = "grey30")

dev.off()

