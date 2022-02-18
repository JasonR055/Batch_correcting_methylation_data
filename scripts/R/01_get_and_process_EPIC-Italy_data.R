rm(list=ls())


#####  Dependencies  #####

library(Harman)
library(limma)
library(sva)
library(minfi)
library(R.utils)
library(EpiDISH)
library(wateRmelon)
library(RColorBrewer)
library(IlluminaHumanMethylation450kmanifest)


#####  Paths  #####

path_base <- "/home/ros259/R_projects/450K_Batch_effects"
path_peapod <- "/home/ros259/R_projects/SGA Infants"
path_rgset <- file.path(path_peapod, "data", "rdata")
path_data <- file.path(path_base, "data")
path_save_rdata <- file.path(path_base, "data", "rdata")
path_idats <- file.path(path_base, "data", "idats")
path_code <- file.path(path_base, "scripts", "R")


#####  Functions  #####

ArrayNumber <- function(pd) {
    arrayorder <- character()
    prev <- ''
    for(a in pd$array) {
        if(a != prev) {
            arrayorder <- c(arrayorder, a)
        }
        prev <- a
    }
    arrayorder
}


#####  Load Atlas data  #####

atlas_folder = "/media/data/R_projects/450K_methylome_atlas/data/GEO_rdata"
load(file.path(atlas_folder, "GSE51032.RData"))


#####  Make phenotype data  #####

pd <- pData(my_geo)
#rm(my_geo)
pd$array_id <- as.character(pd$title)
pd$array <- sub('(\\d+)_R\\d+C\\d+', '\\1', as.character(pd$title))
pd$row <- sub('\\d+_R(\\d+)C\\d+', '\\1', as.character(pd$title))
pd$col <- sub('\\d+_R\\d+C(\\d+)', '\\1', as.character(pd$title))
pd$array_num <- as.integer(factor(pd$array, levels=unique(ArrayNumber(pd))))

pd$harman_batch <- factor(paste(pd$array_num, toupper(substr(pd$gender, 1, 1)), sep='_'),
                          levels=paste(rep(1:length(unique(pd$array_num)),each=2), c('F', 'M'), sep='_'),
                          ordered=TRUE)

batch_cols <- rep(brewer.pal(10, 'Paired'), 15)[1:length(levels(pd$harman_batch))]
pd$batch_col <- batch_cols[pd$harman_batch]
pd$Basename <- paste(pd$geo_accession, pd$array_id, sep="_")


#####  Download Raw IDAT files  #####

path_epic_data = file.path(path_data, "idats_EPIC-Italy")

data_file <- file.path(path_data, "GSE51032_RAW.tar")
if(!file.exists(data_file)) {
    download.file(url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51032/suppl/GSE51032_RAW.tar",
                  destfile=data_file)
    untar(data_file, exdir = path_epic_data)
    
    for(f in list.files(path_epic_data, pattern = "idat.gz", full.names = TRUE)) {
        R.utils::gunzip(f)    
    }
}

RGSet <- read.metharray.exp(base = path_epic_data, extended = TRUE, verbose=TRUE)
stopifnot(identical(1:ncol(RGSet), match(colnames(RGSet), as.character(pd$Basename))))

# record detection rate
detP <- detectionP(RGSet)
pd$n_failed_probes <- colSums(detP > 0.01)
pd$mean_detection_p <- colMeans(detP)

#####  Load Cell references  #####

# EpiDISH cell-types
blood <- as.matrix(read.csv(file.path(path_data, "blood_reference.csv"), row.names = 1))

# Use BMIQ normalisation as the EpiDISH authors used this for the reference
MSet <- preprocessRaw(RGSet)
betas_bmiq <- BMIQ(MSet)
cells_IC <- epidish(betas_bmiq, blood, method = "RPC")

stopifnot(identical(rownames(cells_IC$estF), colnames(MSet)))

pd <- cbind(pd, cells_IC$estF)
pd$immune_group <- cut(pd$Neutro, breaks=c(0.0, mean(pd$Neutro), 1), labels = c("ic_low", "ic_high"))
pData(RGSet) <- DataFrame(pd)
rm(betas_bmiq, MSet, cells_IC)

MSet.noob <- preprocessNoob(RGSet, dyeCorr = TRUE)
GMSet.noob <- mapToGenome(MSet.noob)
array_gender <- getSex(GMSet.noob, cutoff = -2)
pData(RGSet)$predictedSex <- array_gender$predictedSex
pData(MSet.noob)$predictedSex <- array_gender$predictedSex

table(pData(RGSet)$predictedSex == pData(RGSet)$gender.ch1)

my_prefix <- "GSE51032"

save(RGSet, detP,
     file=file.path(path_save_rdata, paste(my_prefix, "RGSet.RData", sep="_")),
     compress=TRUE)
rm(detP, pd)
gc()


#####  Call normalisation script  #####

#my_rndseed <- 42

# Set up model variables
table(expt = paste(pData(RGSet)$gender, pData(RGSet)$immune_group), batch = pData(RGSet)$array_num)
harman_expt <- paste(pData(RGSet)$gender, pData(RGSet)$immune_group)
batch <- pData(RGSet)$array_num
modcombat = model.matrix(~gender + Neutro, data=pData(RGSet))


##### noob  #####

noob_m <- getM(MSet.noob)

MSet.noob_hm <- harman(datamatrix = noob_m, expt = harman_expt, batch = batch)
noob_harman_m <- reconstructData(MSet.noob_hm)
noob_combat_m <- ComBat(dat=noob_m, batch=batch, mod=modcombat)

session <- sessionInfo()

save(MSet.noob, noob_m, MSet.noob_hm, noob_harman_m, noob_combat_m, session,
     file=file.path(path_save_rdata, paste(my_prefix, "MSet.noob.RData", sep="_")),
     compress=TRUE)


#####  PCA plots  #####

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
autosomal_probes <- rownames(anno)[anno$chr != "chrX"]

pca <- list(all=prcomp(t(noob_m))$x,
            auto=prcomp(t(noob_m[autosomal_probes, ]))$x)

PCAPlot(pca$all, 1, 2, labels = pd$plate, cex=0.7, col=pd$batch_col)

PCAPlot(pca$auto, 3, 4, labels = pd$array_num, cex=0.7, col=c("#FF000060", "#00000060")[factor(pd$prop_IC > 0.02)])
PCAPlot(pca$auto, 1, 2, labels = pd$array_num, cex=0.7, col=c("#", "#00000060")[factor(pd$gender)])
PCAPlot(pca$auto, 5, 6, labels = pd$array_num, cex=0.7, col=pd$batch_col)
PCAPlot(pca$all, 1, 2, labels = pd$array_num, cex=0.7, col=pd$batch_col)


diff_m = abs(noob_combat_m - noob_m)
med_diff_m = colMeans(diff_m)
bad_batches <- c(29, 32, 36, 39, 41, 43, 44, 45, 47, 50, 52, 54, 56, 59, 62, 63, 64, 66, 68)
boxplot(med_diff_m ~ pd$array_num, col=c("black", "red")[factor(unique(pd$array_num) %in% bad_batches)], las=2, cex.axis=0.7)
abline(h=0.2, col="grey50")
x=tapply(med_diff_m, pd$array_num, median)
names(x)[x > 0.2]

q("yes")
