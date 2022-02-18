rm(list=ls())


#####  Dependencies  #####

library(Harman)
library(limma)
library(sva)
library(R.utils)
library(GEOquery)
library(RColorBrewer)
library(minfi)
library(wateRmelon)
library(EpiDISH)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


#####  Paths  #####

path_base <- "/home/ros259/R_projects/450K_Batch_effects"
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

genderRecode <- function(x) {
    
    z <- c("1"="male", "2"="female", "M"="male", "F"="female")
    as.character(z[as.character(x)])
}


extract_generic <- function(my_geo, this_col) {
    
    x <- as.character(pData(my_geo)[, this_col])
    generic <- sub(".+: (.+)", "\\1", x)
    generic <- gsub("_", " ", generic)
    generic <- gsub("\\.", " ", generic)
    generic[toupper(generic) == "NA"] <- NA
    generic
}


PCAPlot <- function(pca, dimx=1, dimy=2, labels, main='', ...) {
    
    plot(pca[, dimx], pca[, dimy],
         type='n',
         xlab=paste("Dim", dimx),
         ylab=paste("Dim", dimy),
         main=main,
         ...
    )
    text(pca[, dimx], pca[, dimy],
         labels=labels,
         ...
    )
    
}


#####  Get GEO pheno data  #####

# Neonatal Neurobehavior and Outcomes in Very Preterm Infants (NOVI) Study
# Buccal swab samples collected near term-equivalent age using the Isohelix Buccal Swab kit
# https://www.ncbi.nlm.nih.gov/pubmed/31004082
# GSE128821

my_gse <- "GSE128821"
my_geo <- getGEO(GEO = my_gse, destdir = path_data, getGPL=FALSE)[[1]]


#####  Make phenotype data  #####

pd <- pData(my_geo)
#rm(my_geo)
pd$array_id <- as.character(pd$title)
pd$array <- sub('buccal swab \\d+ \\[(\\d+)_R\\d+C\\d+\\]', '\\1', as.character(pd$title))
pd$row <- sub('buccal swab \\d+ \\[\\d+_R(\\d+)C\\d+\\]', '\\1', as.character(pd$title))
pd$col <- sub('buccal swab \\d+ \\[\\d+_R\\d+C(\\d+)\\]', '\\1', as.character(pd$title))
pd$array_num <- as.integer(factor(pd$array, levels=unique(ArrayNumber(pd))))
pd$array_name <- sub('buccal swab \\d+ \\[(\\d+_R\\d+C\\d+)\\]', '\\1', as.character(pd$title))
pd$Basename <- paste(pd$geo_accession, pd$array_name, sep="_")
pd$gender <- pd$`gender:ch1`

pd$harman_batch <- factor(paste(pd$array_num, toupper(substr(pd$gender, 1, 1)), sep='_'),
                          levels=paste(rep(1:length(unique(pd$array_num)),each=2), c('F', 'M'), sep='_'),
                          ordered=TRUE)

batch_cols <- rep(brewer.pal(10, 'Paired'), 16)[1:length(levels(pd$harman_batch))]
pd$batch_col <- batch_cols[pd$harman_batch]


#####  Download Raw IDAT files  #####

path_epic_data = file.path(path_data, "idats_GSE128821")

data_file <- file.path(path_data, "GSE128821_RAW.tar")
if(!file.exists(data_file)) {
    download.file(url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE128nnn/GSE128821/suppl/GSE128821_RAW.tar",
                  destfile=data_file)
    untar(data_file, exdir = path_epic_data)
    
    for(f in list.files(path_epic_data, pattern = "idat.gz", full.names = TRUE)) {
        R.utils::gunzip(f)    
    }
}


epic_files <- list.files(path_epic_data, pattern = "Grn.idat", full.names = TRUE)
# GSM3686020_201557540048_R07C01 looks to be a bad truncated array file.
bad_array <- "GSM3686020_201557540048_R07C01"
RGSet <- read.metharray(basenames = epic_files[!grepl(bad_array, epic_files)], verbose=TRUE, force = TRUE)
pd <- pd[!grepl(bad_array, pd$Basename), ]
rm(bad_array)
# sanity check that the RGSet and pData are in the same order
stopifnot(identical(1:ncol(RGSet), match(colnames(RGSet), as.character(pd$Basename))))


# record detection rate
detP <- detectionP(RGSet)
pd$n_failed_probes <- colSums(detP > 0.01)
pd$mean_detection_p <- colMeans(detP)

pd$plate <- LETTERS[factor(pd[, "plate:ch1"])]
# Fill in the missing plate data for 4 samples.
# Assume other arrays from the same slide were on the same plate
pd$plate[pd$geo_accession == "GSM3686258"] <- "D"
pd$plate[pd$geo_accession %in% c("GSM3686281", "GSM3686387", "GSM3686432")] <- "E"

pData(RGSet) <- DataFrame(pd)
RGSet <- RGSet[, !pd$n_failed_probes > nrow(detP) * 0.05]
rm(pd)
# Remove 1 array with a probe failure rate > 5%


#####  Load Cell references  #####

# EpiDISH cell-types
epiFibIC <- as.matrix(read.csv(file.path(path_data, "EpiFibIC_reference.csv"), row.names = 1))

# Use BMIQ normalisation as the EpiDISH authors used this for the reference
MSet <- preprocessRaw(RGSet)
MSet.noob <- preprocessNoob(RGSet)
GMSet.noob <- mapToGenome(MSet.noob)
betas_bmiq <- BMIQ(MSet)

array_gender <- getSex(GMSet.noob, cutoff = -2)
pData(RGSet)$predictedGender <- genderRecode(array_gender$predictedSex)

cell_proportions <- epidish(betas_bmiq, epiFibIC, method = "RPC")
# Sanity check
stopifnot(identical(rownames(cell_proportions$estF), colnames(MSet)))
pData(RGSet)$prop_epi <- cell_proportions$estF[, "Epi"]
pData(RGSet)$prop_IC <- cell_proportions$estF[, "IC"]

immune <- cut(pData(RGSet)$prop_IC, breaks=c(0.0, 0.025, 1), include.lowest=TRUE, labels = c("ic_low", "ic_high"))
table(paste(immune, pData(RGSet)$gender), pData(RGSet)$array_num)
pData(RGSet)$immune_group <- immune

# Insert modified phenoData back into the MSet
pData(MSet.noob) <- pData(RGSet)

#####  PCA plots  #####

anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
autosomal_probes <- rownames(anno)[anno$chr != "chrX"]

noob_m <- getM(MSet.noob)
pca <- list(all=prcomp(t(noob_m))$x,
            auto=prcomp(t(noob_m[autosomal_probes, ]))$x)
PCAPlot(pca$all, 3, 4, labels = pData(RGSet)$array_num, cex=0.7, col=c("#FF8000C0", "#000000C0")[factor(pData(RGSet)$prop_IC > 0.025)])
PCAPlot(pca$auto, 1, 2, labels = pData(RGSet)$array_num, cex=0.7, col=pData(RGSet)$batch_col)
PCAPlot(pca$auto, 3, 4, labels = pData(RGSet)$array_num, cex=0.7, col=pData(RGSet)$batch_col)
PCAPlot(pca$auto, 5, 6, labels = pData(RGSet)$array_num, cex=0.7, col=pData(RGSet)$batch_col)
PCAPlot(pca$auto, 7, 8, labels = pData(RGSet)$array_num, cex=0.7, col=pData(RGSet)$batch_col)

rm(betas_bmiq, MSet, cell_proportions, GMSet.noob, array_gender, immune, epiFibIC)

save(RGSet, detP,
     file=file.path(path_save_rdata, paste(my_gse, "RGSet.RData", sep="_")),
     compress=TRUE)

rm(detP)
gc()


##### Batch correct  #####

# Set up model variables
table(expt = paste(pData(RGSet)$gender, pData(RGSet)$immune_group), batch = pData(RGSet)$array_num)
harman_expt <- paste(pData(RGSet)$gender, pData(RGSet)$immune_group)
batch <- pData(RGSet)$array_num
modcombat = model.matrix(~gender + prop_IC, data=pData(RGSet))

MSet.noob_hm <- harman(datamatrix = noob_m, expt = harman_expt, batch = batch)
noob_harman_m <- reconstructData(MSet.noob_hm)
noob_combat_m = ComBat(dat=noob_m, batch=batch, mod=modcombat)

session <- sessionInfo()

save(MSet.noob, noob_m, MSet.noob_hm, noob_harman_m, noob_combat_m, session,
     file=file.path(path_save_rdata, paste(my_gse, "MSet.noob.RData", sep="_")),
     compress=TRUE)


##### bad batches ######

diff_m = abs(noob_combat_m - noob_m)
med_diff_m = colMeans(diff_m)
bad_batches <- c(2, 3, 4, 9, 16, 17, 19, 20, 28, 29, 30, 32, 35, 39, 40, 45, 47, 50, 53, 60, 69, 70, 72, 73, 75)
boxplot(med_diff_m ~ pd$array_num, col=c("black", "red")[factor(unique(pd$array_num) %in% bad_batches)], las=2, cex.axis=0.7)
abline(h=0.2, col="grey50")

x=tapply(med_diff_m, pd$array_num, median)
names(x)[x > 0.2]
bad_batches

q("no")
