rm(list=ls())

# Change timeout download options to 600 min, 10 hrs (36000 secs)
options(timeout = 36000)

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

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132181
# https://www.tandfonline.com/doi/full/10.1080/15592294.2020.1817290
# GSE132181

my_gse <- "GSE132181"
my_geo <- getGEO(GEO = my_gse, destdir = path_data, getGPL=FALSE)[[1]]


#####  Make phenotype data  #####

pd <- pData(my_geo)
#rm(my_geo)
pd$array_id <- pd$title
pd$sample_id <- substr(pd$title, 1, 9)
pd$array_name <- sub('.*/suppl/GSM\\d+_(.*)_Grn.idat.gz$', '\\1', as.character(pd$supplementary_file))
pd$array <- sub('(\\d+)_R\\d+C\\d+', '\\1', pd$array_name)
pd$row <- sub('\\d+_R(\\d+)C\\d+', '\\1', pd$array_name)
pd$col <- sub('\\d+_R\\d+C(\\d+)', '\\1', pd$array_name)
pd$array_num <- as.integer(factor(pd$array, levels=unique(ArrayNumber(pd))))
pd$Basename <- paste(pd$geo_accession, pd$array_name, sep="_")
pd$gender <- substr(toupper(as.character(pd$`Sex:ch1`)), 1, 1)
pd$race <- as.character(pd$`selfreportedrace:ch1`)
pd$plate <- as.integer(pd$`methylationplatenumber:ch1`)
pd$ga <- as.numeric(pd$`gestationalage_birth:ch1`)
pd$collection_site <- as.character(pd$`samplecollectionsite:ch1`)
pd$age <- as.character(pd$`age:ch1`)
pd$tissue <- as.character(pd$`tissue:ch1`)

pd$harman_batch <- factor(paste(pd$array_num, pd$gender, sep='_'),
                          levels=paste(rep(1:length(unique(pd$array_num)),each=2), c('F', 'M'), sep='_'),
                          ordered=TRUE)

batch_cols <- rep(brewer.pal(10, 'Paired'), 16)[1:length(levels(pd$harman_batch))]
pd$batch_col <- batch_cols[pd$harman_batch]


#####  Download Raw IDAT files  #####

path_epic_data = file.path(path_data, "idats_GSE132181")

# ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4787nnn/GSM4787785/suppl/GSM4787785_202242410127_R02C01_Grn.idat.gz

data_file <- file.path(path_data, "GSE132181_RAW.tar")
if(!file.exists(data_file)) {
    download.file(url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132181/suppl/GSE132181_RAW.tar",
                  destfile=data_file)
    untar(data_file, exdir = path_epic_data)
    
    for(f in list.files(path_epic_data, pattern = "idat.gz", full.names = TRUE)) {
        R.utils::gunzip(f)    
    }
}


epic_files <- list.files(path_epic_data, pattern = "Grn.idat", full.names = TRUE)
# GSM3686020_201557540048_R07C01 looks to be a bad truncated array file.
RGSet <- read.metharray.exp(base = path_epic_data, extended = TRUE, verbose=TRUE)
# sanity check that the RGSet and pData are in the same order
stopifnot(identical(1:ncol(RGSet), match(colnames(RGSet), as.character(pd$Basename))))


# record detection rate
detP <- detectionP(RGSet)
pd$n_failed_probes <- colSums(detP > 0.01)
pd$mean_detection_p <- colMeans(detP)

pData(RGSet) <- DataFrame(pd)
# Remove 1 array with a probe failure rate > 5% (GSM3852526)
RGSet <- RGSet[, !pd$n_failed_probes > nrow(detP) * 0.05]

#rm(pd)


#####  Load Cell references  #####

# EpiDISH cell-types
epiFibIC <- as.matrix(read.csv(file.path(path_data, "EpiFibIC_reference.csv"), row.names = 1))
blood <- as.matrix(read.csv(file.path(path_data, "blood_reference.csv"), row.names = 1))

# Use BMIQ normalisation as the EpiDISH authors used this for the reference
MSet <- preprocessRaw(RGSet)
MSet.noob <- preprocessNoob(RGSet)
GMSet.noob <- mapToGenome(MSet.noob)
betas_bmiq <- BMIQ(MSet)

array_gender <- getSex(GMSet.noob, cutoff = -2)
pData(RGSet)$predictedGender <- genderRecode(array_gender$predictedSex)

cell_proportions <- epidish(betas_bmiq, epiFibIC, method = "RPC")
cells_IC <- epidish(betas_bmiq, blood, method = "RPC")

# Sanity check
stopifnot(identical(rownames(cell_proportions$estF), colnames(MSet)))
pData(RGSet)$prop_epi <- cell_proportions$estF[, "Epi"]
pData(RGSet)$prop_IC <- cell_proportions$estF[, "IC"]
pData(RGSet) <- cbind(pData(RGSet), cells_IC$estF)


anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
autosomal_probes <- rownames(anno)[anno$chr != "chrX"]

noob_m <- getM(MSet.noob)
pca <- list(all=prcomp(t(noob_m))$x,
            auto=prcomp(t(noob_m[autosomal_probes, ]))$x)

# Split on dim 1 and 2 by birth/age7, gender and age
PCAPlot(pca$all, 1, 2, labels = pData(RGSet)$age, cex=0.7, col=pData(RGSet)$batch_col)
PCAPlot(pca$auto, 1, 2, labels = substr(pData(RGSet)$tissue, 1, 4), cex=0.6, col=pData(RGSet)$batch_col)
# Split on dim 4 by slide, plate, row
PCAPlot(pca$auto, 3, 4, labels = pData(RGSet)$array_num, cex=0.7, col=pData(RGSet)$batch_col)
PCAPlot(pca$auto, 3, 4, labels = pData(RGSet)$plate, cex=0.7, col=pData(RGSet)$batch_col)
PCAPlot(pca$auto, 3, 4, labels = pData(RGSet)$row, cex=0.7, col=as.integer(pData(RGSet)$row))
# Split on dimension 5 by neutrophil counts
PCAPlot(pca$auto, 5, 6, labels = pData(RGSet)$array_num, cex=0.7, col=brewer.pal(9, "YlGnBu")[cut(pData(RGSet)$Neutro, 9)])
# Dim 7, split by autosomal gender probes
PCAPlot(pca$auto, 7, 8, labels = (1:196)[factor(pData(RGSet)$sample_id)], cex=0.7, col=pData(RGSet)$batch_col)
# Dim 9 and 10, split by race 
PCAPlot(pca$auto, 9, 10, labels = substr(pData(RGSet)$race, 1, 1), cex=0.7, col=pData(RGSet)$batch_col)

immune <- cut(pData(RGSet)$Neutro, breaks=c(0.0, 0.07, 1), include.lowest=TRUE, labels = c("neutro_low", "neutro_high"))
#table(paste(immune, pData(RGSet)$gender), pData(RGSet)$array_num)
table(immune, pData(RGSet)$array_num)
table(paste(immune, pData(RGSet)$gender, pData(RGSet)$tissue), pData(RGSet)$array_num)
table(paste(pData(RGSet)$gender, pData(RGSet)$tissue), pData(RGSet)$array_num)
pData(RGSet)$immune_group <- immune

# Insert modified phenoData back into the MSet
pData(MSet.noob) <- pData(RGSet)

rm(betas_bmiq, MSet, cell_proportions, GMSet.noob, array_gender, immune, epiFibIC)

save(RGSet, detP,
     file=file.path(path_save_rdata, paste(my_gse, "RGSet.RData", sep="_")),
     compress=TRUE)

rm(detP)
gc()


#####  Call normalisation script  #####

# Set up model variables
table(expt = paste(pData(RGSet)$tissue, pData(RGSet)$gender, pData(RGSet)$immune_group), batch = pData(RGSet)$array_num)
table(expt = paste(pData(RGSet)$tissue, pData(RGSet)$gender), batch = pData(RGSet)$array_num)
harman_expt <- paste(pData(RGSet)$tissue, pData(RGSet)$gender)
batch <- pData(RGSet)$array_num
table(harman_expt, batch)
modcombat = model.matrix(~gender + tissue, data=pData(RGSet))

##### noob  #####

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
x=tapply(med_diff_m, pData(MSet.noob)$array_num, median)
names(x)[x > 0.2]
bad_batches <- c(24, 25, 28)
boxplot(med_diff_m ~ pData(MSet.noob)$array_num, col=c("black", "red")[factor(unique(pData(MSet.noob)$array_num) %in% bad_batches)], las=2, cex.axis=0.7)
abline(h=0.2, col="grey50")
boxplot(med_diff_m ~ pData(MSet.noob)$array_num, col=unique(pData(MSet.noob)[, c("array_num", "plate")])$plate, las=2, cex.axis=0.7)

q("no")
