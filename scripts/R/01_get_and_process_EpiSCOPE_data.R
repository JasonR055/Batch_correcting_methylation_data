rm(list=ls())


#####  Dependencies  #####

library(RColorBrewer)
library(GEOquery)
library(limma)
library(minfi)
library(minfiData)
library(EpiDISH)
library(wateRmelon)
library(IlluminaHumanMethylation450kmanifest)
library(Harman)
library(sva)
library(R.utils)


#####  Paths  #####

path_base <- "/home/ros259/R_projects/450K_Batch_effects"
path_data <- file.path(path_base, "data")
path_save_rdata <- file.path(path_base, "data", "rdata")
path_idats <- file.path(path_base, "data", "idats")
path_code <- file.path(path_base, "scripts", "R")
path_ftp <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89278/suppl"


#####  Functions  #####


extract_trait <- function(this_col) {
    
    x <- pData(my_geo)[, this_col]
    x <- tolower(as.character(x))
    sub(".+: (.+)", "\\1", x)
}


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


#####  Load RData  #####

#####  GSE89278: EpiSCOPE Infant PBL   #####

my_prefix <- "GSE89278"

# Download TAR of idat files if not present
raw_file <- paste(my_prefix, "RAW.tar", sep="_")
if(!file.exists(file.path(path_data, raw_file))) {
    download.file(url=file.path(path_ftp, raw_file),
                  destfile=file.path(path_data, raw_file))
    untar(file.path(path_data, raw_file), exdir = path_idats)
    
    for(f in list.files(path_idats, pattern=".idat.gz", full.names = TRUE)) {
        cat("Unzipping", f, "\n")
        gunzip(f, remove=TRUE)    
    }
}

# Get the phenotype data from GEO
my_geo <- getGEO(GEO = my_prefix, destdir = path_data, getGPL=FALSE)[[1]]

# Craft a new phenotype file
sex <- extract_trait("characteristics_ch1")
sex <- toupper(substr(sex, 1, 1))
pheno <- as.character(pData(my_geo)$source_name_ch1)
pheno <- tolower(sub("Neonatal blood spot,(.+) supplement", "\\1", pheno))

pd <- data.frame(array_id=as.character(pData(my_geo)$title),
                 geo_accession=as.character(pData(my_geo)$geo_accession),
                 tissue=extract_trait("characteristics_ch1.1"),
                 pheno=pheno,
                 sex=sex,
                 row.names = as.character(pData(my_geo)$title),
                 stringsAsFactors = FALSE)

pd$array <- sub('(\\d+)_R\\d+C\\d+', '\\1', rownames(pd))
pd$row <- sub('\\d+_R(\\d+)C\\d+', '\\1', rownames(pd))
pd$col <- sub('\\d+_R\\d+C(\\d+)', '\\1', rownames(pd))
pd$array_num <- as.integer(factor(pd$array, levels=ArrayNumber(pd)))
pd$harman_batch <- factor(paste(pd$array_num, substr(pd$sex, 1, 1), sep='_'),
                          levels=paste(rep(1:length(unique(pd$array_num)),each=2), c('F', 'M'), sep='_'),
                          ordered=TRUE)

batch_cols <- rep(brewer.pal(10, 'Paired'), 7)[1:length(levels(pd$harman_batch))]
pd$batch_col <- batch_cols[pd$harman_batch]
pd$Basename <- paste(pd$geo_accession, pd$array_id, sep="_")


#####  Load Cell references  #####

# EpiDISH cell-types
epiFibIC <- as.matrix(read.csv(file.path(path_data, "EpiFibIC_reference.csv"), row.names = 1))
blood <- as.matrix(read.csv(file.path(path_data, "blood_reference.csv"), row.names = 1))


#####  Load and preprocess IDATs  #####

RGSet <- read.metharray.exp(base = path_idats, targets = pd)
detP <- detectionP(RGSet)

# Use BMIQ normalisation as the EpiDISH authors used this for the reference
MSet <- preprocessRaw(RGSet)
betas_bmiq <- BMIQ(MSet)
#epiFibIC_EPIC <- epiFibIC[rownames(epiFibIC) %in% rownames(betas_bmiq), ]
cells_epiFibIC <- epidish(betas_bmiq, epiFibIC, method = "RPC")
cells_IC <- epidish(betas_bmiq, blood, method = "RPC")

stopifnot(identical(rownames(cells_IC$estF), colnames(MSet)))

pd <- cbind(pd, cells_IC$estF)
pData(RGSet) <- DataFrame(pd)
rm(betas_bmiq, MSet)

save(RGSet, detP,
     file=file.path(path_save_rdata, paste(my_prefix, "RGSet.RData", sep="_")),
     compress=TRUE)
rm(detP)


#####  Call normalisation script  #####

my_prefix <- "GSE89278"
my_rndseed <- 42

# Set up model variables
harman_expt <- paste(pData(RGSet)$sex, pData(RGSet)$pheno)
batch <- pData(RGSet)$array_num
modcombat = model.matrix(~sex + pheno, data=pData(RGSet))

source(file.path(path_code, "array_normalisation.R"))
