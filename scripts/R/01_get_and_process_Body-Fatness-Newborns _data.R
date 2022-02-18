rm(list=ls())


#####  Dependencies  #####

library(minfi)


#####  Paths  #####

path_base <- "/home/ros259/R_projects/450K_Batch_effects"
path_peapod <- "/home/ros259/R_projects/SGA Infants"
path_rgset <- file.path(path_peapod, "data", "rdata")
path_data <- file.path(path_base, "data")
path_save_rdata <- file.path(path_base, "data", "rdata")
path_idats <- file.path(path_base, "data", "idats")
path_code <- file.path(path_base, "scripts", "R")


#####  Functions  #####



#####  Load data  #####

load(file=file.path(path_rgset, "SGA_Infants_RGSet.RData"))

RGSet <- rgSet
rm(rgSet, detP)

# EPIC row, col fix
pData(RGSet)$row <- sub('\\d+_R(\\d+)C\\d+', '\\1', as.character(pData(RGSet)$Basename))
pData(RGSet)$col <- sub('\\d+_R\\d+C(\\d+)', '\\1', as.character(pData(RGSet)$Basename))


#####  Call normalisation script  #####

my_prefix <- "EPIC"
my_rndseed <- 42

# Set up model variables
immune <- cut(pData(RGSet)$prop_IC, breaks=c(0.0, 0.06, 1), labels = c("ic_low", "ic_high"))
table(expt = paste(pData(RGSet)$gender, immune), batch = pData(RGSet)$array_num)
harman_expt <- paste(pData(RGSet)$gender, immune)
batch <- pData(RGSet)$array_num
modcombat = model.matrix(~gender + prop_IC, data=pData(RGSet))

source(file.path(path_code, "array_normalisation.R"))
