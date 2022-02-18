rm(list=ls())

library(RColorBrewer)
library(minfi)
library(FDb.InfiniumMethylation.hg19)
library(Ckmeans.1d.dp)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(BSgenome.Hsapiens.UCSC.hg19)


#####  Load files  #####

path_snps <- "/home/ros259/R_projects/utility_scripts"
dbsnp <- "snp151"
load(file=file.path(path_snps, paste(dbsnp, "Common.RData", sep="")))

path_base <- "/home/ros259/R_projects/450K_Batch_effects"
path_data <- file.path(path_base, "data")
xrct_file <- file.path(path_data, "48639-non-specific-probes-Illumina450k.csv")
xrct_probes <- read.csv(file=xrct_file, stringsAsFactors=FALSE)
rm(xrct_file)


for(i in list(c(1,3), 2, 4, 5)) {
  prefix_set <- c("EpiSCOPE", "EPIC-Italy", "BodyFatness", "NOVI", "URECA")
  names(prefix_set) <- c("GSE89278", "GSE51032", "EPIC", "GSE128821", "GSE132181")
  prefix_alts <- prefix_set[i]
  print(prefix_alts)
  path_base <- "/home/ros259/R_projects/450K_Batch_effects"
  path_code <- file.path(path_base, "scripts", "R")
  source(file.path(path_code, "create_probe_clusters.R"))
  rm(list=setdiff(ls(), c("dbsnp", "dbsnp_gr", "xrct_probes")))
  gc()
  Sys.sleep(10)
}

