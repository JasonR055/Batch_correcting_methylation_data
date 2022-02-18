#####  Dependencies  #####

library(IlluminaHumanMethylationEPICmanifest)
library(HELP)

#####  Paths  #####

path_base <- "/home/ros259/R_projects/450K_Batch_effects"
path_rdata <- file.path(path_base, "data", "rdata")

cores <- 20

probes_typeI <- getProbeInfo(IlluminaHumanMethylationEPICmanifest, type = "I")
probes_typeI$Type="typeI"
probes_typeII <- getProbeInfo(IlluminaHumanMethylationEPICmanifest, type = "II")
probes_typeII$Type="typeII"

# Type I

probes_typeI$tm_A <- unlist(mclapply(as.character(probes_typeI$ProbeSeqA),
                                     calcTm,
                                     mc.preschedule=TRUE, mc.cores=cores))
probes_typeI$tm_B <- unlist(mclapply(as.character(probes_typeI$ProbeSeqB),
                                     calcTm,
                                     mc.preschedule=TRUE, mc.cores=cores))
probes_typeI$tm_avg <- (probes_typeI$tm_A + probes_typeI$tm_B) / 2  # mean

# Type II

seq_typeII <- as.character(probes_typeII$ProbeSeqA)

# TypeII probes have probes with the degenerate base R. Need to resolve this
# to A and G, then take the average Tm.

probes_typeII$tm_A <- unlist(mclapply(gsub('R', 'A', seq_typeII),
                                      calcTm,
                                      mc.preschedule=TRUE, mc.cores=cores))
probes_typeII$tm_B <- unlist(mclapply(gsub('R', 'G', seq_typeII),
                                      calcTm,
                                      mc.preschedule=TRUE, mc.cores=cores))
probes_typeII$tm_avg <- (probes_typeII$tm_A + probes_typeII$tm_B) / 2  # mean

rm(seq_typeII)

these_cols <- c("Name", "Type", "nCpG", "tm_A", "tm_B", "tm_avg")
probes <- rbind(probes_typeI[, these_cols], probes_typeII[, these_cols])
row.names(probes) <- probes$Name

#####  Save  #####

session <- sessionInfo()
save(probes, session, file=file.path(path_rdata, "EPIC_probe_tm_info.RData"))

q("no")
