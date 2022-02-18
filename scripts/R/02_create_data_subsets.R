rm(list = ls())

#####  Dependencies  #####

library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

#####  Paths  #####

path_base <- "/home/ros259/R_projects/450K_Batch_effects"
path_save_rdata <- file.path(path_base, "data", "rdata")


#####  Main loop  #####

for(prefix in c("GSE89278", "EPIC")) {
    
    cat(prefix, "loading MSets")
    load(file.path(path_save_rdata, paste(prefix, "MSet.Raw.RData", sep="_")))
    load(file.path(path_save_rdata, paste(prefix, "MSet.Illumina.RData", sep="_")))
    load(file.path(path_save_rdata, paste(prefix, "MSet.noob.RData", sep="_")))
    load(file.path(path_save_rdata, paste(prefix, "MSet.swan.RData", sep="_")))
    load(file.path(path_save_rdata, paste(prefix, "MSet.noob.swan.RData", sep="_")))
    load(file.path(path_save_rdata, paste(prefix, "MSet.bmiq.RData", sep="_")))
    load(file.path(path_save_rdata, paste(prefix, "MSet.noob.bmiq.RData", sep="_")))
    load(file.path(path_save_rdata, paste(prefix, "gset.funnorm.RData", sep="_")))
    load(file.path(path_save_rdata, paste(prefix, "MSet.enmix.RData", sep="_")))
    load(file.path(path_save_rdata, paste(prefix, "MSet.enmix.quantile.RData", sep="_")))
    load(file.path(path_save_rdata, paste(prefix, "MSet.dasen.RData", sep="_")))
    load(file.path(path_save_rdata, paste(prefix, "MSet.noob.dasen.RData", sep="_")))
    
    #  Process data to form lists of M-values
    batch <- MSet_hm$factors$batch
    expt <- MSet_hm$factors$expt
    pd <- as.data.frame(pData(MSet.noob))
    

    cat(", saving")
    hm <- list("raw"=MSet_hm, "illumina"=MSet.illumina_hm,
               "noob"=MSet.noob_hm, "swan"=MSet.swan_hm,
               "noobswan"=MSet.noob.swan_hm, "func"=gset.funnorm_hm,
               "bmiq"=MSet.bmiq_hm, "noobbmiq"=MSet.noob.bmiq_hm,
               "enmix"=MSet.enmix_hm, "enmixquantile"=MSet.enmix.quantile_hm,
               "dasen"=MSet.dasen_hm, "noobdasen"=MSet.noob.dasen_hm)
    save(hm, pd, file=file.path(path_save_rdata, paste(prefix, "harmanresults.RData", sep="_")))
    rm(hm)
    rm(list=ls(pattern = "_hm"))
    gc()

    # BMIQ doesn't return MSets
    save(MSet, MSet.illumina, MSet.noob, MSet.swan,
         MSet.noob.swan, gset.funnorm,
         MSet.enmix, MSet.enmix.quantile,
         MSet.dasen, MSet.noob.dasen,
         file=file.path(path_save_rdata, paste(prefix, "methylsets.RData", sep="_")))
    rm(MSet, MSet.illumina, MSet.noob, MSet.swan, MSet.noob.swan, gset.funnorm,
       MSet.enmix, MSet.enmix.quantile, MSet.dasen, MSet.noob.dasen)
    gc()

    ms <- list("raw"=raw_m, "illumina"=illumina_m,
               "noob"=noob_m, "swan"=swan_m,
               "noobswan"=noob.swan_m, "func"=funnorm_m,
               "bmiq"=bmiq_m,"noobbmiq"=noob.bmiq_m,
               "enmix"=enmix_m, "enmixquantile"=enmix.quantile_m,
               "dasen"=dasen_m, "noobdasen"=noob.dasen_m)
    
    ms_harman <- list("raw_hm"=raw_harman_m, "illumina_hm"=illumina_harman_m,
                      "noob_hm"=noob_harman_m, "swan_hm"=swan_harman_m,
                      "noobswan_hm"=noob.swan_harman_m, "func_hm"=funnorm_harman_m,
                      "bmiq_hm"=bmiq_harman_m,"noobbmiq_hm"=noob.bmiq_harman_m,
                      "enmix_hm"=enmix_harman_m, "enmixquantile_hm"=enmix.quantile_harman_m,
                      "dasen_hm"=dasen_harman_m, "noobdasen_hm"=noob.dasen_harman_m)

    ms_combat <- list("raw_ct"=raw_combat_m, "illumina_ct"=illumina_combat_m,
                      "noob_ct"=noob_combat_m, "swan_ct"=swan_combat_m,
                      "noobswan_ct"=noob.swan_combat_m, "func_ct"=funnorm_combat_m,
                      "bmiq_ct"=bmiq_combat_m,"noobbmiq_ct"=noob.bmiq_combat_m,
                      "enmix_ct"=enmix_combat_m, "enmixquantile_ct"=enmix.quantile_combat_m,
                      "dasen_ct"=dasen_combat_m, "noobdasen_ct"=noob.dasen_combat_m)

    cat(", done. Running PCAs")
    pcas <- list(original=lapply(ms, function(x) prcomp(t(x))$x),
                 harman=lapply(ms_harman, function(x) prcomp(t(x))$x),
                 combat=lapply(ms_combat, function(x) prcomp(t(x))$x))
    
    if(prefix == "EPIC") {
        anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    }
    if(prefix == "GSE89278") {
        anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)    
    }
    
    autosomal_probes <- rownames(anno)[anno$chr != "chrX"]
    cat(", now Autosomal PCAs")

    pcas_notx <- list(original=lapply(ms, function(x) prcomp(t(x[autosomal_probes, ]))$x),
                 harman=lapply(ms_harman, function(x) prcomp(t(x[autosomal_probes, ]))$x),
                 combat=lapply(ms_combat, function(x) prcomp(t(x[autosomal_probes, ]))$x))
    cat(", done.")
    rm(list=ls(pattern="_m$"))
    gc()
    
    #####  Save and quit  #####
    cat(" Saving files")
    save(pcas, pcas_notx, pd, file=file.path(path_save_rdata, paste(prefix, "pcas.RData", sep="_")))
    rm(pcas, pcas_notx)
    gc()
    save(ms, pd, batch, expt,
         file=file.path(path_save_rdata, paste(prefix, "original_data_matrices.RData", sep="_")))
    betas <- lapply(ms, ilogit2)
    save(betas, pd, batch, expt,
         file=file.path(path_save_rdata, paste(prefix, "original_betas.RData", sep="_")))
    rm(ms, betas)
    gc()
    
    save(ms_harman, pd, batch, expt,
         file=file.path(path_save_rdata, paste(prefix, "harman_data_matrices.RData", sep="_")))
        betas_harman <- lapply(ms_harman, ilogit2)
    save(betas_harman, pd, batch, expt,
         file=file.path(path_save_rdata, paste(prefix, "harman_betas.RData", sep="_")))
    rm(ms_harman, betas_harman)
    gc()
    
    save(ms_combat, pd, batch, expt,
         file=file.path(path_save_rdata, paste(prefix, "combat_data_matrices.RData", sep="_")))
    betas_combat <- lapply(ms_combat, ilogit2)
    save(betas_combat, pd, batch, expt,
         file=file.path(path_save_rdata, paste(prefix, "combat_betas.RData", sep="_")))
    rm(ms_combat, betas_combat)
    cat(", done.\n")
    gc()
    Sys.sleep(10)
}


q("no")