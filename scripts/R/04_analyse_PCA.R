rm(list=ls())

library(RColorBrewer)


#####  Paths  #####

path_base <- "/home/ros259/R_projects/450K_Batch_effects"
path_data <- file.path(path_base, "data")
path_save_rdata <- file.path(path_base, "data", "rdata")
path_results <- file.path(path_base, "results")
path_idats <- file.path(path_base, "data", "idats")
path_code <- file.path(path_base, "scripts", "R")


#####  Functions  #####

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



PCA_Montage <- function(method, pcas, pd, x, y) {
    
    plot_names <- c("raw"="Raw", "illumina"="Illumina", "noob"="noob",
                    "swan"="SWAN", "noobswan"="noob+SWAN", "func"="Functional",
                    "bmiq"="BMIQ", "noobbmiq"="noob+BMIQ",
                    "enmix"="ENmix", "enmixquantile"="ENmix-Quantile",
                    "dasen"="Dasen", "noobdasen"="noob+Dasen")
    
    stopifnot(method %in% names(plot_names))
    
    my_correction <- c("Original" = method,
                       "Harman" = paste(method, "hm", sep="_"),
                       "ComBat" = paste(method, "ct", sep="_"))
    

    pca_xlim <- range(sapply(1:3, function(i) pcas[[i]][[my_correction[i]]][, x]))
    pca_ylim <- range(sapply(1:3, function(i) pcas[[i]][[my_correction[i]]][, y]))
    
    
    for(i in 1:length(my_correction)) {
        PCAPlot(pcas[[i]][[as.vector(my_correction[i])]], dimx=x, dimy=y,
                labels=pd$array_num,
                main=paste(plot_names[method], names(my_correction)[i]),
                xlim=pca_xlim,
                ylim=pca_ylim,
                col=pd$batch_col,
                cex=0.8, cex.axis=0.7, cex.lab=0.8, cex.main=1
        )
    }
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

prefixes <- c("GSE89278", "EPIC")
prefix_alts <- c("450K", "EPIC")

pca <- list()

for(i in 1:length(prefixes)) {
    load(file=file.path(path_save_rdata, paste(prefixes[i], "pcas.RData", sep="_")))
    
    if(prefix_alts[i] == "450K") {
        pd$array <- sub('(\\d+)_R\\d+C\\d+', '\\1', rownames(pd))
        pd$row <- sub('\\d+_R(\\d+)C\\d+', '\\1', rownames(pd))
        pd$col <- sub('\\d+_R\\d+C(\\d+)', '\\1', rownames(pd))
        pd$harman_batch <- factor(paste(pd$array_num, pd$sex, sep='_'),
                                  levels=paste(rep(1:length(unique(pd$array_num)),each=2), c('F', 'M'), sep='_'),
                                  ordered=TRUE)
        batch_cols <- rep(brewer.pal(10, 'Paired'), 7)[1:length(levels(pd$harman_batch))]
        pd$batch_col <- batch_cols[pd$harman_batch]
        pd$harman_expt <- paste(pd$sex, pd$pheno)
    }
    if(prefix_alts[i] == "EPIC") {
        
        pd$array <- sub('(\\d+)_R\\d+C\\d+', '\\1', rownames(pd))
        pd$row <- sub('\\d+_R(\\d+)C\\d+', '\\1', rownames(pd))
        pd$col <- sub('\\d+_R\\d+C(\\d+)', '\\1', rownames(pd))
        pd$harman_batch <- factor(paste(pd$array_num, toupper(substr(pd$gender, 1, 1)), sep='_'),
                                  levels=paste(rep(1:length(unique(pd$array_num)),each=2), c('F', 'M'), sep='_'),
                                  ordered=TRUE)
        batch_cols <- rep(brewer.pal(10, 'Paired'), 5)[1:length(levels(pd$harman_batch))]
        pd$batch_col <- batch_cols[pd$harman_batch]
        
        pd$immune <- cut(pd$prop_IC, breaks=c(0.0, 0.06, 1), labels = c("ic_low", "ic_high"))
        pd$harman_expt <- paste(pd$gender, pd$immune)
    }
    pca[[prefix_alts[i]]]$pca <- pcas
    pca[[prefix_alts[i]]]$pca_notx <- pcas_notx
    pca[[prefix_alts[i]]]$pd <- pd
    rm(pcas, pcas_notx, pd)
}


#####  PCA plots  #####

for(prefix_alt in prefix_alts) {
    
    the_methods <- names(pca[[prefix_alt]]$pca$original)

    pdf(file.path(path_results, paste(prefix_alt, "Supp_Figure0X_PCA.pdf", sep="_")), paper = "a4r")
    par(mfrow=c(3, 3), mar=c(4, 4, 2, 2) + 0.1)
    for(method in the_methods) {
        PCA_Montage(method, pcas = pca[[prefix_alt]][["pca"]], pd = pca[[prefix_alt]][["pd"]], 1, 2)
    }
    dev.off()
    
    pdf(file.path(path_results, paste(prefix_alt, "Supp_Figure0X_PCA_noChrX.pdf", sep="_")))
    par(mfrow=c(3, 3), mar=c(4, 4, 2, 2) + 0.1)
    for(method in the_methods) {
        PCA_Montage(method, pcas = pca[[prefix_alt]][["pca_notx"]], pd = pca[[prefix_alt]][["pd"]], 1, 2)
    }
    for(method in the_methods) {
        PCA_Montage(method, pcas = pca[[prefix_alt]][["pca_notx"]], pd = pca[[prefix_alt]][["pd"]], 3, 4)
    }
    dev.off()
}


pdf(file.path(path_results, "Figure04_PCA_panel.pdf"))
par(mfrow=c(4, 4), mar=c(4, 4, 1, 1) + 0.1)

for(prefix_alt in prefix_alts) {

    PCA_Montage("raw", pcas = pca[[prefix_alt]][["pca_notx"]], pd = pca[[prefix_alt]][["pd"]], 1, 2)
    if(prefix_alt == "450K") {
        my_cells = "Neutro"
        colmetabatch <- c("grey80", "red")[factor(pca[[prefix_alt]]$pd$array_num %in% c(1, 5, 9, 17, 25))]
    } else {
        my_cells = "prop_IC"
        colmetabatch <- c("grey80", "red")[factor(pca[[prefix_alt]]$pd$plate)]
    }
    
    cell_factor <- cut(pca[[prefix_alt]]$pd[, my_cells], breaks=5)
    colcells <- brewer.pal(n=length(levels(cell_factor)), "RdPu")
    colcells <- colcells[cell_factor]
    
    PCAPlot(pca = pca[[prefix_alt]][["pca_notx"]]$original$noob, dimx=1, dimy=2,
            labels=pca[[prefix_alt]]$pd$array_num,
            main="noob Original: super",
            col=colmetabatch,
            cex=0.8, cex.axis=0.7, cex.lab=0.8, cex.main=1)
    
    PCA_Montage("raw", pcas = pca[[prefix_alt]][["pca_notx"]], pd = pca[[prefix_alt]][["pd"]], 3, 4)
    
    if(prefix_alt == "450K") { # superbatch for 450K
        PCAPlot(pca = pca[[prefix_alt]][["pca_notx"]]$original$noob, dimx=3, dimy=4,
                labels=pca[[prefix_alt]]$pd$array_num,
                main="noob Original: cells",
                col=colcells,
                cex=0.8, cex.axis=0.7, cex.lab=0.8, cex.main=1)
    } else { # cells for EPIC
        PCAPlot(pca = pca[[prefix_alt]][["pca_notx"]]$original$noob, dimx=1, dimy=2,
                labels=pca[[prefix_alt]]$pd$array_num,
                main="noob Original: cells",
                col=colcells,
                cex=0.8, cex.axis=0.7, cex.lab=0.8, cex.main=1)
    }
}
dev.off()


