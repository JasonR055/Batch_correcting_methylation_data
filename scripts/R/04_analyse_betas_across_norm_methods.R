rm(list=ls())

library(RColorBrewer)


#####  Paths  #####

path_base <- "/home/ros259/R_projects/450K_Batch_effects"
path_data <- file.path(path_base, "data")
path_save_rdata <- file.path(path_base, "data", "rdata")
path_results <- file.path(path_base, "results")
path_idats <- file.path(path_base, "data", "idats")
path_code <- file.path(path_base, "scripts", "R")
path_snps <- "/home/ros259/R_projects/utility_scripts"

#####  Functions  #####

betaDiffs <- function(betas, betas_corrected) {
    
    # Differences of corrected betas - original betas
    beta_diffs <- list()
    for(method in names(betas)) {
        # sanity check
        cat(method, " ", sep="")
        i <- grep(paste("^", method, "_", sep=""), names(betas_corrected))
        stopifnot(identical(dimnames(betas[[method]]),
                            dimnames(betas_corrected[[i]])))
        beta_diffs[[method]] <- betas_corrected[[i]] - betas[[method]]
    }
    cat(".\n")
    beta_diffs
}


plotLogDiffs <- function(log_diffs, array="450K", correction="harman", this=names(log_diffs[[correction]][[1]]), legend_placement="topright", legend_cex=1, ...) {
    
    # Set1 order is: red, blue, green, purple, orange, yellow, brown, pink, grey
    # Add, black, lgt red #fb8072, lgt purple #bebada
    plot_cols <- c("black", "#fb8072", "#bebada", brewer.pal(9, "Set1"))
    
    names(plot_cols) <- c("raw", "bmiq", "swan", "noob", "noobbmiq", "noobswan", "illumina", "func",  "enmix", "enmixquantile", "dasen", "noobdasen")
    # Replace Color Brewer yellow
    plot_cols[plot_cols == "#FFFF33"] <- "gold1"
    plot_cols <- c(plot_cols, c("raw_typeI"="grey50", "raw_typeII"="grey50"))
    
    my_cols <- plot_cols[this]
    my_lwds <- rep(2, length(this))
    names(my_lwds) <- this
    my_lwds["raw"] <- 3
    my_ltys <- rep(1, length(this))
    names(my_ltys) <- this
    my_ltys["raw_typeI"] <- 2
    my_ltys["raw_typeII"] <- 3
    
    plot(NULL,
         xlab="log10 maximal probe-wise beta difference",
         ylab="Density", ...)
    abline(v=c(-1, -2), col="grey90")
    
    for(method in this) {
        lines(x=log_diffs[[array]][[correction]][[method]]$x,
              y=log_diffs[[array]][[correction]][[method]]$y,
              col=my_cols[method], lwd=my_lwds[method], lty=my_ltys[method])
    }
    if(!is.null(legend_placement)) {
        legend(legend_placement, legend=as.vector(nice_names[this]),
               col=my_cols, cex=legend_cex, lwd=my_lwds, lty=my_ltys)
    }
}


#####  Globals  #####

nice_names <- c("raw"="Raw", "illumina"="Illumina", "noob"="noob",
                "swan"="SWAN", "noobswan"="noob+SWAN", "func"="Functional",
                "bmiq"="BMIQ", "noobbmiq"="noob+BMIQ",
                "enmix"="ENmix", "enmixquantile"="ENmix-Quantile",
                "dasen"="Dasen", "noobdasen"="noob+Dasen",
                "raw_typeI"="Raw (I)", "raw_typeII"="Raw (II)",
                "harman"="Harman", "combat"="ComBat")


#####  Main loops  #####

log_diffs <- list()
for(prefix in c("GSE89278", "EPIC")) {
    prefix_alt <- prefix
    if(prefix == "GSE89278") prefix_alt <- "450K"
    # load files
    load(file=file.path(path_save_rdata, paste(prefix_alt, "probe_tm_info.RData", sep="_")))
    load(file=file.path(path_save_rdata, paste(prefix, "original_betas.RData", sep="_")))
    # If using the IlluminaHumanMethylationEPICmanifest library, 745 probes have been removed from the B4 array,
    # however the library was made with B2 annotations.
    Sys.sleep(5)
    probes <- as.data.frame(probes)
    stopifnot(sum(rownames(betas$raw) %in% rownames(probes)) == nrow(betas$raw))
    probes <- probes[rownames(betas$raw), ]
    
    typeI_probe_names <- probes$Name[probes$Type == "typeI"]
    typeII_probe_names <- probes$Name[probes$Type == "typeII"]

    beta_diffs <- list()
    cgdiffs <- list()
    
    for(correction in c("harman", "combat")) {
        # load corrected betas
        betas_corrected_name <- paste("betas", correction, sep="_")
        load(file=file.path(path_save_rdata, paste(prefix, correction, "betas.RData", sep="_")))
        betas_corrected = get(betas_corrected_name)
        
        cat(correction, "\n")
        # Make beta diffs
        beta_diffs[[correction]] <- betaDiffs(betas, betas_corrected)
        rm(list=c("betas_corrected", "betas_corrected_name"))
    }
    rm(betas)
    gc()
    
    for(correction in names(beta_diffs)) {
        cgranges <- mclapply(beta_diffs[[correction]], function(x) apply(x, 1, range), mc.cores=6)
        for(n in names(cgranges)) rownames(cgranges[[n]]) <- c('min', 'max')
        # For each matrix of feature-wise ranges in cg diffs,
        # find the differences of the min/max values in the columns
        # So, cgdiffs is the difference in the beta adjustments for each probe across samples
        cgdiffs[[correction]] <- mclapply(cgranges, function(x) apply(x, 2, diff), mc.cores=6)
        log_diffs[[prefix_alt]][[correction]] <- mclapply(cgdiffs[[correction]], function(x) density(log10(x)), mc.cores=length(cgdiffs[[correction]]))
        log_diffs[[prefix_alt]][[correction]]$raw_typeI <- density(log10(cgdiffs[[correction]]$raw[names(cgdiffs[[correction]]$raw) %in% typeI_probe_names]))
        log_diffs[[prefix_alt]][[correction]]$raw_typeII <- density(log10(cgdiffs[[correction]]$raw[names(cgdiffs[[correction]]$raw) %in% typeII_probe_names]))
        rm(cgranges)
    }
    
    for(correction in names(beta_diffs)) {
        # With functional normalisation, some probes are lost
        df <- data.frame(over_10_percent=sapply(cgdiffs[[correction]], function(x) sum(x > 10^-1)),
                         over_10_typeI=sapply(cgdiffs[[correction]], function(x) sum(x[names(x) %in% typeI_probe_names] > 10^-1)),
                         over_10_typeII=sapply(cgdiffs[[correction]], function(x) sum(x[names(x) %in% typeII_probe_names] > 10^-1)),
                         under_1_percent=sapply(cgdiffs[[correction]], function(x) sum(x < 10^-2)),
                         under_1_typeI=sapply(cgdiffs[[correction]], function(x) sum(x[names(x) %in% typeI_probe_names] < 10^-2)),
                         under_1_typeII=sapply(cgdiffs[[correction]], function(x) sum(x[names(x) %in% typeII_probe_names] < 10^-2)),
                         mean=round(sapply(cgdiffs[[correction]], mean), 4),
                         quantile_15=round(sapply(cgdiffs[[correction]], quantile, 0.15), 4),
                         quantile_25=round(sapply(cgdiffs[[correction]], quantile, 0.25), 4),
                         median=round(sapply(cgdiffs[[correction]], median), 4),
                         quantile_75=round(sapply(cgdiffs[[correction]], quantile, 0.75), 4),
                         quantile_85=round(sapply(cgdiffs[[correction]], quantile, 0.85), 4)
        )
        
        df$over10_propTypeI <- round(df$over_10_typeI/df$over_10_percent, 4)
        df$under1_propTypeI <- round(df$under_1_typeI/df$under_1_percent, 4)
        
        write.csv(df, file=file.path(path_results, paste(prefix_alt, correction, "Table_Betadiff.csv", sep="_")))
        rm(df)
    }
    rm(beta_diffs, cgdiffs, probes, typeI_probe_names, typeII_probe_names)
    gc()
}

# plots

log_diffs_lims <- list()
for(prefix in names(log_diffs)) {
    for(correction in names(log_diffs[[prefix]])) {
    log_diffs_lims[[prefix]][[correction]]$x <- range(sapply(log_diffs[[prefix]][[correction]], function(x) range(x$x)))
    log_diffs_lims[[prefix]][[correction]]$y <- range(sapply(log_diffs[[prefix]][[correction]], function(x) range(x$y)))
    }
}

my_xlim <- c(-3.5, 0.25)
my_ylim <- list("450K"=c(0, 1.7), "EPIC"=c(0, 2.8))

pdf(file=file.path(path_results, "Figure0X_Betadiff_densities_new3.pdf"))
par(mfcol=c(2,2))
for(my_array in names(log_diffs)) {
    for(my_correction in names(log_diffs[[my_array]])) {
        
        my_legend_placement <- "topleft"
        if(my_correction == "combat") my_legend_placement <- NULL
        
        plotLogDiffs(log_diffs, array=my_array, correction=my_correction,
                     this=c("raw", "bmiq", "swan", "noobswan", "noob", "noobbmiq", "enmix"),
                     xlim=my_xlim,
                     ylim=my_ylim[[my_array]],
                     main=paste(nice_names[my_correction], my_array), legend_placement=my_legend_placement, legend_cex=0.7)
        
        plotLogDiffs(log_diffs, array=my_array, correction=my_correction,
                     this=c("raw", "illumina", "func", "dasen", "noobdasen", "enmixquantile"),
                     xlim=my_xlim,
                     ylim=my_ylim[[my_array]],
                     main=paste(nice_names[my_correction], my_array), legend_placement=my_legend_placement, legend_cex=0.7)
    }
}
dev.off()


my_xlim <- c(-2.5, 0.17)
my_ylim <- c(0, 2.1)


pdf(file=file.path(path_results, "Figure0X_Betadiff_densities_raw_typeIandII.pdf"))
par(mfrow=c(2,2))
for(my_array in names(log_diffs)) {
    for(my_correction in names(log_diffs[[my_array]])) {
        my_legend_placement <- "topleft"
        if(my_correction == "combat") my_legend_placement <- NULL

        plotLogDiffs(log_diffs, array=my_array, correction=my_correction,
                     this=c("raw", "raw_typeI", "raw_typeII"),
                     xlim=my_xlim,
                     ylim=my_ylim,
                     main=paste(nice_names[my_correction], my_array), legend_placement=my_legend_placement, legend_cex=0.7)
    }
}
dev.off()


q("no")