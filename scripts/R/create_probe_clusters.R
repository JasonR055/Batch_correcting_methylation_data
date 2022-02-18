#rm(list=ls())

#library(RColorBrewer)
#library(minfi)
#library(FDb.InfiniumMethylation.hg19)
#library(Ckmeans.1d.dp)
#library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
#library(BSgenome.Hsapiens.UCSC.hg19)


#####  Paths  #####

path_base <- "/home/ros259/R_projects/450K_Batch_effects"
path_data <- file.path(path_base, "data")
path_save_rdata <- file.path(path_base, "data", "rdata")
path_results <- file.path(path_base, "results")
path_idats <- file.path(path_base, "data", "idats")
path_code <- file.path(path_base, "scripts", "R")
#path_results <- file.path(path_base, "results", "GEO")
path_snps <- "/home/ros259/R_projects/utility_scripts"


#####  Functions  #####

testCluster <- function(b, ks, min_cluster_size=5, min_cluster_dist=0.1) {
    
    ckset <- list(Ckmeans.1d.dp(b, ks[1]),
                  Ckmeans.1d.dp(b, ks[2])) # y = myweights
    names(ckset) <- ks
    # Are there at least min_cluster_size number of samples in a cluster
    has_min_n <- sapply(ckset, function(x) min(x$size) >= min_cluster_size)
    
    # If more than one cluster, the min gap must be at least min_cluster_dist
    has_min_gap <- sapply(ckset, function(x) {
        if(length(x$centers) > 1) {
            min(diff(x$centers)) >= min_cluster_dist
        } else {
            TRUE
        }
    })
    #ckset[[as.character(max(ks[has_min_n & has_min_gap]))]]
    # cg05230888 returns -Inf, so build edge case for this
    max(1, max(ks[has_min_n & has_min_gap]))
}


#####  Load files  #####

#dbsnp <- "snp151"
#load(file=file.path(path_snps, paste(dbsnp, "Common.RData", sep="")))

#xrct_file <- file.path(path_data, "48639-non-specific-probes-Illumina450k.csv")
#xrct_probes <- read.csv(file=xrct_file, stringsAsFactors=FALSE)
#rm(xrct_file)


#####  Load RData and form diffs  #####

betas <- list()
#prefix_alts <- c("EpiSCOPE", "EPIC-Italy", "BodyFatness", "NOVI", "URECA")
#names(prefix_alts) <- c("GSE89278", "GSE51032", "EPIC", "GSE128821", "GSE132181")

#prefix_alts <- c("450K", "Italy", "EPIC", "UCL")
#names(prefix_alts) <- c("GSE89278", "GSE51032", "EPIC", "GSE111631")


for(i in 1:length(prefix_alts)) {
    prefix_alt <- prefix_alts[i]
    load(file=file.path(path_save_rdata, paste(names(prefix_alt), "MSet.noob.RData", sep="_")))
    
    betas[[prefix_alt]] = list(pd = pData(MSet.noob),
                               original = ilogit2(noob_m),
                               harman = ilogit2(noob_harman_m),
                               combat = ilogit2(noob_combat_m))
    
    rm(noob_combat_m, noob_harman_m, noob_m, MSet.noob, MSet.noob_hm)
    gc()
    betas[[prefix_alt]]$diffs_harman = betas[[prefix_alt]]$harman - betas[[prefix_alt]]$original
    betas[[prefix_alt]]$diffs_combat = betas[[prefix_alt]]$combat - betas[[prefix_alt]]$original
}


#####  Map to SNPs  #####

for(prefix_alt in prefix_alts) {
    slide <- betas[[prefix_alt]]$pd$array_num
    
    if(prefix_alt == "EpiSCOPE") {
        InfiniumMethylation <- features(FDb.InfiniumMethylation.hg19)
        InfiniumMethylation <- InfiniumMethylation[rownames(betas[[prefix_alt]]$original)]
        
        bad_batches <- c(1, 5, 9, 17, 25)
        super <- c(bad_batches, max(betas[[prefix_alt]]$pd$array_num))
        superbatches <- rep(NA_character_, nrow(betas[[prefix_alt]]$pd))
        for(i in 1:(length(super)-1)) {
            superbatches[betas[[prefix_alt]]$pd$array_num %in% super[i]:super[i+1]] <- LETTERS[i]
        }
        rm(super)
        is_bad_sample <- slide %in% bad_batches
        bad_cols <- c("green", "red")[factor(is_bad_sample)]
        # Slide 30 is marginal whether it's a bad batch or not
        is_good_sample <- !is_bad_sample
    }
    
    if(prefix_alt == "EPIC-Italy") {
        InfiniumMethylation <- features(FDb.InfiniumMethylation.hg19)
        InfiniumMethylation <- InfiniumMethylation[rownames(betas[[prefix_alt]]$original)]
                
        #bad_batches <- c(29, 32, 35, 39, 43, 44, 47, 52, 53, 55, 56, 62, 63, 66, 68)
        bad_batches <- c(29, 32, 36, 39, 41, 43, 44, 45, 47, 50, 52, 54, 56, 59, 62, 63, 64, 66, 68)
        # Many slides after batch 28 are far worse in quality
        superbatches <- as.character(factor(LETTERS)[as.integer(betas[[prefix_alt]]$pd$array_num <= 28) + 1])
        is_bad_sample <- slide %in% bad_batches
        bad_cols <- c("green", "red")[factor(is_bad_sample)]
        is_good_sample <- !is_bad_sample
    }
    
    if(prefix_alt == "BodyFatness") {
        data("Locations")
        # Get EPIC location data
        # Make GRanges
        # Width 2 as CpG sites are 2 nts
        InfiniumMethylation = GRanges(seqnames=IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations$chr,
                                      ranges=IRanges(start=IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations$pos, width=2),
                                      strand = IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations$strand,
                                      seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19))
        names(InfiniumMethylation) = rownames(IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations)
        # Remove the probes that went missing from B2 to B4 array
        InfiniumMethylation <- InfiniumMethylation[rownames(betas[[prefix_alt]]$original)]
        bad_batches <- 3
        # The second plate is batchy
        superbatches <- as.character(factor(LETTERS)[betas[[prefix_alt]]$pd$plate])
        is_bad_sample <- slide %in% bad_batches
        bad_cols <- c("green", "red")[factor(is_bad_sample)]
        is_good_sample <- !is_bad_sample
                
        col_factor <- factor(paste(betas[[prefix_alt]]$pd$sentrix_barcode, betas[[prefix_alt]]$pd$gender))
        batch_cols <- rep(brewer.pal(10, 'Paired'), 5)[1:length(levels(col_factor))]
        betas[[prefix_alt]]$pd$batch_col <- batch_cols[col_factor]
        rm(batch_cols, col_factor, bad_cols)
    }
    
    if(prefix_alt == "UCL") {
        data("Locations")
        # Get EPIC location data
        # Make GRanges
        # Width 2 as CpG sites are 2 nts
        InfiniumMethylation = GRanges(seqnames=IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations$chr,
                                      ranges=IRanges(start=IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations$pos, width=2),
                                      strand = IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations$strand,
                                      seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19))
        names(InfiniumMethylation) = rownames(IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations)
        # Remove the probes that went missing from B2 to B4 array
        InfiniumMethylation <- InfiniumMethylation[rownames(betas[[prefix_alt]]$original)]
        bad_batches <- c(3, 7)
        # There are two plates
        superbatches <- as.character(factor(LETTERS)[betas[[prefix_alt]]$pd$sample_plate])
        is_bad_sample <- slide %in% bad_batches
        bad_cols <- c("green", "red")[factor(is_bad_sample)]
        is_good_sample <- !is_bad_sample
    }
    
    
    if(prefix_alt == "NOVI") {
        data("Locations")
        # Get EPIC location data
        # Make GRanges
        # Width 2 as CpG sites are 2 nts
        InfiniumMethylation = GRanges(seqnames=IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations$chr,
                                      ranges=IRanges(start=IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations$pos, width=2),
                                      strand = IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations$strand,
                                      seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19))
        names(InfiniumMethylation) = rownames(IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations)
        # Remove the probes that went missing from B2 to B4 array
        InfiniumMethylation <- InfiniumMethylation[rownames(betas[[prefix_alt]]$original)]
        #bad_batches <- c(9, 22, 33, 36, 37, 40, 41, 47, 53, 60, 62, 64, 70, 75)
        bad_batches <- c(2, 3, 4, 9, 16, 17, 19, 20, 28, 29, 30, 32, 35, 39, 40, 45, 47, 50, 53, 60, 69, 70, 72, 73, 75)
        # There are seven plates
        superbatches <- betas[[prefix_alt]]$pd$plate
        is_bad_sample <- slide %in% bad_batches
        bad_cols <- c("green", "red")[factor(is_bad_sample)]
        is_good_sample <- !is_bad_sample
    }
    
    
    if(prefix_alt == "URECA") {
        data("Locations")
        # Get EPIC location data
        # Make GRanges
        # Width 2 as CpG sites are 2 nts
        InfiniumMethylation = GRanges(seqnames=IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations$chr,
                                      ranges=IRanges(start=IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations$pos, width=2),
                                      strand = IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations$strand,
                                      seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19))
        names(InfiniumMethylation) = rownames(IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations)
        # Remove the probes that went missing from B2 to B4 array
        InfiniumMethylation <- InfiniumMethylation[rownames(betas[[prefix_alt]]$original)]
        #bad_batches <- c(9, 22, 33, 36, 37, 40, 41, 47, 53, 60, 62, 64, 70, 75)
        #bad_batches <- c(2, 3, 4, 9, 16, 17, 19, 20, 28, 29, 30, 32, 35, 39, 40, 45, 47, 50, 53, 60, 69, 70, 72, 73, 75)
        bad_batches <- c(24, 25, 28)
        # There are five plates
        superbatches <- betas[[prefix_alt]]$pd$plate
        is_bad_sample <- slide %in% bad_batches
        bad_cols <- c("green", "red")[factor(is_bad_sample)]
        is_good_sample <- !is_bad_sample
    }
    
    # sanity check
    stopifnot(identical(rownames(betas[[prefix_alt]]$pd), colnames(betas[[prefix_alt]]$original)))
    # Now add
    betas[[prefix_alt]]$pd$is_bad_sample <- is_bad_sample
    betas[[prefix_alt]]$pd$is_good_sample <- is_good_sample
    betas[[prefix_alt]]$pd$superbatches <- superbatches
    
    ol <- findOverlaps(InfiniumMethylation, dbsnp_gr, select = "first", ignore.strand=TRUE)
    InfiniumMethylation$snp151common[!is.na(ol)] <- dbsnp_gr$rsid[ol[!is.na(ol)]]
    InfiniumMethylation$observed[!is.na(ol)] <- dbsnp_gr$observed[ol[!is.na(ol)]]
    InfiniumMethylation$avHet[!is.na(ol)] <- dbsnp_gr$avHet[ol[!is.na(ol)]]
    
    # How far is the nearest SNP from the CpG site under inspection?
    cpg_dist <- as.data.frame(distanceToNearest(InfiniumMethylation, dbsnp_gr, select="arbitrary", ignore.strand=TRUE))
    InfiniumMethylation$distanceSnp[cpg_dist$queryHits] <- cpg_dist$distance
    InfiniumMethylation$distanceRsid[cpg_dist$queryHits] <- dbsnp_gr$rsid[cpg_dist$subjectHits]

    betas[[prefix_alt]]$ranges <- InfiniumMethylation
    
    save(InfiniumMethylation, file=file.path(path_save_rdata, paste(prefix_alt, "InfiniumMethylation", dbsnp, "overlap.RData", sep="_")))
    
    rm(ol, InfiniumMethylation, cpg_dist, slide)
}

rm(superbatches, i, prefix_alt, is_bad_sample, is_good_sample)
gc()


#####  Modal probes - Ckmeans  #####

myks <- 1:10
num_clusters <- list()
num_opt_clusters <- list()

for(prefix_alt in prefix_alts) {
    my_cores = 12
    
    # Find bi- or tri-modal methylation patterns suggestive of SNPs
    cat("First clustering for", prefix_alt)
    bs <- betas[[prefix_alt]]$original[, betas[[prefix_alt]]$pd$is_good_sample]
    
    unbiased_clusters <- mclapply(1:nrow(bs), function(p) {
        Ckmeans.1d.dp(sort(bs[p, ]), k=myks)
    },
    mc.cores = my_cores)
    names(unbiased_clusters) <- rownames(betas[[prefix_alt]]$original)
    num_clusters[[prefix_alt]] <- sapply(unbiased_clusters, function(x) length(x$centers))
    
    # Cycle over clusters and find the best two options by BIC
    ktext_to_int <- myks
    names(ktext_to_int) <- names(unbiased_clusters[[1]]$BIC)
    cat(", find top two k")
    # get the two top ks by BIC
    top_ks <- lapply(unbiased_clusters[num_clusters[[prefix_alt]] > 1], function(x) {
        bic <- x$BIC
        bic <- bic[!is.nan(bic)]
        bic_diffs <- diff(c(0, bic))
        highest_bic_diffs <- sort(bic_diffs, decreasing = TRUE)[1:2]
        as.integer(ktext_to_int[names(highest_bic_diffs)])
    })
    rm(ktext_to_int)
    
    # Recluster for multiple cluster probes and specify a distance and number cutoff
    cat(", done. Reclustering", prefix_alt)
    bs_small <- bs[names(top_ks), ]
    best_k <- mclapply(1:length(top_ks), function(i) {
        testCluster(sort(bs_small[i, ]), top_ks[[i]])
    }, mc.cores=my_cores)
    names(best_k) <- names(top_ks)
    best_k <- do.call("c", best_k)
    cat(", done.\n")
    
    num_best_clusters <- num_clusters[[prefix_alt]]
    num_best_clusters[names(best_k)] <- best_k
    num_opt_clusters[[prefix_alt]] <- num_best_clusters
    rm(bs, bs_small, best_k, top_ks, num_best_clusters)
    gc()
}


table(num_clusters$`450K`)
# Is good samples
#     1      2      3      4      5      6      7      8      9     10 
#386210  54364  35982   5548   1888    788    378    233     90     31 
table(num_opt_clusters$`450K`)
#     1      2      3      4      5      6 
#463231  19816   2379     78      7      1
table(num_clusters$Italy)
#     1      2      3      4      5      6      7      8      9     10 
#265355  90493  94564  26836   5169   1603    696    419    240    137 
table(num_opt_clusters$Italy)
#     1      2      3      4      5      6      7 
#456184  26601   2636     75     13      2      1 
table(num_clusters$`EPIC`)
#     1      2      3      4      5      6      7      8      9     10 
#655542 133504  59469  13269   2821    895    347    160     62     22 
table(num_opt_clusters$EPIC)
#     1      2      3      4      5 
#785339  76041   4588    112     11 
table(num_clusters$Buccal)
#     1      2      3      4      5      6      7      8      9     10 
#360156 298672 162104  35400   6354   1718    684    397    247    127 
table(num_opt_clusters$Buccal)
#     1      2      3      4      5      6 
#698809 162070   4527    388     60      5 

# Set any optimal cluster number greater than 4 to 4.
for(prefix_alt in prefix_alts) {
    num_opt_clusters[[prefix_alt]][num_opt_clusters[[prefix_alt]] > 4] <- 4
    write.csv(as.data.frame(table(num_opt_clusters[[prefix_alt]])),
              file=file.path(path_results, paste(prefix_alt, "Probe_Ckmeans_n5.csv", sep="_")),
              row.names = FALSE)
}


#####  Recluster all samples (not just good samples) with the optimal k  #####

for(prefix_alt in prefix_alts) {
    
    # For probes with greater than 1 cluster, redo clustering with a set k and
    # inclusion of all samples (batchy or not)
    ps <- names(num_opt_clusters[[prefix_alt]][num_opt_clusters[[prefix_alt]] > 1])
    ks <- num_opt_clusters[[prefix_alt]][ps]
    bs <- betas[[prefix_alt]]$original[ps, ]
    
    stopifnot(identical(names(ks), rownames(bs)))
    
    x <- mclapply(1:length(ps), function(i) {
        Ckmeans.1d.dp::Ckmeans.1d.dp(bs[i, ], k=ks[i])[c("cluster", "size")]
    }, mc.cores=12)
    names(x) <- ps
    
    myclusters <- t(sapply(x, "[[", "cluster"))
    colnames(myclusters) <- colnames(bs)
    
    stopifnot(identical(rownames(myclusters), rownames(bs)))
    
    # Pre-fill ones matrix
    clusters <- matrix(1,
                       nrow = nrow(betas[[prefix_alt]]$original),
                       ncol = ncol(betas[[prefix_alt]]$original),
                       dimnames = dimnames(betas[[prefix_alt]]$original))
    # Stuff in values
    clusters[rownames(myclusters), ] <- myclusters
    betas[[prefix_alt]]$clusters <- clusters
    
    mycluster_props <- lapply(x, "[[", "size")
    num <- ncol(clusters)
    # Pref-fill with the number of samples
    cluster_props <- vector("list", nrow(clusters))
    cluster_props <- lapply(cluster_props, function(x) num)
    names(cluster_props) <- rownames(clusters)
    # Stuff in values
    cluster_props[names(mycluster_props)] <- mycluster_props
    betas[[prefix_alt]]$cluster_props <- cluster_props
    
    modal <- data.frame(probe_id=ps,
                        chr=as.character(seqnames(betas[[prefix_alt]]$ranges[ps])),
                        start=start(betas[[prefix_alt]]$ranges[ps]),
                        end=end(betas[[prefix_alt]]$ranges[ps]),
                        strand=strand(betas[[prefix_alt]]$ranges[ps]),
                        num_opt_clusters=ks,
                        commonSNP=betas[[prefix_alt]]$ranges[ps]$snp151common,
                        observed_alleles=betas[[prefix_alt]]$ranges[ps]$observed,
                        avHet=betas[[prefix_alt]]$ranges[ps]$avHet,
                        distanceSnp=betas[[prefix_alt]]$ranges[ps]$distanceSnp,
                        distanceRsid=betas[[prefix_alt]]$ranges[ps]$distanceRsid,
                        stringsAsFactors = FALSE)
    betas[[prefix_alt]]$modal <- modal
    
    rm(ps, ks, bs, x, clusters, myclusters, modal, mycluster_props, cluster_props)
}

rm(num_clusters, num_opt_clusters, unbiased_clusters)


#####  Dispersion measures after batch correction  #####

sum_of_squares <- function(x) {
    x <- x[!is.na(x)]
    # the 'classic' approach
    sum((x - mean(x))^2)
}


for(prefix_alt in prefix_alts) {
    
    ps <- rownames(betas[[prefix_alt]]$clusters)
    
    iqr_clust = mclapply(1:nrow(betas[[prefix_alt]]$original), function(i) {
    #iqr_clust = lapply(1:nrow(betas[[prefix_alt]]$original), function(i) {
        # for each probe, find the interquartile range fold-change for each cluster
        # Find cluster with the largest IQRs after correction
        print(i)
        this_cluster <- betas[[prefix_alt]]$clusters[i, ]
        cluster_table <- betas[[prefix_alt]]$cluster_props[[i]]
        # check if any of the clusters is 1 in size, so it will have no variance and IQR
        is_one <- cluster_table == 1
        if(sum(is_one) > 0) {
            # Stuff into cluster 1
            this_cluster[this_cluster == as.integer(names(is_one)[is_one])] <- 1
            cluster_table <- as.vector(table(this_cluster))
        }
        harman_b <- betas[[prefix_alt]]$harman[i, ]
        combat_b <- betas[[prefix_alt]]$combat[i, ]
        original_b <- betas[[prefix_alt]]$original[i, ]
        harman_diff <- abs(harman_b - original_b)
        combat_diff <- abs(combat_b - original_b)
        
        original <- tapply(original_b, INDEX = this_cluster, stats::IQR,  na.rm=TRUE)
        df <- data.frame("probe"=ps[i],
                         "total_clusters"=length(original),
                         "clust"=as.integer(names(original)),
                         "num"=cluster_table,
                         "original_iqr"=as.vector(original),
                         "harman_iqr"=as.vector(tapply(harman_b, INDEX = this_cluster, stats::IQR,  na.rm=TRUE)),
                         "combat_iqr"=as.vector(tapply(combat_b, INDEX = this_cluster, stats::IQR,  na.rm=TRUE)),
                         "harman_diff"=as.vector(tapply(harman_diff, INDEX = this_cluster, mean,  na.rm=TRUE)),
                         "combat_diff"=as.vector(tapply(combat_diff, INDEX = this_cluster, mean,  na.rm=TRUE)),
                         "harman_diff_iqr"=as.vector(tapply(harman_diff, INDEX = this_cluster, stats::IQR,  na.rm=TRUE)),
                         "combat_diff_iqr"=as.vector(tapply(combat_diff, INDEX = this_cluster, stats::IQR,  na.rm=TRUE)),
                         "original_ssq"=tapply(original_b, INDEX = this_cluster, sum_of_squares),
                         "harman_ssq"=tapply(harman_b, INDEX = this_cluster, sum_of_squares),
                         "combat_ssq"=tapply(combat_b, INDEX = this_cluster, sum_of_squares),
                         "is_biggest"=FALSE,
                         "is_maxiqr_harman"=FALSE,
                         "is_maxiqr_combat"=FALSE,
                         "is_maxdiff_harman"=FALSE,
                         "is_maxdiff_combat"=FALSE,
                         "is_maxdiffiqr_harman"=FALSE,
                         "is_maxdiffiqr_combat"=FALSE,
                       stringsAsFactors = FALSE)
        df$iqr_ratio_harman <- as.vector(df$harman_iqr / df$original_iqr)
        df$iqr_ratio_combat <- as.vector(df$combat_iqr / df$original_iqr)
        # the sum of squares divided by the number of observations, less one, for sample variance
        df$original_var <- sum(df$original_ssq) / (sum(df$num) - 1)
        df$harman_var <- sum(df$harman_ssq) / (sum(df$num) - 1)
        df$combat_var <- sum(df$combat_ssq) / (sum(df$num) - 1)
        df$var_ratio_harman <- as.vector(df$harman_var / df$original_var)
        df$var_ratio_combat <- as.vector(df$combat_var / df$original_var)
        is_biggest <- as.vector(cluster_table == max(cluster_table))
        if(sum(is_biggest) > 1) {
            # case with 2 clusters or more, with at least 2 of them equal to max size.
            # Record 'is_biggest' as the largest cluster with the highest IQR
            is_biggest[is_biggest] <- df$original_iqr[is_biggest] == max(df$original_iqr[is_biggest])
            stopifnot(sum(is_biggest) == 1)
        }
        df$is_biggest <- is_biggest
        df$is_maxiqr_harman[df$iqr_ratio_harman == max(df$iqr_ratio_harman)] <- TRUE
        df$is_maxiqr_combat[df$iqr_ratio_combat == max(df$iqr_ratio_combat)] <- TRUE
        df$is_maxdiff_harman[df$harman_diff == max(df$harman_diff)] <- TRUE
        df$is_maxdiff_combat[df$combat_diff == max(df$combat_diff)] <- TRUE
        df$is_maxdiffiqr_harman[df$harman_diff_iqr == max(df$harman_diff_iqr)] <- TRUE
        df$is_maxdiffiqr_combat[df$combat_diff_iqr == max(df$combat_diff_iqr)] <- TRUE
        df
    }, mc.cores = 8)
    #})
    rm(ps)
    
    iqrs <- do.call(rbind, iqr_clust)
    rownames(iqrs) <- paste(iqrs$probe, iqrs$clust, sep=".")
    #iqrs$iqr_ratio_harman <- as.vector(iqrs$harman_iqr / iqrs$original_iqr)
    #iqrs$iqr_ratio_combat <- as.vector(iqrs$combat_iqr / iqrs$original_iqr)
    #iqrs$ssq_ratio_harman <- as.vector(iqrs$harman_ssq / iqrs$original_ssq)
    #iqrs$ssq_ratio_combat <- as.vector(iqrs$combat_ssq / iqrs$original_ssq)
    #iqrs$ssq_diff_harman <- as.vector(iqrs$harman_ssq - iqrs$original_ssq)
    #iqrs$ssq_diff_combat <- as.vector(iqrs$combat_ssq - iqrs$original_ssq)
    #iqrs$iqr_logFC_harman <- log2(iqrs$harman_iqr) - log2(iqrs$original_iqr)
    #iqrs$iqr_logFC_combat <- log2(iqrs$combat_iqr) - log2(iqrs$original_iqr)
    
    betas[[prefix_alt]]$iqrs <- iqrs
    rm(iqrs, iqr_clust, prefix_alt)
    gc()
}


for(prefix_alt in prefix_alts) {
    for(this in c("harman", "combat")) {
        this_diff <- paste("diffs", this, sep="_")
        betas[[prefix_alt]][[paste("meandiffs", this, sep="_")]] <- rowMeans(abs(betas[[prefix_alt]][[this_diff]]))
        cgranges <- matrixStats::rowRanges(betas[[prefix_alt]][[this_diff]], na.rm = TRUE)
        colnames(cgranges) <- c('min', 'max')
        cgdiffs <- as.vector(matrixStats::rowDiffs(cgranges))
        names(cgdiffs) <- rownames(betas[[prefix_alt]][[this_diff]])
        betas[[prefix_alt]][[paste("maxdiffs", this, sep="_")]] <- cgdiffs
        rm(cgranges, cgdiffs)
    }
}


#####  Fisher tests for clusters  #####

tests <- list()

for(prefix_alt in prefix_alts) {
    # prefix_alt=prefix_alts[1]
    cat("Establish total clusters for", prefix_alt)
    total_clusters <- sapply(betas[[prefix_alt]]$cluster_props, length)
    total_clusters <- total_clusters[total_clusters > 1]
    clust2or3 <- total_clusters[total_clusters %in% c(2, 3)]
    is_autosomal <- as.logical(!(seqnames(betas[[prefix_alt]]$ranges) %in% c("chrX", "chrY")))
    autosomal_clust2or3 <- intersect(names(betas[[prefix_alt]]$ranges)[is_autosomal], names(clust2or3))
    
    cluster_indices <- match(names(total_clusters), rownames(betas[[prefix_alt]]$clusters))
    #stopifnot(identical(names(betas[[prefix_alt]]$ranges), rownames(betas[[prefix_alt]]$clusters)))    
    #autosomal_clust2or3_indices <- match(autosomal_clust2or3, names(betas[[prefix_alt]]$cluster_props))
    
    cat(", now test association with sex, batch, superbatch, cell composition")
    # sex
    if(prefix_alt == "EpiSCOPE") sex <- factor(betas[[prefix_alt]]$pd$sex)
    if(prefix_alt %in% c("EPIC-Italy", "BodyFatness", "NOVI", "URECA")) sex <- factor(betas[[prefix_alt]]$pd$gender)

    batchy <- factor(betas[[prefix_alt]]$pd$is_bad_sample)
    superbatchy <- factor(betas[[prefix_alt]]$pd$superbatches)
    # cell composition
    if(prefix_alt %in% c("EpiSCOPE", "EPIC-Italy", "URECA")) cell_comp <- betas[[prefix_alt]]$pd$Neutro
    if(prefix_alt %in% c("BodyFatness", "NOVI")) cell_comp <- betas[[prefix_alt]]$pd$prop_IC
    
    x <- mclapply(cluster_indices, function(i) {
        cluster_factor <- factor(betas[[prefix_alt]]$clusters[i, ])
        res_cell <- lm(cell_comp ~ betas[[prefix_alt]]$clusters[i, ])
        # handle higher dimension contingency tables with ChiSquare and simulated p-values
        if(prefix_alt %in% c("EpiSCOPE", "NOVI", "URECA")) {
            superbatch <- chisq.test(x=cluster_factor, y=superbatchy, simulate.p.value = TRUE, B=5e3)$p.value
        }
        # traditional Fisher Exact
        if(prefix_alt %in% c("EPIC-Italy", "BodyFatness")) {
            superbatch <- fisher.test(x=cluster_factor, y=superbatchy, hybrid=TRUE)$p.value
        }
        
        c(sex=fisher.test(x=cluster_factor, y=sex, hybrid=TRUE)$p.value,
          batch=fisher.test(x=cluster_factor, y=batchy, hybrid=TRUE)$p.value,
          superbatch=superbatch,
          cells=summary(res_cell)$coefficients[2, "Pr(>|t|)"])
    }, mc.cores=8)
    
    tests[[prefix_alt]] <- data.frame(do.call(rbind, x))
    rownames(tests[[prefix_alt]]) <- names(total_clusters) 
    rm(sex, batchy, cell_comp, x)
    
    tests[[prefix_alt]]$hwe <- NA_real_
    # Hardy-Weinberg
    cat(", Hardy Weinberg, ")
    hwe <- mclapply(autosomal_clust2or3, function(p) {
        prop <- betas[[prefix_alt]]$cluster_props[[p]]
        if(length(prop) == 2) {
            # Ensure that the intermediate methylation cluster stays in the middle
            if(prop[1] < prop[2]) {
                prop = c(0, prop)
            } else {
                prop = c(prop, 0)
            }
        }
        if(length(prop) == 3) {
            allele_order <- prop[c(1, 3)][order(prop[c(1, 3)], decreasing = TRUE)]
            # Keep cluster 2 as the het alleles, reorder clusters 1 and 3 such that
            # the largest proportion cluster is first in order
            prop <- c(allele_order[1], prop[2], allele_order[2])
        }
        names(prop) <- c("AA", "AB", "BB")
        #HWChisq(prop, cc=0, verbose=FALSE)
        HardyWeinberg::HWExact(prop, verbose=FALSE)$pval
    }, mc.cores=8)
    cat("done.\n")
    tests[[prefix_alt]][autosomal_clust2or3, "hwe"] <- unlist(hwe)
    rm(is_autosomal, autosomal_clust2or3, clust2or3)
}


#####  Add to modal  #####

cutoff <- 1e-3
for(prefix_alt in prefix_alts) {
    # sanity check
    stopifnot(identical(rownames(tests[[prefix_alt]]), rownames(betas[[prefix_alt]]$modal)))
    betas[[prefix_alt]]$modal$sex_p <- p.adjust(tests[[prefix_alt]]$sex, method="BH")
    betas[[prefix_alt]]$modal$batch_p <- p.adjust(tests[[prefix_alt]]$batch, method="BH")
    betas[[prefix_alt]]$modal$superbatch_p <- p.adjust(tests[[prefix_alt]]$superbatch, method="BH")
    betas[[prefix_alt]]$modal$cells_p <- p.adjust(tests[[prefix_alt]]$cells, method="BH")
    betas[[prefix_alt]]$modal$hwe_p <- p.adjust(tests[[prefix_alt]]$hwe, method="BH")
    betas[[prefix_alt]]$modal$xhyb_47 <- xrct_probes$X47[match(rownames(betas[[prefix_alt]]$modal), xrct_probes$TargetID)]
    betas[[prefix_alt]]$modal$xhyb_48 <- xrct_probes$X48[match(rownames(betas[[prefix_alt]]$modal), xrct_probes$TargetID)]
    betas[[prefix_alt]]$modal$xhyb_49 <- xrct_probes$X49[match(rownames(betas[[prefix_alt]]$modal), xrct_probes$TargetID)]
    betas[[prefix_alt]]$modal$xhyb_50 <- xrct_probes$X49[match(rownames(betas[[prefix_alt]]$modal), xrct_probes$TargetID)]
    betas[[prefix_alt]]$modal$is_xhyb <- rownames(betas[[prefix_alt]]$modal) %in% xrct_probes$TargetID
    betas[[prefix_alt]]$modal$is_batchy <- betas[[prefix_alt]]$modal$batch_p <= cutoff
    betas[[prefix_alt]]$modal$is_superbatchy <- betas[[prefix_alt]]$modal$superbatch_p <= cutoff
    betas[[prefix_alt]]$modal$is_imprinted <- betas[[prefix_alt]]$modal$sex_p <= cutoff
    betas[[prefix_alt]]$modal$is_snp <- !is.na(betas[[prefix_alt]]$modal$commonSNP)
    betas[[prefix_alt]]$modal$is_chrX <- betas[[prefix_alt]]$modal$chr=="chrX"
    betas[[prefix_alt]]$modal$is_chrY = betas[[prefix_alt]]$modal$chr == "chrY"
    betas[[prefix_alt]]$modal$is_snp_10bp = !is.na(betas[[prefix_alt]]$modal$distanceSnp) & betas[[prefix_alt]]$modal$distanceSnp > 0 & betas[[prefix_alt]]$modal$distanceSnp <=10
    betas[[prefix_alt]]$modal$in_hwe = !is.na(betas[[prefix_alt]]$modal$hwe_p) & betas[[prefix_alt]]$modal$hwe_p > cutoff
    betas[[prefix_alt]]$modal$is_celltype <- betas[[prefix_alt]]$modal$cells_p <= cutoff
}


save(betas, file=file.path(path_save_rdata, paste(paste(prefix_alts, collapse = "_and_"), "noob_beta_diffs.RData", sep="_")))

#q("no")