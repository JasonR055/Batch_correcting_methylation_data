rm(list=ls())

library(RColorBrewer)
library(easyPubMed)

#####  Paths  #####

path_base <- "/home/ros259/R_projects/450K_Batch_effects"
path_data <- file.path(path_base, "data")
path_save_rdata <- file.path(path_base, "data", "rdata")
path_results <- file.path(path_base, "results")
path_idats <- file.path(path_base, "data", "idats")
path_code <- file.path(path_base, "scripts", "R")
path_snps <- "/home/ros259/R_projects/utility_scripts"

#####  Functions  #####

publication_plot <- function(myprobes, these = prefix_alts) {
    
    toplot <- function(p, this, correction, cex=0.5, title="") {
        
        myorder <- order(betas[[this]]$pd$array_num, paste(betas[[this]]$pd$row, betas[[this]]$pd$col))
        
        plot(betas[[this]][[correction]][p, myorder], cex=cex, cex.main=0.7,
             ylim=c(0,1), xaxt="n",
             las=2, cex.axis=0.6,
             xlab="", ylab="",
             col=betas[[this]]$pd$batch_col[myorder], main=title) 
    }
    
    myiqrs <- list()
    mymodal <- list()
    
    for(prefix in these) {
        myiqrs[[prefix]] <- betas[[prefix]]$iqrs[betas[[prefix]]$iqrs$is_biggest == TRUE, ]
        mymodal[[prefix]] <- betas[[prefix]]$modal
    }
    
    for(p in myprobes) {

        for(prefix in these) {
            cex=0.5
            
            if(prefix %in% c("EPIC-Italy", "NOVI")) cex=0.3
            idx <- match(p, myiqrs[[prefix]]$probe)
            toplot(p, prefix, correction="original", cex=cex,
                   title=paste(p, ": ", prefix, ",",
                               " SD=", round(var(betas[[prefix]]$original[idx, ]), 3), sep=""))
            toplot(p, prefix, correction="harman", cex=cex,
                   title=paste("Harman,",
                               " SD=", round(var(betas[[prefix]]$harman[idx, ]), 3),
                               ", LVR=", round(log2(myiqrs[[prefix]]$var_ratio_harman[idx]), 2),
                               ", Shift=", round(betas[[prefix]]$meandiffs_harman[idx], 4), sep=""))  # \U0394x\U0304
            toplot(p, prefix, correction="combat", cex=cex,
                   title=paste("ComBat,",
                               " SD=", round(var(betas[[prefix]]$harman[idx, ]), 3),
                               ", LVR=", round(log2(myiqrs[[prefix]]$var_ratio_combat[idx]), 2),
                               ", Shift=", round(betas[[prefix]]$meandiffs_combat[idx], 4), sep=""))
        }
    }
}


#####  Load files  #####

prefix_alts <- c("EpiSCOPE", "EPIC-Italy", "BodyFatness", "NOVI", "URECA")
names(prefix_alts) <- c("GSE89278", "GSE51032", "EPIC", "GSE128821", "GSE132181")

load(file=file.path(path_save_rdata, "EpiSCOPE_and_BodyFatness_noob_beta_diffs.RData"))
all_betas <- betas
rm(betas)

for(prefix in c("EPIC-Italy", "NOVI", "URECA")) {
    load(file=file.path(path_save_rdata, paste(prefix, "noob_beta_diffs.RData", sep="_")))
    all_betas[[prefix]] <- betas[[prefix]]
    rm(betas)
}

betas <- all_betas
rm(all_betas)
gc()


load(file=file.path(path_results, "bad_ps.RData"))


#####  EWAS Pubmed hits  #####

# Takes hours to run, so keep the running conditional
these_ps <- bad_ps
sink(file = file.path(path_results, "Bad_offenders_pubmed_hits.txt"))
for(i in 1:length(these_ps)) {
    cat(i, these_ps[i], "\n")
    my_query <- paste(these_ps[i], "AND EWAS[Title/Abstract]")
    my_entrez_id <- get_pubmed_ids(my_query)
    if(as.integer(my_entrez_id$Count) > 0) {
        my_abstracts_txt <- fetch_pubmed_data(my_entrez_id, format = "abstract")
        cat(rep("*", 30), "\n    ", these_ps[i], "\n", rep("*", 30), "\n", sep="")
        print(my_abstracts_txt)
    }
}
sink()


#####  Make figures  #####

dodgy_ewas <- c("cg11963436", "cg18368637", "cg22385669")

pdf(file.path(path_results, paste("bad_ewas_figure.pdf", sep="_")))
par(mfrow=c(5, 3),  mar=c(0.1, 1.8, 1.0, 0.1))
publication_plot(dodgy_ewas)
dev.off()

