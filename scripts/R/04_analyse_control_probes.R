rm(list=ls())

library(RColorBrewer)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(reshape)
library(ggplot2)
library(ggpubr)


#####  Paths  #####

path_base <- "/home/ros259/R_projects/450K_Batch_effects"
path_data <- file.path(path_base, "data")
path_save_rdata <- file.path(path_base, "data", "rdata")
path_results <- file.path(path_base, "results")
path_peapod <- "/home/ros259/R_projects/SGA Infants"

# HTML RGB
w3c <- read.csv(file=file.path(path_data, "w3c_cols.csv"),
                    stringsAsFactors=FALSE)
html_cols <- w3c$HEX
names(html_cols) <- tolower(w3c$ColorName)
html_cols["-99"] <- "grey50"  # To account for missing colours in the probe info table
rm(w3c)


#####  Cex settings  #####

# Make consistent cex pattern

mycex <- c("axis"=0.6, "lab"=0.6, "main"=0.7, "cex"=0.1)

prefix_alts <- c("EpiSCOPE", "EPIC-Italy", "BodyFatness", "NOVI", "URECA")
names(prefix_alts) <- c("GSE89278", "GSE51032", "EPIC", "GSE128821", "GSE132181")

rg <- list()

for(prefix in prefix_alts) {
    print(prefix)
    rgset_file = names(prefix_alts)[prefix_alts == prefix]

    # load files
    if(prefix == "BodyFatness") {
        load(file=file.path(path_save_rdata, "SGA_Infants_RGSet.RData"))
    } else {
        load(file=file.path(path_save_rdata, paste(rgset_file, "RGSet.RData", sep="_")))
    }
    
    rg[[prefix]] <- list(r=getRed(RGSet),
                         g=getGreen(RGSet),
                         detp=detP[, pData(RGSet)$Basename],
                         pd=pData(RGSet),
                         control_df=as.data.frame(getProbeInfo(RGSet, type="Control")))
    
    rg[[prefix]]$control_grn <- getGreen(RGSet)[rg[[prefix]]$control_df$Address, ]
    rg[[prefix]]$control_red <- getRed(RGSet)[rg[[prefix]]$control_df$Address, ]
    
    rm(RGSet, detP)
}


sentrix_cols <- brewer.pal(8, "Set2")


for(prefix in prefix_alts) {
    
    if(prefix %in% c("EpiSCOPE", "EPIC-Italy")) {
        mypos <- paste("R", as.integer(rg[[prefix]]$pd$row), "C", as.integer(rg[[prefix]]$pd$col), sep="")
        rg[[prefix]]$pd$sentrix_description <- mypos
    }
    if(prefix == "BodyFatness") {
        mypos <- sub("R0(\\d)C\\d{2}", "R\\1", rg[[prefix]]$pd$sentrix_position)
        rg[[prefix]]$pd$sentrix_description <- mypos
    }
    if(prefix == "UCL") {
        mypos <- sub("0(\\d)", "R\\1", rg[[prefix]]$pd$row)
        rg[[prefix]]$pd$sentrix_description <- mypos
    }
    if(prefix %in% c("NOVI", "URECA")) {
        mypos <- sub("0(\\d)", "R\\1", rg[[prefix]]$pd$row)
        rg[[prefix]]$pd$sentrix_description <- mypos
    }

    rg[[prefix]]$control_df$Web_color <- html_cols[tolower(rg[[prefix]]$control_df$Color)]
    
    # Make summary stats
    rg[[prefix]]$pd$detected_probes_05 <- colSums(rg[[prefix]]$detp <= 0.05)
    rg[[prefix]]$pd$detected_probes_05_pct <- rg[[prefix]]$pd$detected_probes_05 / nrow(rg[[prefix]]$detp) * 100
    rg[[prefix]]$signal_avg_grn <- colMeans(rg[[prefix]]$g)
    rg[[prefix]]$signal_avg_red <- colMeans(rg[[prefix]]$r)
}


#####  Control probes  #####

ctrls <- c("STAINING", "TARGET REMOVAL",
           "HYBRIDIZATION", "EXTENSION",
           "SPECIFICITY I", "SPECIFICITY II",
           "BISULFITE CONVERSION I", "BISULFITE CONVERSION II",
           "NON-POLYMORPHIC", "NEGATIVE")


processChannel <- function(prefix, ctrl, colour) {

    ctrl_df <- rg[[prefix]]$control_df[rg[[prefix]]$control_df$Type == ctrl, ]
    this_channel <- paste("control", colour, sep="_")
    
    ch <- log2(rg[[prefix]][[this_channel]][ctrl_df$Address, ] + 0.5)
    nsamples <- ncol(ch)
    nctrls <- nrow(ch)
    ch <- melt(ch)
    names(ch) <- c("address", "id", "log_signal")
    if(colour == "red" & ctrl != "NEGATIVE") {
        ch_offset <- 0.2
    }
    if(colour == "grn" & ctrl != "NEGATIVE") {
        ch_offset <- -0.2
    }
    if(ctrl == "NEGATIVE") {
        ch_offset <- 0
    }
    ch$x <- rep(rg[[prefix]]$pd$array_num, each=nctrls) + ch_offset
    ch$type <- rep(ctrl_df$ExtendedType, nsamples)
    ch$channel <- colour
    ch$colour <- rep(ctrl_df$Web_color, nsamples)
    ch
}


for(prefix in prefix_alts) {
    
    print(table(rg[[prefix]]$control_df$Type))
    
    xtext = 8
    # make smaller for large datasets
    if(prefix %in% c("EPIC-Italy", "NOVI", "URECA")) {
        xtext = 3
    }
    
    gplots <- list()
    i=1
    for(this_ctrl in ctrls) {
        
        tidy <- rbind(processChannel(prefix, ctrl=this_ctrl, colour="red"),
                      processChannel(prefix, ctrl=this_ctrl, colour="grn"))
        
        if(this_ctrl != "NEGATIVE") {
            myplot <- ggplot(data = tidy, mapping = aes(x=x, y=log_signal, colour=type))
            myplot <- myplot + geom_vline(xintercept = unique(rg[[prefix]]$pd$array_num), colour="grey90", size=0.5)
            myplot <- myplot + geom_point(position = "jitter", size=0.5, show.legend = TRUE)
            myplot <- myplot + scale_x_continuous(breaks=unique(rg[[prefix]]$pd$array_num))
            myplot <- myplot + scale_colour_manual(values = tidy$colour)
            myplot <- myplot + theme_light() + theme(panel.grid = element_blank(), legend.position = "right",
                                                     legend.title = element_blank(), legend.key.size = unit(4, "pt"),
                                                     axis.text.y = element_text(size=8), axis.text.x = element_text(size=xtext)) # top, right, bottom, left  plot.margin = unit(c(1, 2, 1, 2), "cm")
            myplot <- myplot + labs(x="Array number", y="Log2 Intensity", colour="Type", title=tools::toTitleCase(tolower(this_ctrl)), tag=LETTERS[i])
            
            gplots[[this_ctrl]] <- myplot
            i <- i + 1
        }
        
        if(this_ctrl == "NEGATIVE") {
            tidy$x <- factor(tidy$x, levels = 1:max(rg[[prefix]]$pd$array_num), ordered = TRUE)
            myfill <- c("red"="red", "grn"="green")
            for(my_channel in c("red", "grn")) {
                myplot <- ggplot(data = tidy[tidy$channel == my_channel, ], mapping = aes(x=x, y=log_signal, group=x)) + ylim(range(tidy$log_signal)) 
                myplot <- myplot + geom_vline(xintercept = unique(rg[[prefix]]$pd$array_num), colour="grey90", size=0.5) # + xlim(range(rg[[prefix]]$pd$array_num))
                myplot <- myplot + geom_boxplot(show.legend = FALSE, notch = TRUE, fill=myfill[my_channel], outlier.size = 0.5, outlier.colour = "grey50")
                myplot <- myplot + scale_x_discrete(breaks=unique(rg[[prefix]]$pd$array_num))
                myplot <- myplot + theme_light() + theme(panel.grid = element_blank(),
                                                         axis.text.y = element_text(size=8), axis.text.x = element_text(size=xtext))
                myplot <- myplot + labs(x="Array number", y="Log2 Intensity", colour="Type",
                                        title=tools::toTitleCase(tolower(paste(this_ctrl, my_channel, sep=": "))), tag=LETTERS[i])
                gplots[[paste(this_ctrl, my_channel, sep="_")]] <- myplot
                i <- i + 1
            }
        }
    }
    
    
    pdf(file=file.path(path_results, paste(prefix, "test_ggplot_ControlProbes.pdf", sep="_")), paper = "a4")
    print(ggarrange(plotlist=gplots[1:3], ncol=1, nrow=3))
    print(ggarrange(plotlist=gplots[4:6], ncol=1, nrow=3))
    print(ggarrange(plotlist=gplots[7:9], ncol=1, nrow=3))
    print(ggarrange(plotlist=gplots[10:11], ncol=1, nrow=3))
    dev.off()
}

# Delete rg we won't need later.

rg[c("EPIC-Italy", "NOVI", "URECA")] <- NULL
gc()

# Make figure

myfig <- data.frame(set=c(rep("EpiSCOPE", 3), rep("BodyFatness", 2)),
                    ctrls=c("STAINING", "TARGET REMOVAL", "EXTENSION", "BISULFITE CONVERSION I", "BISULFITE CONVERSION II"))


xtext = 8
gplots <- list()
i=1
for(f in 1:nrow(myfig)) {
    prefix = myfig$set[f]
    this_ctrl = myfig$ctrls[f]
    
    tidy <- rbind(processChannel(prefix, ctrl=this_ctrl, colour="red"),
                  processChannel(prefix, ctrl=this_ctrl, colour="grn"))
    
    myplot <- ggplot(data = tidy, mapping = aes(x=x, y=log_signal, colour=type))
    myplot <- myplot + geom_vline(xintercept = unique(rg[[prefix]]$pd$array_num), colour="grey90", size=0.5)
    myplot <- myplot + geom_point(position = "jitter", size=0.5, show.legend = TRUE)
    myplot <- myplot + scale_x_continuous(breaks=unique(rg[[prefix]]$pd$array_num))
    myplot <- myplot + scale_colour_manual(values = tidy$colour)
    myplot <- myplot + theme_light(base_size = 7) + theme(panel.grid = element_blank(), legend.position = "right",
                                                          legend.title = element_blank(), legend.key.size = unit(4, "pt"))
    myplot <- myplot + labs(x="Array number", y="Log2 Intensity", colour="Type", title=tools::toTitleCase(tolower(this_ctrl)), tag=LETTERS[i])
    gplots[[this_ctrl]] <- myplot
    i <- i + 1
}


pdf(file=file.path(path_results, "Figure02_control_probes.pdf"), paper = "a4")
print(ggarrange(plotlist=gplots, ncol=1, nrow=5))
dev.off()

rm(prefix, this_ctrl, gplots)


#####  Make samplemeans meth  #####

meth <- list()
for(prefix in c("EpiSCOPE", "BodyFatness")) {
    
    print(prefix)
    mset_file = names(prefix_alts)[prefix_alts == prefix]
    
    # load files
    load(file=file.path(path_save_rdata, paste(mset_file, "methylsets.RData", sep="_")))
    
    # Can't get meth and unmeth for BMIQ or noob+BMIQ, the method doesn't operate on raw meth/unmeth signal data
    # For functional normalisation, an MSet can be returned by setting ratio=FALSE.
    meth[[prefix]] <- list("Green (Cy3)"=rg[[prefix]]$signal_avg_grn,
                           "Red (Cy5)"=rg[[prefix]]$signal_avg_red,
                           "Raw meth"=colMeans(getMeth(MSet)),
                           "Raw unmeth"=colMeans(getUnmeth(MSet)),
                           "SWAN meth"=colMeans(getMeth(MSet.swan)),
                           "SWAN unmeth"=colMeans(getUnmeth(MSet.swan)),
                           "noob meth"=colMeans(getMeth(MSet.noob)),
                           "noob unmeth"=colMeans(getUnmeth(MSet.noob)),
                           "noob+SWAN meth"=colMeans(getMeth(MSet.noob.swan)),
                           "noob+SWAN unmeth"=colMeans(getUnmeth(MSet.noob.swan)),
                           "ENmix meth"=colMeans(getMeth(MSet.enmix)),
                           "ENmix unmeth"=colMeans(getUnmeth(MSet.enmix)),
                           "Illumina meth"=colMeans(getMeth(MSet.illumina)),
                           "Illumina unmeth"=colMeans(getUnmeth(MSet.illumina)),
                           "Functional meth"=colMeans(getMeth(gset.funnorm)),
                           "Functional unmeth"=colMeans(getUnmeth(gset.funnorm)),
                           "Dasen meth"=colMeans(getMeth(MSet.dasen)),
                           "Dasen unmeth"=colMeans(getUnmeth(MSet.dasen)),
                           "noob+Dasen meth"=colMeans(getMeth(MSet.noob.dasen)),
                           "noob+Dasen unmeth"=colMeans(getUnmeth(MSet.noob.dasen)),
                           "ENmix+Quantile meth"=colMeans(getMeth(MSet.enmix.quantile)),
                           "ENmix+Quantile unmeth"=colMeans(getUnmeth(MSet.enmix.quantile)))
    
    rm(list=ls(pattern = "MSet"))
    rm(gset.funnorm)
    gc()
}



#####  Red-Green by Sentrix  #####

# Now as ggplot facet

for(prefix in c("EpiSCOPE", "BodyFatness")) {
    
    tidy_meth = list()
    for(this in names(meth[[prefix]])) {
        tidy_meth[[prefix]][[this]] = data.frame(sample=names(meth[[prefix]][[this]]),
                                            means=meth[[prefix]][[this]],
                                            sentrix_position=factor(rg[[prefix]]$pd$sentrix_description), stringsAsFactors = FALSE)
    }
    tidy_meth = reshape::melt.list(tidy_meth)
    
    names(meth[[prefix]])
    tidy_meth$L2 = factor(tidy_meth$L2, levels=names(meth[[prefix]]), ordered=TRUE)
    
    if(prefix == "EpiSCOPE") {
        mycols <- rep(sentrix_cols[1:6], each=2)
    }
    if(prefix == "BodyFatness") {
        mycols <- sentrix_cols
    }
    names(mycols) <- sort(unique(rg[[prefix]]$pd$sentrix_description))

    pdf(file=file.path(path_results, paste(prefix, "Figure0X_norm_method_comparison_ggplot.pdf", sep="_")), paper = "a4")
    p = ggplot(tidy_meth, aes(x=sentrix_position, y=value, fill=sentrix_position))
    p = p + geom_boxplot(outlier.size = 0.3, outlier.colour = "grey50") + facet_wrap(facets=vars(L2), ncol = 4, nrow = 6) + scale_fill_manual(values=mycols)
    p = p + theme_bw() + labs(x="Sentrix Position", y = "Fluorescence Intensity")
    p = p + theme(axis.text.x = element_text(size=8, angle=90), axis.text.y = element_text(size=8), legend.position = "bottom")
    p = p + guides(fill = guide_legend(label.position = "right", title = "Position", nrow = 2, keywidth = 0.6, keyheight = 0.6))
    print(p)
    dev.off()
}

