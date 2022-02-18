rm(list=ls())

#####  Paths  #####

path_base <- "/home/ros259/R_projects/450K_Batch_effects"
path_data <- file.path(path_base, "data")
path_save_rdata <- file.path(path_base, "data", "rdata")
path_results <- file.path(path_base, "results")


#####  450K Figure  #####

load(file=file.path(path_save_rdata, "GSE89278_RGSet.RData"))
pd <- as.data.frame(pData(RGSet))

unique(pd$array)

table(pd$array_num, pd$sex)
table(pd$array_num, pd$pheno)

xrange <- c(1,16)
yrange <- c(1,24)

pdf(file.path(path_results, "GSE89278_Figure0X_Slide_design.pdf"))

plot(NA, NA, xlim=xrange, ylim=yrange,
     frame.plot=FALSE, axes=FALSE, xlab="", ylab="")

x_start <- min(xrange)
y_start <- max(yrange)

for(a in split(pd, factor(pd$array_num))) {
    
    x <- x_start + (as.integer(a$col) - 1)
    y <- y_start - (as.integer(a$row) - 1)
    points(x, y,
           bg=c("pink", "powderblue")[factor(a$sex)],
           pch=c(21, 22)[factor(a$pheno)],
           cex=2)
    text(min(x) + 0.5, max(y) - 2.5, labels=unique(a$array_num),
         col="grey30", cex=0.7)
    
    y_start <- min(y) - 1
    if(unique(a$array_num %% 4) == 0) {
        x_start <- x_start + 2
        y_start <- max(yrange)
    }
}
    
abline(v=seq(from=min(xrange), to=max(xrange) + 1, by = 2) - 0.5,
       h=seq(from=min(yrange), to=max(yrange) + 1, by = 6) - 0.5,
       col="grey30")

dev.off()
rm(pd, detP, RGSet)


#####  EPIC Figure  #####

load(file=file.path(path_save_rdata, "SGA_Infants_RGSet.RData"))
pd <- as.data.frame(pData(RGSet))
pd$row <- sub('\\d+_R(\\d+)C\\d+', '\\1', rownames(pd))
pd$col <- sub('\\d+_R\\d+C(\\d+)', '\\1', rownames(pd))

unique(pd$array)

table(pd$array_num, pd$gender)
bf_thirds <- cut(pd$bodyfat, 3)

xrange <- c(1,16)
yrange <- c(1,24)

pdf(file.path(path_results, "EPIC_Figure0X_Slide_design.pdf"))

plot(NA, NA, xlim=xrange, ylim=yrange,
     frame.plot=FALSE, axes=FALSE, xlab="", ylab="")

x_start <- min(xrange)
y_start <- max(yrange)

for(a in split(pd, factor(pd$array_num))) {
    
    x <- x_start + (as.integer(a$col) - 1)
    y <- y_start - (as.integer(a$row) - 1)
    points(x, y,
           bg=c("pink", "powderblue")[factor(a$sex)],
           pch=c(25, 23, 24)[bf_thirds],
           #pch=21,
           cex=1.8)
    text(min(x) + 1, max(y) - 3.5, labels=unique(a$array_num),
         col="grey30", cex=0.7)
    
    y_start <- min(y) - 1
    if(unique(a$array_num %% 3) == 0) {
        x_start <- x_start + 2
        y_start <- max(yrange)
    }
}

abline(v=seq(from=min(xrange), to=max(xrange) + 1, by = 2) - 0.5,
       h=seq(from=min(yrange), to=max(yrange) + 1, by = 8) - 0.5,
       col="grey30")

dev.off()
