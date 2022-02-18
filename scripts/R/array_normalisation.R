library(minfi)
library(minfiData)
library(Harman)
library(ENmix)
library(wateRmelon)
library(limma)
library(sva)


#####  Functions  #####


combat_var_hack <- function(bs, noise_sd=1e-8) {
    is_zero <- is.na(bs)
    bs[is_zero] <- rnorm(n=sum(is_zero), mean=0.5, sd=noise_sd)
    bs
}



#####  Normalise  #####

##### Raw  #####

# Use betas as we can get infinite M-values with exact values of 0 or 1.
# Convert NA beta values to 0.5, upon logit transformation these will become
# M-values of 0, which is consistent with SWAN
MSet <- preprocessRaw(RGSet)
bs <- getBeta(MSet)
# on the EPIC arrays, probe cg06180910 has all M-vals of 0. This crashes combat
# as the variance will be 0. Use the hack function to add a little noise around 0
#raw_m <- combat_var_hack(raw_m)
bs <- combat_var_hack(bs)
raw_m <- logit2(Harman::shiftBetas(bs))
rm(bs)

MSet_hm <- harman(datamatrix = raw_m, expt = harman_expt, batch = batch)
raw_harman_m <- reconstructData(MSet_hm)
raw_combat_m = ComBat(dat=raw_m, batch=batch, mod=modcombat)

meth_zero <- getMeth(MSet) == 0
unmeth_zero <- getUnmeth(MSet) == 0
both_zero <- meth_zero & unmeth_zero
either_zero <- meth_zero | unmeth_zero

zero_calls <- data.frame(meth_zero=colSums(meth_zero),
                         unmeth_zero=colSums(unmeth_zero),
                         both_zero=colSums(both_zero),
                         either_zero=colSums(either_zero))
rm(meth_zero, unmeth_zero, both_zero, either_zero)

save(MSet, raw_m, MSet_hm, raw_harman_m, raw_combat_m, zero_calls,
     file=file.path(path_save_rdata, paste(my_prefix, "MSet.Raw.RData", sep="_")),
     compress=TRUE)

rm(raw_m, MSet_hm, raw_harman_m, raw_combat_m, zero_calls)


#####  Illumina  #####

MSet.illumina <- preprocessIllumina(RGSet, bg.correct = TRUE,
                                    normalize = "controls")
bs <- getBeta(MSet.illumina)
bs <- combat_var_hack(bs)
illumina_m <- logit2(Harman::shiftBetas(bs))
rm(bs)

MSet.illumina_hm <- harman(datamatrix = illumina_m, expt = harman_expt, batch = batch)
illumina_harman_m <- reconstructData(MSet.illumina_hm)
illumina_combat_m = ComBat(dat=illumina_m, batch=batch, mod=modcombat)

save(MSet.illumina, illumina_m, MSet.illumina_hm, illumina_harman_m, illumina_combat_m,
     file=file.path(path_save_rdata, paste(my_prefix, "MSet.Illumina.RData", sep="_")),
     compress=TRUE)

rm(MSet.illumina, illumina_m, MSet.illumina_hm, illumina_harman_m, illumina_combat_m)


##### Swan  #####

# Set random seed as SWAN results will differ slightly by seed
set.seed(my_rndseed)
MSet.swan <- preprocessSWAN(RGSet)
swan_m <- getM(MSet.swan)

MSet.swan_hm <- harman(datamatrix = swan_m, expt = harman_expt, batch = batch)
swan_harman_m <- reconstructData(MSet.swan_hm)
swan_combat_m = ComBat(dat=swan_m, batch=batch, mod=modcombat)

save(MSet.swan, swan_m, MSet.swan_hm, swan_harman_m, swan_combat_m,
     file=file.path(path_save_rdata, paste(my_prefix, "MSet.swan.RData", sep="_")),
     compress=TRUE)

rm(MSet.swan, swan_m, MSet.swan_hm, swan_harman_m, swan_combat_m)


##### noob  #####

MSet.noob <- preprocessNoob(RGSet, dyeCorr = TRUE)
noob_m <- getM(MSet.noob)

MSet.noob_hm <- harman(datamatrix = noob_m, expt = harman_expt, batch = batch)
noob_harman_m <- reconstructData(MSet.noob_hm)
noob_combat_m = ComBat(dat=noob_m, batch=batch, mod=modcombat)

session <- sessionInfo()

save(MSet.noob, noob_m, MSet.noob_hm, noob_harman_m, noob_combat_m, session,
     file=file.path(path_save_rdata, paste(my_prefix, "MSet.noob.RData", sep="_")),
     compress=TRUE)

rm(noob_m, MSet.noob_hm, noob_harman_m, noob_combat_m)
gc()

##### Swan + noob  #####

# Set random seed as SWAN results will differ slightly by seed
set.seed(my_rndseed)
MSet.noob.swan <- preprocessSWAN(RGSet, mSet=MSet.noob)
noob.swan_m <- getM(MSet.noob.swan)

MSet.noob.swan_hm <- harman(datamatrix = noob.swan_m, expt = harman_expt, batch = batch)
noob.swan_harman_m <- reconstructData(MSet.noob.swan_hm)
noob.swan_combat_m = ComBat(dat=noob.swan_m, batch=batch, mod=modcombat)

save(MSet.noob.swan, noob.swan_m, MSet.noob.swan_hm, noob.swan_harman_m, noob.swan_combat_m,
     file=file.path(path_save_rdata, paste(my_prefix, "MSet.noob.swan.RData", sep="_")),
     compress=TRUE)

rm(MSet.noob.swan, noob.swan_m, MSet.noob.swan_hm, noob.swan_harman_m, noob.swan_combat_m)


##### BMIQ  #####

bs <- getBeta(MSet)
bs <- combat_var_hack(bs)
type_ii <- getProbeInfo(MSet, type="II")$Name
design.v <- rep(1, times = nrow(bs))
design.v[rownames(bs) %in% type_ii] <- 2

# Have to build a loop for beta values, no vectorised constructor :(
bs_bmiq <- matrix(nrow = nrow(bs), ncol = ncol(bs),
                  dimnames = dimnames(bs))
for (i in seq_len(ncol(bs))) {
    this_col <- colnames(bs)[i]
    # Same as default signature for MethylLumiSet
    bmiq_res <- BMIQ(beta.v = bs, design.v = design.v,
                     sampleID = this_col,
                     nL=3, doH=TRUE, nfit=5000,
                     th1.v=c(0.2,0.75), th2.v=NULL,
                     niter=5, tol=0.001, plots=FALSE, pri=TRUE)
    bs_bmiq[, this_col] <- bmiq_res$nbeta[, this_col]
}

bmiq_m <- logit2(Harman::shiftBetas(bs_bmiq))
rm(bs_bmiq)

MSet.bmiq_hm <- harman(datamatrix = bmiq_m, expt = harman_expt, batch = batch)
bmiq_harman_m <- reconstructData(MSet.bmiq_hm)
bmiq_combat_m = ComBat(dat=bmiq_m, batch=batch, mod=modcombat)

save(bmiq_m, MSet.bmiq_hm, bmiq_harman_m, bmiq_combat_m,
     file=file.path(path_save_rdata, paste(my_prefix, "MSet.bmiq.RData", sep="_")),
     compress=TRUE)

rm(bmiq_m, MSet.bmiq_hm, bmiq_harman_m, bmiq_combat_m)


##### BMIQ + noob  #####

bs_noob.bmiq <- BMIQ(MSet.noob)
noob.bmiq_m <- logit2(bs_noob.bmiq)
rm(bs_noob.bmiq)

MSet.noob.bmiq_hm <- harman(datamatrix = noob.bmiq_m, expt = harman_expt, batch = batch)
noob.bmiq_harman_m <- reconstructData(MSet.noob.bmiq_hm)
noob.bmiq_combat_m = ComBat(dat=noob.bmiq_m, batch=batch, mod=modcombat)

save(noob.bmiq_m, MSet.noob.bmiq_hm, noob.bmiq_harman_m, noob.bmiq_combat_m,
     file=file.path(path_save_rdata, paste(my_prefix, "MSet.noob.bmiq.RData", sep="_")),
     compress=TRUE)

rm(noob.bmiq_m, MSet.noob.bmiq_hm, noob.bmiq_harman_m, noob.bmiq_combat_m)


##### Dasen  #####

MSet.dasen <- dasen(RGSet)
bs <- getBeta(MSet.dasen)
# Some betas are 0, makes for NaN M values
bs <- combat_var_hack(bs)
dasen_m <- logit2(bs)
rm(bs)

MSet.dasen_hm <- harman(datamatrix = dasen_m, expt = harman_expt, batch = batch)
dasen_harman_m <- reconstructData(MSet.dasen_hm)
dasen_combat_m = ComBat(dat=dasen_m, batch=batch, mod=modcombat)

save(MSet.dasen, dasen_m, MSet.dasen_hm, dasen_harman_m, dasen_combat_m,
     file=file.path(path_save_rdata, paste(my_prefix, "MSet.dasen.RData", sep="_")),
     compress=TRUE)

rm(MSet, MSet.dasen, dasen_m, MSet.dasen_hm, dasen_harman_m, dasen_combat_m)


##### Dasen + noob  #####

MSet.noob.dasen <- dasen(MSet.noob)
noob.dasen_m <- getM(MSet.noob.dasen)

MSet.noob.dasen_hm <- harman(datamatrix = noob.dasen_m, expt = harman_expt, batch = batch)
noob.dasen_harman_m <- reconstructData(MSet.noob.dasen_hm)
noob.dasen_combat_m = ComBat(dat=noob.dasen_m, batch=batch, mod=modcombat)

save(MSet.noob.dasen, noob.dasen_m, MSet.noob.dasen_hm, noob.dasen_harman_m, noob.dasen_combat_m,
     file=file.path(path_save_rdata, paste(my_prefix, "MSet.noob.dasen.RData", sep="_")),
     compress=TRUE)

rm(MSet.noob.dasen, noob.dasen_m, MSet.noob.dasen_hm, noob.dasen_harman_m, noob.dasen_combat_m, MSet.noob)


##### Funnorm  #####

# Have ratioConvert set to FALSE to return a GenomicRatioSet instead, which allows getMeth()
gset.funnorm <- preprocessFunnorm(RGSet, ratioConvert = FALSE)
#gset.funnorm <- preprocessFunnorm(RGSet)
funnorm_beta <- shiftBetas(combat_var_hack(minfi::getBeta(gset.funnorm)))
funnorm_m <- logit2(funnorm_beta)
rm(funnorm_beta)

gset.funnorm_hm <- harman(datamatrix = funnorm_m, expt = harman_expt, batch = batch)
funnorm_harman_m <- reconstructData(gset.funnorm_hm)
funnorm_combat_m = ComBat(dat=funnorm_m, batch=batch, mod=modcombat)

save(gset.funnorm, funnorm_m, gset.funnorm_hm, funnorm_harman_m, funnorm_combat_m,
     file=file.path(path_save_rdata, paste(my_prefix, "gset.funnorm.RData", sep="_")),
     compress=TRUE)

rm(gset.funnorm, funnorm_m, gset.funnorm_hm, funnorm_harman_m, funnorm_combat_m)
gc()

##### ENmix  #####

MSet.enmix <- preprocessENmix(RGSet, nCores=20)
enmix_m <- getM(MSet.enmix)

MSet.enmix_hm <- harman(datamatrix = enmix_m, expt = harman_expt, batch = batch)
enmix_harman_m <- reconstructData(MSet.enmix_hm)
enmix_combat_m = ComBat(dat=enmix_m, batch=batch, mod=modcombat)

save(MSet.enmix, enmix_m, MSet.enmix_hm, enmix_harman_m, enmix_combat_m,
     file=file.path(path_save_rdata, paste(my_prefix, "MSet.enmix.RData", sep="_")),
     compress=TRUE)

rm(enmix_m, MSet.enmix_hm, enmix_harman_m, enmix_combat_m)


##### ENmix Quantile  #####

MSet.enmix.quantile <- norm.quantile(MSet.enmix, method="quantile1")
enmix.quantile_m <- getM(MSet.enmix.quantile)

MSet.enmix.quantile_hm <- harman(datamatrix = enmix.quantile_m, expt = harman_expt, batch = batch)
enmix.quantile_harman_m <- reconstructData(MSet.enmix.quantile_hm)
enmix.quantile_combat_m = ComBat(dat=enmix.quantile_m, batch=batch, mod=modcombat)

save(MSet.enmix.quantile, enmix.quantile_m, MSet.enmix.quantile_hm, enmix.quantile_harman_m, enmix.quantile_combat_m,
     file=file.path(path_save_rdata, paste(my_prefix, "MSet.enmix.quantile.RData", sep="_")),
     compress=TRUE)

rm(MSet.enmix, MSet.enmix.quantile, enmix.quantile_m, MSet.enmix.quantile_hm, enmix.quantile_harman_m, enmix.quantile_combat_m)

