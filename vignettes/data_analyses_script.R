#!/usr/bin/env R

# Script to generate env objects used in data_analyses.Rmd
# Notes:
# Currently running this with hdf5 db
# Rerun with h5se RGChannelSet when this is working (e.g. syncs to server, etc.)

savepath <- "data_analyses/env.RData"
# list of objects to store for package builds
obj.store  <- c("sst1", "sst2", "plot1", "plot2", "st.query") 

suppressMessages(library(rhdf5))
suppressMessages(library(minfi))
suppressMessages(library(recountmethylation))
suppressMessages(library(knitr))

path <- "remethdb2.h5"
mdf <- data_mdpost(path)

# specify filter values
cn <- "tissue" # colname to query
labs <- c("blood", "brain") # labels of interest
regex.patt <- paste0("(^|;)", labs,"(;|$)") # pattern for strict label matches
# define and filter on sgroup variable
mdf$sgroup <- "NA"
sg.blood <- grepl(regex.patt[1], mdf[,cn]) & 
  !grepl(regex.patt[2], mdf[,cn])
sg.brain <- grepl(regex.patt[2], mdf[,cn]) & 
  !grepl(regex.patt[1], mdf[,cn])
mdf[sg.blood,]$sgroup <- labs[1]
mdf[sg.brain,]$sgroup <- labs[2]
mdf <- mdf[sg.blood|sg.brain,]
# select 50 samples per groups
set.seed(2)
ngsm.sg <- 50
gsm.blood <- which(mdf$sgroup == labs[1])
gsm.brain <- which(mdf$sgroup == labs[2])
gsm.indices <- c(sample(gsm.blood, ngsm.sg), sample(gsm.brain, ngsm.sg))
mdf <- mdf[gsm.indices,]

get_sst <- function(labs = c("blood", "brain")){
  # get summary stats by label
  sstat <- lapply(labs, function(x){
    mv <- mdf[mdf$sgroup==x,]
    v1 <- nrow(mv)
    numgsm.gse <- as.data.frame(table(mv$gseid))
    v2 <- round(mean(numgsm.gse[,2]), 3)
    v3 <- round(sd(numgsm.gse[,2]), 3)
    v4 <- nrow(numgsm.gse)
    # predage summaries
    pv <- as.numeric(mv$predage)
    v5 <- min(pv)
    v6 = max(pv)
    v7 <- round(mean(pv, na.rm = TRUE), 3)
    v8 <- round(sd(pv, na.rm = TRUE), 3)
    v9 <- length(pv[is.na(pv)])
    # predsex summaries
    v10 <- round(100*nrow(mv[mv$predsex=="F",])/nrow(mv), 3) # perc. F
    v11 <- nrow(mv[is.na(mv$predsex),])
    return(data.frame("ngsm" = v1, "meangsm.gse" = v2, "sdgsm.gse" = v3,  "numgse" = v4,
                      "min.predage" = v5, "max.predage" = v6, "mean.predage" = v7, 
                      "sd.predage" = v8, "numna.predage" = v9, "percfemale.predsex" = v10, 
                      "numna.predsex" = v11))
  })
  sst <- sapply(sstat, rbind) # bind label groups
  colnames(sst) <- labs
  return(sst)
}
sst1 <- get_sst()
save.image(file = savepath)

library(minfi)

# get the RGChannelSet for samples of interest
# note this can take several minutes
s1 <- Sys.time()
rg <- getrg(gsmv = mdf$gsm)
st.query <- Sys.time() - s1
save.image(file = savepath)

ms.raw <- preprocessRaw(rg)

# get some minfi qc metrics
mdf <- mdf[order(match(mdf$gsm, rg$gsm)),]
identical(mdf$gsm, rg$gsm)

# all samples, aggregate
jpeg("data_analyses/densityplot.jpg")
plot1 <- densityPlot(rg, sampGroups = mdf$sgroup)
dev.off()
jpeg("data_analyses/densitybeanplot.jpg")
densityBeanPlot(rg, sampGroups = mdf$sgroup)
dev.off()
jpeg("data_analyses/mdsplot.jpg")
mdsPlot(rg, sampGroups = mdf$sgroup)
dev.off()

jpeg("data_analyses/qcplot.jpg", 10, 5, units = "in", res = 400)

qc.blood <- getQC(ms.blood)
qc.brain <- getQC(ms.brain)

par(mfrow = c(1, 2))

plotQC(qc.blood)
legend("topright", legend = "BLOOD", bty = "n")
plotQC(qc.brain)
legend("topright", legend = "BRAIN", bty = "n")

dev.off()

# subsets and group-specific plots
is.blood <- mdf$sgroup == "blood"
is.brain <- mdf$sgroup == "brain"
rg.blood <- rg[,is.blood]
rg.brain <- rg[,is.brain]
ms.blood <- ms.raw[,is.blood]
ms.brain <- ms.raw[,is.brain]

# blood-only
jpeg("data_analyses/ctrlstripplot_blood.jpg")
controlStripPlot(rg.blood, controls = c("BISULFITE CONVERSION II"))
dev.off()
jpeg("data_analyses/betasbytype_blood.jpg")
plotBetasByType(ms.blood[,sample(ncol(ms.blood), 1)])
dev.off()

# brain-only
jpeg("data_analyses/ctrlstripplot_brain.jpg")
controlStripPlot(rg.brain, controls = c("BISULFITE CONVERSION II"))
dev.off()
jpeg("data_analyses/betasbytype_brain.jpg")
plotBetasByType(ms.brain[,sample(ncol(ms.brain), 1)])
dev.off()

# do some qc filters

# rerun summary stats table
mdf <- mdf[mdf$gsm %in% rg$gsm,]
sst2 <- get_sst()
save.image(file = savepath)

#-------------------
# diff dnam analyses
#-------------------
identical(mdf$gsm, rg$gsm)
var <- "sgrp"

# ftest
dft <- dmpFinder(mdf$sgrp)

# filtered table
fdr.thresh <- 1e-3
dftf <- dft[dft$fdr <= fdr.thresh,]

# summary plots

#--------------------
# enrichment analyses
#--------------------

