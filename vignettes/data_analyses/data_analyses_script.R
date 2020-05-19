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
suppressMessages(library(limma))
suppressMessages(library(GenomicRanges))

# helper functions
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
get_regionrep <- function(cgid, man, nround = 3){
  manf <- man[unique(cgid),]
  dmp.tot <- nrow(manf)
  bg.tot <- nrow(man)
  # filt vars
  manf.filt.gen <- !manf$UCSC_RefGene_Name == ""
  man.filt.gen <- !man$UCSC_RefGene_Name == ""
  manf.filt.isl <- !manf$Relation_to_Island == "OpenSea"
  man.filt.isl <- !man$Relation_to_Island == "OpenSea"
  # genic
  dmp.genic <- round(100*nrow(manf[manf.filt.gen,])/dmp.tot, nround)
  bg.genic <- round(100*nrow(man[man.filt.gen,])/bg.tot, nround)
  # island
  dmp.isl <- round(100*nrow(manf[manf.filt.isl,])/dmp.tot, nround)
  bg.isl <- round(100*nrow(man[man.filt.isl,])/bg.tot, nround)
  # genic-island
  dmp.genisl <- round(100*nrow(manf[manf.filt.isl & manf.filt.gen,])/dmp.tot, nround)
  bg.genisl <- round(100*nrow(man[man.filt.isl & man.filt.gen,])/bg.tot, nround)
  cgt <- matrix(c(dmp.genic, bg.genic, dmp.isl, 
                  bg.isl, dmp.genisl, bg.genisl),
                ncol = 6)
  rownames(cgt) <- rownames(manf)
  colnames(cgt) <- c("dmp.genic", "bg.genic", "dmp.isl", "bg.isl", 
                     "dmp.genic.isl", "bg.genic.isl")
  return(cgt)
}

# path to full dataset
path <- "remethdb2.h5"
mdf <- data_mdpost(path) # all metadata

#---------------------------------
# analysis 1: blood, young vs. old
#---------------------------------

#----------------------------
# analysis 2: brain vs. blood
#----------------------------
# define sgroup and filter samples
{
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
}
# get summary stats table1
sst1 <- get_sst()
save.image(file = savepath)
# query data and assign samples group
{
  rg <- getrg(gsmv = mdf$gsm)
  ms.raw <- preprocessRaw(rg)
  # match metadata and assign sgroup
  mdf <- mdf[order(match(mdf$gsm, rg$gsm)),]
  identical(mdf$gsm, rg$gsm)
  rg$sgroup <- mdf$sgroup
}
save.image(file = savepath)
# qc plots
{
  # all samples, aggregate summaries
  jpeg("data_analyses/densityplot.jpg")
  plot1 <- densityPlot(rg, sampGroups = mdf$sgroup)
  dev.off()
  jpeg("data_analyses/densitybeanplot.jpg")
  densityBeanPlot(rg, sampGroups = mdf$sgroup)
  dev.off()
  jpeg("data_analyses/mdsplot.jpg")
  mdsPlot(rg, sampGroups = mdf$sgroup)
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
}
# quality assurance on raw data
{
  jpeg("data_analyses/qcplot.jpg", 10, 5, units = "in", res = 400)
  qc.blood <- getQC(ms.blood)
  qc.brain <- getQC(ms.brain)
  par(mfrow = c(1, 2))
  plotQC(qc.blood)
  legend("topright", legend = "BLOOD", bty = "n")
  plotQC(qc.brain)
  legend("topright", legend = "BRAIN", bty = "n")
  meth.raw <- getMeth(ms.raw)
  l2meth <- log2(colMedians(meth.raw))
  unmeth.raw <- getUnmeth(ms.raw)
  l2unmeth <- log2(colMedians(unmeth.raw))
  jpeg("data_analyses/l2med.jpg", 10, 5, units = "in", res = 400)
  par(mfrow = c(1, 2))
  boxplot(l2meth ~ rg$sgroup, xlab = "", ylab = "", 
          main = "Meth. log2(median)",
          ylim = c(9, 14))
  boxplot(l2unmeth ~ rg$sgroup, xlab = "", ylab = "", 
          main = "Unmeth. log2(median)",
          ylim = c(9, 14))
  dev.off()
  # remove low-signal samples
  min.sig <- 11
  filt.cond <- which(l2meth < min.sig & l2unmeth < min.sig)
  filt.gsmv <- colnames(rg)[filt.cond]
  rgf <- rg[,!colnames(rg) %in% filt.gsmv]
  # post-qc summary stats
  mdff <- mdf[mdf$gsm %in% colnames(rgf),]
  sst2 <- get_sst(mdff)
}
save.image(file = savepath)
# normalization and linear adjustment
set.seed(2)
ncg <- 5000
cgidv <- sample(rownames(ms.raw), ncg) # random probe subset
{
  # noob normalization
  ms.noob <- preprocessNoob(rgf)
  ms.noob <- ms.noob[cgidv,]
  ms.noob <- mapToGenome(ms.noob)
  ms.unnorm <- ms.raw[cgidv, colnames(ms.noob)]
  # make adjustment on mval
  var.sgroup <- as.factor(ms.noob$sgroup)
  var.gseid <- as.factor(ms.noob$gseid)
  var.predage <- as.numeric(ms.noob$predage)
  var.predsex <- as.factor(ms.noob$predsex)
  mv <- getM(ms.noob)
  cm <- matrix(c(var.predage, var.predsex), ncol = 2)
  madj <- removeBatchEffect(mv, batch = ms.noob$gseid, 
                            covariates = cm)
  gr.adj <- GenomicRatioSet(gr = granges(ms.noob), M = madj)
}
# anova comparisons
{
  lb <- list("raw" = getBeta(ms.unnorm), 
             "noob" = getBeta(ms.noob), 
             "adj" = getBeta(gr.adj))
  lp <- lv <- list("raw" = matrix(ncol = 4, nrow = 0), 
                   "noob" = matrix(ncol = 4, nrow = 0), 
                   "adj" = matrix(ncol = 4, nrow = 0))
  for(r in 1:nrow(gr.adj)){
    for(d in c("raw", "noob", "adj")){
      datr <- as.numeric(lb[[d]][r,])
      ld <- lm(datr ~ var.sgroup + var.gseid + 
                 var.predage + var.predsex)
      an <- anova(ld)
      ap <- an[c(1:4),5]
      av <- round(100*an[c(1:4),2]/sum(an[,2]), 3)
      lp[[d]] <- rbind(lp[[d]], ap)
      lv[[d]] <- rbind(lv[[d]], av)
    }
    message(r)
  }
}
# anova plots
{
  for(d in c("raw", "noob", "adj")){
    rownames(lp[[d]]) <- rownames(lv[[d]]) <- cgidv
    colnames(lp[[d]]) <- colnames(lv[[d]]) <- c("sgroup", "gseid", "predage", "predsex")
  }
  pairs(data.frame("raw" = lp[["raw"]][,2],
                   "noob" = lp[["noob"]][,2],
                   "adj" = lp[["adj"]][,2]))
  
  pairs(data.frame("raw" = lv[["raw"]][,2],
                   "noob" = lv[["noob"]][,2],
                   "adj" = lv[["adj"]][,2]))
  
  plot(density(lv[["adj"]][,2]), lwd = 2, col = "forestgreen", 
       xlab = "GSE ID Percent Variance", main = "", xlim = c(-20, 120))
  lines(density(lv[["noob"]][,2]), lwd = 2, col = "blue")
  lines(density(lv[["raw"]][,2]), lwd = 2, col = "red")
  legend("top", legend = c("raw", "noob", "adj"), bty = "n",
         lwd = c(2,2,2), col = c("red", "blue", "forestgreen"))
  # check contributions within adj set
  pairs(lp$adj)
  pairs(lv$adj)
}
# analysis and results summaries
{
  # set1: blood > brain 
  cg.grp1 <- rowMeans(gr.adj[,var.sgroup == "blood"]) >
    rowMeans(gr.adj[,var.sgroup == "blood"])
  # set2: brain > blood
  cg.grp2 <- rowMeans(gr.adj[,var.sgroup == "brain"]) >
    rowMeans(gr.adj[,var.sgroup == "blood"])
  # differential variances
  dmp.allcg <- dmpFinder(getBeta(gr.adj), var.sgroup)
  # region-level summaries
  cgt <- get_regionrep(rownames(dmp), 
                       as.data.frame(getAnnotation(rg)))
}
# cg filter
{
  max.pval <- 1e-2
  dat <- lp[["adj"]]
  cond.filt <- dat[,1] <= max.pval & 
    dat[,2] > max.pval & 
    dat[,3] > max.pval & 
    dat[,4] > max.pval
  cgid.analyze <- rownames(gr.adj)[cond.filt]
  gr.analysis <- gr.adj[cgid.analyze,]
}
dmp.filt <- dmpFinder(getBeta(gr.analysis), 
                      var.sgroup)
kable(head(dmp))