#!/usr/bin/env R

# Script to generate env objects used in data_analyses.Rmd
# Notes:
# Currently running this with hdf5 db
# Rerun with h5se RGChannelSet when this is working (e.g. syncs to server, etc.)

savepath <- "data_analyses/env.RData"

# list of objects to store for package builds
library(rhdf5)
library(minfi)
library(recountmethylation)
library(knitr)
library(limma)
library(GenomicRanges)

# path to full dataset
path <- "remethdb2.h5"
mdf <- data_mdpost(path) # all metadata

#------------
# global vars
#------------
# get local metadata
mddir <- system.file("extdata", "metadata", package = "recountmethylation")
mdpath <- list.files(mddir)
md <- get(load(mdpath))

# load MethylSet
gmdn <- "remethdb-h5se_gm_0-0-1_1590090412"
gm <- loadHDF5SummarizedExperiment(gmdn)

# load noob norm GenomicRatioSet
grdn <- "remethdb-h5se_gr_0-0-1_1590090412"
gr <- loadHDF5SummarizedExperiment(grdn)

#-----------------
# helper functions
#-----------------
get_sst <- function(sgroup.labs = c("blood", "brain"), mdf, nround = 2){
  # get summary stats by label
  sstat <- lapply(sgroup.labs, function(x){
    mv <- mdf[mdf$sgroup==x,]
    v1 <- nrow(mv)
    numgsm.gse <- as.data.frame(table(mv$gseid))
    v2 <- round(mean(numgsm.gse[,2]), nround)
    v3 <- round(sd(numgsm.gse[,2]), nround)
    v4 <- nrow(numgsm.gse)
    # predage summaries
    pv <- as.numeric(mv$predage)
    v5 <- min(pv)
    v6 = max(pv)
    v7 <- round(mean(pv, na.rm = TRUE), nround)
    v8 <- round(sd(pv, na.rm = TRUE), nround)
    v9 <- length(pv[is.na(pv)])
    # predsex summaries
    v10 <- round(100*nrow(mv[mv$predsex=="F",])/nrow(mv), nround) # perc. F
    v11 <- nrow(mv[is.na(mv$predsex),])
    return(data.frame("ngsm" = v1, "meangsm.gse" = v2, "sdgsm.gse" = v3,  "numgse" = v4,
                      "min.predage" = v5, "max.predage" = v6, "mean.predage" = v7, 
                      "sd.predage" = v8, "numna.predage" = v9, "percfemale.predsex" = v10, 
                      "numna.predsex" = v11))
  })
  sst <- sapply(sstat, rbind) # bind label groups
  colnames(sst) <- sgroup.labs
  return(sst)
}
getblocks <- function(slength, bsize){
  iv <- list()
  if(slength < bsize){
    iv[[1]] <- seq(1, slength, 1)
  } else{
    sc <- 1; ec <- sc + bsize - 1
    nblocks <- slength %/% bsize
    for(b in 1:nblocks){
      iv[[b]] <- seq(sc, ec, 1)
      sc <- ec + 1; ec <- ec + bsize
    }
    # add final indices
    if(nblocks < (slength/bsize)){
      iv[[length(iv) + 1]] <- seq(sc, slength, 1)
    }
  }
  return(iv)
}
makevp <- function(lfilt, ltxcg){
  bpdf.mean <- bpdf.var <- matrix(nrow = 0, ncol = 2)
  for(t in 1:length(lfilt)){
    tname = names(lfilt)[t]
    bt = as.data.frame(lfilt[[t]])
    btf <- bt[rownames(bt) %in% ltxcg[[t]],]
    dt <- data.frame(btf$mean, btf$var, rep(tname, nrow(btf)), 
                     stringsAsFactors = F)
    bpdf.mean = rbind(bpdf.mean, dt[,c(1, 3)])
    bpdf.var = rbind(bpdf.var, dt[,c(2, 3)])
  }
  bpdf.mean = as.data.frame(bpdf.mean, stringsAsFactors = F)
  bpdf.var = as.data.frame(bpdf.var, stringsAsFactors = F)
  bpdf.mean[,1] = as.numeric(as.character(bpdf.mean[,1]))
  bpdf.var[,1] = as.numeric(as.character(bpdf.var[,1]))
  colnames(bpdf.mean) = c("mean", "tissue")
  colnames(bpdf.var) = c("var", "tissue")
  
  vp1 = ggplot(bpdf.mean, aes(x = tissue, y = mean, fill = tissue)) + 
    geom_violin(trim = F, show.legend = F, draw_quantiles = c(0.5)) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90)) + 
    ggtitle("") + xlab("") + ylab("Mean")
  
  vp2 = ggplot(bpdf.var, aes(x = tissue, y = var, fill = tissue)) + 
    geom_violin(trim = F, show.legend = F, draw_quantiles = c(0.5)) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90)) + 
    ggtitle("") + xlab("") + ylab("Variance")
  
  return(list("vp.mean" = vp1, "vp.var" = vp2))
}
get_cga <- function(anno){
  prom.stat = grepl("TSS|5'", anno$UCSC_RefGene_Group)
  body.stat = grepl("Body|Exon|3'", anno$UCSC_RefGene_Group)
  anno$gene.type = ifelse(anno$UCSC_RefGene_Name=="", "intergenic",
                          ifelse(prom.stat & body.stat, "intragenic_promoter-body",
                                 ifelse(prom.stat & !body.stat, "intragenic_promoter", 
                                        "intragenic_body")))
  
  anno$isl.type = ifelse(anno$Relation_to_Island=="OpenSea", "interisland_opensea", 
                         ifelse(anno$Relation_to_Island=="Island", "intraisland_main", 
                                "intraisland_other"))
  anno$type.composite = paste0(anno$isl.type,";",anno$gene.type)
  return(anno)
}
hmsets <- function(datlf, cga, mincg = 2){
  ngrp <- length(datlf)
  annom = matrix(nrow = 0, ncol = ngrp)
  # summary stats -- annotated region overlaps with high-var cgs
  utype = unique(cga$type.composite)
  hmma.mean = hmma.var = hmma.size = matrix(nrow = 0, ncol = ngrp)
  for(a in utype){
    cgidv = rownames(cga[cga$type.composite==a,])
    rdat.mean = rdat.var = rdat.size = c()
    for(t in 1:ngrp){
      # get intersect with txmvp and type
      datx <- as.data.frame(datlf[[t]])
      mvpt = rownames(datlf[[t]])
      cgint = intersect(mvpt, cgidv)
      cgint.size = length(cgint)
      if(length(cgint) >= mincg){
        rdat.mean = c(rdat.mean, mean(datx[rownames(datx) %in% cgint,]$mean))
        rdat.var = c(rdat.var, var(datx[rownames(datx) %in% cgint,]$mean))
        rdat.size = c(rdat.size, cgint.size)
      } else{
        rdat.mean = c(rdat.mean, "NA")
        rdat.var = c(rdat.var, "NA")
        rdat.size = c(rdat.size, "NA")
      }
    }
    hmma.mean = rbind(hmma.mean, rdat.mean)
    hmma.var = rbind(hmma.var, rdat.var)
    hmma.size = rbind(hmma.size, rdat.size)
    # message(a)
  }
  rownames(hmma.mean) <- rownames(hmma.var) <- rownames(hmma.size) <- utype
  colnames(hmma.mean) <- colnames(hmma.var) <- colnames(hmma.size) <- names(datlf)
  class(hmma.mean) <- class(hmma.var) <- class(hmma.size) <- "numeric"
  return(list("hm.mean" = hmma.mean, "hm.var" = hmma.var, "hm.size" = hmma.size))
}
hmplots <- function(hmma.mean, hmma.var, hmma.size){
  # coerce to tall tables for plots
  hdx <- hdv <- matrix(nrow = 0, ncol = 4)
  for(r in rownames(hmma.mean)){
    for(c in colnames(hmma.mean)){
      hdx = rbind(hdx, matrix(c(hmma.mean[r, c], c, r, hmma.size[r, c]), nrow = 1))
      hdv = rbind(hdv, matrix(c(hmma.var[r, c], c, r, hmma.size[r, c]), nrow = 1))
    }
  }
  hdx = as.data.frame(hdx, stringsAsFactors = F)
  hdv = as.data.frame(hdv, stringsAsFactors = F)
  colnames(hdx) = c("mean", "tissue", "anno", "size")
  colnames(hdv) = c("var", "tissue", "anno", "size")
  hdx$mean <- as.numeric(hdx$mean)
  hdx$size <- as.numeric(hdx$size)
  hdv$var <- as.numeric(hdv$var)
  hdv$size <- as.numeric(hdv$size)
  # get heatmap plot objects
  hm.mean = ggplot(hdx, aes(tissue, anno)) +
    geom_tile(aes(fill = mean)) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5,
                         limits = c(0, 1)) +
    geom_text(aes(label = size), color = "black") +
    theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("Mean") + ylab("") + xlab("")
  
  hm.var = ggplot(hdv, aes(tissue, anno)) +
    geom_tile(aes(fill = var)) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.04, 
                         limits = c(0, round(max(hdv$var), 2))) +
    geom_text(aes(label = size), color = "black") +
    theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
    xlab("") + ylab("") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ggtitle("Variance")
  return(list("hm.mean.plot" = hm.mean, "hm.var.plot" = hm.var))
}

#---------------------------------------
# example 1 -- Mined and predicted ages
#---------------------------------------
# get new vars and do pre-filtering
mdf <- md[!md$age == "valm:NA",]
mdf$chron.age <- as.numeric(gsub(";.*", "", gsub("^valm:", "", mdf$age)))
mdf$predage <- as.numeric(mdf$predage)
mdf <- mdf[!is.na(mdf$chron.age),]
mdf <- mdf[!is.na(mdf$predage),]

mdf$stype <- as.character(gsub(";.*", "", gsub("^msraptype:", "", mdf$sampletype)))
mdf <- mdf[!is.na(mdf$stype),]

mdf$is.cx <- ifelse(grepl(".*cancer.*", mdf$disease), TRUE, FALSE)

# get study-wise errors
xdif <- c()
for(g in unique(mdf$gseid)){
  mdff <- mdf[mdf$gseid==g, ]
  xdif <- c(xdif, mean(abs(mdff$chron.age - as.numeric(mdff$predage))))
}
names(xdif) <- unique(mdf$gseid)

# get post-filtered metadata
filt <- mdf$stype == "tissue" & !mdf$is.cx
filt <- filt & !mdf$gseid %in% names(xdif[xdif > 10])
mdff <- mdf[filt, ]

# do multiple regression
lm1 <- lm(mdf$predage ~ mdf$chron.age + mdf$gseid + mdf$stype + mdf$is.cx)
lm2 <- lm(mdff$predage ~ mdff$chron.age + mdff$gseid)

# anovas
av1 <- anova(lm1)
av2 <- anova(lm2)
# results summaries
sperc1 <- round(100*av1$`Sum Sq`[1:4]/sum(av1$`Sum Sq`), 2)
pval1 <- av1$`Pr(>F)`[1:4]
sperc2 <- round(100*av2$`Sum Sq`[1:2]/sum(av2$`Sum Sq`), 2)
pval2 <- av2$`Pr(>F)`[1:2]
# summary table
dan <- data.frame(varperc.anova1 = c(sperc1), 
                  pval.anova1 = c(pval1),
                  varperc.anova2 = c(sperc2, NA, NA), 
                  pval.anova2 = c(pval2, NA, NA), 
                  stringsAsFactors = FALSE)
rownames(dan) <- c("Chron.Age", "GSEID", "SampleType", "Cancer")
# r-squared
rsq1 <- round(summary(lm1)$r.squared, 2)
rsq2 <- round(summary(lm2)$r.squared, 2)
# cor test estimates
rho1 <- round(cor.test(mdf$predage, mdf$chron.age, 
                       method = "spearman")$estimate, 2)
rho2 <- round(cor.test(mdff$predage, mdff$chron.age, 
                       test = "spearman")$estimate, 2)
# mean abs differences
mad1 <- round(mean(abs(mdf$chron.age - mdf$predage)), 2)
mad2 <- round(mean(abs(mdff$chron.age - mdff$predage)), 2)
# summary stats df
dss <- data.frame(groups = c("all.available", "filtered.gseerr.nocx"),
                  r.squared = c(rsq1, rsq2), rho = c(rho1, rho2),
                  mad = c(mad1, mad2), stringsAsFactors = FALSE)

# barplot of study-wise differences
xdif <- xdif[order(xdif)]
jpeg("barplot_ages.jpg", 9, 3, units = "in", res = 400)
barplot(xdif, las = 2, cex.names = 0.2)
abline(h = 10, col = "red")
dev.off()
# scatterplot of ages among post-filtered samples
jpeg("mainfig1b_agedif-mfilt.jpg", 3, 3)
ggplot(mdff, aes(x = chron.age, y = predage)) +
  geom_point(size = 1.2, alpha = 0.2) + geom_smooth(method = "lm", size = 1.2) +
  theme_bw() + xlab("Chronological Age") + ylab("Epigenetic (DNAm) Age")
dev.off()

#-------------------------------------------------------
# example 2 -- FFPE versus Frozen tissue quality signals
#-------------------------------------------------------
# get samples with storage condition 
md <- pData(gm)
mdf <- md[!md$storage == "NA",]
gmf <- gm[, gm$gsm %in% mdf$gsm]
table(mdf$storage)

# get signal matrices
meth.all <- getMeth(gmf)
unmeth.all <- getUnmeth(gmf)

# get blocks for processing
blocks <- getblocks(slength = ncol(gmf), bsize = 100)

# process data in blocks
ms <- matrix(nrow = 0, ncol = 2)
for(i in 1:length(blocks)){
  b <- blocks[[i]]
  gmff <- gmf[, b]
  methb <- as.matrix(meth.all[, b])
  unmethb <- as.matrix(unmeth.all[, b])
  l2meth <- l2unmeth <- c()
  for(c in 1:ncol(methb)){
    l2meth <- c(l2meth, log2(median(methb[,c])))
    l2unmeth <- c(l2unmeth, log2(median(unmethb[,c])))
  }
  ms <- rbind(ms, matrix(c(l2meth, l2unmeth), ncol = 2))
  message(i)
}
rownames(ms) <- colnames(meth.all)
colnames(ms) <- c("meth.l2med", "unmeth.l2med")
ds <- as.data.frame(ms)
ds$storage <- ifelse(grepl("FFPE", gmf$storage), "ffpe", "frozen")
save(ds, file = "df-l2med-signals.rda")

# signals scatterplot by storage type
jpeg("scatterplot_storage.jpg", 7, 7, units = "in", res = 400)
ggplot(ds, aes(x = meth.l2med, y = unmeth.l2med, color = storage)) + 
  geom_point() + theme_bw() +
  scale_color_manual(values = c("ffpe" = "orange", "frozen" = "purple"))
dev.off()

# violin plots of signals by storage type
vp <- matrix(nrow = 0, ncol = 2)
vp <- rbind(vp, matrix(c(ds$meth.l2med, paste0("meth.", ds$storage)), ncol = 2))
vp <- rbind(vp, matrix(c(ds$meth.l2med, paste0("unmeth.", ds$storage)), ncol = 2))
vp <- as.data.frame(vp, stringsAsFactors = FALSE)
vp[,1] <- as.numeric(vp[,1])
colnames(vp) <- c("signal", "group")
vp$col <- ifelse(grepl("ffpe", vp$group), "orange", "purple")

jpeg("vp_storage.jpg", 5, 3, units = "in", res = 400)
ggplot(vp, aes(x = group, y = signal, color = group)) + 
  scale_color_manual(values = c("meth.ffpe" = "orange", "unmeth.ffpe" = "orange", 
                                "meth.frozen" = "purple", "unmeth.frozen" = "purple")) +
  geom_violin() + theme_bw() + theme(legend.position = "none")
dev.off()

#-----------------------------------------------------------------
# example 3 -- Tissue-specific DNAm variation in liver and adipose
#-----------------------------------------------------------------
# get samples GSM IDs
fn <- "old/nct_7tissues_gsm7k.rda" # GSM IDs of interest
gsmv <- get(load(fn))
tv <- c("liver", "adipose")
tstr <- paste0(".*", tv, ".*")
mdf <- md[md$gsm %in% gsmv & grepl(paste(tstr, collapse = "|"), md$tissue),]
mdf$sgroup <- ifelse(grepl(tstr[1], mdf$tissue), tv[1], tv[2])
# get samples summary
sst.tvar <- get_sst(sgroup.labs = c("liver", "adipose"), mdf)

# get signal data
l2med.meth <- log2(apply(getMeth(gm), 2, median))
l2med.unmeth <- log2(apply(getUnmeth(gm), 2, median))
lqc <- list("l2med.meth" = l2med.meth, "l2med.unmeth" = l2med.unmeth, "gsmv" = gsmv)

# plot signal data
jpeg("qc_2tissues.jpg", 5, 5, units = "in", res = 400)
plot(lqc[["l2med.meth"]], lqc[["l2med.unmeth"]], 
     col = ifelse(names(lqc[["gsmv"]]) == "adipose", "red", "blue"),
     xlab = "Methylated Signal (log2 medians)",
     ylab = "Unmethylated Signal (log2 medians)")
legend("bottomright", legend = c("adipose", "liver"), col = c("red", "blue"), pch = c(1,1))
dev.off()

# do normalization and study ID linear adjustments
lgr <- lmd <- lb <- lan <- list()
tv <- c("adipose", "liver")
# get noob norm data
dn <- grdn
gr <- loadHDF5SummarizedExperiment(dn)
for(t in tv){
  # do noob normalization
  lmv[[t]] <- gr[, gr$sgroup == t]
  msi <- lmv[[t]]
  madj <- limma::removeBatchEffect(getM(msi), batch = msi$gseid)
  # store adjusted data in a new se object
  lgr[[t]] <- GenomicRatioSet(GenomicRanges::granges(msi), M = madj, 
                              annotation = annotation(msi))
  # append samples metadata
  lmd[[t]] <- pData(lgr[[t]]) <- pData(lmv[[t]])
  # append preprocessing metadata
  metadata(lgr[[t]]) <- list("preprocess" = "noobbeta;removeBatchEffect_gseid")
  # make betavals list
  lb[[t]] <- getBeta(lgr[[t]]) # beta values list
}

# ANOVAs
lan <- list() # results list
for(t in c("adipose", "liver")){
  lan[[t]] <- list("pval" = matrix(ncol = 9, nrow = 0),
                   "var.fract" = matrix(ncol = 9, nrow = 0))
}
# define vars
var.gseid.adipose <- as.factor(lmd[[1]]$gseid)
var.gseid.liver <- as.factor(lmd[[2]]$gseid)
var.predsex.adipose <- as.factor(lmd[[1]]$predsex)
var.predsex.liver <- as.factor(lmd[[2]]$predsex)
var.predage.adipose <- as.numeric(lmd[[1]]$predage)
var.predage.liver <- as.numeric(lmd[[2]]$predage)
var.predcellCD8T.adipose <- as.numeric(lmd[[1]]$predcell.CD8T)
var.predcellCD8T.liver <- as.numeric(lmd[[2]]$predcell.CD8T)
var.predcellCD4T.adipose <- as.numeric(lmd[[1]]$predcell.CD4T)
var.predcellCD4T.liver <- as.numeric(lmd[[2]]$predcell.CD4T)
var.predcellNK.adipose <- as.numeric(lmd[[1]]$predcell.NK)
var.predcellNK.liver <- as.numeric(lmd[[2]]$predcell.NK)
var.predcellBcell.adipose <- as.numeric(lmd[[1]]$predcell.Bcell)
var.predcellBcell.liver <- as.numeric(lmd[[2]]$predcell.Bcell)
var.predcellMono.adipose <- as.numeric(lmd[[1]]$predcell.Mono)
var.predcellMono.liver <- as.numeric(lmd[[2]]$predcell.Mono)
var.predcellGran.adipose <- as.numeric(lmd[[1]]$predcell.Gran)
var.predcellGran.liver <- as.numeric(lmd[[2]]$predcell.Gran)
# run anovas
bv <- lb[[1]]
for(r in 1:nrow(bv)){
  for(t in tv){
    datr <- as.numeric(lb[[t]][r,])
    if(t == "adipose"){
      ld <- lm(datr ~ var.gseid.adipose + var.predsex.adipose + var.predage.adipose + var.predcellCD8T.adipose +
                 var.predcellCD4T.adipose + var.predcellNK.adipose + var.predcellBcell.adipose +
                 var.predcellMono.adipose + var.predcellGran.adipose)
    } else{
      ld <- lm(datr ~ var.gseid.liver + var.predsex.liver + var.predage.liver + var.predcellCD8T.liver +
                 var.predcellCD4T.liver + var.predcellNK.liver + var.predcellBcell.liver +
                 var.predcellMono.liver + var.predcellGran.liver)
    }
    an <- anova(ld)
    ap <- an[c(1:9),5]
    av <- round(100*an[c(1:4),2]/sum(an[,2]), 3)
    lan[[t]][["pval"]] <- rbind(lan[[t]][["pval"]], ap)
    lan[[t]][["var.fract"]] <- rbind(lan[[t]][["var.fract"]], av)
  }
  message(r)
}

# assign dimnames
rnv <- rownames(bv)
cnv <- c("gseid", "predsex", "predage", "predcell.CD8T", "predcell.CD4T",
         "predcell.NK", "predcell.Bcell", "predcell.Mono", "predcell.Gran")
colnames(lan[[1]]$pval) <- colnames(lan[[1]]$var.fract) <- colnames(lan[[2]]$pval) <- colnames(lan[[2]]$var.fract) <- cnv
rownames(lan[[1]]$pval) <- rownames(lan[[1]]$var.fract) <- rownames(lan[[2]]$pval) <- rownames(lan[[2]]$var.fract) <- rnv
save(lan, file = "lanova_2tissues.rda")

# ANOVA results probe filters
pfilt <- 1e-3
varfilt <- 10
lcgkeep <- list() # list of filtered probe sets
anno <- getAnnotation(gr)
xy.cg <- rownames(anno[anno$chr %in% c("chrY", "chrX"),]) # sex chr probes
for(t in names(lan)){
  pm <- lan[[t]]$pval
  vm <- lan[[t]]$var.fract
  # retain autosome cgids
  pm <- pm[!rownames(pm) %in% xy.cg,]
  vm <- vm[!rownames(vm) %in% xy.cg,]
  # parse variable thresholds
  cm <- as.data.frame(matrix(nrow = nrow(pm), ncol = ncol(pm)))
  for(c in 1:ncol(pm)){
    pc <- pm[,c]; 
    pc.adj <- as.numeric(p.adjust(pc, method = "BH"))
    pc.filt <- pc.adj < pfilt
    vc.filt <- vm[,c] >= varfilt
    cm[,c] <- (pc.filt & vc.filt)
  }
  cgkeep <- apply(cm, 1, function(x){return((length(x[x == TRUE]) == 0))})
  lcgkeep[[t]] <- rownames(pm)[cgkeep]
}
save(lcgkeep, file = "lcgkeep_2tissues.rda")
lgr.filt <- list("adipose" = lgr[[1]][lcgkeep[[1]],],
                 "liver" = lgr[[2]][lcgkeep[[2]],])
save(lgr.filt, file = "lgratio-filt_2tissues.rda")

# barplots of probes retained and removed
require(ggplot2)
cgtot <- 485512
bpdf <- data.frame(tissue = c(rep("adipose", 2), rep("liver",2)),
                   nprobes = c(cgtot - length(lcgkeep[[1]]), length(lcgkeep[[1]]), 
                               cgtot - length(lcgkeep[[2]]), length(lcgkeep[[2]])),
                   type = rep(c("removed", "retained"), 2), stringsAsFactors = F)
jpeg("barplot_2tissues.jpg", 3, 2, units = "in", res = 400)
ggplot(bpdf, aes(x = tissue, y = nprobes, fill = type)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

# DNAm summary statistics for filtered probe sets
lgr.filt <- get(load("lgratio-filt_2tissues.rda"))
# min, max, mean, median, sd, var
tv <- c("adipose", "liver")
cnv <- c("min", "max", "mean", "median", "sd", "var")
mt <- matrix(nrow = 0, ncol = 6)
lcg.ss <- list(mt, mt); names(lcg.ss) <- tv
for(t in tv){
  bt <- as.matrix(getBeta(lgr.filt[[t]]))
  for(r in 1:nrow(bt)){
    x <- bt[r,]
    newrow <- c(min(x), max(x), mean(x), median(x), sd(x), var(x))
    lcg.ss[[t]] <- rbind(lcg.ss[[t]], newrow)
    message(r)
  }
  #bm <- apply(bt, 1, function(x){
  #  newrow <- c(min(x), max(x), mean(x), median(x), sd(x), var(x))
  #  lcg.ss[[t]] <- rbind(lcg.ss[[t]], newrow)
  #})
  colnames(lcg.ss[[t]]) <- cnv
  rownames(lcg.ss[[t]]) <- rownames(bt)
}
save(lcg.ss, file = "lcg-filt_summarystats_2tissues.rda")

# probe variance analysis -- absolute quantile method
lfilt <- get(load("lcg-filt_summarystats_2tissues.rda"))
qiv = seq(0, 1, 0.01)
qwhich = c(100)
lmvp.abs <- list()
lci <- list()
for(t in 1:length(lfilt)){
  ba = as.data.frame(lfilt[[t]], stringsAsFactors = FALSE)
  q <- qf <- quantile(ba$var, qiv)[qwhich]
  lci[[names(q)]] <- rownames(ba)[ba$var > q]
  lmvp.abs[[names(lfilt)[t]]] = lci
}
save(lmvp.abs, file = "lmvp-absquantfilt_2tissues.rda")

# probe variance analysis -- binned quantiles method
# quantile bin method
lfilt <- get(load("lcg-filt_summarystats_2tissues.rda"))
qiv = seq(0, 1, 0.01) # quantile filter
qwhich = c(100)
bin.xint <- 0.01
binv = seq(0, 1, bin.xint)[1:100] # binned bval mean
# iter on ncts
lmvp.bin = list()
for(t in 1:length(lfilt)){
  ba = as.data.frame(lfilt[[t]], stringsAsFactors = FALSE)
  lcg = list()
  # iterate on betaval bins
  for(b in binv){
    # get probes in bin
    bf = ba[ba$mean >= b & ba$mean < b + bin.xint, ] 
    # do bin-specific quantile filter
    q <- qf <- quantile(bf$var, qiv)[qwhich]
    # append probes list
    lcg[[names(q)]] = c(lcg[[names(q)]], rownames(bf)[bf$var > q])
  }
  names(lcg) = paste0("ci:", names(q))
  lmvp.bin[[names(lfilt)[t]]] = lcg
}
save(lmvp.bin, file = "lmvp-binquantfilt_2tissues.rda")

# categorize high-variance probes
cgav <- c()
for(t in 1:length(labs)){
  txcg <- unique(c(labs[[t]][[1]], 
                   lbin[[t]][[1]]))
  cgav <- c(cgav, txcg)
}
cgdf <- as.data.frame(table(cgav))
cgdf$type <- ifelse(cgdf[,2] > 1, 
                    "non-specific", 
                    "tissue-specific")
table(cgdf$type)
save(cgdf, file = "cgtypes_2tissues.rda")

# select tissue-specific probes
lfilt <- get(load("lcg-filt_summarystats_2tissues.rda"))
cgdf <- get(load("cgtypes_2tissues.rda"))
cgfilt <- cgdf$type == "non-specific"
cgdff <- cgdf[!cgfilt,]
ltxcg <- list()
for(t in 1:length(labs)){
  cgtx <- c()
  cgabs <- labs[[t]][[1]]
  cgbin <- lbin[[t]][[1]]
  bt <- as.data.frame(lfilt[[t]])
  # get t tissue specific probes
  filtbt <- rownames(bt) %in% cgdff[,1]
  bt <- bt[filtbt,]
  # get top 1k t tissue specific abs probes
  filt.bf1 <- rownames(bt) %in% cgabs
  bf1 <- bt[filt.bf1,]
  bf1 <- bf1[rev(order(bf1$var)),]
  cgtx <- rownames(bf1)[1:1000]
  # get top 1k t tissue specific bin probes, after filt
  filt.bf2 <- rownames(bt) %in% cgbin &
    !rownames(bt) %in% rownames(bf1)
  bf2 <- bt[filt.bf2,]
  bf2 <- bf2[rev(order(bf2$var)),]
  cgtx <- c(cgtx, rownames(bf2)[1:1000])
  ltxcg[[names(labs)[t]]] <- cgtx
}

# get filtered cg summaries
lfcg <- lapply(lfilt, function(x){x <- x[rownames(x) %in% unique(unlist(ltxcg)),]})
save(ltxcg, file = "ltxcg_topvar2k_2tissues.rda")
# get annotation subset
anno <- getAnnotation(lgr[[1]]) # save anno for cga
anno <- anno[,c("Name", "UCSC_RefGene_Name", "Relation_to_Island")]
anno <- anno[rownames(anno) %in% unique(unlist(ltxcg)),]

# violin plots of probe set means, vars across tissues
lvp <- makevp(lfcg, ltxcg)
jpeg("violinplots_mean-var_2tissues.jpg", 5, 5, units = "in", res = 400)
grid.arrange(lvp[[1]], lvp[[2]], ncol = 1, bottom = "Tissue")
dev.off()

# heatmaps of region-specific DNAm among probe sets
cga <- get_cga(anno)
lhmset <- hmsets(ltxcg, cga)
lhmplot <- hmplots(hmma.mean = lhmset[[1]], hmma.var = lhmset[[2]], 
                   hmma.size = lhmset[[3]])
jpeg("heatmaps_mean-var_2tissues.jpg", 9, 4, units = "in", res = 400)
grid.arrange(lhmplot$hm.mean.plot, lhmplot$hm.var.plot, 
             layout_matrix = matrix(c(rep(1, 7), rep(2, 4)), nrow = 1),
             bottom = "Tissue", left = "Annotation/Region Type")
dev.off()

#-----------
# rdata file
#-----------
# make the .RData file for vignette to load
