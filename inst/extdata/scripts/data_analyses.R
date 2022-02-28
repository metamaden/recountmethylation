#!/usr/bin/env R

# Script to generate env objects used in data_analyses.Rmd
# Notes:
# Currently running this with hdf5 db
# Rerun with h5se RGChannelSet when this is working (e.g. syncs to server, etc.)

library(recountmethylation)
library(rhdf5)
library(HDF5Array)
library(minfi)
library(limma)
library(GenomicRanges)
library(ggplot2)
library(gridExtra)

dfp <- "data_analyses"
env.name <- "data_analyses.RData"
dir.create(dfp)
savepath <- paste(dfp, env.name, sep = "/")

#------------
# global vars
#------------
# get local metadata
path <- system.file("extdata", "metadata", package = "recountmethylation")
mdpath <- paste(path, list.files(path)[1], sep = "/")
md <- get(load(mdpath))
dim(md) # [1] 35360    19

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
                     stringsAsFactors = FALSE)
    bpdf.mean = rbind(bpdf.mean, dt[,c(1, 3)])
    bpdf.var = rbind(bpdf.var, dt[,c(2, 3)])
  }
  bpdf.mean = as.data.frame(bpdf.mean, stringsAsFactors = FALSE)
  bpdf.var = as.data.frame(bpdf.var, stringsAsFactors = FALSE)
  bpdf.mean[,1] = as.numeric(as.character(bpdf.mean[,1]))
  bpdf.var[,1] = as.numeric(as.character(bpdf.var[,1]))
  colnames(bpdf.mean) = c("mean", "tissue")
  colnames(bpdf.var) = c("var", "tissue")
  
  vp1 = ggplot(bpdf.mean, aes(x = tissue, y = mean, fill = tissue)) + 
    geom_violin(trim = FALSE, show.legend = FALSE, draw_quantiles = c(0.5)) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90)) + 
    ggtitle("") + xlab("") + ylab("Mean")
  
  vp2 = ggplot(bpdf.var, aes(x = tissue, y = var, fill = tissue)) + 
    geom_violin(trim = FALSE, show.legend = FALSE, draw_quantiles = c(0.5)) + 
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
hmsets <- function(ltxcg, lcg.ss, cga, mincg = 2){
  ngrp <- length(ltxcg)
  annom = matrix(nrow = 0, ncol = ngrp)
  # summary stats -- annotated region overlaps with high-var cgs
  utype = unique(cga$type.composite)
  hmma.mean = hmma.var = hmma.size = matrix(nrow = 0, ncol = ngrp)
  for(a in utype){
    cgidv = rownames(cga[cga$type.composite==a,])
    rdat.mean = rdat.var = rdat.size = c()
    for(t in 1:ngrp){
      # get intersect with txmvp and type
      cgint = intersect(ltxcg[[t]], cgidv)
      cgint.size = length(cgint)
      cgsst <- as.data.frame(lcg.ss[[t]])
      if(cgint.size >= mincg){
        cg.select <- rownames(cgsst) %in% cgint
        rdat.mean = c(rdat.mean, mean(cgsst[cg.select,]$mean))
        rdat.var = c(rdat.var, var(cgsst[cg.select,]$mean))
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
    message(a)
  }
  rownames(hmma.mean) <- rownames(hmma.var) <- rownames(hmma.size) <- utype
  colnames(hmma.mean) <- colnames(hmma.var) <- colnames(hmma.size) <- names(lcg.ss)
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
  hdx = as.data.frame(hdx, stringsAsFactors = FALSE)
  hdv = as.data.frame(hdv, stringsAsFactors = FALSE)
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
                         limits = c(0, max(hdv$var))) +
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
{
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
  jpeg("mainfig1b_agedif-mfilt.jpg", 3, 3, units = "in", res = 400)
  ggplot(mdff, aes(x = chron.age, y = predage)) +
    geom_point(size = 1.2, alpha = 0.2) + 
    geom_smooth(method = "lm", size = 1.2) +
    theme_bw() + xlab("Chronological Age") + ylab("Epigenetic (DNAm) Age")
  dev.off()
}


#-------------------------------------------------------
# example 2 -- FFPE versus Frozen tissue quality signals
#-------------------------------------------------------
{
  # get samples with storage condition 
  mdf <- md[!md$storage == "NA",]
  gmf <- gm[, gm$gsm %in% mdf$gsm]
  mdf <- mdf[order(match(mdf$gsm, gmf$gsm)),]
  identical(gmf$gsm, mdf$gsm)
  # add complete available storage info
  gmf$storage <- mdf$storage
  table(gmf$storage)
  
  # get signal matrices
  meth.all <- getMeth(gmf)
  unmeth.all <- getUnmeth(gmf)
  
  # get blocks for processing
  blocks <- getblocks(slength = ncol(gmf), bsize = 200)
  
  # process data in blocks
  ms <- matrix(nrow = 0, ncol = 2)
  l2meth <- l2unmeth <- c()
  for(i in 1:length(blocks)){
    b <- blocks[[i]]
    gmff <- gmf[, b]
    methb <- as.matrix(meth.all[, b])
    unmethb <- as.matrix(unmeth.all[, b])
    l2meth <- apply(methb, 2, function(x){
      log2(median(as.numeric(x)))
    })
    l2unmeth <- apply(unmethb, 2, function(x){
      log2(median(as.numeric(x)))
    })
    ms <- rbind(ms, matrix(c(l2meth, l2unmeth), ncol = 2))
    message(i)
  }
  rownames(ms) <- colnames(meth.all)
  colnames(ms) <- c("meth.l2med", "unmeth.l2med")
  ds <- as.data.frame(ms)
  ds$storage <- ifelse(grepl("FFPE", gmf$storage), "ffpe", "frozen")
  save(ds, file = file.path("data_analyses", "df-l2med-signals.rda"))
  
  # signals scatterplot by storage type
  jpeg("scatterplot_storage.jpg", 7, 7, units = "in", res = 400)
  ggplot(ds, aes(x = meth.l2med, y = unmeth.l2med, color = storage)) + 
    geom_point() + theme_bw() +
    scale_color_manual(values = c("ffpe" = "orange", "frozen" = "purple"))
  dev.off()
  
  # violin plots of signals by storage type
  vp <- matrix(nrow = 0, ncol = 2)
  vp <- rbind(vp, matrix(c(ds$meth.l2med, paste0("meth.", ds$storage)), ncol = 2))
  vp <- rbind(vp, matrix(c(ds$unmeth.l2med, paste0("unmeth.", ds$storage)), ncol = 2))
  vp <- as.data.frame(vp, stringsAsFactors = FALSE)
  vp[,1] <- as.numeric(vp[,1])
  colnames(vp) <- c("signal", "group")
  vp$col <- ifelse(grepl("ffpe", vp$group), "orange", "purple")
  
  jpeg("vp_storage.jpg", 5, 3, units = "in", res = 400)
  ggplot(vp, aes(x = group, y = signal, color = group)) + 
    scale_color_manual(values = c("meth.ffpe" = "orange", "unmeth.ffpe" = "orange", 
                                  "meth.frozen" = "purple", "unmeth.frozen" = "purple")) +
    geom_violin(draw_quantiles = c(0.5)) + theme_bw() + theme(legend.position = "none")
  dev.off()
}


#-----------------------------------------------------------------
# example 3 -- Tissue-specific DNAm variation in liver and adipose
#-----------------------------------------------------------------
# get samples GSM IDs
adipose.gsmv <- c('GSM1505062','GSM1505031','GSM1505051','GSM2781519','GSM3455770','GSM1505058','GSM1505035','GSM1505207','GSM1505169','GSM1505061','GSM1505055','GSM2781515','GSM1505168','GSM1505187','GSM1505189','GSM1505180','GSM1505166','GSM1505185','GSM1505172','GSM1505016','GSM2781541','GSM1505079','GSM1505153','GSM1505167','GSM1505149','GSM1505069','GSM1505148','GSM1505151','GSM1505082','GSM1505056','GSM1505203','GSM1505152','GSM1505032','GSM1505191','GSM1505176','GSM2781501','GSM3455772','GSM1505059','GSM1505182','GSM1505177','GSM1505076','GSM1505196','GSM1505043','GSM1505027','GSM1505019','GSM2781518','GSM1505078','GSM1505178','GSM1505164','GSM1505184','GSM1505050','GSM1505018','GSM1505048','GSM2781506','GSM1505150','GSM1505197','GSM1505204','GSM1505213','GSM1505026','GSM1505192','GSM1505210','GSM2781502','GSM2781512','GSM2781529','GSM1505044','GSM1505072','GSM1505147','GSM1505206','GSM1505060','GSM1505080','GSM1505024','GSM1505075','GSM1505188','GSM2781497','GSM1505049','GSM1505053','GSM1505209','GSM1505064','GSM1505183','GSM1505181','GSM1505022','GSM2781504','GSM2781546','GSM2781526','GSM2781514','GSM1505190','GSM1505173','GSM2781536','GSM1505195','GSM1505174','GSM2781543','GSM1505023','GSM1505034','GSM1505170','GSM1505054','GSM1505038','GSM1505199','GSM1505162','GSM1505161','GSM2781492','GSM2781513','GSM1505057','GSM2781509','GSM1505171')
liver.gsmv <- c('GSM2859944','GSM1504960','GSM1504940','GSM1504980','GSM2859948','GSM3455774','GSM2859954','GSM1586522','GSM2859937','GSM1586536','GSM2859962','GSM2859957','GSM1504966','GSM2859941','GSM1504988','GSM2859959','GSM1504946','GSM1504977','GSM2859942','GSM1586530','GSM1504949','GSM1504948','GSM2084821','GSM1504937','GSM2770836','GSM2859955','GSM2859972','GSM1504994','GSM1504972','GSM1586520','GSM2859960','GSM1504959','GSM1586535','GSM2859952','GSM2859949','GSM2859958','GSM1586525','GSM1504933','GSM2859974','GSM1504965','GSM2859950','GSM1504978','GSM1504962','GSM2859961','GSM2859966','GSM1504931','GSM2859970','GSM1586534','GSM2770838','GSM1504934','GSM2770841','GSM1586515','GSM1504989','GSM1586528','GSM2859946','GSM1504955','GSM2859968','GSM2859973','GSM1586529','GSM2770832','GSM2859953','GSM2859945','GSM1586526','GSM2859964','GSM1504969','GSM1504945','GSM2770840','GSM2859963','GSM1504928','GSM2770833','GSM1586531','GSM1504973','GSM1504967','GSM1586532','GSM1504964','GSM1504954','GSM2859967','GSM2859965','GSM1586524','GSM1586519','GSM1504947','GSM2859943','GSM2859976','GSM1504970','GSM2859956','GSM2859940','GSM1504957','GSM1586521','GSM2770837','GSM1504961','GSM1586514','GSM1504968','GSM2859971','GSM2770839','GSM1586523','GSM1504938','GSM1504990','GSM2770834','GSM1504952','GSM1504956','GSM1647847','GSM2859975','GSM1504971','GSM2859947','GSM1586527','GSM1504985','GSM1586518','GSM2770835','GSM1586513','GSM1504935','GSM2859969','GSM2859939')
gsmv <- c(adipose.gsmv, liver.gsmv)
mdf <- md[md$gsm %in% gsmv,]
mdf$sgroup <- ifelse(mdf$gsm %in% adipose.gsmv, "adipose", "liver")
table(mdf$sgroup)
# adipose   liver 
# 104     112 

# get metadata summaries
sst.tvar <- get_sst(sgroup.labs = c("liver", "adipose"), mdf)

# subset gm, append sgroup, map to genome
ms <- gm[, gm$gsm %in% mdf$gsm]
ms <- ms[, order(match(ms$gsm, mdf$gsm))]
identical(ms$gsm, mdf$gsm)
ms$sgroup <- mdf$sgroup
ms <- mapToGenome(ms)
dim(ms)

# get log2 medians
meth.tx <- getMeth(ms)
unmeth.tx <- getUnmeth(ms)
blocks <- getblocks(slength = ncol(ms), bsize = 10)
# process data in blocks
l2m <- matrix(nrow = 0, ncol = 2)
for(i in 1:length(blocks)){
  b <- blocks[[i]]
  gmff <- ms[, b]
  methb <- as.matrix(meth.tx[, b])
  unmethb <- as.matrix(unmeth.tx[, b])
  l2meth <- l2unmeth <- c()
  l2meth <- c(l2meth, apply(methb, 2, function(x){
    log2(median(as.numeric(x)))
  }))
  l2unmeth <- c(l2unmeth, apply(unmethb, 2, function(x){
    log2(median(as.numeric(x)))
  }))
  l2m <- rbind(l2m, matrix(c(l2meth, l2unmeth), ncol = 2))
  message(i)
}
lqc <- list("l2med.meth" = l2m[,1], "l2med.unmeth" = l2m[,2], 
            "gsmv" = ms$gsm, "sgroup" = ms$sgroup)

# plot signal data
plot(lqc[["l2med.meth"]], lqc[["l2med.unmeth"]], 
     col = ifelse(lqc[["sgroup"]] == "adipose", "red", "blue"),
     xlab = "Methylated Signal (log2 medians)",
     ylab = "Unmethylated Signal (log2 medians)")
legend("bottomright", legend = c("adipose", "liver"), col = c("red", "blue"), pch = c(1,1))

# ggplot signal data
ds <- as.data.frame(do.call(cbind, lqc))
ds$tissue <- as.factor(ds$sgroup)
ds$l2med.meth <- as.numeric(ds$l2med.meth)
ds$l2med.unmeth <- as.numeric(ds$l2med.unmeth)
ggplot(ds2, aes(x = l2med.meth, y = l2med.unmeth, color = tissue)) + 
  geom_point(alpha = 0.3, cex = 3) + theme_bw()

# do study ID adjustment
lmv <- lgr <- lb <- list()
tv <- c("adipose", "liver")
# get noob norm data
gr <- gr[,colnames(gr) %in% colnames(ms)]
gr <- gr[,order(match(colnames(gr), colnames(ms)))]
identical(colnames(gr), colnames(ms))
gr$sgroup <- ms$sgroup
# do study ID adj
for(t in tv){
  msi <- gr[, gr$sgroup == t]
  madj <- limma::removeBatchEffect(getM(msi), batch = msi$gseid)
  lgr[[t]] <- GenomicRatioSet(GenomicRanges::granges(msi), M = madj, 
                              annotation = annotation(msi))
  metadata(lgr[[t]]) <- list("preprocess" = "noobbeta;removeBatchEffect_gseid")
  pData(lgr[[t]]) <- pData(gr[, gr$sgroup == t])
  lb[[t]] <- getBeta(lgr[[t]])
}

# prepare ANOVAs
# get autosomal probes
anno <- getAnnotation(gr)
chr.xy <-c("chrY", "chrX")
cg.xy <- rownames(anno[anno$chr %in% chr.xy,])
lbf <- list()
for(t in tv){
  bval <- lb[[t]]
  lbf[[t]] <- bval[!rownames(bval) %in% cg.xy,]
}
bv <- lbf[[1]]
# define ANOVA vars
lvar <- list()
cnf <- c("gseid", "predsex", "predage", "predcell.CD8T",
         "predcell.CD4T", "predcell.NK", "predcell.Bcell",
         "predcell.Mono", "predcell.Gran")
for(t in tv){
  for(c in cnf){
    if(c %in% c("gseid", "predsex")){
      lvar[[t]][[c]] <- as.factor(pData(lgr[[t]])[,c])
    } else{
      lvar[[t]][[c]] <- as.numeric(pData(lgr[[t]])[,c])
    }
  }
}
# set up ANOVAs
bv <- lbf[[1]]
blocks <- getblocks(slength = nrow(bv), bsize = 100000)
lan <- list("adipose" = matrix(nrow = 0, ncol = 18),
            "liver" = matrix(nrow = 0, ncol = 18))
# run ANOVAs with vectorization
t1 <- Sys.time()
for(bi in 1:length(blocks)){
  for(t in tv){
    datr <- lbf[[t]][blocks[[bi]],]
    tvar <- lvar[[t]]
    newchunk <- t(apply(datr, 1, function(x){
      # do multiple regression and anova
      x <- as.numeric(x)
      ld <- lm(x ~ tvar[[1]] + tvar[[2]] + tvar[[3]] + tvar[[4]] +
                 tvar[[5]] + tvar[[6]] + tvar[[7]] + tvar[[8]] + tvar[[9]])
      an <- anova(ld)
      # get results
      ap <- an[c(1:9),5] # pval
      av <- round(100*an[c(1:9),2]/sum(an[,2]), 3) # percent var
      return(as.numeric(c(ap, av)))
    }))
    # append new results
    lan[[t]] <- rbind(lan[[t]], newchunk)
  }
  message(bi, "tdif: ", Sys.time() - t1)
}
# append colnames
for(t in tv){colnames(lan[[t]]) <- rep(cnf, 2)}
save(lan, file = "lanova_2tissues.rda")

# ANOVA results probe filters
pfilt <- 1e-3
varfilt <- 10
lcgkeep <- list() # list of filtered probe sets
for(t in tv){
  pm <- lan[[t]][,c(1:9)]
  vm <- lan[[t]][,c(10:18)]
  # parse variable thresholds
  cm <- as.data.frame(matrix(nrow = nrow(pm), ncol = ncol(pm)))
  for(c in 1:ncol(pm)){
    pc <- pm[,c]; 
    pc.adj <- as.numeric(p.adjust(pc))
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

# barplots of probes retained and removed
cgtot <- 485512
bpdf <- data.frame(tissue = c(rep("adipose", 2), rep("liver",2)),
                   nprobes = c(cgtot - length(lcgkeep[[1]]), length(lcgkeep[[1]]), 
                               cgtot - length(lcgkeep[[2]]), length(lcgkeep[[2]])),
                   type = rep(c("removed", "retained"), 2), stringsAsFactors = FALSE)
ggplot(bpdf, aes(x = tissue, y = nprobes, fill = type)) +
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

# DNAm summary statistics for filtered probe sets (new)
tv <- c("adipose", "liver")
cnv <- c("min", "max", "mean", "median", "sd", "var")
bv <- getBeta(lgr.filt[[t]])
lbt <- lcg.ss <- list()
bsize = 100000
for(t in tv){
  lcg.ss[[t]] <- matrix(nrow = 0, ncol = 6)
  lbt[[t]] <- bt <- as.matrix(getBeta(lgr.filt[[t]]))
  blockst <- getblocks(slength = nrow(bt), bsize = bsize)
  for(bi in 1:length(blockst)){
    bc <- bt[blockst[[bi]],]
    newchunk <- t(apply(bc, 1, function(x){
      newrow <- c(min(x), max(x), mean(x), median(x), sd(x), var(x))
      return(as.numeric(newrow))
    }))
    lcg.ss[[t]] <- rbind(lcg.ss[[t]], newchunk)
    message(t, ";", bi)
  }
  colnames(lcg.ss[[t]]) <- cnv
}

# variance analysis method1 -- absolute quantiles
qiv = seq(0, 1, 0.01)
qwhich = c(100)
lmvp.abs <- list()
lci <- list()
for(t in tv){
  cgv <- c()
  sa <- lcg.ss[[t]]
  sa <- as.data.frame(sa, stringsAsFactors = FALSE)
  q <- quantile(sa$var, qiv)[qwhich]
  lmvp.abs[[t]] <- rownames(sa[sa$var > q,])
}

# variance analysis method2 -- binned quantiles
qiv = seq(0, 1, 0.01) # quantile filter
qwhich = c(100)
bin.xint <- 0.1
binv = seq(0, 1, bin.xint)[1:10] # binned bval mean
# iter on ncts
lmvp.bin = list()
for(t in tv){
  sa <- as.data.frame(lcg.ss[[t]])
  cgv <- c()
  # iterate on betaval bins
  for(b in binv){
    bf <- sa[sa$mean >= b & sa$mean < b + bin.xint, ] # get probes in bin
    q <- qf <- quantile(bf$var, qiv)[qwhich] # do bin filter
    cgv <- c(cgv, rownames(bf)[bf$var > q]) # append probes list
  }
  lmvp.bin[[t]] <- cgv
}

# get cg specificity
cgav <- c()
for(t in tv){
  txcg <- unique(c(lmvp.abs[[t]], lmvp.bin[[t]]))
  cgav <- c(cgav, txcg)
}
cgdf <- as.data.frame(table(cgav))
cgdf$type <- ifelse(cgdf[,2] > 1, "non-specific", "tissue-specific")
table(cgdf$type)

# cg specificity filter
cgfilt <- cgdf$type == "non-specific"
cgdff <- cgdf[!cgfilt,]
ltxcg <- list()
for(t in tv){
  cgtx <- c()
  cgabs <- lmvp.abs[[t]]
  cgbin <- lmvp.bin[[t]]
  st <- as.data.frame(lcg.ss[[t]])
  # get t tissue specific probes
  filtbt <- rownames(st) %in% cgdff[,1]
  st <- st[filtbt,]
  # get top 1k t tissue specific abs probes
  filt.bf1 <- rownames(st) %in% cgabs
  sf1 <- st[filt.bf1,]
  sf1 <- sf1[rev(order(sf1$var)),]
  cgtx <- rownames(sf1)[1:1000]
  # get top 1k t tissue specific bin probes, after filt
  filt.bf2 <- rownames(st) %in% cgbin &
    !rownames(st) %in% rownames(sf1)
  sf2 <- st[filt.bf2,]
  sf2 <- sf2[rev(order(sf2$var)),]
  cgtx <- c(cgtx, rownames(sf2)[1:1000])
  ltxcg[[t]] <- cgtx
}

# probe set summaries, annotations
# filtered cg summaries
lfcg <- lapply(lcg.ss, function(x){x <- x[rownames(x) %in% unique(unlist(ltxcg)),]})
# annotation subset
anno <- getAnnotation(gr) # save anno for cga
anno <- anno[,c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island")]
anno <- anno[rownames(anno) %in% unique(unlist(ltxcg)),]
# filtered beta values
lcgssf <- list()
for(t in tv){
  bv <- lcg.ss[[t]]
  bvf <- bv[rownames(bv) %in% ltxcg[[t]],]
  lcgssf[[t]] <- bvf
}

# means table of probe set statistics by tissue group
tcgss <- matrix(nrow = 0, ncol = 6)
for(t in tv){
  datt <- apply(lcgssf[[t]], 2, function(x){
    round(mean(x), digits = 2)
  })
  mt <- matrix(datt, nrow = 1)
  tcgss <- rbind(tcgss, mt)
}
colnames(tcgss) <- colnames(lcgssf$adipose)
rownames(tcgss) <- tv
kable(t(tcgss))

# violin plots
lvp <- makevp(lfcg, ltxcg)
grid.arrange(lvp[[1]], lvp[[2]], ncol = 1, bottom = "Tissue")

# heatmaps
cga <- get_cga(anno)
lhmset <- hmsets(ltxcg, lfcg, cga)
lhmplots <- hmplots(lhmset$hm.mean, lhmset$hm.var, lhmset$hm.size)
grid.arrange(lhmplots$hm.mean.plot, lhmplots$hm.var.plot, 
             layout_matrix = matrix(c(rep(1, 7), rep(2, 4)), nrow = 1),
             bottom = "Tissue", left = "Annotation/Region Type")

#-----------
# rdata file
#-----------
# make the .RData file for vignette to load
save(list = c("tv", "ds", "ds2", "lfcg", "ltxcg", "anno", 
              "lcgssf", "adipose.gsmv", "liver.gsmv", "cgdf",
              "lmvp.abs", "lmvp.bin", "get_sst", 
              "getblocks", "makevp", "get_cga", "hmsets", "hmplots"), 
     file = "data_analyses.RData")
