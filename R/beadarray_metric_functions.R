#!/usr/bin/env R

# Functions for BeadArray metric calculations.

#' ba.background
#'
#' Get background signal for BeadArray metric calculations.
#' 
#' @param cdf Control probe data frame.
#' @param ct Column name/metric title.
#' @rdname ba
ba.background <- function(cdf, ct = "Extension"){
  # Notes: single metric, uses extension grn channel
  message("Getting background signal...")
  ci <- cdf[grepl(ct, cdf$ExtendedType),]
  if(nrow(ci) == 0){return("NA")}
  addr.bkg <- ci[grepl("(A)|(T)", ci$ExtendedType),]$Address
  return(addr.bkg)
}

#' ba.biotinstaining.red
#' 
#' Get the Biotin staining Red BeadArray metric.
#' 
#' @param rs Red signal matrix (required).
#' @param biotin.baseline Baseline signal for the biotin stain 
#' assay (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @param ct Column name/metric title (character, required).
#' @rdname ba
ba.biotinstaining.red <- function(rs, biotin.baseline, rm, cdf, 
                                  mnum = 1, mtot = 17,
                                  cnamei = "biotin.stain.red", 
                                  ct = "Biotin|DNP"){
  message("(",mnum,"/",mtot,") Calculating Biotin Staining Red metric...")
  ci <- cdf[grepl(ct, cdf$ExtendedType),]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  which.stain <- which(rownames(rs) %in% 
                         ci[grepl("DNP \\(High", ci$ExtendedType),]$Address)
  which.bkg <- which(rownames(rs) %in% 
                       ci[grepl("DNP \\(Bkg", ci$ExtendedType),]$Address)
  m1 <- apply(rs, 2, function(x){x[which.stain]/(x[which.bkg] + biotin.baseline)})
  rm <- cbind(rm, m1); colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}

#' ba.biotinstaining.grn
#' 
#' Get the Biotin staining Green BeadArray metric.
#' 
#' 
#' @param gs Green signal matrix (required).
#' @param biotin.baseline Baseline signal for the biotin stain 
#' assay (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @param ct Column name/metric title (character, required).
#' @rdname ba
ba.biotinstaining.grn <- function(gs, biotin.baseline, 
                                  rm, cdf, mnum = 2, 
                                  mtot = 17, 
                                  cnamei = "biotin.stain.grn", 
                                  ct = "Biotin|DNP"){
  message("(",mnum,"/",mtot,") Calculating Biotin Staining Green metric...")
  ci <- cdf[grepl(ct, cdf$ExtendedType),]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  which.stain <- which(rownames(gs) %in% 
                         ci[grepl("Biotin \\(High", ci$ExtendedType),]$Address)
  which.bkg <- which(rownames(gs) %in% 
                       ci[grepl("Biotin \\(Bkg", ci$ExtendedType),]$Address)
  m2 <- apply(gs, 2, function(x){
    x[which.stain]/(x[which.bkg] + biotin.baseline)})
  rm <- cbind(rm, m2)
  colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}

#' ba.nonpolymorphic.red
#' 
#' Get the Non-polymorphic Red BeadArray metric.
#' 
#' @param rs Red signal matrix (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @param ct Column name/metric title (character, required).
#' @rdname ba
ba.nonpolymorphic.red <- function(rs, rm, cdf, mnum = 3, 
                                  mtot = 17, 
                                  cnamei = "nonpolymorphic.red", 
                                  ct = "NP"){
  message("(",mnum,"/",mtot,") Calculating Non-polymorphic Red metric...")
  ci <- cdf[grepl(ct, cdf$ExtendedType),]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  which.cg <- which(rownames(rs) %in% ci$Address[grepl("(C)|(G)", ci$ExtendedType)])
  which.at <- which(rownames(rs) %in% ci$Address[grepl("(A)|(T)", ci$ExtendedType)])
  m1 <- apply(rs, 2, function(x){min(x[which.at])/max(x[which.cg])})
  rm <- cbind(rm, m1); colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}

#' ba.nonpolymorphic.grn
#' 
#' Get the Non-polymorphic Green BeadArray metric.
#' 
#' @param gs Green signal matrix (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @param ct Column name/metric title (character, required).
#' @rdname ba
ba.nonpolymorphic.grn <- function(gs, rm, cdf, mnum = 4, mtot = 17, 
  cnamei = "nonpolymorphic.grn", ct = "NP"){
  message("(",mnum,"/",mtot,") Calculating Non-polymorphic Grn metric...")
  ci <- cdf[grepl(ct, cdf$ExtendedType),]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  which.cg <- which(rownames(gs) %in% 
                      ci$Address[grepl("(C)|(G)", ci$ExtendedType)])
  which.at <- which(rownames(gs) %in% 
                      ci$Address[grepl("(A)|(T)", ci$ExtendedType)])
  m2 <- apply(gs, 2, function(x){
    min(x[which.cg])/max(x[which.at])})
  rm <- cbind(rm, m2)
  colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}

#' ba.bisulfiteconv1.red
#' 
#' Get the Bisulfite Conversion I Red BeadArray metric.
#' 
#' @param rs Red signal matrix (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @param ct Column name/metric title (character, required).
#' @rdname ba
ba.bisulfiteconv1.red <- function(rs, rm, cdf, mnum = 5, mtot = 17, 
  cnamei = "bisulfite.conv1.red", ct = "Conversion I-"){
  message("(",mnum,"/",mtot,") Calculating Bisulfite Conversion I metric...")
  ci <- cdf[grepl(ct, cdf$ExtendedType),]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  which.c456 <- which(rownames(rs) %in% 
                        ci$Address[grepl("C4|C5|C6", ci$ExtendedType)])
  which.u456 <- which(rownames(rs) %in% 
                        ci$Address[grepl("U4|U5|U6", ci$ExtendedType)])
  m1 <- apply(rs, 2, function(x){
    min(x[which.c456])/max(x[which.u456])})
  rm <- cbind(rm, m1)
  colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}

#' ba.restoration
#' 
#' Get the Restoration BeadArray metric.
#' 
#' @param gs Green signal matrix (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param addr.bkg Background signal probe address (required).
#' @param baseline Baseline measure for signals (integer, 3000).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @param ct Column name/metric title (character, required).
#' @rdname ba
ba.restoration <- function(gs, rm, cdf, addr.bkg, baseline = 3000, mnum = 6, 
  mtot = 17, cnamei = "restoration.grn", ct = "RESTORATION"){
  # note: single metric, uses only grn channel
  message("(",mnum,"/",mtot,") Calculating Restoration metric...")
  ci <- cdf[cdf$Type == ct,]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  which.rest <- which(rownames(gs) %in% ci$Address)
  m1 <- apply(gs, 2, function(x){
    x[which.rest]/(max(x[addr.bkg]) + baseline)})
  rm <- cbind(rm, m1)
  colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}

#' ba.specificity1.red
#' 
#' Get the Specificity I Red BeadArray metric.
#' 
#' @param rs Red signal matrix (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @rdname ba
ba.specificity1.red <- function(rs, rm, cdf, mnum = 7, mtot = 17, 
  cnamei = "specificityI.red"){
  message("(",mnum,"/",mtot,") Calculating Specificity I Red metric...")
  ci <- cdf[grepl("Mismatch (1|2|3)", cdf$ExtendedType),]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  addr.mm.index <- which(rownames(rs) %in% 
                           ci[grepl("MM", ci$ExtendedType),]$Address)
  addr.pm.index <- which(rownames(rs) %in% 
                           ci[grepl("PM", ci$ExtendedType),]$Address)
  m1 <- apply(rs, 2, function(x){
    min(x[addr.pm.index])/max(x[addr.mm.index])})
  rm <- cbind(rm, m1)
  colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}

#' ba.specificity1.grn
#' 
#' Get the Specificity I Green BeadArray metric.
#' 
#' @param gs Green signal matrix (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @rdname ba
ba.specificity1.grn <- function(gs, rm, cdf, mnum = 8, mtot = 17, 
  cnamei = "specificityI.grn"){
  message("(",mnum,"/",mtot,") Calculating Specificity I...")
  ci <- cdf[grepl("Mismatch (1|2|3)", cdf$ExtendedType),]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  addr.mm.index <- which(rownames(gs) %in% 
                           ci[grepl("MM", ci$ExtendedType),]$Address)
  addr.pm.index <- which(rownames(gs) %in% 
                           ci[grepl("PM", ci$ExtendedType),]$Address)
  m2 <- apply(gs, 2, function(x){
    min(x[addr.pm.index])/max(x[addr.mm.index])})
  rm <- cbind(rm, m2)
  colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}

#' ba.specificity2
#' 
#' Get the Specificity II Green BeadArray metric.
#' 
#' @param rs Red signal matrix (required).
#' @param gs Green signal matrix (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @param ct Column name/metric title (character, required).
#' @rdname ba
ba.specificity2 <- function(rs, gs, rm, cdf, mnum = 9, mtot = 17, 
  cnamei = "specificityII", ct = "Specificity 2"){
  message("(",mnum,"/",mtot,") Calculating Specificity II metric...")
  ci <- cdf[cdf$ExtendedType == ct,]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  which.addr.red <- which(rownames(rs) %in% ci$Address)
  which.addr.grn <- which(rownames(gs) %in% ci$Address)
  m0.1 <- apply(rs, 2, function(x){
    min(x[which.addr.red])})
  m0.2 <- apply(gs, 2, function(x){
    max(x[which.addr.grn])
  })
  m1 <- m0.1/m0.2
  rm <- cbind(rm, m1)
  colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}

#' ba.extension.red
#' 
#' Get the Extension Red BeadArray metric.
#' 
#' @param rs Red signal matrix (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @param ct Column name/metric title (character, required).
#' @rdname ba
ba.extension.red <- function(rs, rm, cdf, mnum = 10, mtot = 17,
  cnamei = "extension.red", ct = "Extension.*"){
  message("(",mnum,"/",mtot,") Calculating Extension Red metric...")
  ci <- cdf[grepl(ct, cdf$ExtendedType),]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  addr.ext.cg <- ci$Address[grepl("(C)|(G)", ci$ExtendedType)]
  addr.ext.at <- ci$Address[grepl("(A)|(T)", ci$ExtendedType)]
  which.cg <- which(rownames(rs) %in% addr.ext.cg)
  which.at <- which(rownames(rs) %in% addr.ext.at)
  m1 <- apply(rs, 2, function(x){
    min(x[which.at])/max(x[which.cg])})
  rm <- cbind(rm, m1)
  colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}

#' ba.extension.grn
#' 
#' Get the Extension Green BeadArray metric.
#' 
#' @param gs Green signal matrix (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @param ct Column name/metric title (character, required).
#' @rdname ba
ba.extension.grn <- function(gs, rm, cdf, mnum = 11, mtot = 17, 
  cnamei = "extension.grn", ct = "Extension.*"){
  message("(",mnum,"/",mtot,") Calculating Extension Green metric...")
  ci <- cdf[grepl(ct, cdf$ExtendedType),]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  addr.ext.cg <- ci$Address[grepl("(C)|(G)", ci$ExtendedType)]
  addr.ext.at <- ci$Address[grepl("(A)|(T)", ci$ExtendedType)]
  which.cg <- which(rownames(gs) %in% addr.ext.cg)
  which.at <- which(rownames(gs) %in% addr.ext.at)
  m2 <- apply(gs, 2, function(x){
    min(x[which.cg])/max(x[which.at])})
  rm <- cbind(rm, m2)
  colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}

#' ba.hybridization.hi.vs.med
#' 
#' Get the Hybridication (high vs. medium) BeadArray metric.
#' 
#' @param gs Green signal matrix (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @param ct Column name/metric title (character, required).
#' @rdname ba
ba.hybridization.hi.vs.med <- function(gs, rm, cdf, mnum = 12, mtot = 17,
  cnamei = "hyb.hi.med", ct = "High|Medium"){
  message("(",mnum,"/",mtot,") Calculating Hybridization High vs. Medium metric...")
  ci <- cdf[grepl(ct, cdf$ExtendedType),]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  which.hi <- which(rownames(gs) %in% 
                      ci$Address[grepl("(High)", ci$ExtendedType)])
  which.med <- which(rownames(gs) %in% 
                       ci$Address[grepl("(Medium)", ci$ExtendedType)])
  m1 <- apply(gs, 2, function(x){
    x[which.hi]/x[which.med]})
  rm <- cbind(rm, m1)
  colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}

#' ba.hybridization.med.vs.low
#' 
#' Get the Hybridication (medium vs. low) BeadArray metric.
#' 
#' @param gs Green signal matrix (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @param ct Column name/metric title (character, required).
#' @rdname ba
ba.hybridization.med.vs.low <- function(gs, rm, cdf, mnum = 13, mtot = 17, 
  cnamei = "hyb.med.low", ct = "Low|Medium"){
  message("(",mnum,"/",mtot,") Calculating Hybridization Medium vs. Low metric...")
  ci <- cdf[grepl(ct, cdf$ExtendedType),]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  which.med <- which(rownames(gs) %in% 
                       ci$Address[grepl("(Medium)", ci$ExtendedType)])
  which.low <- which(rownames(gs) %in% 
                       ci$Address[grepl("(Low)", ci$ExtendedType)])
  m2 <- apply(gs, 2, function(x){
    x[which.med]/x[which.low]})
  rm <- cbind(rm, m2)
  colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}

#' ba.targetremoval1
#' 
#' Get the Target Removal 1 BeadArray metric.
#' 
#' @param gs Green signal matrix (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param baseline Baseline measure for signals (integer, 3000).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @param ct Column name/metric title (character, required).
#' @rdname ba
ba.targetremoval1 <- function(gs, rm, cdf, baseline = 3000, mnum = 14, 
  mtot = 17, cnamei = "target.removal1", ct = "Extension|Target Removal 1"){
  message("(",mnum,"/",mtot,") Calculating Target Removal 1 metric...")
  ci <- cdf[grepl("Extension|Target Removal 1", cdf$ExtendedType),]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  which.bkg <- which(rownames(gs) %in% 
                       ci[grepl("\\(A\\)|\\(T\\)", 
                                ci$ExtendedType), ]$Address)
  which.t1 = which(rownames(gs) %in% 
                     ci[grepl("Removal 1", 
                              ci$ExtendedType),]$Address)
  m1 <- apply(gs, 2, function(x){
    (max(x[which.bkg]) + baseline)/x[which.t1]})
  rm <- cbind(rm, m1)
  colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}

#' ba.targetremoval2
#' 
#' Get the Target Removal 2 BeadArray metric.
#' 
#' @param gs Green signal matrix (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param baseline Baseline measure for signals (integer, 3000).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @param ct Column name/metric title (character, required).
#' @rdname ba
ba.targetremoval2 <- function(gs, rm, cdf, baseline = 3000, mnum = 15, mtot = 17, 
  cnamei = "target.removal2", ct = "Extension|Target Removal 2"){
  message("(",mnum,"/",mtot,") Calculating Target Removal 2 metric...")
  ci <- cdf[grepl(ct, cdf$ExtendedType),]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  which.bkg <- which(rownames(gs) %in% 
                       ci[grepl("(A)|(T)", ci$ExtendedType), ]$Address)
  which.t2 <- which(rownames(gs) %in% 
                      ci[grepl("Target Removal 2", ci$ExtendedType),]$Address)
  m2 <- apply(gs, 2, function(x){
    (max(x[which.bkg]) + baseline)/x[which.t2]})
  rm <- cbind(rm, m2)
  colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}

#' ba.bisulfiteconv1.grn
#' 
#' Get the Bisulfite Conversion I Green BeadArray metric.
#' 
#' @param gs Green signal matrix (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @param ct Column name/metric title (character, required).
#' @rdname ba
ba.bisulfiteconv1.grn <- function(gs, rm, cdf, mnum = 16, mtot = 17, 
  cnamei = "bisulfite.conv1.grn", ct = "Conversion I-"){
  message("(",mnum,"/",mtot,") Calculating Bisulfite Conversion I metric...")
  ci <- cdf[grepl(ct, cdf$ExtendedType),]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  which.c123 <- which(rownames(gs) %in% 
                        ci$Address[grepl("C1|C2|C3", ci$ExtendedType)])
  which.u123 <- which(rownames(gs) %in% 
                        ci$Address[grepl("U1|U2|U3", ci$ExtendedType)])
  m2 <- apply(gs, 2, function(x){
    min(x[which.c123])/max(x[which.u123])})
  rm <- cbind(rm, m2)
  colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}

#' ba.bisulfiteconv2
#' 
#' Get the Bisulfite Conversion II BeadArray metric.
#' 
#' @param rs Red signal matrix (required).
#' @param gs Green signal matrix (required).
#' @param rm Results matrix containing control metrics (required).
#' @param cdf Control probe data frame (required).
#' @param mnum Metric number out of total (numeric).
#' @param mtot Total metrics to be calculated (numeric, 17).
#' @param cnamei Column name (character, required).
#' @param ct Column name/metric title (character, required).
#' @rdname ba
ba.bisulfiteconv2 <- function(rs, gs, rm, cdf, mnum = 17, mtot = 17, 
  cnamei = "bisulfite.conv2", ct = "Conversion II-"){
  message("(",mnum,"/",mtot,") Calculating Bisulfite Conversion II metric...")
  ci <- cdf[grepl(ct, cdf$ExtendedType),]
  if(nrow(ci) == 0){
    message("Warning: couldn't find control probe addresses for metric ", 
            cnamei, ". Continuing..."); return(rm)}
  which.ci.red <- which(rownames(rs) %in% ci$Address)
  which.ci.grn <- which(rownames(gs) %in% ci$Address)
  m0.1 <- apply(rs, 2, function(x){
    min(x[which.ci.red])})
  m0.2 <- apply(gs, 2, function(x){
    max(x[which.ci.grn])})
  m1 <- m0.1/m0.2
  rm <- cbind(rm, m1)
  colnames(rm)[ncol(rm)] <- cnamei # update rm
  return(rm)
}