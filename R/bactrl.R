#!/usr/bin/env R

# Functions to generate and assess Illumina BeadArray probe metrics

# Control probe data

#' Get control probe signals
#'
#' Returns signal tables for control probes in an RGChannelSet object.
#' @param rgset An RGChannelSet object.
#' @return list of signal tables.
#' @examples 
#' library(minfiData)
#' rg <- get(data(RGsetEx))
#' lctrl <- ctrl_signal(rg)
#' names(lctrl)
#' dim(lctrl[[1]])
ctrl_signal <- function(rgset){
  require(IlluminaHumanMethylation450kmanifest)
  if(!class(rgset) == "RGChannelSet"){
      stop("Object must be of class RGChannelSet")
    }
  data(IlluminaHumanMethylation450kmanifest)
  ctrladdr <- minfi::getControlAddress(IlluminaHumanMethylation450kmanifest)
  ca.int <- intersect(rownames(rgset), ctrladdr)
  return(list("red" = getRed(rgset[ca.int,]), 
              "green" = getGreen(rgset[ca.int,])))
}


# Control metrics

c("restoration.grn","biotin.stain.red", "biotin.stain.grn",
  "specificityI.red", "specificityI.grn", "specificityII",
  "extension.red", "extension.grn", "hyb.hi.med", "hyb.med.low",
  "target.removal.1", "target.removal.2", "bisulfite.conv.I.red",
  "bisulfite.conv.I.grn", "bisulfite.conv.II", "nonpolymorphic.red",
  "nonpolymorphic.grn")


# restoration.grn
#'
#' @param
#' @returns
rest.grn <- function(cdf){
  # Background addr, # 1 metric using extension grn channel
  ct = "EXTENSION" 
  ci = cdf[cdf$Type==ct,]
  addr.bkg = ci[grepl("(A)|(T)", ci$ExtendedType),]$Address
  return()
}

#' background
#' 
#' @param
#' @returns Background address
get_bgaddr <- function(ct = "EXTENSION"){
  # Background addr, # 1 metric using extension grn channel
  ci = cdf[cdf$Type==ct,]
  addr.bkg = ci[grepl("(A)|(T)", ci$ExtendedType),]$Address
  return(addr.bkg)
}

# biotin.stain
#'
#' @param
#' @returns

# specificityI
#'
#' @param
#' @returns

# specificityII
#'
#' @param
#' @returns

# extension
#'
#' @param
#' @returns

# hyb.hi.med
#'
#' @param
#' @returns

# hyb.med.low
#'
#' @param
#' @returns

# target.removal.1
#'
#' @param
#' @returns

# target.removal.2
#'
#' @param
#' @returns 

# bisulfite.conv.I
#'
#' @param
#' @returns

# bisulfite.conv.II
#'
#' @param
#' @returns

# nonpolymorphic
#'
#' @param
#' @returns



#' Get Illumina control metrics from red and green signal
#'
#' Uses bead signals from control probes to calculate Illumina's control metrics.
#' 
#' @param baseline Baseline offset for Extension green signal based background 
#' (see source 2)
#' @param biotin.baseline Special baseline offset for biotin staining background 
#' probes.
#' @param rs Red signal data frame (cols = probes, rows = samples, colnames = 
#' addresses).
#' @param gs Grn signal data frame (cols = probes, rows = samples, colnames = 
#' addresses).
#' @param cdf Control probe annotations data frame (cols = properties, 
#' rows = probes).
#' @param verbose Whether to return verbose messages.
#' @details 
#' This resource takes guidance from 3 key sources: the Illumina 
#' GenomeStudio Methylation 
#' Module (v1.8, source 1), the BeadArray Controls Reporter Software Guide 
#' (v00, source 2), 
#' and the ewastools resource (v1.5, source 3). Notes on how function 
#' relates to these sources follows:
#' 
#' 1. Use C and U 1-3 and 4-6 for Bisulfite Conversion I per source 2, 
#'   source 3 only uses 1-2/4-5, and may have string matching error.
#' 
#' 2. Use probe address "34648333" and "43603326" (DNP, Biotin subtype) for Biotin Staining 
#'   Background where appropriate, per sources 2 and 3.
#' 
#' 3. To offset many Inf/-Inf values due to 0 background signal, set a special Biotin
#'   specific baseline offset (default: 1) in denom of Biotin red/grn Bkg value. 
#'   No offset specified in source 2.
#' 
#' 4. For Specificity I, use PM, MM 1-3 for grn, 4-6 for red (per source 3, source 2 ambiguous)
#' 
#' 5. For Specificity II, use just probes S1-3, as probe S4 unavailable in control probe 
#'   annotation (per source 2, corroborated in source 3).
#' 
#' 6. Note definition for background is ambiguous (can be either Red or Grn channel, from extension),
#'   per source 2; source 3 uses extension Grn A/T probes for system background.
#' 
#' 7. Note additional source 2 metrics using background (e.g. specificity background/U1-3) not 
#'   calculated here (yet).
#'   
#' @return results matrix of control metrics for samples
bactrl = function(baseline = 3000, biotin.baseline = 1, rs = rs, gs = gs, cdf = cdf){
  rm = data.frame(gsm = rownames(rs)) # return matrix
  cnl = c("restoration.grn","biotin.stain.red", "biotin.stain.grn",
          "specificityI.red", "specificityI.grn", "specificityII",
          "extension.red", "extension.grn", "hyb.hi.med", "hyb.med.low",
          "target.removal.1", "target.removal.2", "bisulfite.conv.I.red",
          "bisulfite.conv.I.grn", "bisulfite.conv.II", "nonpolymorphic.red",
          "nonpolymorphic.grn")
  # Background addr, # 1 metric using extension grn channel
  {
    ct = "EXTENSION" 
    ci = cdf[cdf$Type==ct,]
    addr.bkg = ci[grepl("(A)|(T)", ci$ExtendedType),]$Address
  }
  
  # Restoration, 1 metric, uses just grn channel
  {
    message("calculating restore metric...")
    ct = "RESTORATION"; ci = cdf[cdf$Type == ct,]
    which.rest = which(colnames(gs) %in% ci$Address)
    m1 = apply(gs, 1, function(x){x[which.rest]/(max(x[addr.bkg]) + baseline)})
    mi = 1
    rm = cbind(rm, m1)
    colnames(rm)[ncol(rm)] = cnl[mi]
  }
  
  # BIOTIN STAINING, 2 metrics, 1 per chan
  {
    message("calculating biotin staining metrics...")
    ci = cdf[grepl("Biotin|DNP", cdf$ExtendedType),]
    # red 
    which.stain = which(colnames(rs) %in% ci[grepl("DNP \\(High", ci$ExtendedType),]$Address)
    which.bkg = which(colnames(rs) %in% ci[grepl("DNP \\(Bkg", ci$ExtendedType),]$Address)
    m1 = apply(rs, 1, function(x){x[which.stain]/(x[which.bkg] + biotin.baseline)})
    mi = 2
    rm = cbind(rm, m1)
    colnames(rm)[ncol(rm)] = cnl[mi]
    # grn
    which.stain = which(colnames(gs) %in% ci[grepl("Biotin \\(High", ci$ExtendedType),]$Address)
    which.bkg = which(colnames(gs) %in% ci[grepl("Biotin \\(Bkg", ci$ExtendedType),]$Address)
    m2 = apply(gs, 1, function(x){x[which.stain]/(x[which.bkg] + biotin.baseline)})
    mi = 3
    rm = cbind(rm, m2)
    colnames(rm)[ncol(rm)] = cnl[mi]
  }
  
  # SPECIFICITY I, 2 metrics, 1 per channel
  {
    message("calculating specificity I...")
    ct = "SPECIFICITY I"
    ci1 = cdf[cdf$Type == ct,]
    # red
    ci = ci1[grepl("Mismatch (4|5|6)", ci1$ExtendedType),]
    addr.mm.index = which(colnames(rs) %in% ci[grepl("MM", ci$ExtendedType),]$Address)
    addr.pm.index = which(colnames(rs) %in% ci[grepl("PM", ci$ExtendedType),]$Address)
    m1 = apply(rs, 1, function(x){
      min(x[addr.pm.index])/max(x[addr.mm.index])
    })
    mi = 4
    rm = cbind(rm, m1)
    colnames(rm)[ncol(rm)] = cnl[mi]
    # grn
    ci = ci1[grepl("Mismatch (1|2|3)", ci1$ExtendedType),]
    addr.mm.index = which(colnames(gs) %in% ci[grepl("MM", ci$ExtendedType),]$Address)
    addr.pm.index = which(colnames(gs) %in% ci[grepl("PM", ci$ExtendedType),]$Address)
    m2 = apply(gs, 1, function(x){
      min(x[addr.pm.index])/max(x[addr.mm.index])
    })
    mi = 5
    rm = cbind(rm, m2)
    colnames(rm)[ncol(rm)] = cnl[mi]
  }
  
  # SPECIFICITY II, 1 metric, uses both channels
  {
    message("calculating specificityII metric...")
    ct = "SPECIFICITY II"
    ci = cdf[cdf$Type == ct,]
    which.addr.red = which(colnames(rs) %in% ci$Address)
    which.addr.grn = which(colnames(gs) %in% ci$Address)
    m0.1 = apply(rs, 1, function(x){
      min(x[which.addr.red])
    })
    m0.2 = apply(gs, 1, function(x){
      max(x[which.addr.grn])
    })
    m1 = m0.1/m0.2
    mi = 6
    rm = cbind(rm, m1)
    colnames(rm)[ncol(rm)] = cnl[mi]
  }
  
  # EXTENSION, 2 metrics, 1 per channel
  {
    message("calculating extension metrics...")
    ct = "EXTENSION"
    ci = cdf[cdf$Type == ct,]
    addr.ext.cg = ci$Address[grepl("(C)|(G)", ci$ExtendedType)]
    addr.ext.at = ci$Address[grepl("(A)|(T)", ci$ExtendedType)]
    # red
    which.cg = which(colnames(rs) %in% addr.ext.cg)
    which.at = which(colnames(rs) %in% addr.ext.at)
    m1 = apply(rs, 1, function(x){
      min(x[which.at])/max(x[which.cg])
    })
    # grn
    which.cg = which(colnames(gs) %in% addr.ext.cg)
    which.at = which(colnames(gs) %in% addr.ext.at)
    m2 = apply(gs, 1, function(x){
      min(x[which.cg])/max(x[which.at])
    })
    mi = 7
    rm = cbind(rm, m1)
    colnames(rm)[ncol(rm)] = cnl[mi]
    mi = 8
    rm = cbind(rm, m2)
    colnames(rm)[ncol(rm)] = cnl[mi]
  }
  
  # HYBRIDIZATION, 2 metrics, 1 per channel
  {
    message("calculating hybridization metrics...")
    ct = "HYBRIDIZATION"
    ci = cdf[cdf$Type == ct,]
    which.hi = which(colnames(gs) %in% ci$Address[grepl("(High)", ci$ExtendedType)])
    which.med = which(colnames(gs) %in% ci$Address[grepl("(Medium)", ci$ExtendedType)])
    which.low = which(colnames(gs) %in% ci$Address[grepl("(Low)", ci$ExtendedType)])
    m1 = apply(gs, 1, function(x){
      x[which.hi]/x[which.med]
    })
    m2 = apply(gs, 1, function(x){
      x[which.med]/x[which.low]
    })
    mi = 9
    rm = cbind(rm, m1)
    colnames(rm)[ncol(rm)] = cnl[mi]
    mi = 10
    rm = cbind(rm, m2)
    colnames(rm)[ncol(rm)] = cnl[mi]
  }
  
  # TARGET REMOVAL, 2 metrics, both from grn channel
  {
    message("calculating target removal metics...")
    ct = "EXTENSION"
    ci = cdf[cdf$Type == ct,]
    which.bkg = which(colnames(gs) %in% ci[grepl("(A)|(T)", ci$ExtendedType), ]$Address)
    ct = "TARGET REMOVAL"
    ci = cdf[cdf$Type == ct,]
    which.t1 = which(colnames(gs) %in% ci[grepl("Removal 1", ci$ExtendedType),]$Address)
    which.t2 = which(colnames(gs) %in% ci[grepl("Removal 2", ci$ExtendedType),]$Address)
    # rem 1
    m1 = apply(gs, 1, function(x){
      (max(x[which.bkg]) + baseline)/x[which.t1]
    })
    # rem 2
    m2 = apply(gs, 1, function(x){
      (max(x[which.bkg]) + baseline)/x[which.t2]
    })
    mi = 11
    rm = cbind(rm, m1)
    colnames(rm)[ncol(rm)] = cnl[mi]
    mi = 12
    rm = cbind(rm, m2)
    colnames(rm)[ncol(rm)] = cnl[mi]
  }
  
  # BISULFITE CONVERSION I, 2 metrics, 1 per channel
  {
    message("calculating bisulfite conversion I metrics...")
    ct = "BISULFITE CONVERSION I"
    ci = cdf[cdf$Type == ct,]
    which.c123 = which(colnames(gs) %in% ci$Address[grepl("C1|C2|C3", ci$ExtendedType)])
    which.u123 = which(colnames(gs) %in% ci$Address[grepl("U1|U2|U3", ci$ExtendedType)])
    which.c456 = which(colnames(rs) %in% ci$Address[grepl("C4|C5|C6", ci$ExtendedType)])
    which.u456 = which(colnames(rs) %in% ci$Address[grepl("U4|U5|U6", ci$ExtendedType)])
    # red
    m1 = apply(rs, 1, function(x){
      min(x[which.c456])/max(x[which.u456])
    })
    # grn
    m2 = apply(gs, 1, function(x){
      min(x[which.c123])/max(x[which.u123])
    })
    mi = 13
    rm = cbind(rm, m1)
    colnames(rm)[ncol(rm)] = cnl[mi]
    mi = 14
    rm = cbind(rm, m2)
    colnames(rm)[ncol(rm)] = cnl[mi]
  }
  
  # BISULFITE CONVERSION II, 2 metrics, 1 per channel
  {
    message("calculating bisulfite conversion II metric...")
    ct = "BISULFITE CONVERSION II"
    ci = cdf[cdf$Type == ct,]
    which.ci.red = which(colnames(rs) %in% ci$Address)
    which.ci.grn = which(colnames(gs) %in% ci$Address)
    m0.1 = apply(rs, 1, function(x){
      min(x[which.ci.red])
    })
    m0.2 = apply(gs, 1, function(x){
      max(x[which.ci.grn])
    })
    m1 = m0.1/m0.2
    mi = 15
    rm = cbind(rm, m1)
    colnames(rm)[ncol(rm)] = cnl[mi]
  }
  
  # NON-POLYMORPHIC, 2 metrics, 1 per channel
  {
    message("calculating non-polymorphic metrics...")
    ct = "NON-POLYMORPHIC"
    ci = cdf[cdf$Type == ct,]
    # red
    which.cg = which(colnames(rs) %in% ci$Address[grepl("(C)|(G)", ci$ExtendedType)])
    which.at = which(colnames(rs) %in% ci$Address[grepl("(A)|(T)", ci$ExtendedType)])
    m1 = apply(rs, 1, function(x){
      min(x[which.at])/max(x[which.cg])
    })
    # grn
    which.cg = which(colnames(gs) %in% ci$Address[grepl("(C)|(G)", ci$ExtendedType)])
    which.at = which(colnames(gs) %in% ci$Address[grepl("(A)|(T)", ci$ExtendedType)])
    m2 = apply(gs, 1, function(x){
      min(x[which.cg])/max(x[which.at])
    })
    mi = 16
    rm = cbind(rm, m1)
    colnames(rm)[ncol(rm)] = cnl[mi]
    mi = 17
    rm = cbind(rm, m2)
    colnames(rm)[ncol(rm)] = cnl[mi]
  }
  
  return(rm)
}

# get control quality assessments from an object returned by `bactrl`
bathresh = function(rm = rm){
  # bathresh
  # Apply BeadArray thresholds for quality flagging. Thresholds obtained from 
  #   the BeadArray Controls Reporter Software Guide (v00, source 2) and ewastools
  #   resource (v1.5, source 3).
  #   
  # Arguments:
  #   rm (matrix): Matrix of control probe metrics returned by `bactrl()` 
  #     (rows = samples, cols = metrics).
  # Returns:
  #   rm (matrix): Matrix of threshold assessments, either 'FAIL' or 'PASS', of same 
  #     dimensions as input matrix.
  dft = data.frame(restoration.grn = 0, 
                            biotin.stain.red = 5, 
                            biotin.stain.grn = 5,
                            specificityI.red = 1, 
                            specificityI.grn = 1, 
                            specificityII = 1,
                            extension.red = 5, 
                            extension.grn = 5, 
                            hyb.hi.med = 1, 
                            hyb.med.low = 1,
                            target.removal.1 = 1, 
                            target.removal.2 = 1, 
                            bisulfite.conv.I.red = 1,
                            bisulfite.conv.I.grn = 1, 
                            bisulfite.conv.II = 1, 
                            nonpolymorphic.red = 5, 
                            nonpolymorphic.grn = 5, 
                            stringsAsFactors = F)
  for(c in colnames(rm)[2:ncol(rm)]){
    rm[,c] = ifelse(rm[,c] < dft[, c], "FAIL", "PASS")
  }
  return(rm)
}