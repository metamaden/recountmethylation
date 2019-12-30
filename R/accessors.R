#!/usr/bin/env R

# Access recount-methylation HDF5 database.

# HDF5 connection objects
#

#' Connection to HDF5 database file
#'
#' Get '.h5' file connection (db connection).
#'
#' @param dbn.path Path to HDF5 database file.
#' @return An HDF5 file connection object.
#' @export
hfcon = function(dbn.path = "remethdb.h5"){
  hfile = rhdf5::H5Fopen(dbn.path)
  return(hfile)
}

#' Connect to an HDF5 dataset file.
#'
#' Get a dataset connection object from an HDF5 database ('.h5') file.
#'
#' @param dsn Name of dataset or group of dataset to connect with.
#' @param hfile Connection to an HDF5 database file.
#' @return HDF5 database connection object.
#' @export
hdcon = function(dsn = "noobbeta", hfile = hfcon()){
  hdopen = rhdf5::H5Dopen(hfile, dsn)
  return(hdopen)
}

# Access and query the metadata

#---------------------------------
# Query and retrieve HDF5 datasets
#---------------------------------

#' Retrieve samples metadata.
#'
#' Retrieves sample postprocessed metadata from an HDF5 database.
#'
#' @param dbn.path Path to HDF5 database file.
#' @param dsn Name or group path to HDF5 dataset containing postprocessed metadata.
#' @return Postprocessed metadata as a `data.frame`.
#' @export
data.mdpost = function(dbn.path = "remethdb.h5", dsn = "mdpost"){
  return(rhdf5::h5read(file = dbn.path, name = dsn))
}

# get noobbeta data

#' Get noob-normalized Beta-values.
#'
#' Retrieves noob-normalized Beta-values dataset from HDF5 database file.
#'
#' @param dbn Name of HDF5 database ('.h5') file.
#' @param dsn Name of HDF5 dataset containing noob-normalized Beta-values.
#' @param which.gsm Numeric indices of samples whose data to retrieve.
#' @param which.cg Numeric indices of CpG probes whose data to retrieve.
#' @return Returns subset of noob-normalized Beta-values HDF5 dataset matching query.
#' @export
data.nb = function(dbn = "remethdb.h5", dsn = 'noobbeta',
                  which.gsm = sample(35700, 5), which.cg = c(1:485512)){
  message("Connecting to hdf5 entities...")
  hfile = hfcon(dbn.path = dbn)
  dsr = hdcon(dsn = dsn, hfile = hfile)
  # filter rownames as necessary
  message("subsetting and coercing data...")
  # validate provided indices for rows and cols
  message("Validating provided row and col indices...")
  rn.nb = rhdf5::h5read(dbn, paste0(dsn, ".rownames")) # rows/cpgs
  cn.nb = rhdf5::h5read(dbn, paste0(dsn, ".colnames")) # cols/samples
  which.gsm = which.gsm[which.gsm %in% seq(1, length(cn.nb), 1)]
  which.cg = which.cg[which.cg %in% seq(1, length(rn.nb), 1)]
  message("Subsetting the signal matrix...")
  sd = dsr[which.cg, which.gsm]
  rownames(sd) = as.character(rn.nb[which.cg])
  colnames(sd) = as.character(cn.nb[which.gsm])
  return(sd)
}

#' Get raw probe signal plots
#'
#' Retrieves raw signal from either red or green array color channel HDF5 dataset.
#'
#' @param dbn Name of HDF5 database ('.h5') file.
#' @param dsn Name or group path to a raw signal dataset (either 'redsignal' or 'greensignal').
#' @param which.gsm Numeric indices of samples (signal dataset rows) whose data to retrieve.
#' @param which.cg Numeric indices of CpG probes (signal dataset columns) whose data to retrieve.
#' @return Returns the subset of the raw signal dataset matching the query.
#' @export
data.rgsignal = function(which.gsm, which.cg, 
                         dbn = "remethdb.h5", dsn = c("redsignal", "greensignal"),
                         rnrg = rhdf5::h5read(dbn, paste0(dsn, ".rownames")),
                         cnrg = rhdf5::h5read(dbn, paste0(dsn, ".colnames"))){
  message("Connecting to hdf5 entities...")
  hfile = hfcon(dbn.path = dbn)
  dsr = hdcon(dsn = dsn, hfile = hfile)
  # validate provided indices for rows and cols
  message("Validating provided row and col indices...")
  message("Subsetting the signal matrix...")
  # sdv = dsr[eval(which.gsm), eval(which.cg)]
  wgsmv = eval(which.gsm); wcgv = eval(which.cgv)
  sdv = dsr[wgsmv, wcgv]
  rnv = eval(as.character(rnrg[eval(which.gsm)]))
  message(length(rnv))
  message(dim(sdv))
  rownames(sdv) = rnv
  # colnames(sdv) = as.character(cnrg[which.cg])
  return(sdv)
}

#-------------------------------------------------------
# Get SummarizedExperiment objects from dataset queries
#-------------------------------------------------------

#' Form an object of class `RGChannelSet` from a raw signal dataset query
#'
#' Forms a `RGChannelSet` object from an HDF5 databse file query to the red and green raw signal tables. See `getrg()` function for implementation.
#'
#' @param ldat List of raw signal data query results. Must include 2 `data.frame` objects named 'redsignal' and 'greensignal'.
#' @return Returns a `RGChannelSet` object from raw signal dataset queries.
#' @export
rgse = function(ldat){
  if(!("greensignal" %in% names(ldat) & "redsignal" %in% names(ldat))){
    message("Error: invalid ldat object passed. Verify the necessary data types are available.")
    return()
  }
  # match and check CpG addresses for signal datasets
  message("Matching CpG addresses for signal matrices...")
  rga = ldat[["redsignal"]]; gga = ldat[["greensignal"]]
  addrv = unique(c(colnames(rga), colnames(gga))) # get the unique CpG addresses
  rgs = rga[, colnames(rga) %in% addrv]; rgs = rgs[, order(match(colnames(rgs), addrv))]
  ggs = gga[, colnames(gga) %in% addrv]; ggs = ggs[, order(match(colnames(ggs), addrv))]
  if(!identical(colnames(rgs), colnames(ggs))){
    message("Error: couldn't match CpG addresses for signal data! Returning...")
    return()
  }
  # match and check GMS IDs for signal datasets
  message("Matching GSM IDs for signal matrices...")
  gsmidv = unique(c(rownames(rgs), rownames(ggs)))
  rgf = rgs[rownames(rgs) %in% gsmidv,]
  rgf = rgf[order(match(rownames(rgs), gsmidv)),]
  ggf = ggs[rownames(ggs) %in% gsmidv,]
  ggf = ggf[order(match(rownames(ggf), gsmidv)),]
  if(!identical(rownames(rgf), rownames(ggf))){
    message("Error: couldn't match GSM IDs for signal data! Returning...")
    return()
  }
  # match and check GSM IDs for metadata
  if("mdpost" %in% names(ldat)){
    message("Checking provided postprocessed metadata...")
    mdp = ldat[["mdpost"]]
    mdp$gsm = as.character(mdp$gsm)
    # append blank rows for missing GSM IDs
    gsmov = gsmidv[!gsmidv %in% mdp$gsm]
    if(length(gsmov) > 0){
      message("Appending data for ", length(gsmov), " GSM IDs lacking metadata...")
      numo = length(gsmov)
      nmm = matrix(c(gsmov, rep(rep("NA", numo), ncol(mdp) - 1)), nrow = numo)
      colnames(nmm) = colnames(mdp)
      mdp = rbind(mdp, nmm)
    }
    message("Checking metadata match...")
    mdf = mdp[mdp$gsm %in% gsmidv,]
    mdf = mdf[order(match(mdf$gsm, gsmidv)),]
    if(!(identical(mdf$gsm, rownames(rgf)) |
         identical(mdf$gsm, rownames(ggf)))){
      message("Error: couldn't match metadata GSM IDs with signal ID matrices...")
      return()
    }
    rownames(mdf) = mdf$gsm
  }
  message("forming the RGset...")
  anno = c("IlluminaHumanMethylation450k", "ilmn12.hg19")
  names(anno) = c("array", "annotation")
  rgi = minfi::RGChannelSet(Green = t(ggf), Red = t(ggf),
                     annotation = anno)
  if("mdpost" %in% names(ldat)){
    message("Adding postprocessed metadata as pheno data to SE set...")
    minfi::pData(rgi) = S4Vectors::DataFrame(mdf)
  }
  return(rgi)
}

#' Form an object of class `GenomicRatioSet` from a raw signal dataset query
#'
#' Forms a `GenomicRatioSet` object from an HDF5 databse file query to the red and green raw signal tables. See `getgr()` function for implementation.
#'
#' @param ldat List of raw signal data query results. Must include 2 `data.frame` objects named 'redsignal' and 'greensignal'.
#' @return Returns a `GenomicRatioSet` object from raw signal dataset queries.
#' @export
grse = function(ldat, granges.obj = data("granges_minfi")){
  if(!"noobbeta" %in% names(ldat)){
    message("Error: invalid ldat object passed. Verify the necessary data type(s) are available.")
    return()
  }
  nb = ldat[["noobbeta"]]
  # check provided metadata
  if("mdpost" %in% names(ldat)){
    gsmidv = colnames(nb)
    message("Checking provided postprocessed metadata...")
    mdp = ldat[["mdpost"]]; mdp$gsm = as.character(mdp$gsm)
    # append blank rows for missing GSM IDs
    gsmov = gsmidv[!gsmidv %in% mdp$gsm]
    if(length(gsmov) > 0){
      message("Appending rows for ", length(gsmov), " GSM IDs lacking metadata...")
      numo = length(gsmov)
      nmm = matrix(c(gsmov, rep(rep("NA", numo), ncol(mdp) - 1)), nrow = numo)
      colnames(nmm) = colnames(mdp)
      mdp = rbind(mdp, nmm)
    }
    message("Checking metadata match...")
    mdf = mdp[mdp$gsm %in% gsmidv,]
    mdf = mdf[order(match(mdf$gsm, gsmidv)),]
    if(!identical(as.character(mdf$gsm),
                  as.character(colnames(nb)))){
      message("Error: couldn't match metadata GSM IDs with noob Beta-value ID matrices...")
      return()
    }
    rownames(mdf) = mdf$gsm
  }
  # form the SE object
  message("Forming the GenomicRatioSet...")
  # data("granges_minfi")
  anno = c("IlluminaHumanMethylation450k", "ilmn12.hg19")
  names(anno) = c("array", "annotation")
  int.cg = intersect(names(granges.obj), rownames(nb))
  grmf = granges.obj[names(granges.obj) %in% int.cg]
  nb = nb[rownames(nb) %in% int.cg,]; nb = nb[order(match(rownames(nb), names(grmf))),]
  if(!identical(rownames(nb), names(grmf))){
    message("Error: couldn't match CpG names from noobbeta with GRanges names...")
    return()
  }
  gri = minfi::GenomicRatioSet(gr = grmf, Beta = nb, annotation = anno)

  if("mdpost" %in% names(ldat)){
    message("Adding postprocessed metadata as pheno data to SE set...")
    minfi::pData(gri) = S4Vectors::DataFrame(mdf)
  }
  return(gri)
}

#' Get raw signal data as either a list of `data.frame`'s or an `RGChannelSet` object
#'
#' Retrieves query matches from raw signal HDF5 datasets. Handles identity queries to rows (GSM IDs) or columns (CpG probe addresses). Returns query matches either as a list of 2 `data.frame`s or a signle `RGChannelSet` object.
#'
#' @param dbn Name of the HDF5 database file.
#' @param gsmv Vector valid GSM IDs (rows) to query in the raw signal datasets.
#' @param cgv Vector of valid CpG probe addresses (columns) to query in the raw signal datasets.
#' @param data.type Format for returned query matches, either as datasets 'df' or `RGChannelSet` 'se' object.
#' @param dsv Vector of raw signal datasets or group paths to query, including both the red channel 'redsignal' and green channel 'greensignal' datasets.
#' @param metadata Whether to access available postprocessed metadata for queries samples.
#' @return Returns either an `RGChannelSet` or list of `data.frame` objects from dataset query matches.
#' @export
getrg = function(gsmv, cgv, dbn = "remethdb.h5", dat.type = c("se", "df"),
                 dsv = c("redsignal", "greensignal"), metadata = TRUE){
  # form the datasets list
  if(length(gsmv) == 0 | length(cgv) == 0){
    message("Error: invalid GSM or CpG IDs. Check arguments for 'gsmv' and 'cgv'...")
    return()
  }
  ldat = list()
  
  for(d in dsv){
    message("Working on ", d, "...")
    if(d %in% c("redsignal", "greensignal")){
      rnd = rhdf5::h5read(dbn, paste0(d, ".rownames")) # rownames, GSM IDs
      cnd = rhdf5::h5read(dbn, paste0(d, ".colnames")) # colnames, CpG addr
      ddat = hdcon(d, hfcon(dbn))[which(rnd %in% gsmv), which(cnd %in% cgv)]
      rownames(ddat) = rnd[which(rnd %in% gsmv)]
      colnames(ddat) = cnd[which(cnd %in% cgv)]
      ldat[[d]] <- ddat
    } else{
      message("Warning: Invalid data type detected: ", d, ", continuing...")
    }
  }
  # append metadata
  if(metadata){
    mdpost = data.mdpost(); mdpost$gsm = as.character(mdpost$gsm)
    mdf = mdpost[mdpost$gsm %in% gsmv,]
    ldat[["mdpost"]] = mdf
  }
  # return desired data type
  if(dat.type == "df"){
    message("Returning the datasets list...")
    robj = ldat
  }
  if(dat.type == "se"){
    message("Forming the RGChannelSet...")
    robj = rgse(ldat = ldat)
  }
  return(robj)
}

#' Get noob-normalized Beta-value data as either a list of `data.frame` or an `GenomicRatioSet` object
#'
#' Retrieves query matches from raw signal HDF5 datasets. Handles identity queries to rows (GSM IDs) or columns (CpG probe addresses). Returns query matches either as a list of 2 `data.frame`s or a signle `RGChannelSet` object.
#'
#' @param dbn Name of the HDF5 database file.
#' @param gsmv Vector valid GSM IDs (rows) to query in the raw signal datasets.
#' @param cgv Vector of valid CpG probe addresses (columns) to query in the raw signal datasets.
#' @param data.type Format for returned query matches, either as datasets 'df' or `GenomicRatioSet` 'se' object.
#' @param dsv Vector of raw signal datasets or group paths to query, including both the red channel 'redsignal' and green channel 'greensignal' datasets.
#' @param metadata Whether to access available postprocessed metadata for queries samples.
#' @return Returns either an `GenomicRatioSet` or list of `data.frame` objects from dataset query matches.
#' @export
getgs = function(gsmv, cgv, dbn = "remethdb.h5", 
                 dat.type = c("se", "df"),
                 dsv = c("noobbeta"), metadata = TRUE, 
                 granges.obj = data("granges_minfi")){
  # form the datasets list
  if(length(gsmv) == 0 | length(cgv) == 0){
    message("Error: invalid GSM or CpG IDs. Check arguments for 'gsmv' and 'cgv'...")
    return()
  }
  ldat = list()
  for(d in dsv){
    message("working on ", d, "...")
    if(d %in% c("noobbeta")){
      rnd = rhdf5::h5read(dbn, paste0(d, ".rownames")) # rownames, CpG ID
      cnd = rhdf5::h5read(dbn, paste0(d, ".colnames")) # colnames, GSM ID
      ddat = hdcon(d, hfcon(dbn))[which(rnd %in% cgv), which(cnd %in% gsmv)]
      colnames(ddat) = cnd[which(cnd %in% gsmv)]
      message("cn")
      rownames(ddat) = rnd[which(rnd %in% cgv)]
      message("rn")
      ldat[[d]] <- ddat
    } else{
      message("Invalid data type detected: ", d, ", continuing...")
    }
  }
  # append metadata
  if(metadata){
    mdpost = data.mdpost(); mdpost$gsm = as.character(mdpost$gsm)
    mdf = mdpost[mdpost$gsm %in% gsmv,]
    ldat[["mdpost"]] = mdf
  }
  # return desired data type
  if(dat.type == "df"){
    message("Returning the datasets list...")
    robj = ldat
  }
  if(dat.type == "se"){
    message("Forming the RGChannelSet...")
    robj = grse(ldat = ldat, granges.obj = granges.obj)
  }
  return(robj)
}
