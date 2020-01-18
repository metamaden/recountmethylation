#!/usr/bin/env R

# Access recount-methylation HDF5 database.

#--------------------------
# HDF5 connection utilities
#--------------------------

#' Connect to an HDF5 dataset file.
#'
#' Get a dataset connection object from an HDF5 database ('.h5') file.
#' @param ri rows indices in dataset.
#' @param ci columns indices in dataset.
#' @param dsn Name of dataset or group of dataset to connect with.
#' @param dbn Path to h5 database file.
#' @return HDF5 database connection object.
#' @export
hread = function(ri, ci, dsn = "redsignal", dbn = "remethdb.h5"){
  return(rhdf5::h5read(dbn, dsn, index = list(ri, ci)))
}

#---------------------------------
# Query and retrieve HDF5 datasets
#---------------------------------

#' Retrieve samples metadata.
#'
#' Retrieves sample postprocessed metadata from an HDF5 database.
#'
#' @param dbn Path to HDF5 database file.
#' @param dsn Name or group path to HDF5 dataset containing postprocessed metadata.
#' @return Postprocessed metadata as a `data.frame`.
#' @export
data.mdpost = function(dbn = "remethdb.h5", dsn = "metadata"){
  mdp = as.data.frame(rhdf5::h5read(file = dbn, name = dsn), stringsAsFactors = F)
  colnames(mdp) = rhdf5::h5read(file = dbn, name = paste(dsn, "colnames", sep = "."))
  return(mdp)
}

#-------------------------------------------------------
# Get SummarizedExperiment objects from dataset queries
#-------------------------------------------------------

#' Form an object of class `RGChannelSet` from a raw signal dataset query
#'
#' Forms a `RGChannelSet` object from an HDF5 databse file query to the red and green raw signal tables. See `getrg()` function for implementation.
#'
#' @param ldat List of raw signal data query results. Must include 2 `data.frame` objects named 'redsignal' and 'greensignal'.
#' @param verbose Whether to post status messages.
#' @return Returns a `RGChannelSet` object from raw signal dataset queries.
#' @export
rgse = function(ldat, verbose = FALSE){
  if(!("greensignal" %in% names(ldat) & "redsignal" %in% names(ldat))){
    stop("Invalid ldat object passed. Verify the necessary data types are available.")
  }
  # match and check CpG addresses for signal datasets
  if(verbose){
    message("Matching CpG addresses for signal matrices...")
  }
  rga = ldat[["redsignal"]]; gga = ldat[["greensignal"]]

  # match probe IDs
  addrv = unique(c(rownames(rga), rownames(gga))) # get the unique CpG addresses
  rgs = rga[rownames(rga) %in% addrv, ]; rgs = rgs[order(match(rownames(rgs), addrv)), ]
  ggs = gga[rownames(gga) %in% addrv, ]; ggs = ggs[order(match(rownames(ggs), addrv)), ]
  if(!identical(rownames(rgs), rownames(ggs))){
    stop("Couldn't match CpG addresses for signal data.")
  }

  # match and check GMS IDs for signal datasets
  message("Matching GSM IDs for signal matrices...")
  gsmidv = unique(c(colnames(rgs), colnames(ggs)))
  rgf = rgs[, colnames(rgs) %in% gsmidv]
  ggf = ggs[, colnames(ggs) %in% gsmidv]
  rgf = rgf[, order(match(colnames(rgf), gsmidv))]
  ggf = ggf[, order(match(colnames(ggf), gsmidv))]
  if(!identical(colnames(rgf), colnames(ggf))){
    stop("Couldn't match GSM IDs for signal data.")
  }

  # match and check GSM IDs for metadata
  if("metadata" %in% names(ldat)){
    if(verbose){
      message("Checking provided postprocessed metadata...")
    }
    mdp = ldat[["metadata"]]
    mdp$gsm = as.character(mdp$gsm)
    # append blank rows for missing GSM IDs
    gsmov = gsmidv[!gsmidv %in% mdp$gsm]
    if(length(gsmov) > 0){
      if(verbose){
        message("Appending data for ", length(gsmov), " GSM IDs lacking metadata...")
      }
      numo = length(gsmov)
      nmm = matrix(c(gsmov, rep(rep("NA", numo), ncol(mdp) - 1)), nrow = numo)
      colnames(nmm) = colnames(mdp)
      mdp = rbind(mdp, nmm)
    }
    if(verbose){
      message("Checking metadata match...")
    }
    mdf = mdp[mdp$gsm %in% gsmidv,]
    mdf = mdf[order(match(mdf$gsm, gsmidv)),]
    checkid = identical(mdf$gsm, colnames(rgf)) & identical(mdf$gsm, colnames(ggf))
    if(!checkid){
      stop("Couldn't match metadata GSM IDs with signal ID matrices.")
    }
    rownames(mdf) = mdf$gsm
  }
  if(verbose){
    message("forming the RGset...")
  }
  anno = c("IlluminaHumanMethylation450k", "ilmn12.hg19")
  names(anno) = c("array", "annotation")
  rgi = minfi::RGChannelSetExtended(Green = ggf, Red = rgf,
                                    annotation = anno)
  if("metadata" %in% names(ldat)){
    if(verbose){
      message("Adding postprocessed metadata as pheno data to SE set...")
    }
    pData(rgi) = S4Vectors::DataTable(mdf)
  return(rgi)
  }
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
#' @param md.dsn Name of metadata dataset in h5 file.
#' @param verbose Whether to post status messages.
#' @return Returns either an `RGChannelSet` or list of `data.frame` objects from dataset query matches.
#' @export
getrg = function(gsmv = "random", cgv = "all",
                 dbn = "remethdb.h5", data.type = c("se", "df"),
                 dsv = c("redsignal", "greensignal"), metadata = TRUE,
                 md.dsn = "metadata", verbose = FALSE){
  # form the datasets list
  if(length(gsmv) == 0 | length(cgv) == 0){
    stop("Invalid GSM or CpG IDs. Check arguments for 'gsmv' and 'cgv'.")
  }
  ldat = list()
  for(d in dsv){
    if(verbose){
      message("Working on ", d, "...")
    }
    if(d %in% c("redsignal", "greensignal")){
      rnd = rhdf5::h5read(dbn, paste(d, "rownames", sep = ".")) # rownames, GSM IDs
      rnd = gsub("\\..*", "", rnd) # clean GSM IDs
      cnd = rhdf5::h5read(dbn, paste(d, "colnames", sep = ".")) # colnames, CpG addr
      # parse the index values
      if(cgv == "all"){
        cgvp = seq(1, length(cnd), 1)
      } else{
        cgvp = which(cnd %in% cgv)
      }
      if(gsmv == "random"){
        gsmvp = sample(length(rnd), 5)
      } else{
        gsmvp = which(rnd %in% gsmv)
      }
      # get data matrix
      ddat = hread(ri = gsmvp, ci = cgvp, d, dbn)
      rownames(ddat) = rnd[gsmvp]
      colnames(ddat) = cnd[cgvp]
      ldat[[d]] <- t(ddat) # append transpose of data
    } else{
      if(verbose){
        message("Invalid data type detected: ", d, ", continuing...")
      }
    }
  }
  # append metadata
  if(metadata){
    mdpost = data.mdpost(dbn = dbn, dsn = md.dsn)
    mdpost$gsm = as.character(mdpost$gsm)
    mdf = mdpost[mdpost$gsm %in% gsmv,]
    ldat[["metadata"]] = mdf
  }
  # return desired data type
  if(data.type == "df"){
    if(verbose){
      message("Returning the datasets list...")
    }
    robj = ldat
  }
  if(data.type == "se"){
    if(verbose){
      message("Forming the RGChannelSet...")
    }
    robj = rgse(ldat = ldat, verbose = verbose)
  }
  rhdf5::h5closeAll() # close all open connections
  return(robj)
}
