#!/usr/bin/env R

# Access recount-methylation HDF5 database.

#--------------------------
# HDF5 connection utilities
#--------------------------

#' Query and store an HDF5 dataset on row and column indices.
#'
#' Get a dataset connection object from an HDF5 database 
#' ('.h5') file and return the indexed table subset.
#' @param ri rows indices in dataset.
#' @param ci columns indices in dataset.
#' @param dsn Name of dataset or group of dataset to connect with.
#' @param dbn Path to h5 database file.
#' @return HDF5 database connection object.
#' @examples 
#' # get red signal for first 2 probe addresses, first 3 samples
#' st <- hread(1:3, 1:2, d = "redsignal", dbn = "remethdb2.h5")
#' @seealso data_mdpost()
#' @export
hread = function(ri, ci, dsn = "redsignal", dbn = "remethdb2.h5"){
    return(rhdf5::h5read(dbn, dsn, index = list(ri, ci)))
}

#--------------------------
# Query the sample metadata
#--------------------------

#' Query and store sample metadata.
#'
#' Retrieves sample postprocessed metadata from an HDF5 database.
#'
#' @param dbn Path to HDF5 database file.
#' @param dsn Name or group path to HDF5 dataset containing 
#' postprocessed metadata.
#' @return Postprocessed metadata as a `data.frame`.
#' @examples 
#' # get all available sample metadata
#' mdp <- data_mdpost(dbn = "remethdb2.h5", dsn = "mdpost")
#' @seealso hread()
#' @export
data_mdpost = function(dbn = "remethdb2.h5", dsn = "mdpost"){
    mdp <- as.data.frame(rhdf5::h5read(file = dbn, name = dsn),
        stringsAsFactors = FALSE)
    mcn <- paste(dsn, "colnames", sep = ".")
    colnames(mdp) <- rhdf5::h5read(file = dbn, name = mcn)
    return(mdp)
}

#-------------------------------------------------------
# Get SummarizedExperiment objects from dataset queries
#-------------------------------------------------------

#' Form an `RGChannelSet` from a signal table query
#'
#' Forms an object of `RGChannelSet` class from a query to the red 
#' and green signal tables.
#'
#' @param ldat List of raw signal data query results. Must include 2 
#' `data.frame` objects named 'redsignal' and 'greensignal'.
#' @param verbose Whether to post status messages.
#' @return Returns a `RGChannelSet` object from raw signal dataset queries.
#' @examples 
#' # get the list of datasets for all probe addresses, 3 samples
#' ldat = getrg(gsmv = c("GSM1235984", "GSM1236090", "GSM1506278"), 
#' data.type = "df", metadata = FALSE)
#' 
#' # get the rg set object
#' rg = rgse(ldat)
#' @seealso getrg()
#' @export
rgse = function(ldat, verbose = FALSE){
    if(!("greensignal" %in% names(ldat) & "redsignal" %in% names(ldat))){
        stop(paste0("Invalid datasets list passed."))
    }
    if(verbose){message("Matching probe IDs in signal matrices...")}
    rga <- ldat[["redsignal"]]; gga <- ldat[["greensignal"]]
    addrv <- unique(c(rownames(rga), rownames(gga))) # get CpG addresses
    rgs <- rga[rownames(rga) %in% addrv, ]; 
    rgs <- rgs[order(match(rownames(rgs), addrv)), ]
    ggs <- gga[rownames(gga) %in% addrv, ]; 
    ggs <- ggs[order(match(rownames(ggs), addrv)), ]
    cgidmatch = identical(rownames(rgs), rownames(ggs))
    if(!cgidmatch){stop("Couldn't probe IDs.")}
    if(verbose){message("Matching GSM IDs in signal matrices...")}
    gsmidv = unique(c(colnames(rgs), colnames(ggs)))
    rgf <- rgs[, colnames(rgs) %in% gsmidv]
    ggf <- ggs[, colnames(ggs) %in% gsmidv]
    rgf <- rgf[, order(match(colnames(rgf), gsmidv))]
    ggf <- ggf[, order(match(colnames(ggf), gsmidv))]
    gsmidmatch = identical(colnames(rgf), colnames(ggf))
    if(!gsmidmatch){stop("Couldn't match GSM IDs for signal data.")}
    if("metadata" %in% names(ldat)){
        if(verbose){message("Checking metadata...")}
        mdp <- ldat[["metadata"]]; mdp$gsm <- as.character(mdp$gsm)
        gsmov <- gsmidv[!gsmidv %in% mdp$gsm]; numo <- length(gsmov)
        if(numo > 0){
            if(verbose){message("Adding md NAs for ", length(gsmov), " GSMs.")}
            nmdat <- c(gsmov, rep(rep("NA", numo), ncol(mdp) - 1))
            nmm <- matrix(nmdat, nrow = numo); colnames(nmm) <- colnames(mdp)
            mdp <- rbind(mdp, nmm)
        }
        mdf <- mdp[mdp$gsm %in% gsmidv,]
        mdf <- mdf[order(match(mdf$gsm, gsmidv)),]
        mdf$gsm <- as.character(mdf$gsm)
        mdmatchid <- identical(mdf$gsm, colnames(rgf)) & 
            identical(mdf$gsm, colnames(ggf))
        if(!mdmatchid){stop("Couldn't match GSM IDs for md and signal data.")}
        rownames(mdf) <- mdf$gsm
    }
    anno <- c("IlluminaHumanMethylation450k", "ilmn12.hg19")
    names(anno) <- c("array", "annotation")
    rgi = minfi::RGChannelSet(Green = ggf, Red = rgf, annotation = anno)
    if("metadata" %in% names(ldat)){
        minfi::pData(rgi) <- S4Vectors::DataFrame(mdf)
    }
    return(rgi)
}

#' Query and store data from the signal tables
#'
#' Retrieves query matches from raw signal HDF5 datasets. 
#' Handles identity queries to rows (GSM IDs) or columns 
#' (CpG probe addresses). Returns query matches either 
#' as a list of 2 `data.frame`s or a single `RGChannelSet` 
#' object.
#'
#' @param dbn Name of the HDF5 database file.
#' @param gsmv Vector valid GSM IDs (rows) to query, 
#' either NULL or vector of length > 2 valid GSM IDs, 
#' or `all.gsm` should be TRUE.
#' @param cgv Vector of valid CpG probe addresses (columns) query,
#' either NULL or a vector of valid probe addresses, or `all.cg` 
#' should be TRUE.
#' @param data.type Format for returned query matches, either as 
#' datasets 
#' 'df' or `RGChannelSet` 'se' object.
#' @param dsv Vector of raw signal datasets or group paths to query, 
#' including both the red channel 'redsignal' and green channel 
#' 'greensignal' datasets.
#' @param all.gsm Whether to query all available GSM IDs.
#' @param all.cgv Whether to query all available CpG probe addresses.
#' @param metadata Whether to access available postprocessed 
#' metadata for queries samples.
#' @param md.dsn Name of metadata dataset in h5 file.
#' @param verbose Whether to post status messages.
#' @return Returns either an `RGChannelSet` or list of 
#' `data.frame` objects from dataset query matches.
#' @examples
#' gsml = c("GSM1235984", "GSM1236090", "GSM1506278")
#' # get list of data tables for a query
#' ldat = getrg(gsmv = gsml, data.type = "df")
#' 
#' # get the RGChannel set object for a query
#' rgset = getrg(gsmv = gsml, data.type = "se")
#' @seealso rgse()
#' @export
getrg = function(gsmv = NULL, cgv = NULL,
    dbn = "remethdb2.h5", data.type = c("se", "df"),
    dsv = c("redsignal", "greensignal"), all.gsm = FALSE, 
    all.cg = TRUE, metadata = TRUE, md.dsn = "mdpost", 
    verbose = FALSE){
    if(length(gsmv) < 2 | length(cgv) == 0 & !all.cg | all.gsm & all.cg){
        stop("Invalid query indices. Review GSM and probe ID args.")
    }
    ldat <- list() # datasets list
    for(d in dsv){
        if(verbose){message("Working on ", d, "...")}
        if(d %in% c("redsignal", "greensignal")){
            rnd <- rhdf5::h5read(dbn, paste(d, "rownames", sep = ".")) 
            rnd <- gsub("\\..*", "", rnd) # clean GSM IDs
            cnd <- rhdf5::h5read(dbn, paste(d, "colnames", sep = ".")) 
            cgvp = ifelse(all.cg, seq(1, length(cnd), 1),
                which(cnd %in% cgv))  
            gsmvp = ifelse(all.gsm, seq(1, length(rnd), 1), 
                which(rnd %in% gsmv))
            if(length(gsmvp) < 2){stop("Not enough valid GSM IDs found.")}
        }
        ddat <- hread(ri = gsmvp, ci = cgvp, d, dbn)
        rownames(ddat) <- rnd[gsmvp]; colnames(ddat) <- cnd[cgvp]
        ldat[[d]] <- t(ddat) # append transpose of data
    }
    if(metadata){
        mdpost <- data_mdpost(dbn = dbn, dsn = md.dsn)
        mdpost$gsm <- as.character(mdpost$gsm)
        ldat[["metadata"]] <- mdpost[mdpost$gsm %in% gsmv,]
    }
    if(data.type == "df"){
        if(verbose){message("Returning the datasets list...")}
        robj <- ldat
    }
    if(data.type == "se"){
        if(verbose){message("Forming the RGChannelSet...")}
        robj <- rgse(ldat = ldat, verbose = verbose)
    }
    rhdf5::h5closeAll() # close all open connections
    return(robj)
}
