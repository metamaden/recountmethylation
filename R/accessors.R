#!/usr/bin/env R

#' Get DNAm assay data.
#'
#' Uses RCurl to recursively download latest H5SE and HDF5 data objects the from server.
#' @param which.dn  Type of data dir to be downloaded.
#' @param url Server URL containing assay data.
#' @param dfp Target local directory for downloaded files.
#' @param verbose Whether to return verbose messages.
#' @return New filepath to dir with downloaded data.
#' @examples 
#' get_rmdl("h5se-test_gr", verbose = TRUE)
#' @export
get_rmdl <- function(which.dn = c("h5se-test_gr", "h5se_gr", 
                                  "h5se_gm", "h5se_rg", "\\.h5"),
                     url = "https://recount.bio/data/", 
                     dfp = "data", verbose = TRUE){
  if(verbose){message("Retrieving data dirnames from server...")}
  dn <- RCurl::getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  dn <- unlist(strsplit(dn, "\n"))
  catch.str <- paste0(".*", which.dn,".*")
  dn.catch <- grepl(catch.str, dn)
  dn <- unlist(dn)[dn.catch]
  dn.clean <- gsub('<.*', "", gsub('.*">', "", dn))
  if(!length(dn.clean) == 1){stop("There was a problem parsing the file string.")}
  dfp.dn <- paste(c(dfp, dn.clean), collapse = "/")
  if(!dir.exists(dfp.dn)){
    if(verbose){message("Making new dl dir ", dfp.dn)}
    dct1 <- try(dir.create(dfp))
    dct2 <- try(dir.create(dfp.dn))
    if(!(dct1 & dct2)){
      stop("Unable to make new dir path for download.",
           "Check your permissions.")
      }
  }
  if(verbose){message("Retrieving filenames from server...")}
  dn.url <- paste0(url, dn.clean)
  fl = RCurl::getURL(dn.url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  fl <- unlist(strsplit(fl, "\n"))
  fl.str <- paste0(c("assays.h5", "se.rds"), collapse = "|")
  fl.catch.str <- paste0(".*", fl.str,".*")
  fl.catch <- grepl(fl.catch.str, fl)
  fl <- unlist(fl)[fl.catch]
  fl.clean <- gsub('<.*', "", gsub('.*">', "", fl))
  if(verbose){message("Downloading files...")}
  dll <- list()
  for(i in 1:length(fl.clean)){
    f <- fl.clean[i]
    fpath <- paste0(c(dn.url, f), collapse = "")
    dll[[i]] <- try(download.file(fpath, paste0(dfp.dn, f)))
  }
  if(dll[[1]] == 0 & dll[[2]] == 0){
    if(verbose){message("Download completed successfully. ",
                        "Returning new path...")}
    return(dfp.dn)
  } else{
    if(verbose){message("Download incomplete for file ", fl.clean[which(dll!=0)])}
    return(NULL)
  }
  return(NULL)
}

#' IDATs query
#'
#' Queries GDS IDATs and downloads, with exception handling.
#' Files are downloaded 
#' @param gsmvi Param description
#' @param ext Filename extension
#' @param expand Whether to expand compressed files.
#' @param sys.cmd System command to expand compressed files (if expand TRUE).
#' @param verbose Whether to show verbose messages (TRUE/FALSE)
#' @param dfp Download directory
#' @param burl Base URL string for RCurl query
#' @return Describes returned object
#' @examples
#' gsmvi <- c("GSM2465267", "GSM2814572")
#' gds_idatquery(gsmvi)
#' @export
gds_idatquery <- function(gsmvi, ext = "gz", expand = TRUE, 
                          sys.cmd = "gunzip ", verbose = FALSE, dfp = "./idats/",
                          burl = paste0("ftp://ftp.ncbi.nlm.nih.gov/",
                                        "geo/samples/")){
  bnv <- fnv <- c()
  for(gsmi in gsmvi){
    # format URL for query to GDS
    url = paste0(burl, substr(gsmi, 1, nchar(gsmi)-3), 
                 paste(rep("n", 3), collapse = ""), 
                 "/", gsmi, "/suppl/")
    # get urls to idats
    fn = RCurl::getURL(url, ftp.use.epsv = FALSE, 
                       dirlistonly = TRUE)
    fn <- gsub("\\\r", "", unlist(strsplit(fn, "\n")))
    idat.str <- paste0("\\.idat\\.", ext)
    idat.catch <- grepl(idat.str, fn)
    fn <- unlist(fn)[idat.catch]
    # eval conditions
    check.cond <- length(fn) == 2
    fn.grn <- fn[grepl(paste0(".*Grn.idat\\.", ext, "($)"), fn)]
    fn.red <- fn[grepl(paste0(".*Red.idat\\.", ext, "($)"), fn)]
    check.cond <- c(check.cond, 
                    length(fn.grn) > 0, 
                    length(fn.red) > 0)
    if(check.cond[1] & check.cond[2] & 
       check.cond[3]){
      # retain idat basenames
      rg.patt <- "_Red.*|_Grn.*"
      match.rg <- gsub(rg.patt, "", fn)
      idatl <- unique(match.rg)
      bnv = c(bnv, idatl)
      for(f in fn){
        fnv <- c(fnv, f)
        url.dlpath <- paste(url, f, sep = "")
        dest.fpath <- paste(dfp,f, sep = "")
        download.file(url.dlpath, dest.fpath)
        if(expand){
          if(verbose){message("Expanding compressed file...")}
          system(paste0(sys.cmd, dest.fpath))
        }
        message(f)
      }
    } else{
      if(verbose){message("Query didn't identify",
                          " red and grn IDATs for ", 
                          gsmi)}
    }
    if(verbose){message("Finished query for: ", gsmi)}
  }
  return(list("basenames" = bnv, "filenames" = fnv))
}

#' Handles IDAT queries to GDS and returns assay data as `RGChannelSet`.
#'
#' Queries and downloads GSM IDAT files in GDS, and returns
#' assay data as an `RGChannelSet`.
#' ('.h5') file and return the indexed table subset.
#' @param gsmvi A vector of GSM IDs (alphanumeric character strings).
#' @param rmdl Whether to remove downloaded .*idat files (default TRUE).
#' @param ext Filename extension (default "gz").
#' @param verbose Whether to show verbose messages (default FALSE).
#' @param dfp Destination file dir for idats.
#' @param burl Base URL string for idat query (default "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/").
#' @return object of class `RGChannelSet`
#' @examples
#' gsmvi <- c("GSM2465267", "GSM2814572")
#' rg <- gds_idat2rg(gsmvi)
#' @export
gds_idat2rg <- function(gsmvi, rmdl = TRUE, ext = "gz", 
                        verbose = FALSE, dfp = "./idats/", 
                        burl = paste0("ftp://ftp.ncbi.nlm.nih.gov/",
                                      "geo/samples/")){
  dn = "" # download idats to cwd
  bnv = c() # store the idat basenames
  rt = try(gds_idatquery(gsmvi, burl, ext, verbose, dfp)) # idat query and download
  if(!is(rt) == "RGChannelSet"){
    stop("Process ended with the following message: ",
         rt[1])
  }
  rgdl = minfi::read.metharray(basenames = rt[["basenames"]])
  # file cleanup
  if(rmdl){
    if(verbose){message("Removing downloaded files...")}
    for(f in rt[["filenames"]]){
      file.remove(f)
    }
  }
  return(rgdl)
}

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
#' @export
hread <- function(ri, ci, dsn = "redsignal", dbn = "remethdb2.h5"){
    return(rhdf5::h5read(dbn, dsn, index = list(ri, ci)))
}

#' Retrieve available sample metadata
#'
#' Retrieves all sample metadata from an HDF5 database.
#'
#' @param dbn Path to HDF5 database file.
#' @param dsn Name or group path to HDF5 dataset containing 
#' the sample metadata and learned annotations.
#' @return A `data.frame` of all the sample metadata
#' @examples 
#' # get all available sample metadata
#' mdp <- data_mdpost(dbn = "remethdb2.h5", dsn = "mdpost")
#' @seealso hread()
#' @export
data_mdpost <- function(dbn = "remethdb2.h5", dsn = "mdpost"){
    mdp <- as.data.frame(rhdf5::h5read(file = dbn, name = dsn),
        stringsAsFactors = FALSE)
    mcn <- paste(dsn, "colnames", sep = ".")
    colnames(mdp) <- rhdf5::h5read(file = dbn, name = mcn)
    return(mdp)
}

#' Match two datasets
#'
#' Match the character vectors of row or column names 
#' for 2 datasets.
#' 
#' @param ds1 First dataset to match
#' @param ds2 Second dataset to match
#' @param mi1 Match index of ds1 (either "rows" or "columns")
#' @param mi2 Match index of ds2 (either "rows" or "columns")
#' @param subset.match If index lengths don't match, match on the 
#' common subset instead
#' @return A list of the matched datasets.
#' @examples
#' # make 2 mismatched datasets
#' ds1 <- matrix(seq(1, 10, 1), nrow = 5)
#' ds2 <- matrix(seq(11, 20, 1), nrow = 5)
#' rownames(ds1) <- rownames(ds2) <- paste0("row", seq(1, 5, 1))
#' colnames(ds1) <- colnames(ds2) <- paste0("col", c(1, 2))
#' ds2 <- ds2[rev(seq(1, 5, 1)), c(2, 1)]
#' 
#' # match row and column names
#' lmatched = matchds(d1, d2, mi1 = "rows", mi2 = "rows")
#' lmatched = matchds(lmatched[[1]], lmatched[[2]], 
#'     mi1 = "columns", mi2 = "columns)
#' @export
matchds_1to2 <- function(ds1, ds2, mi1 = c("rows", "columns"), 
    mi2 = c("rows", "columns"), subset.match = FALSE){
  if(mi1 == "rows"){ii1 = as.character(rownames(ds1))}
  if(mi1 == "columns"){ii1 = as.character(colnames(ds1))}
  if(mi2 == "rows"){ii2 = as.character(rownames(ds2))}
  if(mi2 == "columns"){ii2 = as.character(colnames(ds2))}
  if(subset.match){
      ii1 = ii1[ii1 %in% intersect(ii1, ii2)]
      ii2 = ii2[ii2 %in% intersect(ii1, ii2)]
  }
  ii1 = ii1[order(match(ii1, ii2))] # match 1 to 2
  if(!identical(ii1, ii2)){
      stop(paste0("Couldn't match provided indices. ",
          "Are they of the same length?"))
  }
  if(mi1 == "rows"){
      ds1m = ds1[rownames(ds1) %in% ii1,]
      ds1m = ds1m[order(match(as.character(rownames(ds1m)), ii1)),]
  }
  if(mi1 == "columns"){
      ds1m = ds1[, colnames(ds1) %in% ii1]
      ds1m = ds1m[, order(match(as.character(colnames(ds1m)), ii1))]
  }
  if(mi2 == "rows"){
      ds2m = ds2[rownames(ds2) %in% ii2,]
      ds2m = ds2m[order(match(as.character(rownames(ds2m)), ii2)),]
  }
  if(mi2 == "columns"){
      ds2m = ds2[, colnames(ds2) %in% ii2]
      ds2m = ds2m[, order(match(as.character(colnames(ds2m)), ii2))]
  }
  return(list(ds1 = ds1m, ds2 = ds2m))
}

#' Form an `RGChannelSet` from signal data
#'
#' Forms an object of `RGChannelSet` class from a signal data object 
#' returned from a query to the red and green signal tables.
#'
#' @param ldat List of raw signal data query results. Must include 2 
#' `data.frame` objects named 'redsignal' and 'greensignal'.
#' @param verbose Whether to post status messages.
#' @return Returns a `RGChannelSet` object from raw signal dataset queries.
#' @examples 
#' # get the list of datasets for all probe addresses, 3 samples
#' gsml = c("GSM1235984", "GSM1236090", "GSM1506278")
#' ldat = getrg(gsmv = gsml, data.type = "df", metadata = FALSE)
#' 
#' # get the rg set object
#' rg = rgse(ldat)
#' @seealso getrg()
#' @export
rgse <- function(ldat, verbose = FALSE){
    if(!("greensignal" %in% names(ldat) & "redsignal" %in% names(ldat))){
        stop(paste0("Invalid datasets list passed."))
    }
    if(verbose){
        message("Matching probe IDs in signal matrices...")
    }
    rga <- ldat[["redsignal"]]; gga <- ldat[["greensignal"]]
    lm.rg <- matchds_1to2(rga, gga, "rows", "rows")
    lm.rg <- matchds_1to2(lm.rg[[1]], lm.rg[[2]], "columns", "columns")
    cgidmatch <- identical(rownames(lm.rg[[1]]), rownames(lm.rg[[2]]))
    gsmidmatch <- identical(colnames(lm.rg[[1]]), colnames(lm.rg[[2]]))
    if(!cgidmatch){
        stop("Couldn't probe IDs.")
    }
    if(!gsmidmatch){
        stop("Couldn't match GSM IDs for signal data.")
    }
    if("metadata" %in% names(ldat)){
        gsmidv <- unique(c(colnames(lm.rg[[1]]), 
                           colnames(lm.rg[[2]])))
        if(verbose){
            message("Checking metadata...")
        }
        mdp <- ldat[["metadata"]]; mdp$gsm <- as.character(mdp$gsm)
        gsmov <- gsmidv[!gsmidv %in% mdp$gsm]; numo <- length(gsmov)
        if(numo > 0){
            if(verbose){message("Adding md NAs for ", length(gsmov), " GSMs.")}
            nmdat <- c(gsmov, rep(rep("NA", numo), ncol(mdp) - 1))
            nmm <- matrix(nmdat, nrow = numo); colnames(nmm) <- colnames(mdp)
            mdp <- rbind(mdp, nmm)
        }
        rownames(mdp) <- mdp$gsm
        lm.pr <- matchds_1to2(mdp, lm.rg[[1]], "rows", "columns")
        mdmatchid <- identical(rownames(lm.pr[[1]]), colnames(lm.pr[[2]]))
        if(!mdmatchid){
            stop("Couldn't match GSM IDs for md and signal data.")
        }
    }
    anno <- c("IlluminaHumanMethylation450k", "ilmn12.hg19")
    names(anno) <- c("array", "annotation")
    rgi <- minfi::RGChannelSet(Red = lm.rg[[1]], 
                               Green = lm.rg[[2]], annotation = anno)
    if("metadata" %in% names(ldat)){
        minfi::pData(rgi) <- S4Vectors::DataFrame(lm.pr[[1]])
    }
    return(rgi)
}

#' Query and store data from the signal tables
#'
#' Retrieves query matches from raw signal HDF5 datasets. 
#' Handles identity queries to rows (GSM IDs) or columns 
#' (CpG probe addresses). Returns query matches either 
#' as a list of datasets or a single `RGChannelSet` 
#' object, with the option of including sample metadata.
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
#' # make samples list
#' gsml = c("GSM1235984", "GSM1236090", "GSM1506278")
#' 
#' # get list of data tables for a query
#' ldat = getrg(gsmv = gsml, data.type = "df")
#' 
#' # get the RGChannel set object for a query
#' rgset = getrg(gsmv = gsml, data.type = "se")
#' @seealso rgse()
#' @export
getrg <- function(gsmv = NULL, cgv = NULL,
    dbn = "remethdb2.h5", data.type = c("se", "df"),
    dsv = c("redsignal", "greensignal"), all.gsm = FALSE, 
    all.cg = TRUE, metadata = TRUE, md.dsn = "mdpost", 
    verbose = FALSE){
    if((length(gsmv) < 2 & !all.gsm) | (length(cgv) == 0 & !all.cg) | 
       (all.gsm & all.cg)){
        stop("Invalid query Review GSM and probe ID args.")
    }
    ldat <- list() # datasets list
    for(d in dsv){
        if(verbose){message("Working on ", d, "...")}
        if(d %in% c("redsignal", "greensignal")){
            rnd <- rhdf5::h5read(dbn, paste(d, "rownames", sep = ".")) 
            rnd <- gsub("\\..*", "", rnd) # clean GSM IDs
            cnd <- rhdf5::h5read(dbn, paste(d, "colnames", sep = ".")) 
            if(all.cg){
              cgvp <- seq(1, length(cnd), 1)
            } else{
              cgvp <- which(cnd %in% cgv)
            }
            if(all.gsm){
              gsmvp <- seq(1, length(rnd), 1)
            } else{
              gsmvp <- which(rnd %in% gsmv)
            }
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
