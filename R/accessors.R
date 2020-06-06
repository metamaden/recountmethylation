#!/usr/bin/env R

#' Query and download IDATs from GEO Data Sets
#'
#' Queries GEO Data Sets for IDATs, and downloads available IDATs. This uses 
#' anticipated string pattern to construct the URL path for the query. 
#' IDATs are detected from the supplement for a GSE record.
#' 
#' @param gsmvi Vector of valid GSM IDs.
#' @param ext Filename extension.
#' @param expand Whether to expand compressed files.
#' @param sys.cmd System command to expand compressed files if expand TRUE 
#' (default "gunzip").
#' @param verbose Whether to show verbose messages (default FALSE).
#' @param dfp Destination directory for downloads.
#' @param burl Base URL string for RCurl query.
#' @return Lists the basename paths and filenames of IDATs downloaded.
#' @examples
#' gsmvi <- c("GSM2465267", "GSM2814572")
#' gds_idatquery(gsmvi)
#' @export
gds_idatquery <- function(gsmvi, ext = "gz", expand = TRUE, 
                          sys.cmd = "gunzip ", verbose = FALSE, dfp = "./idats/",
                          burl = paste0("ftp://ftp.ncbi.nlm.nih.gov/",
                                        "geo/samples/")){
  bnv <- fnv <- c()
  if(verbose){message("Checking dest dir dfp.")}
  if(!dir.exists(dfp)){
    message("Making new dest dir dfp.")
    tdir <- try(dir.create(dfp))
    if(!tdir){stop("There was an issue making the new dest dir dfp.")}
  }
  for(gsmi in gsmvi){
    url = paste0(burl, substr(gsmi, 1, nchar(gsmi)-3), 
                 paste(rep("n", 3), collapse = ""), "/", gsmi, "/suppl/")
    # get urls to idats
    fn = RCurl::getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
    fn <- gsub("\\\r", "", unlist(strsplit(fn, "\n")))
    idat.str <- paste0("\\.idat\\.", ext)
    idat.catch <- grepl(idat.str, fn)
    fn <- unlist(fn)[idat.catch]
    check.cond <- length(fn) == 2
    fn.grn <- fn[grepl(paste0(".*Grn.idat\\.", ext, "($)"), fn)]
    fn.red <- fn[grepl(paste0(".*Red.idat\\.", ext, "($)"), fn)]
    check.cond <- c(check.cond, length(fn.grn) > 0, length(fn.red) > 0)
    if(check.cond[1] & check.cond[2] & 
       check.cond[3]){
      idatl <- unique(gsub("_Red.*|_Grn.*", "", fn))
      bnv = c(bnv, paste(dfp, idatl, sep = "")) # gsm basenames
      for(f in fn){
        url.dlpath <- paste(url, f, sep = "")
        dest.fpath <- paste(dfp, f, sep = "")
        utils::download.file(url.dlpath, dest.fpath)
        fnv <- c(fnv, dest.fpath)
        if(expand){
          if(verbose){message("Expanding compressed file...")}
          system(paste0(sys.cmd, dest.fpath))
        }
        message(f)
      }
    } else{
      if(verbose){message("Query didn't identify", 
                          " red and grn IDATs for ", gsmi)}
    }
    if(verbose){message("Finished query for: ", gsmi)}
  }
  return(list("basenames" = bnv, "filenames" = fnv))
}

#' Get IDATs as an RGChannelSet from GEO/GDS
#'
#' Queries and downloads GSM IDAT files in GEO Data Sets db, then 
#' returns the assay data as an "RGChannelSet", calling gds_idatquery()
#'  then minfi::read.metharray().
#' 
#' @param gsmvi A vector of GSM IDs (alphanumeric character strings).
#' @param rmdl Whether to remove downloaded IDAT files when finished
#'  (default TRUE).
#' @param ext Extension for downloaded files (default "gz").
#' @param dfp Destination for IDAT downloads.
#' @param burl Base URL string for the IDAT query (default "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/").
#' @return An RGChannelSet object
#' @examples
#' gsmvi <- c("GSM2465267", "GSM2814572")
#' rg <- gds_idat2rg(gsmvi)
#' @seealso gds_idatquery(), read.metharray()
#' @export
gds_idat2rg <- function(gsmvi, rmdl = TRUE, ext = "gz", dfp = "./idats/", 
                        burl = paste0("ftp://ftp.ncbi.nlm.nih.gov/",
                                      "geo/samples/")){
  dn = "" # download idats to cwd
  bnv = c() # store the idat basenames
  rt <- try(gds_idatquery(gsmvi = gsmvi, ext = ext, 
                          dfp = dfp, burl = burl))
  if(is(rt)[1] == "try-error"){stop("Process ended with message: ", rt[1])}
  rgdl = minfi::read.metharray(basenames = rt[["basenames"]])
  if(rmdl){
    message("Removing downloaded files...")
    for(f in rt[["filenames"]]){
      file.remove(f)
      file.remove(gsub(paste0("\\.", ext), "", f))
    }
  }
  return(rgdl)
}

#' Query and store an HDF5 dataset on row and column indices.
#'
#' Connect to an HDF5 database h5 file with rhdf5::h5read(). 
#' Returns the subsetted data.
#' 
#' @param ri Row indices in dataset.
#' @param ci Column indices in dataset.
#' @param dsn Name of dataset or group of dataset to connect with.
#' @param dbn Path to h5 database file.
#' @return HDF5 database connection object.
#' @examples
#' # Get tests data pointer
#' path <- system.file("extdata", "h5test", package = "recountmethylation")
#' fn <- list.files(path)
#' dbpath = paste0(path, "/", fn)
#' # red signal, first 2 assay addr, 3 samples
#' reds <- hread(1:2, 1:3, d = "redsignal", dbn = dbpath)
#' dim(reds) # [1] 2 3
#' @seealso h5read()
#' @export
hread <- function(ri, ci, dsn = "redsignal", dbn = "remethdb2.h5"){
    return(rhdf5::h5read(dbn, dsn, index = list(ri, ci)))
}

#' Retrieve all available sample metadata from an HDF5 database.
#'
#' Retrieve all available sample metadata in a dataset from an 
#' HDF5 database. Returns data in metadata dataset "dsn" contained in 
#' an h5 file located at path "dbn."
#'
#' @param dbn Path to h5 HDF5 database file.
#' @param dsn Name or group path to HDF5 dataset containing 
#' the sample metadata and learned annotations.
#' @return data.frame of available sample metadata.
#' @examples 
#' path <- system.file("extdata", "h5test", package = "recountmethylation")
#' fn = list.files(path)
#' dbpath = paste0(path, "/", fn)
#' mdp <- data_mdpost(dbn = dbpath, dsn = "mdpost")
#' dim(mdp) # [1]  2 19
#' @seealso hread()
#' @export
data_mdpost <- function(dbn = "remethdb2.h5", dsn = "mdpost"){
    mdp <- as.data.frame(rhdf5::h5read(file = dbn, name = dsn),
        stringsAsFactors = FALSE)
    mcn <- paste(dsn, "colnames", sep = ".")
    colnames(mdp) <- rhdf5::h5read(file = dbn, name = mcn)
    return(mdp)
}

#' Match two datasets on rows and columns
#'
#' Match 2 datasets using the character vectors of row or column 
#' names. This is used to assemble an "RGChannelSet" from a query 
#' to an h5 dataset.
#' 
#' @param ds1 First dataset to match
#' @param ds2 Second dataset to match
#' @param mi1 Match index of ds1 (either "rows" or "columns")
#' @param mi2 Match index of ds2 (either "rows" or "columns")
#' @param subset.match If index lengths don't match, match on the 
#' common subset instead
#' @return A list of the matched datasets.
#' @examples
#' # get 2 data matrices
#' ds1 <- matrix(seq(1, 10, 1), nrow = 5)
#' ds2 <- matrix(seq(11, 20, 1), nrow = 5)
#' rownames(ds1) <- rownames(ds2) <- paste0("row", seq(1, 5, 1))
#' colnames(ds1) <- colnames(ds2) <- paste0("col", c(1, 2))
#' ds2 <- ds2[rev(seq(1, 5, 1)), c(2, 1)]
#' # match row and column names
#' lmatched = matchds_1to2(ds1, ds2, mi1 = "rows", mi2 = "rows")
#' lmatched = matchds_1to2(lmatched[[1]], lmatched[[2]], mi1 = "columns", mi2 = "columns")
#' # check matches
#' ds1m <- lmatched[[1]]
#' ds2m <- lmatched[[2]]
#' identical(rownames(ds1m), rownames(ds2m))
#' identical(colnames(ds1m), colnames(ds2m))
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

#' Form an RGChannelSet from a list containing signal data matrices
#'
#' Forms an RGChannelSet from signal data list. This is called by 
#' certain queries to h5 files.
#' 
#' @param ldat List of raw signal data query results. Must include 2 
#' data.frame objects named "redsignal" and "greensignal."
#' @param verbose Whether to post status messages.
#' @return Returns a RGChannelSet object from raw signal dataset queries.
#' @examples 
#' path <- system.file("extdata", "h5test", package = "recountmethylation")
#' fn <- list.files(path)
#' dbpath = paste(path, fn, sep = "/")
#' rg = getrg(dbn = dbpath, all.gsm = TRUE, metadata = FALSE)
#' dim(rg) # [1] 11162     2
#' class(rg)
#' # [1] "RGChannelSet"
#' # attr(,"package")
#' # [1] "minfi"
#' @seealso getrg(), RGChannelSet()
#' @export
rgse <- function(ldat, verbose = FALSE){
    if(!("greensignal" %in% names(ldat) & "redsignal" %in% names(ldat))){
        stop(paste0("Invalid datasets list passed."))
    }
    if(verbose){message("Matching probe IDs in signal matrices...")}
    rga <- ldat[["redsignal"]]; gga <- ldat[["greensignal"]]
    lm.rg <- matchds_1to2(rga, gga, "rows", "rows")
    lm.rg <- matchds_1to2(lm.rg[[1]], lm.rg[[2]], "columns", "columns")
    cgidmatch <- identical(rownames(lm.rg[[1]]), rownames(lm.rg[[2]]))
    gsmidmatch <- identical(colnames(lm.rg[[1]]), colnames(lm.rg[[2]]))
    if(!cgidmatch){stop("Couldn't probe IDs.")}
    if(!gsmidmatch){stop("Couldn't match GSM IDs for signal data.")}
    if("metadata" %in% names(ldat)){
        gsmidv <- unique(c(colnames(lm.rg[[1]]), colnames(lm.rg[[2]])))
        if(verbose){message("Checking metadata...")}
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
        if(!mdmatchid){stop("Couldn't match GSM IDs for md and signal data.")}
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

#' Query and store data from h5 file signal tables
#'
#' Queries signal datasets in an h5 HDF5 database file. 
#' Handles identity queries to rows (GSM IDs) or columns 
#' (CpG probe addresses). Returns query matches either 
#' as a list of datasets or a single RGChannelSet, with 
#' option to include sample metadata.
#'
#' @param dbn Name of the HDF5 database file.
#' @param gsmv Vector valid GSM IDs (rows) to query, 
#' either NULL or vector of length > 2 valid GSM IDs, 
#' or "all.gsm" should be TRUE.
#' @param cgv Vector of valid CpG probe addresses (columns) query,
#' either NULL or a vector of valid probe addresses, or "all.cg" 
#' should be TRUE.
#' @param data.type Format for returned query matches, either as 
#' datasets "df" or RGChannelSet "se" object.
#' @param dsv Vector of raw signal datasets or group paths to query, 
#' including both the red channel 'redsignal' and green channel 
#' 'greensignal' datasets.
#' @param all.gsm Whether to query all available GSM IDs.
#' @param all.cg Whether to query all available CpG probe addresses.
#' @param metadata Whether to access available postprocessed 
#' metadata for queried samples.
#' @param md.dsn Name of metadata dataset in h5 file.
#' @param verbose Whether to post status messages.
#' @return Returns either an RGChannelSet or list of 
#' data.frame objects from dataset query matches.
#' @examples
#' path <- system.file("extdata", "h5test", package = "recountmethylation")
#' fn <- list.files(path)
#' dbpath = paste(path, fn, sep = "/")
#' rg = getrg(dbn = dbpath, all.gsm = TRUE, metadata = FALSE)
#' dim(rg) # [1] 11162     2
#' class(rg)
#' # [1] "RGChannelSet"
#' # attr(,"package")
#' # [1] "minfi"
#' @seealso rgse()
#' @export
getrg <- function(dbn, gsmv = NULL, cgv = NULL, data.type = c("se"),
    dsv = c("redsignal", "greensignal"), all.gsm = FALSE, 
    all.cg = TRUE, metadata = TRUE, md.dsn = "mdpost", 
    verbose = FALSE){
    if((length(gsmv) < 2 & !all.gsm) | (length(cgv) == 0 & !all.cg)){
        stop("Invalid query, review GSM and probe ID args.")
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
              if(verbose){message("Retrieving all available CpG IDs.")}
            } else{
              cgvp <- which(cnd %in% cgv)
              if(verbose){message("Found ", length(cgvp)," of ", length(cgv), " CpG addresses.")}
            }
            if(all.gsm){
              gsmvp <- seq(1, length(rnd), 1)
              if(verbose){message("Retrieving all available GSM IDs.")}
            } else{
              gsmvp <- which(rnd %in% gsmv)
              if(verbose){message("Found ", length(gsmvp)," of ", length(gsmv), " GSM IDs.")}
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
        ldat[["metadata"]] <- mdpost[mdpost$gsm %in% rownames(ddat),]
    }
    if(data.type == "df"){
        if(verbose){message("Returning the datasets list...")}
        robj <- ldat
    } else if(data.type == "se"){
        if(verbose){message("Forming the RGChannelSet...")}
        robj <- rgse(ldat = ldat, verbose = verbose)
    }
    rhdf5::h5closeAll() # close all open connections
    return(robj)
}
