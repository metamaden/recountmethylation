#!/usr/bin/env R


#' servermatrix
#'
#' Get a matrix of database files from the recount server. Called by get_rmdl.
#' @param dn Server data returned from RCurl.
#' @param printmatrix Whether to print the data matrix to console (default TRUE).
#' @param verbose Whether to show verbose messages (default FALSE).
#' @param recursive Whether to recursively grab file sizes for h5se object (default TRUE).
#' @returns dm matrix of server files and file metadata
#' @examples 
#' dn <- RCurl::getURL("https://recount.bio/data/", .opts = list(ssl.verifypeer = FALSE))
#' sm <- servermatrix(dn)
#' @seealso get_rmdl
#' @export
servermatrix <- function(dn, url = "https://recount.bio/data/", 
                         printmatrix = TRUE, verbose = FALSE, recursive = TRUE){
  dt <- unlist(strsplit(dn, "\r\n"))
  dt <- gsub('(.*\">|/</a>|</a>)', "", dt)
  dt <- dt[grepl("remethdb", dt)]
  drows <- lapply(as.list(dt), function(x){
    return(unlist(strsplit(gsub("[ ]+", ";", x), ";")))
  })
  dm <- do.call(rbind, drows)
  colnames(dm) <- c("filename", "date", "time", "size (bytes)")
  if(recursive){
    sv <- c() # file sizes vector
    fnv <- dm[grepl("h5se", dm[,1]), 1]
    fnexclude <- c()
    for(f in fnv){
      fv <- RCurl::getURL(paste0(url, f, "/"), dirlistonly = TRUE, 
                          .opts = list(ssl.verifypeer = sslver))
      fvv <- unlist(strsplit(fv, "\r\n"))
      which.start <- which(grepl("Index", fvv))[2] + 1
      which.end <- which(grepl("/pre", fvv)) - 1
      fvf <- fvv[which.start:which.end]
      fniv <- c()
      for(fni in fvf){
        name <- gsub('.*\">', '', gsub("</a>.*", "", fni))
        size <- gsub(".* ", "", fni)
        fniv <- c(fniv, paste0("`", name, "`", " = ", size))
      }
      # check for h5se completeness
      cond.assays <- length(fniv[grepl("assays", fniv)]) == 1
      cond.se <- length(fniv[grepl("se", fniv)]) == 1
      sv <- c(sv, paste(fniv, collapse = ";"))
      if(!(cond.assays & cond.se)){fnexclude <- c(fnexclude, f)}
    }
  }
  dm[grepl("h5se", dm[,1]), 4] <- sv
  dm <- dm[!dm[,1] %in% fnexclude,] # filter incomplete h5se files
  return(dm)
}

#' Get DNAm assay data.
#'
#' Uses RCurl to recursively download latest H5SE and HDF5 data objects the from server.
#' 
#' @param which.class  Class of file to download (either `rg`, `gm`, `gr`, or `test`).
#' @param which.type Type of file (either `h5` for HDF5 file or `h5se` for HDF5-SummarizedExperiment).
#' @param url Server URL containing assay data.
#' @param dfp Target local directory for downloaded files (default "downloads").
#' @param download Whether to download (TRUE) or return queried filename (FALSE).
#' @param verbose Whether to return verbose messages.
#' @param sslver Whether to use server certificate check (default FALSE).
#' @return New filepath to dir with downloaded data.
#' @examples 
#' get_rmdl("h5se-test_gr", verbose = TRUE)
#' @export
get_rmdl <- function(fn = NULL, show.files = FALSE, 
                     which.class = c("rg", "gm", "gr", "test"),
                     which.type = c("h5se", "h5"),
                     url = "https://recount.bio/data/", 
                     dfp = "downloads", download = TRUE,
                     verbose = TRUE, sslver = FALSE){
  if(verbose){message("Retrieving data dirnames from server...")}
  # set up rcurl call
  ftpuseopt <- ifelse(show.files, FALSE, TRUE)
  dirlistopt <- ifelse(show.files, FALSE, TRUE)
  dn <- RCurl::getURL(url, ftp.use.epsv = ftpuseopt, dirlistonly = dirlistopt,
                      .opts = list(ssl.verifypeer = sslver))
  if(verbose){message("Getting file data from server.")}
  sm <- servermatrix(dn)
  if(show.files & !download){
    if(verbose){message("Showing server file data")}
    prmatrix(sm)
  }
  if(is.null(fn)){
    # clean query results
    str1 <- ifelse(which.type == "h5", "\\.", ".*")
    str2 <- ifelse(which.type == "h5", "$", ".*")
    typestr <- paste0(str1, which.type, str2)
    filt.type <- grepl(paste0(".*", which.type,".*"), sm[,1])
    filt.all <- filt.type & grepl(paste0(".*", which.class,".*"), sm[,1])
    dnc <- sm[filt.all, 1]
    if(length(dnc) == 0){
      stop("No files of class and type found.")
    } else if(length(dnc) > 1){
      tsv <- as.numeric(gsub(".*_", "", dnc)) # timestamps
      tsf <- which(tsv == max(tsv))[1] # first instance
      dnc <- dnc[tsf]
    }
  } else{
    dnc <- fn
    check.cond1 <- grepl("(\\.h5$|.*h5se.*)", dnc)
    check.cond2 <- dnc %in% sm[,1]
    condpass <- check.cond1 & check.cond2
    if(!condpass){stop("Provided fn not found on server.")}
  }
  if(!download){
    message("File confirmed on server. Returning.")
    return(dnc)
  }
  # manage download loc
  dct1 <- ifelse(!dir.exists(dfp) & !dfp == "", try(dir.create(dfp)), TRUE)
  dfp.dn <- paste(c(dfp, dnc), collapse = "/")
  if(which.type == "h5"){
    dct2 <- try(file.create(dfp.dn))
  } else{dct2 <- try(dir.create(dfp.dn))}
  if(!(dct1 & dct2)){
    stop("There was a problem with the download dest. Do you have write access?")
  }
  dn.url <- paste0(url, dnc)
  if(which.type == "h5"){fl.clean <- ""} else{
    fnv <- unlist(strsplit(as.character(sm[2,4]), ";"))
    fl.clean <- gsub("( .*|`)", "", fnv)
  }
  if(verbose){message("Downloading file(s)...")}
  dll <- list()
  for(i in 1:length(fl.clean)){
    fpath <- paste0(c(dn.url, fl.clean[i]), collapse = "")
    cf = RCurl::CFILE(paste0(dfp.dn, fl.clean[i]), mode="wb")
    dll[[i]] <- try(RCurl::curlPerform(url = fpath, writedata = cf@ref,
                       .opts = list(ssl.verifypeer = sslver)))
  }
  if(length(dll[dll==0]) == length(dll)){
    if(verbose){message("Finished download, returning file path.")}
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
#' 
#' @param gsmvi Vector of GSM IDs
#' @param ext Filename extension
#' @param expand Whether to expand compressed files.
#' @param sys.cmd System command to expand compressed files (if expand TRUE).
#' @param verbose Whether to show verbose messages (TRUE/FALSE)
#' @param dfp Download directory
#' @param burl Base URL string for RCurl query
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
        download.file(url.dlpath, dest.fpath)
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

#' Handles IDAT queries to GDS and returns assay data as `RGChannelSet`.
#'
#' Queries and downloads GSM IDAT files in GDS, and returns
#' assay data as an `RGChannelSet`.
#' 
#' @param gsmvi A vector of GSM IDs (alphanumeric character strings).
#' @param rmdl Whether to remove downloaded .*idat files (default TRUE).
#' @param ext Filename extension (default "gz").
#' @param dfp Destination file dir for idats.
#' @param burl Base URL string for idat query (default "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/").
#' @return object of class `RGChannelSet`
#' @examples
#' gsmvi <- c("GSM2465267", "GSM2814572")
#' rg <- gds_idat2rg(gsmvi)
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
#' Get a dataset connection object from an HDF5 database 
#' ('.h5') file and return the indexed table subset.
#' 
#' @param ri rows indices in dataset.
#' @param ci columns indices in dataset.
#' @param dsn Name of dataset or group of dataset to connect with.
#' @param dbn Path to h5 database file.
#' @return HDF5 database connection object.
#' @examples
#' # Get tests data pointer
#' path = system.file("extdata", "testh5", package = "recountmethylation")
#' fn = list.files(path)
#' dbpath = paste0(path, "/", fn)
#' # red signal, first 2 assay addr, 3 samples
#' reds <- hread(1:3, 1:2, d = "redsignal", dbn = dbpath)
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
#' path = system.file("extdata", "testh5", package = "recountmethylation")
#' fn = list.files(path)
#' dbpath = paste0(path, "/", fn)
#' mdp <- data_mdpost(dbn = dbpath, dsn = "mdpost")
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
#' path <- system.file("extdata", "testh5", package = "recountmethylation")
#' dbpath <- paste0(path, "/remethdbtest.h5")
#' gsml = c("GSM1235984", "GSM1236090", "GSM1506278")
#' ldat = getrg(gsmv = gsml, dbn = dbpath, data.type = "df", metadata = FALSE)
#' rg = rgse(ldat) # get the rg set object
#' class(rg)
#' dim(minfi::getBeta(rg))
#' @seealso getrg()
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
#' @param all.cg Whether to query all available CpG probe addresses.
#' @param metadata Whether to access available postprocessed 
#' metadata for queried samples.
#' @param md.dsn Name of metadata dataset in h5 file.
#' @param verbose Whether to post status messages.
#' @return Returns either an `RGChannelSet` or list of 
#' `data.frame` objects from dataset query matches.
#' @examples
#' # make samples list
#' path <- system.file("extdata", "testh5", package = "recountmethylation")
#' dbpath <- paste0(path, "/remethdbtest.h5")
#' gsml = c("GSM1235984", "GSM1236090", "GSM1506278")
#' # get list of data tables for a query
#' ldat = getrg(gsmv = gsml, dbn = dbpath, data.type = "df")
#' names(ldat)
#' dim(ldat[[1]])
#' dim(ldat[[2]])
#' dim(ldat[[3]])
#' @seealso rgse()
#' @export
getrg <- function(gsmv = NULL, cgv = NULL,
    dbn = "remethdb2.h5", data.type = c("se", "df"),
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
    }
    if(data.type == "se"){
        if(verbose){message("Forming the RGChannelSet...")}
        robj <- rgse(ldat = ldat, verbose = verbose)
    }
    rhdf5::h5closeAll() # close all open connections
    return(robj)
}
