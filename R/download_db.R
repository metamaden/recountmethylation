#!/usr/bin/env R

# Functions for downloading DNAm datasets/cross-study compilations from the server.

#' servermatrix
#'
#' Called by get_rmdl() to get a matrix of database files and file info from the server. 
#' Verifies valid versions and timestamps in filenames, and that h5se directories contain
#'  both an assays and an se.rds file.
#' @param dn Server data returned from RCurl.
#' @param sslver Whether to use SSL certificate authentication for server connection (default FALSE).
#' @param url Server website url.
#' @param printmatrix Whether to print the data matrix to console (default TRUE).
#' @param verbose Whether to show verbose messages (default FALSE).
#' @param recursive Whether to recursively grab file sizes for h5se objects (default TRUE).
#' @returns dm matrix of server files and file metadata
#' @examples 
#' dn <- RCurl::getURL("https://recount.bio/data/", .opts = list(ssl.verifypeer = FALSE))
#' sm <- servermatrix(dn)
#' @seealso get_rmdl
#' @export
servermatrix <- function(dn, sslver = FALSE, url = "https://recount.bio/data/", 
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
#' This is currently wrapped in the getdb() functions. If tryload = TRUE, successful
#' download completion is tested with either HDF5Array::loadHDF5SummarizedExperiment()
#'  for h5se files or rhdf5::h5ls() for h5 files.
#' 
#' @param which.class  Either "rg", "gm", "gr", or "test" for RGChannelSet, MethylSet, 
#' GenomicRatioSet, or 2-sample subset (default "test").
#' @param which.type Either "h5se" for an HDF5-SummarizedExperiment (default) or 
#' "h5" for an HDF5 database.
#' @param fn Name of file on server to download (optional, default NULL).
#' @param dfp Download destination directory (default "downloads").
#' @param url The server URL to locate files for download.
#' @param show.files Whether to print server file data to console (default FALSE).
#' @param download Whether to download (TRUE) or return queried filename (FALSE).
#' @param tryload Whether to try loading downloaded data (default TRUE).
#' @param verbose Whether to return verbose messages.
#' @param sslver Whether to use server certificate check (default FALSE).
#' @return New filepath to dir with downloaded data.
#' @examples 
#' path <- get_rmdl(which.class = "test", which.type = "h5se", tryload = FALSE)
#' base::unlink("downloads", recursive = TRUE)
#' @seealso servermatrix(), getURL(), loadHDF5SummarizedExperiment(), h5ls()
#' @export
get_rmdl <- function(which.class = c("rg", "gm", "gr", "test"), 
                     which.type = c("h5se", "h5"), fn = NULL, dfp = "downloads", 
                     url = "https://recount.bio/data/", show.files = FALSE, download = TRUE, 
                     tryload = TRUE, verbose = TRUE, sslver = FALSE){
  if(verbose){message("Retrieving data dirnames from server...")}
  # set up rcurl call
  ftpuseopt <- ifelse(show.files, FALSE, TRUE)
  dirlistopt <- ifelse(show.files, FALSE, TRUE)
  dn <- RCurl::getURL(url, ftp.use.epsv = ftpuseopt, dirlistonly = dirlistopt,
                      .opts = list(ssl.verifypeer = sslver))
  if(verbose){message("Getting file data from server.")}
  sm <- servermatrix(dn = dn, sslver = sslver)
  if(show.files){prmatrix(sm)}
  if(is.null(fn)){
    # clean query results
    str1 <- ifelse(which.type == "h5", "\\.", ".*")
    str2 <- ifelse(which.type == "h5", "$", ".*")
    typestr <- paste0(str1, which.type, str2)
    filt.type <- grepl(typestr, sm[,1])
    filt.all <- filt.type & grepl(paste0(".*", which.class,".*"), sm[,1])
    dnc <- sm[filt.all, 1]
    if(!which.class == "test"){dnc <- dnc[!grepl("test", dnc)]}
    if(length(dnc) > 1){
      tsstr <- gsub("(.*_|\\.h5)", "", dnc)
      tsv <- suppressWarnings(as.numeric(tsstr)) # timestamps
      tsv <- tsv[!is.na(tsv)] # rm files without timestamp
      tsf <- which(tsv == max(tsv))[1] # first instance
      dnc <- dnc[tsf]
    }
    if(length(dnc) == 0){stop("No files of class and type found.")}
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
  dfp.dn <- paste(dfp, dnc, sep = "/")
  if(which.type == "h5"){
    dct2 <- try(file.create(dfp.dn))
  } else{dct2 <- try(dir.create(dfp.dn))}
  if(!(dct1 & dct2)){stop("Problem handling download destination.")}
  dn.url <- paste0(url, dnc)
  if(which.type == "h5"){fl.clean <- ""} else{fl.clean <- c("assays.h5", "se.rds")}
  if(verbose){message("Downloading file(s)...")}
  dll <- list() # download statuses list
  for(i in 1:length(fl.clean)){
    fpath <- ifelse(fl.clean[[i]] == "", dn.url, paste(dn.url, fl.clean[i], sep = "/"))
    destpath <- ifelse(fl.clean[[i]] == "", dfp.dn, paste(dfp.dn, fl.clean[i], sep = "/"))
    trydl = try(utils::download.file(url = fpath, destfile = destpath, 
                              .opts = list(ssl.verifypeer = sslver)))
  }
  if(length(dll[dll==0]) == length(dll)){
    if(verbose){message("Completed download.")}
    if(tryload){ message("Testing file load.")
      if(which.type == "h5se"){
        tdl <- try(HDF5Array::loadHDF5SummarizedExperiment(dfp.dn))
      } else{tdl <- try(rhdf5::h5ls(dfp.dn))}
      if(is(tdl)=="try-error"){message("Problem loading, download may be corrupt.")}
    }
    return(dfp.dn)
  } else{
    if(verbose){message("Download incomplete for file ", fl.clean[which(dll!=0)])}
    return(NULL)
  }
  return(NULL)
}

#' @name getdb
#' @rdname getdb
#'
#' @title Access database files.
#'
#' @description Combines download and load functions for databases. 
#' If "name" argument not provided, the latest available file is downloaded.
#' All files include metadata for the available samples.
#' 
#' There are 6 functions. Functions with "h5se" access HDF5-SummarizedExperiment 
#' files, and "h5" functions access HDF5 databases. The 4 h5se functions are 
#' "rg" (RGChannelSet), "gm" (MethylSet), "gr" (GenomicRatioSet), and "test" 
#' (data for 2 samples from "gr"). The 2 h5 functions are "rg" (red and green 
#' signal datasets), and "test" (data for 2 samples from "rg"). See vignette 
#' for details about file types and classes. 
#' 
#' @param name Database file name (optional, default NULL).
#' @param dfp Folder to search for database file specified by "name" 
#' (optional, default "downloads").
#' @param verbose Whether to return verbose messages (default FALSE).
#' @seealso dldb(), get_rmdl()
#' @return Either a SummarizedExperiment object for h5se functions, or a file path
#' for h5 functions.
NULL
#' @rdname getdb
#' @examples
#' path <- getdb_h5se_test()
#' base::unlink("downloads", recursive = TRUE)
#' @export
getdb_h5se_test <- function(name = NULL, dfp = "downloads", verbose = FALSE){
  dbpath <- FALSE
  if(!is.null(name)){
    dbpath <- paste(dfp, name, sep = "/")
    if(!file.exists(dbpath)){
      message("Dataset file not found at dfp.")
    } else{
      message("Dataset file found at dfp.")
    }
  } else{
    message("Downloading database...")
    dbpath <- try(get_rmdl(which.class = "test", dfp = dfp, which.type = "h5se", 
                           tryload = FALSE, verbose = verbose))
    if(!is(dbpath) == "try-errror"){
      message("Download completed.")
    } else{stop("Problem with download.")}
  }
  # parse load
  if(is(dbpath) == "try-error"){stop("Problem with dbpath.")} else{
    message("Loading database file.")
    dbf <- try(HDF5Array::loadHDF5SummarizedExperiment(dbpath))
  }
  if(is(dbf) == "try-error"){stop("Problem loading file.")} else{
    message("Database file loaded.")
    return(dbf)
  }
  return(NULL)
}
#' @rdname getdb
#' @export
getdb_h5_test <- function(name = NULL, dfp = "downloads", verbose = FALSE){
  dbpath <- FALSE
  if(!is.null(name)){
    dbpath <- paste(dfp, name, sep = "/")
    if(!file.exists(dbpath)){
      message("Dataset file not found at dfp.")
    } else{
      message("Dataset file found at dfp.")
    }
  } else{
    message("Downloading database...")
    dbpath <- try(get_rmdl(which.class = "test", dfp = dfp, which.type = "h5", 
                           tryload = FALSE, verbose = verbose))
    if(!is(dbpath) == "try-errror"){
      message("Download completed.")
    } else{stop("Problem with download.")}
  }
  # parse load
  if(is(dbpath) == "try-error"){stop("Problem with dbpath.")} else{
    message("Loading database file.")
    dbf <- try(suppressMessages(rhdf5::h5ls(dbpath)))
  }
  if(is(dbf) == "try-error"){stop("Problem loading file.")} else{
    message("Database file loaded.")
    return(dbpath)
  }
  return(NULL)
}
#' @rdname getdb
#' @export
getdb_h5se_gr <- function(name = NULL, dfp = "downloads", verbose = FALSE){
  dbpath <- FALSE
  if(!is.null(name)){
    dbpath <- paste(dfp, name, sep = "/")
    if(!file.exists(dbpath)){
      message("Dataset file not found at dfp.")
    } else{
      message("Dataset file found at dfp.")
    }
  } else{
    message("Downloading database...")
    dbpath <- try(get_rmdl(which.class = "gr", dfp = dfp, which.type = "h5se", 
                           tryload = FALSE, verbose = verbose))
    if(!is(dbpath) == "try-errror"){
      message("Download completed.")
    } else{stop("Problem with download.")}
  }
  # parse load
  if(is(dbpath) == "try-error"){stop("Problem with dbpath.")} else{
    message("Loading database file.")
    dbf <- try(HDF5Array::loadHDF5SummarizedExperiment(dbpath))
  }
  if(is(dbf) == "try-error"){stop("Problem loading file.")} else{
    message("Database file loaded.")
    return(dbf)
  }
  return(NULL)
}
#' @rdname getdb
#' @export
getdb_h5se_gm <- function(name = NULL, dfp = "downloads", verbose = FALSE){
  dbpath <- FALSE
  if(!is.null(name)){
    dbpath <- paste(dfp, name, sep = "/")
    if(!file.exists(dbpath)){
      message("Dataset file not found at dfp.")
    } else{
      message("Dataset file found at dfp.")
    }
  } else{
    message("Downloading database...")
    dbpath <- try(get_rmdl(which.class = "gm", dfp = dfp, which.type = "h5se", 
                           tryload = FALSE, verbose = verbose))
    if(!is(dbpath) == "try-errror"){
      message("Download completed.")
    } else{stop("Problem with download.")}
  }
  # parse load
  if(is(dbpath) == "try-error"){stop("Problem with dbpath.")} else{
    message("Loading database file.")
    dbf <- try(HDF5Array::loadHDF5SummarizedExperiment(dbpath))
  }
  if(is(dbf) == "try-error"){stop("Problem loading file.")} else{
    message("Database file loaded.")
    return(dbf)
  }
  return(NULL)
}
#' @rdname getdb
#' @export
getdb_h5se_rg <- function(name = NULL, dfp = "downloads", verbose = FALSE){
  dbpath <- FALSE
  if(!is.null(name)){
    dbpath <- paste(dfp, name, sep = "/")
    if(!file.exists(dbpath)){
      message("Dataset file not found at dfp.")
    } else{
      message("Dataset file found at dfp.")
    }
  } else{
    message("Downloading database...")
    dbpath <- try(get_rmdl(which.class = "rg", dfp = dfp, which.type = "h5se", 
                           tryload = FALSE, verbose = verbose))
    if(!is(dbpath) == "try-errror"){
      message("Download completed.")
    } else{stop("Problem with download.")}
  }
  # parse load
  if(is(dbpath) == "try-error"){stop("Problem with dbpath.")} else{
    message("Loading database file.")
    dbf <- try(HDF5Array::loadHDF5SummarizedExperiment(dbpath))
  }
  if(is(dbf) == "try-error"){stop("Problem loading file.")} else{
    message("Database file loaded.")
    return(dbf)
  }
  return(NULL)
}
#' @rdname getdb
#' @export
getdb_h5_rg <- function(name = NULL, dfp = "downloads", verbose = FALSE){
  dbpath <- FALSE
  if(!is.null(name)){
    dbpath <- paste(dfp, name, sep = "/")
    if(!file.exists(dbpath)){
      message("Dataset file not found at dfp.")
    } else{
      message("Dataset file found at dfp.")
    }
  } else{
    message("Downloading database...")
    dbpath <- try(get_rmdl(which.class = "rg", dfp = dfp, which.type = "h5", 
                           tryload = FALSE, verbose = verbose))
    if(!is(dbpath) == "try-errror"){
      message("Download completed.")
    } else{stop("Problem with download.")}
  }
  # parse load
  if(is(dbpath) == "try-error"){stop("Problem with dbpath.")} else{
    message("Loading database file.")
    dbf <- try(suppressMessages(rhdf5::h5ls(dbpath)))
  }
  if(is(dbf) == "try-error"){stop("Problem loading file.")} else{
    message("Database file loaded.")
    return(dbpath)
  }
  return(NULL)
}


