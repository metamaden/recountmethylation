#!/usr/bin/env R

# Functions for downloading DNAm datasets/cross-study compilations from the server.

#' servermatrix
#'
#' Get a matrix of database files from the recount server. Called by get_rmdl.
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
#' 
#' @param which.class  Class of file to download (either "rg" for `RGChannelSet`, 
#' "gm" for `MethylSet`, "gr" for `GenomicRatioSet`, or "test" for the test dataset).
#' @param which.type Type of file. Either "h5se" for HDF5-SummarizedExperiment (default) or 
#' `h5` for HDF5 database.
#' @param fn Name of server file to download (default NULL).
#' @param dfp Target local directory for downloaded files (default "downloads").
#' @param url Server URL containing assay data.
#' @param show.files Whether to print server file data to console (default FALSE).
#' @param download Whether to download (TRUE) or return queried filename (FALSE).
#' @param tryload Whether to try loading downloaded data (default TRUE).
#' @param verbose Whether to return verbose messages.
#' @param sslver Whether to use server certificate check (default FALSE).
#' @return New filepath to dir with downloaded data.
#' @examples 
#' path <- get_rmdl(which.class = "test", which.type = "h5se", tryload = FALSE)
#' base::unlink("downloads", recursive = TRUE)
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
#' @title Functions to access database files.
#'
#' @description Combines download and load for database files. Latest files are downloaded to dfp directory if name not provided.
#'  `HDF5-SummarizedExperiment` files (rg is `RGChannelSet`, gm is `MethylSet`, gr is `GenomicRatioSet`) 
#'  can be downloaded with `h5se` functions, and `HDF5` database can be downloaded with h5 function. 
#'  Files include sample metadata. See vignette for details about file types and classes.
#' @param name Name of database file to load (optional, default NULL). If null, a new database is download.
#' @param dfp Folder containing database to load if name not NULL, otherwise the target folder for 
#' new database download (default "downloads").
#' @param verbose Whether to return verbose messages (default FALSE).
#' @seealso dldb, get_rmdl
#' @return SummarizedExperiment object (for h5se) or file path for (h5).
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


