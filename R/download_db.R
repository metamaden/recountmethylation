#!/usr/bin/env R

# Functions for downloading DNAm datasets/cross-study compilations from 
# the server.

#' servermatrix
#'
#' Called by get_rmdl() to get a matrix of database files and file info from 
#' the server. Verifies valid versions and timestamps in filenames, and that 
#' h5se directories contain both an assays and an se.rds file.
#' 
#' @param dn Server data returned from RCurl (default NULL).
#' @param sslver Whether to use SSL certificate authentication for server 
#' connection (default FALSE).
#' @param printmatrix Whether to print the data matrix to console (default 
#' TRUE).
#' @param url Server website url (default "https://methylation.recount.bio/").
#' @param verbose Whether to show verbose messages (default FALSE).
#' @returns Matrix of server files and file metadata
#' @examples 
#' dn <- "remethdb-h5se_gr-test_0-0-1_1590090412 29-May-2020 07:28 -"
#' sm <- try(servermatrix(dn))
#' @seealso get_rmdl, smfilt
#' @export
servermatrix <- function(dn = NULL, sslver = FALSE, printmatrix = TRUE, 
                         url = "https://methylation.recount.bio/", verbose = FALSE){
  if(verbose){message("Getting server data...")}
  if(is.null(dn)){dn <- RCurl::getURL(url,ftp.use.epsv=TRUE,dirlistonly=TRUE,
                                      .opts = list(ssl.verifypeer = sslver))}
  dt<-unlist(strsplit(dn,"\r\n"));dt <- gsub('(.*\">|/</a>|</a>)', "", dt)
  dt <- dt[grepl("remethdb", dt)]
  sm <- do.call(rbind, lapply(as.list(dt), function(x){
    return(unlist(strsplit(gsub("[ ]+",";",x),";")))
  }));colnames(sm) <- c("filename", "date", "time", "size (bytes)")
  sv <- c(); fnv <- sm[grepl("h5se", sm[,1]), 1];fnexclude <- c()
  for(f in fnv){fniv <- c()
  fv <- RCurl::getURL(paste0(url, f, "/"), dirlistonly = TRUE, 
                      .opts = list(ssl.verifypeer = sslver))
  fvv<-unlist(strsplit(fv,"\r\n"));which.start<-which(grepl("Index",fvv))[2]+1
  which.end <- which(grepl("/pre", fvv)) - 1;fvf <- fvv[which.start:which.end]
  for(fni in fvf){name <- gsub('.*\">', '', gsub("</a>.*", "", fni))
  size<-gsub(".* ","",fni);fniv<-c(fniv,paste0("`",name,"`"," = ",size))}
  cond.assays <- length(fniv[grepl("assays", fniv)]) == 1
  cond.se<-length(fniv[grepl("se",fniv)])==1;sv<-c(sv,paste(fniv,collapse=";"))
  if(!(cond.assays & cond.se)){fnexclude <- c(fnexclude, f)}}
  sm[grepl("h5se",sm[,1]),4]<-sv;sm<-sm[!sm[,1] %in% fnexclude,];return(sm)
}

#' smfilt
#' 
#' Filters the data matrix returned from servermatrix().
#' 
#' @param sm Data matrix returned from servermatrix().
#' @param typesdf Data.frame containing database file info for dm filters.
#' @returns Filtered data matrix of server file info.
#' @examples 
#' dm <- matrix(c("remethdb_h5-rg_epic_0-0-2_1589820348.h5","08-Jan-2021",
#' "09:46","66751358297"), nrow = 1)
#' smfilt(dm)
#' @seealso get_rmdl, servermatrix
#' @export
smfilt <- function(sm, typesdf = NULL){
  if(is.null(typesdf)){typesdf <- data.frame(platform = c(rep("hm450k", 4), 
                                                          rep("epic", 4)),
                                             dbtype=rep(c(paste0("h5se-",c("rg","gr","gm")), 
                                                          "h5-rg"), 2), 
                                             stringsAsFactors = FALSE)};
  smf <- sm;smff <- matrix(nrow = 0, ncol = ncol(smf))
  test.files <- which(grepl("test", smf[,1]))
  if(length(test.files) > 0){smff<-smf[test.files,]}
  for(r in seq(nrow(typesdf))){
    tr <- typesdf[r,,drop = FALSE]
    which.db <- which(grepl(tr$platform, sm[,1]) & grepl(tr$dbtype, sm[,1]))
    db.select <- as.character(sm[which.db, 1])
    if(length(db.select) > 1){
      tsv <- as.numeric(gsub(".*_|\\.h5", "", db.select))
      max.ts <- which(tsv == max(tsv));db.select <- db.select[max.ts][1]}
    smff <- rbind(smf[smf[,1] == db.select,], smff)}
  colnames(smff) <- colnames(smf);return(smff)
}

#' Get DNAm assay data.
#'
#' Uses RCurl to download the latest HDF5-SummarizedExperiment or HDF5 
#' database compilation files objects from the server. Calls servermatrix 
#' and performs various quality checks to validate files and downloads. 
#' This function is wrapped in the getdb() set of functions (type `?getdb` for 
#' details).
#' 
#' @param which.class  Either "rg", "gm", "gr", or "test" for RGChannelSet, 
#' MethylSet, GenomicRatioSet, or 2-sample subset.
#' @param which.type Either "h5se" for an HDF5-SummarizedExperiment or 
#' "h5" for an HDF5 database.
#' @param which.platform Supported DNAm array platform type. Currently 
#' supports either "epic" for EPIC/HM850K, or "hm450k" for HM450K.
#' @param fn Name of file on server to download (optional, default NULL).
#' @param dfp Download destination directory (default "downloads").
#' @param url The server URL to locate files for download.
#' @param show.files Whether to print server file data to console (default 
#' FALSE).
#' @param download Whether to download (TRUE) or return queried filename 
#' (FALSE).
#' @param sslver Whether to use server certificate check (default FALSE).
#' @param verbose Whether to return verbose messages (default TRUE).
#' @returns New filepath to dir containing the downloaded data.
#' @examples 
#' # prints file info from server:
#' path <- try(get_rmdl(which.class = "test", which.type = "h5se", 
#' show.files = TRUE, download = FALSE))
#' @seealso servermatrix(), getURL(), loadHDF5SummarizedExperiment(), h5ls()
#' @export
get_rmdl <- function(which.class = c("rg", "gm", "gr", "test"),
                     which.type = c("h5se", "h5"),
                     which.platform = c("hm450k", "epic"), fn = NULL, 
                     dfp="downloads", url = "https://methylation.recount.bio/", 
                     show.files = FALSE, download = TRUE, sslver = FALSE, 
                     verbose = TRUE){
  if(verbose){message("Retrieving data dirnames from server...")}
  sm <- servermatrix(dn = NULL);smf <- smfilt(sm)
  if(show.files){message("Printing server matrix: ");print(smf)}
  if(is.null(fn)){ # clean query results
    str1 <- ifelse(which.type == "h5", "\\.", ".*")
    str2 <- ifelse(which.type == "h5", "$", ".*")
    filt.type <- grepl(paste0(str1, which.type, str2), smf[,1])
    filt.class <- grepl(paste0(".*", which.class,".*"), smf[,1])
    which.fn<-which(filt.type&filt.class);dnc<-smf[which.fn, 1]
    if(!which.class == "test"){
      dnc <- dnc[grepl(which.platform, dnc) & !grepl("test", dnc)]}
    if(length(dnc) > 1){
      tsv <- suppressWarnings(as.numeric(gsub("(.*_|\\.h5)", "", dnc)))
      tsv <- tsv[!is.na(tsv)];dnc <- dnc[which(tsv == max(tsv))[1]]
    };if(length(dnc) == 0){stop("Error, no files of class and type found.")}
  } else{condpass <- grepl("(\\.h5$|.*h5se.*)", fn) & fn %in% smf[,1]
  if(!condpass){stop("Error, provided fn not found on server.")}}
  if(!download){return(dnc)}
  dct1 <- ifelse(!dir.exists(dfp) & !dfp == "", try(dir.create(dfp)), TRUE)
  dfp.dn <- paste(dfp, dnc, sep = "/") # download loc
  if(file.exists(dfp.dn)){
    ostr<-paste0("Ok to overwrite existing file:\n",dfp.dn,"?\n(yes/no)")
    opt<-readline(ostr);if(!opt%in%c("yes","no")){
      stop("Error, unsupported input")}
    if(opt == "no"){stop("Error, stopping download...")}}
  if(which.type == "h5"){dct2 <- try(file.create(dfp.dn))} else{
    dct2 <- ifelse(!dir.exists(dfp.dn), try(dir.create(dfp.dn)), TRUE)}
  if(!(dct1 & dct2)){stop("Error, problem handling download destination.")}
  dn.url <- paste0(url, dnc);if(which.type=="h5"){fl.clean<-""} else{
    fl.clean<-c("assays.h5","se.rds")};dll <- list()
  for(fi in fl.clean){
    fpath <- ifelse(fi == "", dn.url, paste(dn.url, fi, sep = "/"))
    destpath <- ifelse(fi == "", dfp.dn, paste(dfp.dn, fi, sep="/"))
    trydl <- try(utils::download.file(url = fpath, destfile = destpath,
                                      method = "curl", 
                                      .opts = list(ssl.verifypeer = sslver)))}
  if(is(trydl)[1] == "try-error" | length(dll[dll==0]) < length(dll)){
    message("Download incomplete for ", fl.clean[which(dll!=0)])} else{
      dfp.dn <- gsub("\\\\", "/", dfp.dn)
      return(dfp.dn)};return(NULL)
}

#' @name getdb
#' @rdname getdb
#'
#' @title Access database files.
#'
#' @description Combines download and load functions for databases. 
#' If the "namematch" argument isn't provided, the latest available 
#' file is downloaded. All files include metadata for the available samples.
#' 
#' There are 6 functions. Functions with "h5se" access 
#' HDF5-SummarizedExperiment files, and "h5" functions access HDF5 databases. 
#' The 4 h5se functions are "rg" (RGChannelSet), "gm" (MethylSet), "gr" 
#' (GenomicRatioSet), and "test" (data for 2 samples from "gr"). The 2 h5 
#' functions are "rg" (red and green signal datasets), and "test" (data for 2 
#' samples from "rg"). See vignette for details about file types and classes. 
#' 
#' @param platform Valid supported DNAm array platform type. Currently either
#' "epic" for EPIC/HM850K, or "hm450k" for HM450K.
#' @param namematch Filename pattern to match when searching for database 
#' (see defaults).
#' @param dfp Folder to search for database file 
#' (optional, if NULL then searches cache dir specified by BiocFileCache).
#' @param verbose Whether to return verbose messages (default FALSE).
#' @seealso get_rmdl()
#' @returns Either a SummarizedExperiment object for h5se functions, or a file 
#' path for h5 functions.
NULL
#' @rdname getdb
#' @examples
#' \donttest{
#' h5 <- getdb_h5_test(dfp = tempdir())
#' }
#' @export
getdb_h5se_test <- function(platform = NULL, dfp = NULL,
                            namematch = "remethdb-h5se_gr-test.*",
                            verbose = FALSE){
  download<-FALSE;if(is.null(dfp)){dfp<-BiocFileCache::BiocFileCache()@cache}
  clf <- list.files(dfp);fmatch <- clf[grepl(namematch, clf)]
  if(!is.null(namematch) & length(fmatch) > 0){
    fn1 <- fmatch[1];fpath <- gsub("\\\\", "/", file.path(dfp, fn1))
    ostr <- paste0("Use file:\n", fpath, "?\n(yes/no)");opt <- readline(ostr)
    if(!opt %in% c("yes", "no")){stop("Error, unsupported input.")}
    if(opt == "no"){download <- TRUE}} else{download <- TRUE}
  if(download){
    message("Downloading database...")
    dbpath <- try(get_rmdl(which.class = "test", dfp = dfp, 
                           which.type = "h5se", verbose = verbose))
    if(!is(dbpath)[1] == "try-error"){message("Download completed.")} else{
      stop("Error, problem with download.")}
  } else{dbpath <- fpath}
  if(is(dbpath)[1] == "try-error"){stop("Error, problem with dbpath.")} else{
    message("Loading database...");dbf <- try(
      HDF5Array::loadHDF5SummarizedExperiment(dbpath))
    if(is(dbf)[1] == "try-error"){stop("Error, problem loading file.")} else{
      message("Database file loaded.");return(dbf)}};return(NULL)
}
#' @rdname getdb
#' @export
getdb_h5_test <- function(platform = NULL, namematch = "remethdb-h5_rg-test_.*", 
                          dfp = NULL, verbose = FALSE){
  download<-FALSE;if(is.null(dfp)){dfp<-BiocFileCache::BiocFileCache()@cache}
  clf <- list.files(dfp);fmatch <- clf[grepl(namematch, clf)]
  if(!is.null(namematch) & length(fmatch) > 0){
    fn1 <- fmatch[1];fpath <- gsub("\\\\", "/", file.path(dfp, fn1))
    ostr <- paste0("Use file:\n", fpath, "?\n(yes/no)");opt <- readline(ostr)
    if(!opt %in% c("yes", "no")){stop("Error, unsupported input")}
    if(opt == "no"){download <- TRUE}
  } else{download <- TRUE}
  if(download){
    message("Downloading database...")
    dbpath <- try(get_rmdl(which.class = "test", dfp = dfp, which.type = "h5", 
                           verbose = verbose));if(!is(dbpath)[1]=="try-error"){
                             message("Download completed.")} else{
                               stop("Error, problem with download.")}
  } else{dbpath <- fpath}
  if(is(dbpath)[1] == "try-error"){stop("Error, problem with dbpath.")} else{
    message("Loading database file.")
    dbf <- try(suppressMessages(rhdf5::h5ls(dbpath)))
    if(is(dbf)[1] == "try-error"){stop("Error, problem loading file.")} else{
      message("Database file loaded.");return(dbpath)}};return(NULL)
}
#' @rdname getdb
#' @export
getdb_h5se_gr <- function(platform = c("hm450k", "epic"), dfp = NULL,
                          namematch="remethdb_h5se-gr_.*", verbose=FALSE){
  download<-FALSE;if(is.null(dfp)){dfp<-BiocFileCache::BiocFileCache()@cache}
  clf <- list.files(dfp);fmatch <- clf[grepl(namematch, clf)]
  if(!is.null(namematch) & length(fmatch) > 0){
    fn1 <- fmatch[1];fpath <- gsub("\\\\", "/", file.path(dfp, fn1))
    ostr <- paste0("Use file:\n", fpath, "?\n(yes/no)");opt <- readline(ostr)
    if(!opt %in% c("yes", "no")){stop("Error, unsupported input")}
    if(opt == "no"){download <- TRUE}} else{download <- TRUE}
  if(download){
    message(paste0("Downloading ",platform," database..."))
    dbpath <- try(get_rmdl(which.class = "gr", dfp = dfp, which.type = "h5se",
                           which.platform = platform, verbose = verbose))
    if(!is(dbpath)[1] == "try-error"){message("Download completed.")} else{
      stop("Error, problem with download.")}
  } else{dbpath <- fpath}
  if(is(dbpath)[1] == "try-error"){stop("Error, problem with dbpath.")} else{
    message("Loading database file.")
    dbf <- try(HDF5Array::loadHDF5SummarizedExperiment(dbpath))
    if(is(dbf)[1] == "try-error"){stop("Error, problem loading file.")} else{
      message("Database file loaded.");return(dbf)}};return(NULL)
}
#' @rdname getdb
#' @export
getdb_h5se_gm <- function(platform = c("hm450k", "epic"), dfp = NULL, 
                          namematch = "remethdb_h5se-gm_.*", verbose = FALSE){
  download<-FALSE;if(is.null(dfp)){dfp<-BiocFileCache::BiocFileCache()@cache}
  clf <- list.files(dfp);fmatch <- clf[grepl(namematch, clf)]
  if(!is.null(namematch) & length(fmatch) > 0){
    fn1 <- fmatch[1];fpath <- gsub("\\\\", "/", file.path(dfp, fn1))
    ostr <- paste0("Use file:\n", fpath, "?\n(yes/no)");opt<-readline(ostr)
    if(!opt %in% c("yes", "no")){stop("Error, unsupported input")}
    if(opt == "no"){download <- TRUE}
  } else{download <- TRUE}
  if(download){
    message("Downloading database...")
    dbpath <- try(get_rmdl(which.class = "gm", dfp = dfp, which.type = "h5se",
                           which.platform = platform, verbose = verbose))
    if(!is(dbpath)[1] == "try-error"){message("Download completed.")} else{
      stop("Error, problem with download.")}
  } else{dbpath <- fpath}
  if(is(dbpath)[1] == "try-error"){stop("Error, problem with dbpath.")} else{
    message("Loading database file.")
    dbf <- try(HDF5Array::loadHDF5SummarizedExperiment(dbpath))
    if(is(dbf)[1] == "try-error"){stop("Error, problem loading file.")} else{
      message("Database file loaded.");return(dbf)}};return(NULL)
}
#' @rdname getdb
#' @export
getdb_h5se_rg <- function(platform = c("hm450k", "epic"), dfp = NULL, 
                          namematch = "remethdb-h5se_rg_.*", verbose = FALSE){
  download<-FALSE;if(is.null(dfp)){dfp<-BiocFileCache::BiocFileCache()@cache}
  clf <- list.files(dfp);fmatch <- clf[grepl(namematch, clf)]
  if(!is.null(namematch) & length(fmatch) > 0){
    fn1 <- fmatch[1];fpath <- gsub("\\\\", "/", file.path(dfp, fn1))
    ostr <- paste0("Use file:\n", fpath, "?\n(yes/no)");opt <- readline(ostr)
    if(!opt %in% c("yes", "no")){stop("Error, unsupported input")}
    if(opt == "no"){download <- TRUE}
  } else{download <- TRUE}
  if(download){
    message("Downloading database...")
    dbpath <- try(get_rmdl(which.class = "rg", dfp = dfp, which.type = "h5se",
                           which.platform = platform, verbose = verbose))
    if(!is(dbpath)[1] == "try-error"){message("Download completed.")} else{
      stop("Error, problem with download.")}
  } else{dbpath <- fpath}
  if(is(dbpath)[1] == "try-error"){stop("Error, problem with dbpath.")} else{
    message("Loading database file.")
    dbf <- try(HDF5Array::loadHDF5SummarizedExperiment(dbpath))
    if(is(dbf)[1] == "try-error"){stop("Error, problem loading file.")} else{
      message("Database file loaded.");return(dbf)}};return(NULL)
}
#' @rdname getdb
#' @export
getdb_h5_rg <- function(platform = c("hm450k", "epic"), dfp = NULL, 
                        namematch = "remethdb-h5_rg_.*", verbose = FALSE){
  download<-FALSE;if(is.null(dfp)){dfp<-BiocFileCache::BiocFileCache()@cache}
  clf <- list.files(dfp);fmatch <- clf[grepl(namematch, clf)]
  if(!is.null(namematch) & length(fmatch) > 0){
    fn1 <- fmatch[1];fpath <- gsub("\\\\", "/", file.path(dfp, fn1))
    ostr <- paste0("Use file:\n", fpath, "?\n(yes/no)");opt <- readline(ostr)
    if(!opt %in% c("yes", "no")){stop("Error, unsupported input")}
    if(opt == "no"){download <- TRUE}} else{download <- TRUE}
  if(download){
    message("Downloading database...")
    dbpath <- try(get_rmdl(which.class = "rg", dfp = dfp, which.type = "h5", 
               which.platform = platform, verbose = verbose))
    if(!is(dbpath)[1] == "try-error"){message("Download completed.")} else{
      stop("Error, problem with download.")}
  } else{dbpath <- fpath}
  if(is(dbpath)[1] == "try-error"){stop("Error, problem with dbpath.")} else{
    message("Loading database file.")
    dbf <- try(suppressMessages(rhdf5::h5ls(dbpath)))
    if(is(dbf)[1] == "try-error"){stop("Error, problem loading file.")} else{
      message("Database file loaded.");return(dbpath)}};return(NULL)
}


