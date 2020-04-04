

get_rmdl <- function(which.dn = c("h5se-test_gr", "h5se_gr", "h5se_gm", "h5se_rg", "\\.h5"),
                     url = "https://recount.bio/data/", dfp = "data", verbose = TRUE){
  # use RCurl to download latest data files
  # which.file  --  type of data dirname to be downloaded
  # dfp -- dl dest filepath
  # url -- https://recount.bio/data
  # dfp -- Data file path
  # verbose -- whether to return verbose messages
  if(verbose){message("Retrieving data dirnames from server...")}
  dn <- RCurl::getURL(url, ftp.use.epsv = FALSE, 
                     dirlistonly = TRUE)
  dn <- unlist(strsplit(dn, "\n"))
  catch.str <- paste0(".*", which.dn,".*")
  dn.catch <- grepl(catch.str, dn)
  dn <- unlist(dn)[dn.catch]
  dn.clean <- gsub('<.*', "", gsub('.*">', "", dn))
  if(!length(dn.clean) == 1){
    stop("There was a problem parsing the file string.")
    }
  # recursively get filenames in target dir
  if(verbose){message("Retrieving filenames from server...")}
  dn.url <- paste0(url, dn.clean)
  fl = RCurl::getURL(dn.url, ftp.use.epsv = FALSE,
                       dirlistonly = TRUE)
  fl <- unlist(strsplit(fl, "\n"))
  fl.str <- paste0(c("assays.h5", "se.rds"), collapse = "|")
  fl.catch.str <- paste0(".*", fl.str,".*")
  fl.catch <- grepl(fl.catch.str, fl)
  fl <- unlist(fl)[fl.catch]
  fl.clean <- gsub('<.*', "", gsub('.*">', "", fl))
  # check dest dir
  if(!dir.exists(dfp)){
    dct <- try(dir.create(dfp))
    if(!dft){stop("Attempt to make new dir ", dfp,
                  "failed with error: ", dct[2])}s
  }
  # make data object dir
  new.data.dn <- gsub("/", "", dn.clean)
  path.data.dn <- c(dfp, new.data.dn)
  path.data.dn <- paste(path.data.dn, collapse = "/")
  # dl files to data dir
  if(verbose){message("Downloading files...")}
  dll <- list()
  for(i in 1:length(fl.clean)){
    f <- fl.clean[i]
    fpath <- paste0(c(dn.url, f), collapse = "")
    dll[[i]] <- try(download.file(fpath, path.data.dn))
  }
  if(dll[[1]] & dll[[2]]){
    if(verbose){message("Download completed successfully.")}
    return(TRUE)
  } else{
    if(verbose){message("Download incomplete. ",
                        "Returning with following outcomes: ",
                        dll[[1]], " ", dll[2])}
    return(FALSE)
  }
  return(NULL)
}

# unit tests
# test that catch.str match returns exactly 1 file
# test that dll[[1]] & dll[[2]] both TRUE 

















get_rmdl_curl <- function(which.file = c("h5se-test_gr", "h5se_gr", "h5se_gm", "h5se_rg", "\\.h5"),
                          url = "https://recount.bio/data/", dfp = "."){
  # use curl to download latest data files
  # which.file  --  type of file to download
  # dfp -- dl dest filepath
  # url -- https://recount.bio/data
  fn = RCurl::getURL(url, ftp.use.epsv = FALSE, 
                     dirlistonly = TRUE)
  fn <- unlist(strsplit(fn, "\n"))
  catch.str <- paste0(".*", which.file,".*")
  fn.catch <- grepl(catch.str, fn)
  fn <- unlist(fn)[fn.catch]
  fn.clean <- gsub('<.*', "", gsub('.*">', "", fn))
  if(!length(fn.clean) == 1){stop("There was a problem parsing the file string.")}
  
  # download
  download.file(paste0(url, fn.clean), dfp)
  
}
