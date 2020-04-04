

get_rmdl <- function(which.dn = c("h5se-test_gr", "h5se_gr", "h5se_gm", "h5se_rg", "\\.h5"),
                     url = "https://recount.bio/data/", dfp = "data"){
  # use RCurl to download latest data files
  # which.file  --  type of data dir for dl
  # dfp -- dl dest filepath
  # url -- https://recount.bio/data
  dn <- RCurl::getURL(url, ftp.use.epsv = FALSE, 
                     dirlistonly = TRUE)
  dn <- unlist(strsplit(dn, "\n"))
  catch.str <- paste0(".*", which.dn,".*")
  dn.catch <- grepl(catch.str, dn)
  dn <- unlist(dn)[dn.catch]
  dn.clean <- gsub('<.*', "", gsub('.*">', "", dn))
  if(!length(dn.clean) == 1){stop("There was a problem parsing the file string.")}
  
  # recursively get filenames in target dir
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
  
  # get files from dir
  for(f in fl.clean){
    # download
    download.file(paste0(url, fn.clean), dfp)
  }
  
  
  
}

# unit tests
# test that catch.str match returns exactly 1 file


















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
