get_rmdl <- function(which.file = c("test", "gr", "gm", "rg", "hdf5"),
                     url = "https://recount.bio/data/"){
  # use RCurl to download latest data files
  # which.file  --  type of file to download
  # url -- https://recount.bio/data
  
  fn = RCurl::getURL(url, ftp.use.epsv = FALSE, 
                     dirlistonly = TRUE)
  fn <- unlist(strsplit(fn, "\n"))
  
  catch.str <- ifelse(which.file == "test", "")
  if(which.file == "test"){
    catch.str <- ""
  } 
  
  idat.str <- paste0("\\.idat\\.", ext)
  idat.catch <- grepl(idat.str, fn)
  fn <- unlist(fn)[idat.catch]
  
}