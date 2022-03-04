#!/usr/bin/env R

# Author: Sean Maden
#
# Functions to manage search index construction from DNAm array data, including
# dimension reduction with feature hashing and support for k nearest neighbors
# search lookup. The HNSW implementation in the hnswlib Python library is used
# for search index construction.
#

#' setup_sienv
#'
#' Set up a new virtual environment for search index construction using 
#' the basilisk package.
#' 
#' @param env.name Name of the new virtual environment (required, 
#' "dnam_si-hnswlib")
#' @param pkgv Vector of the dependencies and their versions for the new
#' virtual environment (required, format: "packagename==versionnum").
#' @returns New basilisk environment object.
#' @example
setup_sienv <- function(env.name = "dnam_si_hnswlib", 
                        pkgv = c("hnswlib==0.5.1", "pandas==1.2.2", 
                          "numpy==1.20.1", "mmh3==3.0.0", "h5py==3.2.1")){
  message("Defining the virtual env dependencies...")
  my_env <- basilisk::BasiliskEnvironment(envname = env.name, 
                                          pkgname = "recountmethylation", 
                                          packages = pkgv)
  message("Running virtual environment setup...")
  proc <- basilisk::basiliskStart(my_env) # define env process
  on.exit(basilisk::basiliskStop(proc)) # set exit process
  return(proc)
}

#' get_fh
#'
#' Get the hashed features for a data table. Uses reticulate package to
#' call the Python script to do feature hashing on a table of data. It is
#' assumed the input table has sample data in rows, with probe data in 
#' columns. The input data table should have row names but not column names.
#' 
#' @param csv_savepath Name/path of hashed features table to write (required,
#' string, writes new csv where rows = samples, cols = hashed features).
#' @param csv_openpath Name/path of table to hash (required, string, assumes 
#' a csv where rows = samples, cols = probes).
#' @param ndim Number of hashed features (integer, 1000).
#' @param lstart Line index to start on (0-based for Python, required, int, 0).
#' @returns Path to new hashed featuers table.
#' @examples
#' # get example bval csv
#' of_fpath <- system.file("extdata", "fhtest", 
#'          package = "recountmethylation")
#' of_fpath <- file.path(of_fpath, "tbval_test.csv")
#' # write new hashed features results
#' get_fh(csv_savepath = "bval_fn.csv", csv_openpath = of_fpath, ndim = 100)
#' @export
get_fh <- function(csv_savepath, csv_openpath, ndim = 1000, lstart = 1){
  message("Starting basilisk process...")
  proc <- recountmethylation:::setup_sienv()
  pyscript.path <- system.file("python", package = "recountmethylation")
  pyscript.path <- file.path(pyscript.path, "dnam_search_index.py")
  message("Doing feature hashing...")
  ncol <- length(unlist(strsplit(readLines(csv_openpath, 1), ",")))
  message("Target features will reduce dimensions from ", ncol, 
          " to ", ndim, "...")
  basilisk::basiliskRun(proc, function(pyscript.path, csv_savepath, 
                                       csv_openpath, ndim, lstart){
    message("Sourcing Python functions...")
    reticulate::source_python(pyscript.path)
    make_fhmatrix_autolabel(wf_name = csv_savepath, of_name = csv_openpath, 
                            ndim = as.integer(ndim), lstart = lstart)
  }, pyscript.path = pyscript.path, csv_savepath = csv_savepath, 
  csv_openpath = csv_openpath, ndim = ndim, lstart = lstart)
  if(file.exists(csv_savepath)){return(csv_savepath)} else{return(FALSE)}
  return(NULL)
}

#' make_si
#' 
#' Make search index from table of hashed features. Additional details about 
#' the hnswlib search index parameters (e..g `space_val`, `efc_val`, `m_val`, 
#' and `ef_val`) can be found in the Python package docstrings and ReadMe.
#' 
#' @param fh_csv_fpath Name/path of csv (e.g. a table of hashed features) 
#' containing data for the index (required, string, "bvaltest.csv", where 
#' rows = samples, cols = features).
#' @param si_fname Name of new search index file to save (required, string, 
#' "new_search_index.pickle")
#' @param si_dict_fname Name of new index dictionary, with string labels, to 
#' save (required, string, "new_index_dict.pickle").
#' @param threads Number of threads for processing new index (required, int, 4).
#' @param space_val Space value for new search index (required, valid string, 
#' l2').
#' @param efc_val EFC value for the index (required, int, 2000).
#' @param m_val M value for the index (required, int, 1000).
#' @param ef_val EF value for the index (required, int, 2000).
#' @returns Boolean, TRUE if new search index and dictionary created, FALSE if 
#' creating the new search index and dictionary files failed, otherwise NULL.
#' @example 
#' fh_csv_fpath <- system.file("extdata", "fhtest", 
#' package = "recountmethylation")
#' fh_csv_fpath <- file.path(fh_csv_fpath, "bval_fn.csv")
#' make_si(fh_csv_fpath)
#' @export
make_si <- function(fh_csv_fpath, si_fname = "new_search_index.pickle", 
                    si_dict_fname = "new_index_dict.pickle", threads = 4, 
                    space_val = 'l2', efc_val = 2000, m_val = 1000, ef_val = 2000){
  message("Starting basilisk process...")
  proc <- recountmethylation:::setup_sienv()
  pyscript.path <- system.file("python", package = "recountmethylation")
  pyscript.path <- file.path(pyscript.path, "dnam_search_index.py")
  basilisk::basiliskRun(proc, function(pyscript.path, fh_csv_fpath, si_fname, 
                                       si_dict_fname, threads, space_val, 
                                       efc_val, m_val, ef_val){
    message("Sourcing Python functions...")
    reticulate::source_python(pyscript.path)
    message("Making search index...")
    make_hnsw_si(fname = fh_csv_fpath, index_name = si_fname, 
            dindex_name = si_dict_fname, threads = as.integer(threads), 
            efc_val = as.integer(efc_val), m_val = as.integer(m_val), 
            ef_val = as.integer(ef_val))
  }, 
  pyscript.path = pyscript.path, fh_csv_fpath = fh_csv_fpath,
  si_fname = si_fname, si_dict_fname = si_dict_fname, threads = threads,
  space_val = space_val, efc_val = efc_val, m_val = m_val, ef_val = ef_val
  )
  if(file.exists(si_fname) & file.exists(si_dict_fname)){
    message("Made new search index '",si_fname,
            "' and search index dict '", si_dict_fname,"'")
    return(TRUE)
  } else{
    return(FALSE)
  }
  return(NULL)
}

#' query_si
#'
#' Query an HNSW search index. Does K Nearest Neighbors lookup on a previously 
#' saved search index object, returning the K nearest neighbors of the queried 
#' sample(s). The `query_si()` function returns verbose output, which can be 
#' silenced with suppressMessages()`.
#' 
#' @param sample_idv Vector of valid sample IDs, or GSM IDs, which are included
#' in the rownames of the hashed features table at fh_csv_fpath (requried, 
#' vector of char strings).
#' @param fh_csv_fpath Path to the hashed features table, which includes rownames
#' corresponding to sample ID strings in the sample_idv vector (required, char).
#' @param si_fname Base filename of the search index object, used to find the 
#' search index and index dict files, which are expected to be located at
#' si_fapth (required, char).
#' @param si_fpath Path to the directory containing the search index and index
#' dict files (required, char).
#' @param lkval Vector of K nearest neighbors to return per query (optional, 
#' int, c(1,2)).
#' @returns
#' @example
#' # file paths
#' # fh table
#' fh_csv_fname <- system.file("extdata", "fhtest", 
#' package = "recountmethylation")
#' fh_csv_fname <- file.path(fh_csv_fname, "bval_fh10.csv")
#' # si dict
#' index_dict_fname <- system.file("extdata", "sitest", 
#' package = "recountmethylation")
#' index_dict_fname <- file.path(index_dict_fname, "new_index_dict.pickle")
#' 
#' # set sample ids to query
#' sample_idv <- c("GSM1038308.1548799666.hlink.GSM1038308_5958154021_R01C01",
#'                 "GSM1038309.1548799666.hlink.GSM1038309_5958154021_R02C01")
#' # set a list of k nearest neighbors to query
#' lkval <- c(1,2,3)
#' 
#' # get query results as a data frame (with verbose results messaging)
#' dfk <- query_si(sample_idv = sample_idv, lkval = lkval, 
#'                 fh_csv_fname = "bval_fn.csv", 
#'                 index_dict_fname = "new_index_dict.pickle")
#' # returns:
#' # Starting basilisk process...
#' # Defining the virtual env dependencies...
#' # Running virtual environment setup...
#' # Sourcing Python functions...
#' # Querying the search index...
#' # Getting hashed features data for samples...
#' # Getting index data for sample: 
#' # GSM1038308.1548799666.hlink.GSM1038308_5958154021_R01C01'
#' # Getting index data for sample: 
#' # GSM1038309.1548799666.hlink.GSM1038309_5958154021_R02C01'
#' # Beginning queries of k neighbors from lk...
#' # ii =  0 , ki =  1
#' # Loading search index...
#' # Querying 2 elements in data with k = 1 nearest neighbors...
#' # Query completed, time: 0.0007359981536865234
#' # Applying labels to query results...
#' # Returning data (sample id, k index, and distance)...
#' # ii =  1 , ki =  2
#' # Loading search index...
#' # Querying 2 elements in data with k = 2 nearest neighbors...
#' # Query completed, time: 0.0006208419799804688
#' # Applying labels to query results...
#' # Returning data (sample id, k index, and distance)...
#' # ii =  2 , ki =  3
#' # Provided k '3' > n si samples, skipping...
#' # Returning query results...
#' @export
query_si <- function(sample_idv, fh_csv_fpath, 
                     si_fname = "new_search_index", 
                     si_fpath = ".", lkval = c(1,2)){
  message("Checking index and table locations...")
  if(!file.exists(fh_csv_fpath)){
    stop("Error: didn't find fh table at location:\n'",fh_csv_fpath,"'")}
  si_hnsw_fpath <- file.path(si_fpath, paste0(si_fname, ".pickle"))
  if(!file.exists(si_hnsw_fpath)){
    stop("Error: didn't find search index at location:\n'",v,"'")}
  si_dict_fpath <- file.path(si_fpath, paste0(si_fname, "_dict.pickle"))
  if(!file.exists(si_dict_fpath)){
    stop("Error: didn't find index dict at location:\n'",si_dict_fpath,"'")}
  message("Starting basilisk process...")
  proc <- recountmethylation:::setup_sienv()
  lkval <- lapply(lkval, as.integer) # format query k values list
  pyscript.path <- system.file("python", package = "recountmethylation")
  pyscript.path <- file.path(pyscript.path, "dnam_search_index.py")
  query_result <- basilisk::basiliskRun(proc, function(pyscript.path, 
    sample_idv, lkval, si_dict_fpath, fh_csv_fpath){
    message("Sourcing Python functions...")
    reticulate::source_python(pyscript.path)
    message("Querying the search index...")
    dfk <- make_dfk_sampleid(sample_idv = sample_idv, lk = lkval, 
      index_dict_fname = si_dict_fpath, fh_csv_fname = fh_csv_fpath)
    return(dfk)}, pyscript.path = pyscript.path, sample_idv = sample_idv, 
    lkval = lkval, si_dict_fpath = si_dict_fpath, fh_csv_fpath = fh_csv_fpath)
  return(query_result)
}
