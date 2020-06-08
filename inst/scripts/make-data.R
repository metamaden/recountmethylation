#!/usr/bin/env R

require(rmpipeline)
require(recountmethylation)

# Steps to make recountmethylation databases
#
# This script outlines how database files were made. For details, 
# see the manuscript Maden et al 2020.

#----------
# GEO files
#----------
# This section outlines how files were obtained and extracted using the 
# recount-methylation-server software. Files were obtained from GEO DataSets 
# in March 2019. GSM/sample IDATs and GSE/study SOFT files were obtained, 
# versioned, and managed programmatically using Python in a CentOS 7 remote 
# server environment. For details, consult the recount-methylation-server repo 
# and ReadMe (https://github.com/metamaden/recount-methylation-server).

# * Install recount-methylation-server and dependencies
# * From command line or a new attached screen, run:
# > python3 ./recount-methylation-server/src/server.py
# * From a Python3 session, run the following:
# > import os
# > import sys
# > sys.path.insert(0, os.path.join("recount-methylation-server","src"))
# > from process_soft import expand_soft, extract_gsm_soft, gsm_soft2json, msrap_screens
# > expand_soft()
# > extract_gsm_soft()
# > gsm_soft2json()
# * map GSM JSON-formatted metadata with MetaSRA-pipeline as follows:
# > msrap_screens(nscreensi=200,nmaxscreens=30,qcprint=True)

#----------------
# Sample metadata
#----------------
# Metadata for GSMs/samples was derived from SOFT data and predicted from IDAT
# red/grn signals (Maden et al 2020). 

# Metadata were preprocessed by manual annotation, then postprocessed to format 
# term labels and add additional information. Final or "postprocessed" metadata 
# includes 3 principal data types: 
# * GEO record information (e.g. sample title, GSM ID, GSE ID, etc.)
# * Learned or mapped labels (e.g. tissue, disease, sample type, etc.)
# * Predicted labels (e.g. age, sex, and cell fraction estimates) from red/grn signals 
#
# The variables disease, tissue, age, sex, and storage, were preprocessed by manual 
# annotation, then postprocessed as detailed in the script script md_postprocessing.R. 
# Sample type predictions were generated using MetaSRA-pipeline with GSM JSON metadata files 
# (see above). The functions wateRmelon::agep, minfi::getSex, and minfi::estimateCellCounts
# were used to predict age, sex, and blood cell type estimates, respectively.

#-----------------
# Database objects
#-----------------
# Database objects were created from sample IDATs and metadata. After 
# files were obtained and postprocessed metadata was generated, a pipeline
# was run to obtain paired DNAm assays and metadata as either HDF5 "h5" 
# or HDF5-SummarizedExperiment "h5se" files. Code to generate these files is 
# in the script inst/scripts/runpipeline.R for the "rmpipeline" package
# (https://github.com/metamaden/rmpipeline). This code is reproduced below.

# datasets metadata for pipeline run
newversion <- "0.0.1" # set new version
md <- get_metadata(title = "newrun", version = newversion)
versionfn <- md[["version"]]
timestamp <- md[["timestamp"]]
# red/grn signal data
# navigate to main recount-methylation dir/base dir.
# e.g. 
# > cd recount-methylation
dtables_rg(versionfn, timestamp)
# make the h5 file
# navigate to compilations dir
# e.g. 
# > cd recount-methylation/recount-methylation-analysis/files/mdata/compilations
makeh5db_rg(dbfnstem = "remethdb", version = versionfn, ts = timestamp, 
            mdpath = "mdpost_all-gsm-md.rda", fnpath = ".",
            fnl = c("redsignal_1589820348_0-0-1.mdat.compilation",
                    "greensignal_1589820348_0-0-1.mdat.compilation"))
# make the h5se file
make_h5se(dbn = dbn, newfnstem = fnstem, version = versionfn, ts = timestamp)
# meth/unmeth and betavals data
# make new h5 files
h5name.gm <- make_h5_gm()
h5name.gr <- make_h5_gr()
# make new h5se objects
make_h5se_gm(h5name.gm)
make_h5se_gr(h5name.gr)
