#!/usr/bin/env R

# Make test h5se object

# filename
# testfn <- "remethdb_h5se-test_gr"
# lf <- "remethdb_h5se_gr_00-00-01_1583780004"
# v.ts <- gsub(".*gr\\_", "", lf) # timestamp & version

# subset data
# testfn <- paste0(c(testfn, v.ts), 
#                  collapse = "_")
# h5se <- HDF5Array::loadHDF5SummarizedExperiment(lf)
# h5se.test.gr <- h5se[,c(1:2)]

# colnames(h5se.test.gr)
# 1] "GSM1038308.1548799666.hlink.GSM1038308_5958154021_R01C01"
# [2] "GSM1038309.1548799666.hlink.GSM1038309_5958154021_R02C01"

# HDF5Array::saveHDF5SummarizedExperiment(h5se.test.gr, testfn)