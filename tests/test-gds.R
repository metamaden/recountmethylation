testthat::context("GEO validate SummarizedExperiment")

testthat::test_that("Download and object properties", {

    dldn <- "./temp/"
    gsmv <- c("GSM2465267", "GSM2814572")
    gdq <- recountmethylation::gds_idatquery(gsmv, dfp = dldn)
    # ids in query
    idmatch <- grepl(paste(gsmv, collapse = "|"), gdq$basenames)
    testthat::expect_equal(length(idmatch[idmatch]), 2)
    # idat downloads
    gds.se <- recountmethylation::gds_idat2rg(gsmv, dfp = dldn)
    testthat::expect_that(gds.se, testthat::is_a("RGChannelSet"))
    unlink(dldn, recursive = TRUE)
    
})