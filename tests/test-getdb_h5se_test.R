testthat::context("Test HDF5-SummarizedExperiment download")

testthat::test_that("Database download and object properties", {

    dldn <- "temp"
    h5se.dat <- getdb_h5se_test(dfp = dldn)
    testthat::expect_that(h5se.dat, testthat::is_a("GenomicRatioSet"))
    testthat::expect_that(minfi::pData(h5se.dat), 
        testthat::is_a("DFrame"))
    testthat::expect_that(minfi::getAnnotation(h5se.dat), 
        testthat::is_a("DFrame"))
    testthat::expect_equal(nrow(minfi::pData(h5se.dat)), ncol(h5se.dat))
    testthat::expect_equal(nrow(minfi::getAnnotation(h5se.dat)), 
        nrow(h5se.dat))
    unlink(dldn, recursive = TRUE)

})

