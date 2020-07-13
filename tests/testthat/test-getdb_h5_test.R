context("Test HDF5 dataset download")

testthat::test_that("HDF5 database properties", {

    # database download
    dldn <- "temp"
    h5path <- getdb_h5_test(dfp = dldn)
    testthat::expect_true(file.exists(h5path))
    # database properties
    h5dat <- rhdf5::h5ls(h5path)
    testthat::expect_true(is.data.frame(h5dat))
    testthat::expect_equal(nrow(h5dat[h5dat$name == "greensignal",]), 1)
    testthat::expect_equal(nrow(h5dat[h5dat$name == "greensignal",]), 1)
    testthat::expect_equal(nrow(h5dat[h5dat$name == "redsignal",]), 1)
    testthat::expect_equal(nrow(h5dat[h5dat$name == "mdpost",]), 1)
    testthat::expect_equal(nrow(h5dat[h5dat$name == "greensignal.colnames",]), 1)
    testthat::expect_equal(nrow(h5dat[h5dat$name == "greensignal.rownames",]), 1)
    testthat::expect_equal(nrow(h5dat[h5dat$name == "redsignal.colnames",]), 1)
    testthat::expect_equal(nrow(h5dat[h5dat$name == "redsignal.rownames",]), 1)
    testthat::expect_equal(nrow(h5dat[h5dat$name == "mdpost.colnames",]), 1)
    # remove temp dir
    unlink(dldn, recursive = TRUE)

})

