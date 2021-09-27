testthat::context("Test the server matrix properties")

testthat::test_that("Properties of server matrix returned from get_rmdl 
    and servermatrix", {
    sm <- recountmethylation::get_servermatrix()
    testthat::expect_true(is(sm, "matrix"))
    testthat::expect_equal(ncol(sm), 4)
    testthat::expect_true(is(sm[,1], "character"))
})