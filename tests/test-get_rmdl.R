testthat::context("Database file information from get_rmdl")

testthat::test_that("File info as character string", {

    fileinfo <- get_rmdl(download = FALSE)
    testthat::expect_true(is.character(fileinfo))
    
})