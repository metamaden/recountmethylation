testthat::context("Database file information from get_rmdl")

testthat::test_that("File info as character string", {

    fileinfo <- recountmethylation::get_rmdl(download = FALSE)
    testthat::expect_that(fileinfo, testthat::is_a("character"))
    
})