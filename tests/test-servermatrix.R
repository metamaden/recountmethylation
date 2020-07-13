testthat::context("Servermatrix properties")

testthat::test_that("Properties of server matrix returned from get_rmdl 
    and servermatrix", {
    url = "https://recount.bio/data/"
    show.files = FALSE
    sslver = FALSE
    ftpuseopt <- dirlistopt <- ifelse(show.files, FALSE, TRUE) # rcurl setup
    dn <- RCurl::getURL(url, ftp.use.epsv = ftpuseopt, dirlistonly = dirlistopt,
                      .opts = list(ssl.verifypeer = sslver))
    sm <- recountmethylation::servermatrix(dn = dn, sslver = sslver)
    testthat::expect_true(is.matrix(sm))
    testthat::expect_equal(ncol(sm), 4)
    testthat::expect_true(is.character(sm[,1]))
})