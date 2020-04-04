#!/usr/bin/env R

require("recountmethylation")

# Unit tests for `recountmethylation`

# create tests directory
# utpath <- paste0(c("tests", "testthat"), collapse = "/")
# dir.create("tests/testthat")

gsmvi <- c("GSM2465267", "GSM2814572")

expect_equal(gds_idatquery(gsmvi)[["basenames"]],
             gsmvi)