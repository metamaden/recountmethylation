% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_db.R
\name{smfilt}
\alias{smfilt}
\title{smfilt}
\usage{
smfilt(sm, typesdf = NULL)
}
\arguments{
\item{sm}{Data matrix returned from servermatrix().}

\item{typesdf}{Data.frame containing database file info for dm filters.}
}
\value{
Filtered data matrix of server file info.
}
\description{
Filters the data matrix returned from servermatrix().
}
\examples{
dm <- matrix(c("remethdb_h5-rg_epic_0-0-2_1589820348.h5","08-Jan-2021",
"09:46","66751358297"), nrow = 1)
smfilt(dm)
}
\seealso{
get_rmdl, servermatrix
}
