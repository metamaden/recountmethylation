\name{recountmethylationnews}
\title{Recountmethylation News}
\encoding{UTF-8}

\section{Version 1.1.3}{
  \itemize{
    \item Added v.0.0.2 database compilation files to server (recount.bio/data).
    These reflect IDAT downloads completed in Nov 2020 from GEO/GDS, including 
    the first compilations of EPIC/HM850K arrays.
    \item Added cross-platform support for HM450K and EPIC array files. Added 
    `platform` argument in relevant `getdb()` functions. Updated the 
    `servermatrix` function and added the `smfilt` function to improve server
    file identification across platforms.
    \item Updated the User's Guide to fix typos, reflect v.0.0.2 samples, and 
    add a download troubleshoot section. Revised citations format, and removed 
    evaluation of validation section due to periodic build failures from 
    connectivity issues.
    \item Updated ExperimentHub file metadata script and table to add new 
    v.0.0.2 compilation files. Renamed included sample metadata directory to 
    "gsm_metadata" to avoid confusion with "metadata.csv" file table for hubs.
  }
}

\section{Version 0.99.0}{
  \itemize{
    \item Added getdb functions for database file download and load 
    \item Added servermatrix function to get latest database file metadata
    \item Added User's Guide and Data Analyses vignettes
    \item Added metadata and data_analyses.RData to /inst/extdata
  }
}

\section{Version 0.01}{
    \itemize{
        \item Added key query/accessor functions.
        \item Added package vignette.
    }
}