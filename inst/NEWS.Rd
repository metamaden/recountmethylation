\name{recountmethylationnews}
\title{Recountmethylation News}
\encoding{UTF-8}

\section{Version 1.5.1}{
  \itemize{
    \item Adds QC functions for BeadArray metrics and log M/U signals
    \item Adds data and accessor function for cross-reactive CpGs
    \item Adds vignette showing how to do power analysis with pwrEWAS
    \item Adds vignette showing how to infer genetic ancestry using GLINT/EPISTRUCTURE
    \item Adds vignette showing how to do nearest neighbors search using a search index
    \item Adds functions for feature hashing, search index construction, and KNN search
  }
}

\section{Version 1.1.4}{
  \itemize{
    \item Improves User's Guide, fixes typos, new citation, adds disclaimer text, updates servermatrix() chunk.
    \item Updates examples to further limit downloads. Function get_rmdl() uses download = FALSE, 
    \item Updates Data Analyses vignette. Uses reduced dpi for images in Data Analyses vignette to limit package size. Uses updated metadata file name.
    \item Compresses new metadata v.0.0.2 files to limit package size.
    \item Uses uniform metadata file label for v.0.0.1 file.
  }
}

\section{Version 1.1.3}{
  \itemize{
    \item Added v.0.0.2 database compilation files to server (recount.bio/data)
    and revised recountmethylation functions for cross-platform support.
    The new files reflect IDAT downloads completed in Nov 2020 from GEO/GDS, 
    including the first compilations of EPIC/HM850K arrays.
    \item Added `platform` argument in relevant `getdb` functions.
    \item Added `which.platform` argument to `get_rmdl`
    \item Added new function `smfilt` to filter server data table for newest 
    compilation files, accounting for platform in the filename.
    \item Cleaned up and shoretened the `servermatrix` function. This now 
    handles RCurl call for "dn" (originally from `get_rmdl`) when handling 
    condition `dn = NULL`.
    \item Updated the User's Guide to fix typos, reflect v.0.0.2 samples, and 
    add a download troubleshoot section. Added numeric citations format, removed 
    evaluation of validation section due to possible package build failure 
    from bad internet connection. 
    \item Updated ExperimentHub file metadata script and table to add new 
    v.0.0.2 compilation files. 
    \item Renamed sample metadata directory to "gsm_metadata" to avoid 
    confusion with "metadata.csv" file table for hubs.
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