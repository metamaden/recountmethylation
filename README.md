
[<img style="float: right;" src = "inst/figures/remeth_hexsticker.png" height="180"/>](https://recount.bio/data)
# Authors

Sean K. Maden, Reid F. Thompson, Kasper D. Hansen, Abhinav Nellore

# Overview

The `recountmethylation` package provides access to databases of DNA 
methylation (DNAm) data from over 62,000 cumulative sample records with 
IDATs in the Gene Expression Omnibus (GEO, available by November, 2020). 
Samples were run using either of 2 Illumina BeadArray platform types, 
either the older HM450K platform or the newer EPIC/HM850K platform. The
database compilation files include mined, mapped, and model-based sample 
metadata, and DNAm data in the form of either raw/unnormalized red and 
green signals, raw/unnormalized methylated and unmethylated signals, or 
normalized DNAm fractions (a.k.a. "Beta-values").

# Installation and data access

From an active R session, install the latest available package version from 
Bioconductor with:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("recountmethylation")
```

You may also install the package from GitHub using: 
 
`require(devtools); install_github("metamaden/recountmethylation")`.

See the package User's Guide vignette for more details and background about the 
datasets.

# Manuscript

If you use `recountmethylation` in published work, please cite the following 
paper:

Sean K Maden, Reid F Thompson, Kasper D Hansen, 
Abhinav Nellore, 'Human methylome variation across 
Infinium 450K data on the Gene Expression Omnibus', 
NAR Genomics and Bioinformatics, Volume 3, Issue 2, 
June 2021, lqab025, https://doi.org/10.1093/nargab/lqab025

# Package disclaimer

Please note the following disclaimer regarding the data contained in the 
compilation files.

```
Databases accessed with recountmethylation 
contain data from GEO (ncbi.nlm.nih.gov/geo/), 
a live public database where alterations to 
online records can cause discrepancies with 
stored data over time. We cannot guarantee 
the accuracy of stored data, and advise users 
cross-check their findings with latest available 
records.
```
