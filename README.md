
[<img style="float: right;" src = "inst/figures/remeth_hexsticker.png" height="180"/>](https://recount.bio/data)
# Authors

Sean K. Maden, Reid F. Thompson, Kasper D. Hansen, Abhinav Nellore

# Manuscript
Sean K. Maden, Reid F. Thompson, Kasper D. Hansen, Abhinav Nellore. "Human 
methylome variation across Infinium 450K raw data on the Gene Expression 
Omnibus." 2019 (preprint in preparation for submission).

# Installation and data access

From an active R session, install this package from Bioconductor using:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("recountmethylation")
```

You may alternatively install this package from GitHub using: 
 
`require(devtools); install_github("metamaden/recountmethylation")`.

# Package disclaimer

Databases accessed with recountmethylation 
contain data from GEO (ncbi.nlm.nih.gov/geo/), 
a live public database where alterations to 
online records can cause discrepancies with 
stored data over time. We cannot guarantee 
the accuracy of stored data, and advise users 
cross-check their findings with latest available 
records.
