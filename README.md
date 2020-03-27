
# `recountmethylation` [<img style="float: right;" src = "inst/figures/remeth_hexsticker.png" height="180"/>](https://recount.bio/data)
This package provides methods to access, query, and store DNAm signal and learned sample metadata 
annotations from the Recount Methylation database files.

# Authors

Sean K. Maden, Ph.D. candidate, Computational Biology Program, 
Dept. of Biomedical Engineering, OHSU

Kasper D. Hansen, Ph.D., 
Prof. Dept. of Biostatistics, Johns Hopkins Bloomberg 
School of Public Health

Abhi Nellore, Ph.D.,
Prof. Computational Biology Program, Dept. of Biomedical 
Engineering & Dept. of Surgery, OHSU

# Installation and data access

Install this library from an R session using 
`require(devtools); install_github("metamaden/recountmethylation")`.

# Introduction

# Disclaimer

The `recountmethylation` R package can aid and expand research capabilities 
in epigenetics, especially for initial discovery and exploration phases of an 
experiment utilizing published data. The data contained in Recount Methylation, 
a cross-study compilation of DNAm signal matrices and learned metadata, is 
extensively characterized and described in the accompanying manuscript, which 
also vitally includes transparent methods for how we obtained samples and 
learned sample metadata annotations. 

Recount Methylation contains data from GEO, a live public database where 
alterations (e.g. updates, revisions, additions, and removals) to online 
records can cause discrepancies with stored data over time. Thus, we cannot 
universally guarantee the data contained in Recount Methylation will reflect 
the current state of latest available corresponding records in GEO. We 
showed how to check retrieved data against a fresh download of sample IDATs 
from GEO Data Sets (see vignette). Through this and other means, 
we advise `recountmethylation` users be vigilant in cross-checking their 
findings with latest available records and study data.
