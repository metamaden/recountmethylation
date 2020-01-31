# The `recountmethylation` package
Library to query, access, and store DNAm signal and learned sample metadata annotations from the Recount Methylation HDF5 database.

# Installation

Access this library from GitHub with 
`require(devtools); install_github("metamaden/recountmethylation")`

# Introduction

This package provides accessor functions to query Recount Methylation 
database file, `remethdb2.h5`. The database stores relatively large red- and 
green-channel raw signal tables in 
[hierarchical data format](https://www.hdfgroup.org/), or "HDF5". These DNAm 
signal tables have been extracted directly from IDAT files generated using the 
Illumina HM450k array DNAm array platform. DNAm arrays were accessed from 
the [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/) via the 
GEO Data Sets database.

To help inform sample queries, the HDF5 database includes a sample metadata 
table learned from GEO GSE 
[SOFT](https://www.ncbi.nlm.nih.gov/geo/info/soft.html) files. 
When specifying samples whose DNAm data to return, several arguments assist 
with querying and accessing the database files. For instance, returned objects 
can be either a list of datasets or a single object of class `RGChannelSet`, 
which inherits properties from the `SummarizedExperiment` class.

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
the current state of latest available corresponding records in GEO. In the 
vignette, we showed how to check retrieved data against a fresh download of 
sample IDATs from GEO Datat Sets (see below). Through this and other means, 
we advise `recountmethylation` users be vigilant in cross-checking their 
findings with latest available records and study data.

# Data access

The Recount Methylation [HDF5](https://www.hdfgroup.org/) database can be 
downloaded from (https://recount.bio/data/)[https://recount.bio/data/].
