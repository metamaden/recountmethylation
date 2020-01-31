
# The `recountmethylation` package [<img style="float: right;" src = "inst/figures/remeth_hexsticker.png" height="180"/>](https://recount.bio/data)
Library to query, access, and store DNAm signal and learned sample metadata 
annotations from the Recount Methylation HDF5 database.

# Authors

Sean K. Maden, Ph.D. candidate, Computational Biology Program, 
Dept. of Biomedical Engineering, OHSU

Kasper D. Hansen, Ph.D., 
Professor, Dept. of Biostatistics, Johns Hopkins Bloomberg 
School of Public Health

Abhi Nellore, Ph.D.,
Professor, Computational Biology Program, Dept. of Biomedical 
Engineering & Dept. of Surgery, OHSU

# Installation and data access

Install this library from GitHub with 
`require(devtools); install_github("metamaden/recountmethylation")`.

The latest version of the Recount Methylation HDF5 database is contained
in the `remethdb2.h5` file, which can be downloaded from 
[https://recount.bio/data/](https://recount.bio/data/). We presently 
recommend users download the database file (about 120Gb in size) and
use `recountmethylation` for local access (e.g. by setting 
`dbn = "remethdb2.h5"`).

# Introduction

This package provides accessor functions to query the Recount Methylation HDF5 
database file, `remethdb2.h5`. The database stores relatively large red- and 
green-channel raw signal tables in 
[hierarchical data format](https://portal.hdfgroup.org/display/HDF5/HDF5), 
or "HDF5". These DNAm signal tables have been extracted directly from IDAT 
files generated using the Illumina HM450k array DNAm array platform. 
DNAm arrays were accessed from the 
[Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/) via the 
GEO Data Sets database and 
[Entrez programming utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/) 
software.

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
the current state of latest available corresponding records in GEO. We 
showed how to check retrieved data against a fresh download of sample IDATs 
from GEO Datat Sets (see vigentte). Through this and other means, 
we advise `recountmethylation` users be vigilant in cross-checking their 
findings with latest available records and study data.
