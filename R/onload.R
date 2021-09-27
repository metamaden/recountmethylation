.onAttach = function(libname, pkgname) {
  version = utils::packageDescription(pkgname, fields = "Version")
  msg = paste0("",
"\n***************************************
\n", pkgname, 
"\n\n* VERSION:\n", version,
"\n\n* DISCLAIMER:\n
Databases accessed with recountmethylation 
contain data from GEO (ncbi.nlm.nih.gov/geo/), 
a live public database where alterations to 
online records can cause discrepancies with 
stored data over time. We cannot guarantee 
the accuracy of stored data, and advise users 
cross-check their findings with latest available 
records.
\n\n* CITATION:\n
If you use this package in published research, 
please cite the accompanying manuscript:
\n
Sean K Maden, Reid F Thompson, Kasper D Hansen, 
Abhinav Nellore, 'Human methylome variation across 
Infinium 450K data on the Gene Expression Omnibus', 
NAR Genomics and Bioinformatics, Volume 3, Issue 2, 
June 2021, lqab025, https://doi.org/10.1093/nargab
/lqab025
\n***************************************")
packageStartupMessage(msg)
}
