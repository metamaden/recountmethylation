.onAttach = function(libname, pkgname) {
  version = packageDescription(pkgname, fields = "Version")
  msg = paste0("",
"\n***************************************
\n", pkgname, 
"\n\n* VERSION:\n", version,
"\n\n* DISCLAIMER:\n
Recount Methylation contains data from GEO 
(ncbi.nlm.nih.gov/geo/), a live public 
database where alterations to online records 
can cause discrepancies with stored data 
over time. We cannot gaurantee the accuracy 
of data in Recount Methylation and advise 
users cross-check their findings with 
latest available records.
\n\n* CITATION:\n
If you use this package in published research, 
please cite the accompanying manuscript.
\n***************************************")	
  packageStartupMessage(msg)
}
