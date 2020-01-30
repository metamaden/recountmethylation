.onAttach = function(libname, pkgname) {
  version = packageDescription(pkgname, fields = "Version")
  
  msg = paste0("***
               ", pkgname, " version ", version, "
               #-----------
               # DISCLAIMER
               #-----------
               Recount Methylation contains data from GEO, a live public database where 
               alterations (e.g. updates, revisions, additions, and removals) to online 
               records can cause discrepancies with stored data over time. We advise 
               users to be vigilant in cross-checking their findings with latest available 
               records and study data.
               
               #---------
               # CITATION
               #---------
               If you use this package in published research, please cite the accompanying manuscript.
               ***
               ")	
  
  packageStartupMessage(msg)
}