.onAttach = function(libname, pkgname) {
  version = packageDescription(pkgname, fields = "Version")
  msg = paste0("",
"#---------------------------------------------#
* ", pkgname, " version ", version, 
"* DISCLAIMER
Recount Methylation contains data from GEO 
(https://www.ncbi.nlm.nih.gov/geo/), a live 
public database where alterations to online 
records can cause discrepancies with stored 
data over time. We make no garauntees to the 
accuracy of data in Recount Methylation and 
advise users cross-check their findings with 
latest available records,

* CITATION
If you use this package in published research, 
please cite the accompanying manuscript.
#---------------------------------------------#")	
  packageStartupMessage(msg)
}