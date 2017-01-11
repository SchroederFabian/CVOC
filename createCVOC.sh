cd /home/fabian/Desktop/CVOCcreate
R -e 'library("Rcpp"); Rcpp.package.skeleton( name="CVOC",author="Fabian Schroeder",email="fabian.schroeder.fl@ait.ac.at",license="GPL-2",path="/home/fabian/Desktop/CVOCcreate/",code_files = c("/home/fabian/Desktop/CVOCcreate/ebc.R","/home/fabian/Desktop/CVOCcreate/etc.R","/home/fabian/Desktop/CVOCcreate/roots.R","/home/fabian/Desktop/CVOCcreate/etcND.R","/home/fabian/Desktop/CVOCcreate/countPerm.R"),cpp_files="/home/fabian/Desktop/CVOCcreate/etcND.cpp",example_code=FALSE)'
rm CVOC/Read-and-delete-me
R -e 'library("roxygen2"); roxygenise(package.dir="CVOC", roclets="rd", clean=TRUE)'
cp DESCRIPTION CVOC/DESCRIPTION
cp NAMESPACE CVOC/NAMESPACE
cp CVOC-package.Rd CVOC/man/CVOC-package.Rd
cp Makevars CVOC/Makevars
R CMD check "CVOC"
#R CMD build "CVOC"
