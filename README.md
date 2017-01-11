
This git overviews the creation and further development of the R package CVOC. The creation of the package is controlled by the bash script CVOCcreate.sh. This script calls all the necessary R functions in the following order.

1) using Rcpp::Rcpp.package.skeleton to create the package skeleton.
2) copying the required files, such as DESCRIPTION, NAMESPACE, Makevars, etc. , into the correct folders
3) creating the .Rd documentation files using roxygen2::roxygenise()
4) checking the R-package using R CMD check
5) building the R-package using R CMD build
