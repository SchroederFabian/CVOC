
This git overviews the creation and further development of the R package CVOC. 

In order to run CVOC the C libraries gmp, mfpr and mfprc++ must be installed

for Linux (Debian) do:
 sudo apt-get install libgmp3-dev
 sudo apt-get install libmpfr-dev
 sudo apt-get install libmpfrc++-dev
 
then in the R console:
 library("devtools")
 install_github("Thell/RcppMP")
 install.packages("gmp")
 install.packages("Rmpfr")
 
 
 
 
 
 