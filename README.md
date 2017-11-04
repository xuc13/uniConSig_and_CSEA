# uniConSig_and_CSEA
development of uniConSig and CSEA algorithm
1. Installation of uniConSig in R
1.1 Install development tools for local package installations
UniConSig has been implemented in an R package. Currently it can only be installed by local package installation. To do local package installation in R, you need to install the package "devtools" first. The following R codes can be directly pasted into R for the installation:

install.packages("devtools")
library("devtools")

1.2 Install uniConSig package in R
After installing "devtools", the following R codes can be directly pasted into R console for uniConSig's installation:

setwd("Path_to_Rmodule\\R_module_construction")
install("uniConSigV1.3")
library(uniConSigV1.3)

1.3 Test run
The training gene set and the concept database used in our paper are included in the package. The test run can be performed simply using the following code:

uniConSig(geneListFile="./uniConSigV1.3/tests/CancerGene_CGC20150525.db",output.file="test1")

For details of the other options, please view the results of typing ?uniConSig in R.
