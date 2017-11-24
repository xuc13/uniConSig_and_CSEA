# uniConSig_and_CSEA
This code of uniConSig is for research use only<br>
development of uniConSig and CSEA algorithm<br>

1. Installation of uniConSig in R<br>
1.1 Install development tools for local package installations<br>
UniConSig has been implemented in an R package. Currently it can only be installed by local package installation. To do local package installation in R, you need to install the package "devtools" first. The following R codes can be directly pasted into R for the installation:<br>

install.packages("devtools")<br>
library("devtools")<br>

1.2 Install uniConSig package in R<br>
After installing "devtools", the following R codes can be directly pasted into R console for uniConSig's installation:<br>

setwd("Path_to_Rmodule")<br>
install("uniConSigV1.3")<br>
library(uniConSigV1.3)<br>

Remember to change "Path_to_Rmodule" to the path where you put the folder "uniConSigV1.3" (the folder that contains "uniConSigV1.3"). <br> 

1.3 Test run<br>
The training gene set and the concept database used in our paper are included in the package. The test run can be performed simply using the following code:<br>

uniConSig(geneListFile="./uniConSigV1.3/tests/CancerGene_CGC20150525.db",output.file="test1")<br>

For details of the other options, please type ?uniConSig in R.<br>
