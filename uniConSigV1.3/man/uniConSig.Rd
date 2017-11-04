\name{uniConSig}
\alias{uniConSig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{calculation of uniConSig
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
uniConSig(geneListFile, 
          DB.data.file="./uniConSigV1.3/data/uniConSig_annotationFile_v3.0_rmConceptFilter_conceptName.Rdata", 
		  intMatrix.file="./uniConSigV1.3/data/intersect_uniConSig_JMatrix_v2.2_rmConceptFilter_conceptName.RData", 
		  entrez2sym.file="./uniConSigV1.3/data/entrez2sym_v1.1", 
		  ji.threshold=2,size.min=5,size.max=12000,filter.strings="NA", use.databases=NA, 
		  filtered.file="./filteredConcepts", 
		  output.file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneListFile}{
  gene list file is consisted of 3 columns. First 2 columns specify gene list, 3rd specifies entrezID, seperated by tab.
%%     ~~Describe \code{geneListFile} here~~
}
  \item{DB.data.file}{
  file containing an annotation database .RData file, not necessary to specify if using prebuilt concept database.
  Must have the following variables: 
  'go2eg', 'eg2go', both lists of lists
  'n.genes', a list of integers the same length as go2eg,
  'goHumanLabels', a character vector the same length as go2eg containing human readable annotations the corresponding annotation
%%     ~~Describe \code{DB.data.file} here~~
}
  \item{intMatrix.file}{
  file containing an intersection matrix file, not necessary to specify if using prebuilt concept database.
  Must have the following variable:
  'intMatrix', an integer matrix corresponding to the number of genes shared by the ith and jth annotation terms
  This file can be computed using preProcess_db function of this package provided with the value of argument 2.
%%     ~~Describe \code{intMatrix.file} here~~
}
  \item{entrez2sym.file}{
  file containing mapping from Entrez Gene IDs to Symbols, not necessary to specify if using prebuilt concept database.
%%     ~~Describe \code{entrez2sym.file} here~~
}
  \item{ji.threshold}{
  A Jaccard Index threshold for removing concepts from consideration, default is 2, which means no filtering.
  A Jaccard Index will be calculated between the provied training gene set and every concept. If this JI is 
  greater than the input value, the concept will not be considered in the calculations. If no filtering is 
  requested, use a value larger than 1 for this, which will not filter out any concepts
%%     ~~Describe \code{ji.threshold} here~~
}
  \item{size.min}{
  A minimum concept size to impose, default is 5.
  Any concepts with fewer genes than this will not be considered for analysis.
  Precomputation has already imposed a minimum of 5, so any value of 5 or smaller will have no effect.
%%     ~~Describe \code{size.min} here~~
}
  \item{size.max}{
  A maximum concept size to impose. Currently this function was disabled.
%%     ~~Describe \code{size.max} here~~
}
  \item{filter.strings}{
  A string of concepts to filter out. Default is 'NA'. To filter multiple words, use a long string, separated by commas, i.e. 'CANCER,CARCINOMA,TUMOR'. 
  Any concept that contains any of these strings will be excluded from analysis. 
  This is case insensitive. This does not support using different filters for different training gene sets. 
  If you want to do this, run each training gene set separately.
%%     ~~Describe \code{filter.strings} here~~
}
  \item{use.databases}{
  A string of databases to use. All other databases will be filtered out. Default is NA. 
  You need to know the name of the databases used to make the database from arguments 2 and 3.
  This is case insensitive.
%%     ~~Describe \code{use.databases} here~~
}
  \item{filtered.file}{
  File for filtered concepts, containing the names of all concepts that have been filtered out due to arguments 4 to 7.
%%     ~~Describe \code{filtered.file} here~~
}
  \item{output.file}{
  Output file, no extension. This function will directly make both .RData and .txt versions of results. Need to be specified.
%%     ~~Describe \code{output.file} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
Calculate the uniConSig scores using the cancer gene list from Cancer Gene Census, 
based on preBuilt concept database:

uniConSig(geneListFile="./uniConSigV1.3/tests/CancerGene_CGC20150525.db",output.file="test1")

The result will have a prefix of 'test1'
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
