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
  file containing an annotation database .RData file, not necessary to specify if using pre-built concept database.
  Must have the following variables: 
  'go2eg', 'eg2go' (both lists of lists)
  'n.genes' (a list of integers the same length as go2eg, where n.genes[[i]] <- length(go2eg[[i]])
  'goHumanLabels' (a character vector the same length as go2eg containing human readable annotations the corresponding annotation)
%%     ~~Describe \code{DB.data.file} here~~
}
  \item{intMatrix.file}{
  file containing an intersection matrix (.RData file), not necessary to specify if using pre-built concept database.
  Must have the following variable:
  'intMatrix', a [length(go2eg) x length(go2eg)] matrix. The element accessed by intMatrix[i,j]
  should be an integer corresponding to the number of genes shared by the ith and jth annotation terms
  (i.e. length(intersect(go2eg[[i]], go2eg[[j]]))  )
  This file can be computed using preProcess_db function of this package provided with argument (2) from this file.
%%     ~~Describe \code{intMatrix.file} here~~
}
  \item{entrez2sym.file}{
  file containing mapping from Entrez Gene IDs to Symbols, not necessary to specify if using pre-built concept database.
%%     ~~Describe \code{entrez2sym.file} here~~
}
  \item{ji.threshold}{
  A Jaccard Index threshold for removing concepts from consideration (mainly for testing), default = 2 (not filtering).
  A Jaccard Index will be calculated between the provied gene list and every concept. If this JI is 
  greater than the input value, the concept will not be considered in the calculations. If no filtering is 
  requested, use a value >1 for this, which will not filter out any concepts
%%     ~~Describe \code{ji.threshold} here~~
}
  \item{size.min}{
  A minimum concept size to impose, default =5.
  Any concepts with fewer genes than this will not be considered for analysis.
  Precomputation has already imposed a minimum of 5, so any value of 5 or smaller will have no effect.
%%     ~~Describe \code{size.min} here~~
}
  \item{size.max}{
  A maximum concept size to impose. Currently this function was disabled.
  Any concepts with more genes than this will not be considered for analysis.
  Precomputation has already imposed a maximum of 1500, so any value greater than this will have no effect.
%%     ~~Describe \code{size.max} here~~
}
  \item{filter.strings}{
  A string of concepts to filter out. Default is 'NA' (not filter anything out). 
  To filter multiple words, use a long string, separated by commas, i.e. 'CANCER,CARCINOMA,TUMOR'. 
  Any concept that contains any of these strings will be excluded from analysis. 
  This is case insensitive (i.e. an input of 'cancer', 'Cancer', and 'CANCER' will all filter the same concepts)
  This does not support using different filters for different gene lists. 
  If you want to do this, run each gene list separately.
%%     ~~Describe \code{filter.strings} here~~
}
  \item{use.databases}{
  A string of databases to use.All other databases will be filtered out. Default is NA. 
  You need to know the name of the databases used to make the database (from arguments (2) and (3), the first column###NEED TO FINISH###
  This is case insensitive.
%%     ~~Describe \code{use.databases} here~~
}
  \item{filtered.file}{
  File for filtered concepts (the names of all concepts that have been filtered out due to arguments 4-7)
%%     ~~Describe \code{filtered.file} here~~
}
  \item{output.file}{
  Output file (without an extension, will make both .RData and .txt versions). Need to be specified.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (geneListFile, DB.data.file, intMatrix.file, entrez2sym.file, 
    ji.threshold, size.min, size.max, filter.strings, use.databases, 
    filtered.file, output.file) 
{
    library("Matrix")
    library("methods")
    print("Reading Gene List")
    allGeneLists <- read.table(geneListFile, sep = "\t")
    geneList <- list()
    geneList[[1]] <- list(as.character(allGeneLists[1, 3]))
    first <- c(as.character(allGeneLists[1, 1]))
    second <- c(as.character(allGeneLists[1, 2]))
    print("List-ifying Gene List")
    for (i in 2:dim(allGeneLists)[1]) {
        flag <- 0
        for (j in 1:length(first)) {
            if ((first[j] == allGeneLists[i, 1]) && (second[j] == 
                allGeneLists[i, 2])) {
                flag <- 1
                geneList[[j]] <- append(geneList[[j]], as.character(allGeneLists[i, 
                  3]))
                break
            }
        }
        if (flag == 0) {
            first <- append(first, as.character(allGeneLists[i, 
                1]))
            second <- append(second, as.character(allGeneLists[i, 
                2]))
            geneList[[length(first)]] <- list(as.character(allGeneLists[i, 
                3]))
        }
    }
    geneListLabels <- data.frame(source = first, name = second)
    newGeneList <- list()
    for (i in length(geneList)) {
        newGeneList[[i]] <- unique(geneList[[i]])
        if (length(newGeneList[[i]] != length(geneList[[i]]))) {
            print(paste("WARNING: Duplicate genes detected in ", 
                paste(geneListLabels[i, 1], geneListLabels[i, 
                  2], sep = ":"), ". Duplicates have been removed.", 
                sep = ""))
        }
        geneList[[i]] <- newGeneList[[i]]
    }
    rm(newGeneList)
    print(paste("Loading", DB.data.file))
    load(DB.data.file)
    n.Terms <- length(go2eg)
    conceptGeneListInt <- matrix(0, nrow = n.Terms, ncol = length(geneList))
    conceptGeneListUnion <- matrix(0, nrow = n.Terms, ncol = length(geneList))
    geneListFilter <- matrix(T, nrow = n.Terms, ncol = length(geneList))
    filtered.concepts <- c()
    ji.filtered <- 0
    size.filtered <- 0
    string.filtered <- 0
    print("Computing Intersections and Unions between Database and Gene List(s)")
    for (i in 1:n.Terms) {
        for (j in 1:length(geneList)) {
            conceptGeneListInt[i, j] <- length(intersect(go2eg[[i]], 
                geneList[[j]]))
            conceptGeneListUnion[i, j] <- length(union(go2eg[[i]], 
                geneList[[j]]))
            if ((conceptGeneListInt[i, j]/conceptGeneListUnion[i, 
                j]) > ji.threshold) {
                filtered.concepts <- append(filtered.concepts, 
                  goHumanLabels[i])
                ji.filtered <- ji.filtered + 1
                geneListFilter[i, j] = F
            }
        }
    }
    if (ji.filtered != 0) {
        print(paste("Filtered", ji.filtered, "concept(s) from the Jiccard Index threshold (each gene lists counts separately)"))
    }
    if (size.min > 4) {
        for (i in 1:length(n.genes)) {
            if (n.genes[i] > size.min) {
                filtered.concepts <- append(filtered.concepts, 
                  goHumanLabels[i])
                geneListFilter[i, ] = F
                size.filtered <- size.filtered + 1
            }
        }
        print(paste("Filtered", size.filtered, "concept(s) based on size restrictions"))
    }
    filter.strings <- toupper(filter.strings)
    if (filter.strings == "NA") {
        print("Skipping concept-name based filtering")
    }
    else {
        print(paste("Filtering based on", filter.strings))
        each.str <- strsplit(filter.strings, ",")[[1]]
        for (i in 1:length(each.str)) {
            for (j in 1:length(goHumanLabels)) {
                if (length(grep(each.str[i], goHumanLabels[j])) > 
                  0) {
                  string.filtered <- string.filtered + 1
                  filtered.concepts <- append(filtered.concepts, 
                    goHumanLabels[j])
                  geneListFilter[j, ] = F
                }
            }
        }
        if (string.filtered > 0) {
            print(paste("Filtered", string.filtered, "concept(s) based on string-matching."))
        }
        else {
            print("No concepts were filtered from the string-based filtering")
        }
    }
    if ((ji.filtered + string.filtered + size.filtered) > 0) {
        print(paste("Writing all filtered concepts to", filtered.file))
        write.table(filtered.concepts, quote = F, row.names = F, 
            col.names = F, file = filtered.file)
    }
    print("Reading entrezID -> GeneSymbol Converter")
    entrez2sym <- read.table(entrez2sym.file, sep = "\t")
    output.genelist <- vector("character", 0)
    output.genelist.size <- vector("numeric", 0)
    output.geneSym <- vector("numeric", 0)
    output.geneID <- vector("character", 0)
    output.gene.numConcepts <- vector("numeric", 0)
    output.gene.time <- vector("numeric", 0)
    output.inputStatus <- vector("numeric", 0)
    output.concept.score.sum <- vector("numeric", 0)
    output.penal <- vector("numeric", 0)
    output.score <- vector("numeric", 0)
    output.conceptsWithScore <- vector("character", 0)
    print("Loading Interaction Matrix")
    load(intMatrix.file)
    start.time <- proc.time()
    last.time <- start.time
    for (geneListInd in 1:length(geneList)) {
        print(paste("Starting to process Gene List #", geneListInd, 
            ": ", paste(geneListLabels[geneListInd, 1], geneListLabels[geneListInd, 
                2], sep = ":"), sep = ""))
        valid.GO.inds <- which(geneListFilter[, geneListInd])
        for (geneID in entrez2sym[, 1]) {
            if (is.na(match(geneID, labels(eg2go)))) {
                output.genelist <- append(output.genelist, paste(geneListLabels[geneListInd, 
                  1], geneListLabels[geneListInd, 2], sep = ":"))
                output.genelist.size <- append(output.genelist.size, 
                  length(geneList[[geneListInd]]))
                output.geneSym <- append(output.geneSym, as.character(entrez2sym[match(geneID, 
                  entrez2sym[, 1]), 2]))
                output.geneID <- append(output.geneID, geneID)
                output.gene.numConcepts <- append(output.gene.numConcepts, 
                  0)
                output.gene.time <- append(output.gene.time, 
                  0)
                output.inputStatus <- append(output.inputStatus, 
                  sum(geneID == geneList[[geneListInd]]))
                output.concept.score.sum <- append(output.concept.score.sum, 
                  0)
                output.penal <- append(output.penal, 0)
                output.score <- append(output.score, 0)
                output.conceptsWithScore <- append(output.conceptsWithScore, 
                  0)
            }
            else {
                geneIDInd <- match(geneID, labels(eg2go))
                each.time <- proc.time()
                if ((geneIDInd\%\%100) == 0) {
                  print(paste("In Gene List #", geneListInd, 
                    ";  Done with gene #", geneIDInd - 1, " (of ", 
                    length(eg2go), "). Seconds for this 100: ", 
                    as.character((each.time - last.time)[3]), 
                    sep = ""))
                  last.time <- each.time
                }
                goIDLookup <- vector("numeric", length(eg2go[[geneIDInd]]))
                for (goTermInd in 1:length(eg2go[[geneIDInd]])) {
                  goIDLookup[goTermInd] <- which(eg2go[[geneIDInd]][[goTermInd]] == 
                    labels(go2eg))
                }
                goIDs <- intersect(goIDLookup, valid.GO.inds)
                thisConListInt <- conceptGeneListInt[goIDs, geneListInd]
                if (sum(geneList[[geneListInd]] == labels(eg2go)[geneIDInd]) == 
                  1) {
                  thisConListInt <- thisConListInt - 1
                }
                jiCancer <- thisConListInt/conceptGeneListUnion[goIDs, 
                  geneListInd]
                goLengths <- as.vector(n.genes)[goIDs]
                thisUnionMatrix <- matrix(goLengths, nrow = length(goLengths), 
                  ncol = length(goLengths), byrow = T) + matrix(goLengths, 
                  nrow = length(goLengths), ncol = length(goLengths), 
                  byrow = F) - intMatrix[goIDs, goIDs]
                cc.ji <- intMatrix[goIDs, goIDs]/thisUnionMatrix
                myArray <- as.numeric(cc.ji)
                myArray[myArray <= 0.05] <- 0
                cc.ji <- matrix(myArray, nrow = length(goLengths), 
                  ncol = length(goLengths))
                row.cc.ji <- rowSums(cc.ji) + colSums(cc.ji) - 
                  1
                conceptScore <- jiCancer/row.cc.ji
                denom <- 1/row.cc.ji
                conceptScore[jiCancer == 0] <- 0
                sumDenomScore <- sum(denom)
                if ((length(eg2go[[geneIDInd]]) > 0) && (sumDenomScore > 
                  0)) {
                  sumConceptScore <- sum(conceptScore)
                  score <- sumConceptScore/sqrt(sumDenomScore)
                }
                else {
                  score <- 0
                  sumDenomScore <- 0
                  sumConceptScore <- 0
                }
                conceptOrder <- order(conceptScore, decreasing = T)
                conceptOutputSoFar <- paste(conceptScore[conceptOrder[1]], 
                  labels(go2eg)[goIDs[conceptOrder[1]]], sep = ":")
                for (i in 2:length(conceptScore)) {
                  conceptOutputSoFar <- paste(conceptOutputSoFar, 
                    paste(conceptScore[conceptOrder[i]], labels(go2eg)[goIDs[conceptOrder[i]]], 
                      sep = ":"), sep = ";")
                }
                test <- as.character(entrez2sym[which(as.numeric(entrez2sym[, 
                  1]) == as.numeric(labels(eg2go)[geneIDInd])), 
                  2])
                if (length(test) == 0) {
                  print(paste("No geneID found in entrez2sym for geneID", 
                    labels(eg2go)[geneIDInd]))
                }
                else {
                  output.genelist <- append(output.genelist, 
                    paste(geneListLabels[geneListInd, 1], geneListLabels[geneListInd, 
                      2], sep = ":"))
                  output.genelist.size <- append(output.genelist.size, 
                    length(geneList[[geneListInd]]))
                  output.geneSym <- append(output.geneSym, as.character(entrez2sym[which(as.numeric(entrez2sym[, 
                    1]) == as.numeric(labels(eg2go)[geneIDInd])), 
                    2]))
                  output.geneID <- append(output.geneID, labels(eg2go)[geneIDInd])
                  output.gene.numConcepts <- append(output.gene.numConcepts, 
                    length(goIDs))
                  output.gene.time <- append(output.gene.time, 
                    as.character((proc.time() - each.time)[3]))
                  output.inputStatus <- append(output.inputStatus, 
                    sum(labels(eg2go)[geneIDInd] == geneList[[geneListInd]]))
                  output.concept.score.sum <- append(output.concept.score.sum, 
                    sumConceptScore)
                  output.penal <- append(output.penal, sumDenomScore)
                  output.score <- append(output.score, score)
                  output.conceptsWithScore <- append(output.conceptsWithScore, 
                    conceptOutputSoFar)
                }
            }
        }
    }
    ourOrder <- order(output.score, decreasing = TRUE)
    ourOutput <- data.frame(genelist = output.genelist[ourOrder], 
        genelist.size = output.genelist.size[ourOrder], geneSym = output.geneSym[ourOrder], 
        entrezID = output.geneID[ourOrder], inGenelist = output.inputStatus[ourOrder], 
        raw.score = output.concept.score.sum[ourOrder], penal.score = output.penal[ourOrder], 
        score = output.score[ourOrder], numConcepts = output.gene.numConcepts[ourOrder])
    save(ourOutput, file = paste(output.file, ".RData", sep = ""))
    write.table(ourOutput, file = paste(output.file, ".txt", 
        sep = ""), quote = F, sep = "\t", row.names = F)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
