# Changes for v1.1: Starting to vectorize things: starting with for (goTermInd)

# inputs: (1) gene list file (3 columns, first 2 specify gene list, 3rd specifies entrezID)
# (2) file containing an annotation database (.RData file), must have the following variables: 
#  'go2eg', 'eg2go' (both lists of lists)
#  'n.genes' (a list of integers the same length as go2eg, where n.genes[[i]] <- length(go2eg[[i]])
#  'goHumanLabels' (a character vector the same length as go2eg containing human readable annotations the corresponding annotation)
# (3) file containing an intersection matrix (.RData file), must have the following variable:
#  'intMatrix', a [length(go2eg) x length(go2eg)] matrix. The element accessed by intMatrix[i,j]
#     should be an integer corresponding to the number of genes shared by the ith and jth annotation terms
#     (i.e. length(intersect(go2eg[[i]], go2eg[[j]]))  )
# !!! This file can be computed using preProcess.R provided with argument (2) from this file.
#     Since this matrix can be huge, by default we assume intMatrix is a sparse matrix (from library('Matrix'))
#       in triangular form (i.e. intMatrix[i,j]==0 if (j<i)), since the 'full' intMatrix is symmetric
#     If intMatrix is truly symmetric (i.e. intMatrix[i,j]==intMatrix[j,i] for all integers i,j),
#       then you need to modify the code slightly (the correct line is commented out, but I'd rather 
#	not put in a conditional that is evaluated GL*GN times (GL=#geneLists, GN=length(eg2go)))
# (4) file containing mapping from Entrez Gene IDs to Symbols
# (5) A Jaccard Index threshold for removing concepts from consideration (mainly for testing)
#   A Jaccard Index will be calculated between the provied gene list and every concept. If this JI is 
#   greater than the input value, the concept will not be considered in the calculations. If no filtering is 
#   requested, use a value >1 for this, which will not filter out any concepts
# (6) A minimum concept size to impose. Any concepts with fewer genes than this will not be considered for analysis.
#     Precomputation has already imposed a minimum of 5, so any value of 5 or smaller will have no effect.
# (7) A maximum concept size to impose. Any concepts with more genes than this will not be considered for analysis.
#     Precomputation has already imposed a maximum of 1500, so any value greater than this will have no effect.
# (8) A string of concepts to filter out. Use 'NA' to not filter anything out. To filter multiple words, use
#   a long string, separated by commas, i.e. 'CANCER,CARCINOMA,TUMOR'. Any concept that contains any of these
#   strings will be excluded from analysis. 
#   This is case insensitive (i.e. an input of 'cancer', 'Cancer', and 'CANCER' will all filter the same concepts)
#   This does not support using different filters for different gene lists. 
#     If you want to do this, run each gene list separately.
# (9) A string of databases to use. All other databases will be filtered out.
#   You need to know the name of the databases used to make the database (from arguments (2) and (3), the first column###NEED TO FINISH###
#   This is case insensitive.
# (10) File for filtered concepts (the names of all concepts that have been filtered out due to arguments 4-7)
# (11) Output file (without an extension, will make both .RData and .txt versions)

uniConSig <- function(geneListFile,DB.data.file="./uniConSigV1.3/data/uniConSig_annotationFile_v3.0_rmConceptFilter_conceptName.Rdata",
                      intMatrix.file="./uniConSigV1.3/data/intersect_uniConSig_JMatrix_v2.2_rmConceptFilter_conceptName.RData",
                      entrez2sym.file="./uniConSigV1.3/data/entrez2sym_v1.1",ji.threshold=2,
                      size.min=5,size.max=12000,filter.strings="NA",use.databases=NA,filtered.file="./filteredConcepts",output.file){
library('Matrix')  # for a sparse matrix format
library('methods')

print ('Reading Gene List')
allGeneLists <- read.table(geneListFile, sep='\t')

geneList <- list()	#Chi: training gene list
geneList[[1]] <- list(as.character(allGeneLists[1,3])) #Chi:initialize the first element
first <- c(as.character(allGeneLists[1,1]))
second <- c(as.character(allGeneLists[1,2]))

print ('List-ifying Gene List')
# transform the gene list data into a more R-friendly (and maybe more memory efficient) format
for (i in 2:dim(allGeneLists)[1]) {#Chi: dim(table)[1] is the row number of "table"
  flag <- 0
  # check all previous gene list labels to see if this line is new #Chi: can we do this in a better way? 
  for (j in 1:length(first)) {
    # if we've seen the gene list before, append this gene to the gene list and move on
    if ((first[j] == allGeneLists[i,1]) && (second[j] == allGeneLists[i,2])) {
      flag <- 1
      geneList[[j]] <- append(geneList[[j]], as.character(allGeneLists[i,3]))
      break
    }
  }
  if (flag == 0) { # then we haven't broken and have a new gene list
    first <- append(first, as.character(allGeneLists[i,1]))
    second <- append(second, as.character(allGeneLists[i,2]))
    geneList[[length(first)]] <- list(as.character(allGeneLists[i,3]))
  }
}

# the 'source' and 'name' labels are quite arbitrary
geneListLabels <- data.frame(source=first, name=second)

newGeneList <- list()

# it doesn't make sense to have duplicates in a gene list, so get rid of them
for (i in length(geneList)) {
    newGeneList[[i]] <- unique(geneList[[i]])
    if (length(newGeneList[[i]] != length(geneList[[i]]))) {
        print (paste('WARNING: Duplicate genes detected in ', paste(geneListLabels[i,1], geneListLabels[i,2], sep=':'), 
	             '. Duplicates have been removed.', sep=''))
    }
    geneList[[i]] <- newGeneList[[i]]
}

# Clean up
rm(newGeneList)

print(paste('Loading', DB.data.file))
load(DB.data.file)

n.Terms <- length(go2eg) #Chi: guess go2eg is "concepts to geneID"

# after testing, these are not sparse enough to be worth using Matrix(..., sparse=T). 
############################################
### Re-investigate with larger database! ###
############################################

conceptGeneListInt <- matrix(0, nrow=n.Terms, ncol=length(geneList)) 	#Chi: intersected
conceptGeneListUnion <- matrix(0, nrow=n.Terms, ncol=length(geneList)) 	#Chi: United

geneListFilter <- matrix(T, nrow=n.Terms, ncol=length(geneList))		#Chi: What's "T" here?

# a list of concepts we have filtered so the user knows what they did
filtered.concepts <- c()
# how many are removed due to each category
ji.filtered <- 0
size.filtered <- 0
string.filtered <- 0

print('Computing Intersections and Unions between Database and Gene List(s)')
for (i in 1:n.Terms) {
  for (j in 1:length(geneList)) {
    conceptGeneListInt[i,j] <- length(intersect(go2eg[[i]], geneList[[j]]))
    conceptGeneListUnion[i,j] <- length(union(go2eg[[i]], geneList[[j]]))

    # if the JI between the gene list and this concept is above the threshold,
    #   filter it out (see documentation for argument 4)
    if ((conceptGeneListInt[i,j]/conceptGeneListUnion[i,j]) > ji.threshold) {
      filtered.concepts <- append(filtered.concepts, goHumanLabels[i])
      ji.filtered <- ji.filtered + 1
      geneListFilter[i,j] = F
    }
  }
}
if (ji.filtered != 0) {
  print(paste('Filtered', ji.filtered, 'concept(s) from the Jiccard Index threshold (each gene lists counts separately)'))
}

# If the user wanted to do size-based filtering, do things
if (size.min > 4) { # these values are specified in the preProcess script. If they are changed there, they should be changed here

  for (i in 1:length(n.genes)) {
    if (n.genes[i] > size.min) { #Chi: guess n.genes came from DB.data.file? || n.genes[i] > size.max
      filtered.concepts <- append(filtered.concepts, goHumanLabels[i])
      geneListFilter[i,] = F
      size.filtered <- size.filtered + 1
    }
  }
  print(paste('Filtered', size.filtered, 'concept(s) based on size restrictions'))
}

# If the user wanted to filter out concepts based on string-matching
filter.strings <- toupper(filter.strings)
if (filter.strings == 'NA') {
  print('Skipping concept-name based filtering')
} else {
  print(paste('Filtering based on', filter.strings));
  each.str <- strsplit(filter.strings, ',')[[1]]
  for (i in 1:length(each.str)) {
    for (j in 1:length(goHumanLabels)) {
      # if this string matches any part of this concept, filter it
      if (length(grep(each.str[i], goHumanLabels[j])) > 0) {
        string.filtered <- string.filtered + 1
	filtered.concepts <- append(filtered.concepts, goHumanLabels[j])
	geneListFilter[j,] = F
      }
    }
  }
  if (string.filtered > 0) {
    print (paste('Filtered', string.filtered, 'concept(s) based on string-matching.'))
  } else {
    print ('No concepts were filtered from the string-based filtering')
  }
}

if ((ji.filtered + string.filtered + size.filtered) > 0) {
  print (paste('Writing all filtered concepts to', filtered.file))
  write.table(filtered.concepts, quote=F, row.names=F, col.names=F, file=filtered.file)
}


print ('Reading entrezID -> GeneSymbol Converter')
entrez2sym <- read.table(entrez2sym.file, sep='\t')

output.genelist <- vector('character', 0)
output.genelist.size <- vector('numeric', 0)
output.geneSym <- vector('numeric',0)
output.geneID <- vector('character',0)
output.gene.numConcepts <- vector('numeric',0)
output.gene.time <- vector('numeric',0)
output.inputStatus <- vector('numeric',0)
output.concept.score.sum <- vector('numeric',0)
output.penal <- vector('numeric',0)
output.score <- vector('numeric',0)
output.conceptsWithScore <- vector('character',0)

print ('Loading Interaction Matrix')
load(intMatrix.file) # interactions between known GO terms. May have different versions corresponding to different versions of GO

start.time <- proc.time()
last.time <- start.time

for (geneListInd in 1:length(geneList)) { # for each gene list
  print (paste('Starting to process Gene List #', geneListInd, ': ', paste(geneListLabels[geneListInd,1],geneListLabels[geneListInd,2],sep=':'), sep=''))
  valid.GO.inds <- which(geneListFilter[,geneListInd])
  #for (geneIDInd in 1:length(eg2go)) { # for each gene (with GO annotations)
  myCount <- 0
  for (geneID in entrez2sym[,1]){
    myCount <- myCount+1
    if(is.na(match(geneID,labels(eg2go)))){
    output.genelist <- append(output.genelist, paste(geneListLabels[geneListInd,1],geneListLabels[geneListInd,2],sep=':'))
    output.genelist.size <- append(output.genelist.size, length(geneList[[geneListInd]]))
    output.geneSym <- append(output.geneSym, as.character(entrez2sym[match(geneID,entrez2sym[,1]),2]))
    output.geneID <- append(output.geneID, geneID)
    output.gene.numConcepts <- append(output.gene.numConcepts, 0)
    output.gene.time <- append(output.gene.time, 0)
    output.inputStatus <- append(output.inputStatus, sum(geneID == geneList[[geneListInd]]))
    output.concept.score.sum <- append(output.concept.score.sum, 0)
    output.penal <- append(output.penal, 0)
    output.score <- append(output.score, 0)
    output.conceptsWithScore <- append(output.conceptsWithScore, 0)
    }else{
    geneIDInd<-match(geneID,labels(eg2go))
    each.time <- proc.time()
    if ((myCount %% 100) == 0) {
      print(paste('In Gene List #', geneListInd, ';  Done with gene #', myCount-1, ' (of ', length(entrez2sym[,1]),
        	  '). Seconds for this 100: ', as.character((each.time-last.time)[3]), sep=''))
      last.time <- each.time
    }

    # find the 'original' GO index; precomputing because otherwise we do it n^2 times (n = # concepts assigned to gene)
    goIDLookup <- vector('numeric', length(eg2go[[geneIDInd]]))
    for (goTermInd in 1:length(eg2go[[geneIDInd]])) { # for each GO Term assigned to this gene
      goIDLookup[goTermInd] <- which(eg2go[[geneIDInd]][[goTermInd]] == labels(go2eg))
    }

    goIDs <- intersect(goIDLookup, valid.GO.inds)
    thisConListInt <- conceptGeneListInt[goIDs,geneListInd]

    # if this gene is in the gene list, subtract 1 from each intersection (statistical reasons)
    if (sum(geneList[[geneListInd]] == labels(eg2go)[geneIDInd]) == 1) {
      thisConListInt <- thisConListInt - 1
    }
    jiCancer <- thisConListInt/conceptGeneListUnion[goIDs,geneListInd]

    goLengths <- as.vector(n.genes)[goIDs]

    # |A union B| = |A| + |B| - |A intersect B|
    thisUnionMatrix <- matrix(goLengths,nrow=length(goLengths),ncol=length(goLengths),byrow=T) +
    		       matrix(goLengths,nrow=length(goLengths),ncol=length(goLengths),byrow=F) -
		       intMatrix[goIDs,goIDs]

    cc.ji <- intMatrix[goIDs,goIDs] / thisUnionMatrix ##!!!Chi: this is where the new parameter should come in; please put a "if" condition here, 
	# so that if intMatrix[goIDs,goIDs] / thisUnionMatrix is lower than 0.05, it will not be added to cc.ji
	myArray<-as.numeric(cc.ji)
	myArray[myArray<=0.05]<-0
	cc.ji<-matrix(myArray,nrow=length(goLengths),ncol=length(goLengths))

    # this is for when intMatrix is triangular (either upper or lower)
    # we always double count exactly the diagonal, and the diagonal is always exactly 1, so subtract 1 from our score
    # cc.ji is (intersection/union), and the diagonal compares a GO Term with itself, so |intersection| = |union| = length
    row.cc.ji <- rowSums(cc.ji) + colSums(cc.ji) - 1 #Chi: This is the epsilon.

    # this is for when intMatrix is symmetric
    # row.cc.ji <- rowSums(cc.ji)

    conceptScore <- jiCancer / row.cc.ji
    denom <- 1 / row.cc.ji
    conceptScore[jiCancer==0] <- 0

    # if there's concepts in this gene (may be fiddled with if the user filters GO)
    sumDenomScore <- sum(denom)
    if ((length(eg2go[[geneIDInd]]) > 0) && (sumDenomScore > 0)) {
      sumConceptScore <- sum(conceptScore)
      score <- sumConceptScore / sqrt(sumDenomScore)
    } else {
      score <- 0; sumDenomScore <- 0; sumConceptScore <- 0;
    }

    # record the contribution that each concept is giving, in order of largest contribution to smallest
    conceptOrder <- order(conceptScore,decreasing=T)
    conceptOutputSoFar <- paste(conceptScore[conceptOrder[1]], labels(go2eg)[goIDs[conceptOrder[1]]], sep=':')
    for (i in 2:length(conceptScore)) {
      conceptOutputSoFar <- paste(conceptOutputSoFar, paste(conceptScore[conceptOrder[i]], labels(go2eg)[goIDs[conceptOrder[i]]], sep=':'), sep=';')
    }
	test<-as.character(entrez2sym[which(as.numeric(entrez2sym[,1])==as.numeric(labels(eg2go)[geneIDInd])),2])
    if(length(test)==0){
		print(paste('No geneID found in entrez2sym for geneID',labels(eg2go)[geneIDInd]))
	}else{
		output.genelist <- append(output.genelist, paste(geneListLabels[geneListInd,1],geneListLabels[geneListInd,2],sep=':'))
		output.genelist.size <- append(output.genelist.size, length(geneList[[geneListInd]]))
		output.geneSym <- append(output.geneSym, as.character(entrez2sym[which(as.numeric(entrez2sym[,1])==as.numeric(labels(eg2go)[geneIDInd])),2]))
		output.geneID <- append(output.geneID, labels(eg2go)[geneIDInd])
		output.gene.numConcepts <- append(output.gene.numConcepts, length(goIDs))
		output.gene.time <- append(output.gene.time, as.character((proc.time() - each.time)[3]))
		output.inputStatus <- append(output.inputStatus, sum(labels(eg2go)[geneIDInd] == geneList[[geneListInd]]))
		output.concept.score.sum <- append(output.concept.score.sum, sumConceptScore)
		output.penal <- append(output.penal, sumDenomScore)
		output.score <- append(output.score, score)
		output.conceptsWithScore <- append(output.conceptsWithScore, conceptOutputSoFar)
	}
    }#if...else
  } # for (geneIDInd in eg2go)
} # for (geneListInd)

ourOrder <- order(output.score,decreasing=TRUE)

ourOutput <- data.frame(genelist=output.genelist[ourOrder],
		genelist.size=output.genelist.size[ourOrder],
		geneSym=output.geneSym[ourOrder],
		entrezID=output.geneID[ourOrder],
		inGenelist=output.inputStatus[ourOrder],
		raw.score=output.concept.score.sum[ourOrder],
		penal.score=output.penal[ourOrder],
		score=output.score[ourOrder],
		#concepts=output.conceptsWithScore[ourOrder],
		numConcepts=output.gene.numConcepts[ourOrder]#,
		#time.secs=output.gene.time[ourOrder]
		)

save(ourOutput, file=paste(output.file, '.RData', sep=''))

write.table(ourOutput, file=paste(output.file, '.txt', sep=''), quote=F, sep='\t', row.names=F)
}


# This script translates text-based gene annotation files to a R-friendly format
# This script takes 2 or 3 command line arguments (the 3rd one is optional)
#  (1) The name of a text file to read in, representing all databases to be used.
#  Should have 3 columns, separated by tabs.
#  The first column should be the database this interaction is drawn from. 
#    This is not used, but is nice to help maintaining the database going forward.
#  The second column should have the (Human readable) name of the concept. Two rows with
#    the same name in this column will be treated as the same concept, even if they are from
#    different databases.
#  The third column should have a Entrez Gene ID (integer).
#  To represent a concept with multiple genes, use multiple rows with the same first 
#    two columns, and a different EGID in each one
#  (2) The name of the output file. Should have a .RData extension. This file will have 4 variables:
#  go2eg: a list of lists. Each element of the outer list represents a concept, and the inner list
#    contains all the EGIDs (as strings) that are annotated with that concept.
#  eg2go: a list of lists. Each element of the outer list represents a gene, and the inner list 
#    contains all the concepts annotating that gene. Note that the EGID is a character label of 
#    the outer list, not the index of eg2go. Thus, we may see something like:
#      labels(eg2go)[13] == "1465", which corresponds to the gene symbol 'CSRP1'.
#  n.genes: a vector of integers representing how many genes are in each concept in go2eg
#    length(n.genes) == length(go2eg)
#  goHumanLabels: a vector of strings, where each string corresponds to the concept in go2eg.
#    length(goHumanLabels) == length(go2eg). These strings are from the second column of (1)
#  go2eg.database: a vector of strings, where each string corresponds to the database that the 
#    concept in go2eg came from.
#  (3) The output of a previous run of this program. If this is used, the contents of (1) will be 
#    appended to the variables found in (3)
#    Aside: This will add to go2eg, but recompute eg2go, so don't use this option frivolously 
#    (i.e. adding terms 1 at a time is a bad idea)

#setwd("D:\\Projects\\ConSig\\R_module_construction")

construct_db <- function(database.file,output.file){
# initialize things (use previous work if supplied, otherwise make an empty list)
if (length(args) == 3) {
  old.file <- args[3]
  print (paste('Loading data from', old.file))
  # we only need the old file for it's go2eg, store it in go2eg.all then chuck the rest
  load(old.file)
  rm(eg2go, n.genes, goHumanLabels)
  
  go2eg.all <- go2eg
} else {
  print('Initializing data structures')
  go2eg.all <- list() # a temporary list
}

print (paste('Reading data from', database.file))
X <- read.table(database.file, sep='\t',quote="")

go2eg <- list()  
go2eg.database.list <- list()
go2eg.all.database <- list()

print ('Converting the concepts to an R-friendly format')
# for each line, add the gene to the concept, and keep track of the database
for (i in 1:dim(X)[1]) {
  if(as.character(X[i,3])==c("-")){
    print(paste("Find illigal geneID in ",i,"th gene"))
    next
  }
  conceptName<-paste(as.character(X[i,1]),"^",as.character(X[i,2]))
  go2eg.all[[conceptName]] <- append(go2eg.all[[conceptName]], as.character(X[i,3]))
  go2eg.all.database[[conceptName]] <- as.character(X[i,1])
}

# compute the number of genes corresponding to each concept
n.genes.big <- sapply(go2eg.all, function(x){ length(unlist(x)) }) #calc # of genes per GO term
goBigLabels <- labels(go2eg.all)

print ('Filtering out concepts with too many or few genes')
# only keep concepts with 5-1500 genes (helps computation, doesn't hurt biological results)
for (i in 1:length(n.genes.big)) {
  thisGO <- goBigLabels[i];
  if (n.genes.big[i] > 4) { #this is where the filter "<1500" was removed
    go2eg[[thisGO]] <- go2eg.all[[i]]
    go2eg.database.list[[thisGO]] <- go2eg.all.database[[i]]
  }
}

# clean up and recompute the number of genes per concept
rm(n.genes.big, go2eg.all, go2eg.all.database, goBigLabels);
n.genes <- sapply(go2eg, function(x){ length(unlist(x)) }) #calc # of genes per GO term

eg2go <- list()
goLabels <- labels(go2eg)

# goHumanLabels <- vector('character', length(go2eg))

print ('Constructing the reverse map for EG->GO')
# now build the reverse map (for each gene, which concepts are annotating it?)
for (i in 1:length(go2eg)) {
  thisGO <- goLabels[i]
#  goHumanLabels[i] <- str_replace_all(toupper(as.character(Term(labels(go2eg)[i]))),"\\s+",'_')
  for (j in 1:length(go2eg[[i]])) {
    thisEG <- go2eg[[i]][[j]]
    eg2go[[as.character(thisEG)]] <- append(eg2go[[as.character(thisEG)]],list(thisGO))
  }
}

go2eg.database <- unlist(go2eg.database.list)

print ('Writing output')
# output results
goHumanLabels <- goLabels
save(eg2go, go2eg, n.genes, goHumanLabels, go2eg.database, file=output.file)
}


# This script takes an annotation file and computes the pairwise intersection between each of the terms,
#   based on the genes each term is annotating.
#
# This script requires 3 command line arguments:
# (1) DB.file is a .RData file containing an annotation database, which must have:
#  'go2eg' (a list of lists)
# Additionally, if it is going to be used for ConSig, it will also need to have:
#  'eg2go' (a list of lists)
#  'n.genes' (a list of integers the same length as go2eg, where n.genes[[i]] <- length(go2eg[[i]])
#  'goHumanLabels' (a character vector the same length as go2eg containing human readable annotations the corresponding annotation)
# (2) an output file where the output will be saved. The final output will be in (2).RData, and
#     a checkpoint file of (2).chk.RData will be created to resume your script if it stops prematurely.
# (3) an integer number of nodes to run on for parallelization

# v2 is trying to parallelize the inner loop


preProcess_db <- function(DB.file,output.file,nodes){
load(DB.file)
library('Matrix')
library('doParallel')
library('methods')

output.intersect<-paste("intersect_",output.file,sep="")
#output.union<-paste("union_",output.file,sep="")
n.Terms <- length(go2eg)

if (file.exists(paste(output.intersect, '.chk.RData', sep=''))) {
  load(paste(output.intersect, '.chk.RData', sep=''));
  #load(paste(output.union, '.chk.RData', sep=''));
  fullMatrix <- as(intMatrix, 'matrix');
  #myUnionMatrix <- as(myUnionIntMatrix, 'matrix');
  startIter <- i+1;
} else {
  fullMatrix <- matrix(0, nrow=n.Terms, ncol=n.Terms)
  #myUnionMatrix <- matrix(1, nrow=n.Terms, ncol=n.Terms)
  totalTime <- 0;
  startIter <- 1;
}

hundreth <- floor(n.Terms / 100)

registerDoParallel(cores=nodes)

for (i in startIter:n.Terms) {
  ptm <- proc.time()
  hold <- vector('numeric', n.Terms - i + 1)
  # this takes MUCH longer if you don't store in an intermediate vector. I think it has to do with 
  hold <- foreach (j=i:n.Terms, .combine=cbind) %dopar% {
    hold[i] <- length(intersect(go2eg[[i]], go2eg[[j]]))
#    fullMatrix[j,i] <- fullMatrix[i,j] # uncomment this line to build a symmetric, as opposed to upper triangular matrix
# If the above line is uncommented, changes will need to be made to ConSig.R to ensure that its assumptions about its inputs are correct
  }
  #myUnion <- vector('numeric', n.Terms - i + 1)
  #myUnion <- foreach (j=i:n.Terms, .combine=cbind) %dopar% {
  #  myUnion[i] <- length(union(go2eg[[i]], go2eg[[j]]))
  #}
  crunch.time <- as.numeric(proc.time() - ptm)[3]

  ptm <- proc.time()
  fullMatrix[i,i:n.Terms] <- hold;
  #myUnionMatrix[i,i:n.Terms] <- myUnion;
  copy.time <- as.numeric(proc.time() - ptm)[3]

#  if ((i %% hundreth) == 0) {
    print(paste(i, 'of', n.Terms, 'terms:', n.genes[i], 'genes. Took', round(crunch.time,2), 'seconds to compute,', round(copy.time), 'seconds to copy.'))
#  }
  totalTime <- totalTime + crunch.time + copy.time
  if ((i %% 50) == 0) {
    print('Converting to sparse and saving')
    ptm <- proc.time()
    intMatrix <- as(fullMatrix, 'sparseMatrix')
    #myUnionIntMatrix <- as(myUnionMatrix, 'sparseMatrix')
    print(paste('Converting took', (proc.time()-ptm)[3], 'seconds'))
    ptm <- proc.time()
    save(intMatrix, i, totalTime, file=paste(output.intersect, '.chk.RData', sep=''))
    #save(myUnionIntMatrix, i, totalTime, file=paste(output.union, '.chk.RData', sep=''))
    print(paste('Saving took', (proc.time()-ptm)[3], 'seconds'))
  }
}

print ('Done with all intersections, saving and exiting')
intMatrix <- as(fullMatrix, 'sparseMatrix')
#myUnionIntMatrix <- as(myUnionMatrix, 'sparseMatrix')
save(intMatrix, file=paste(output.intersect, '.RData', sep=''))
#save(myUnionIntMatrix, file=paste(output.union, '.RData', sep=''))
}
