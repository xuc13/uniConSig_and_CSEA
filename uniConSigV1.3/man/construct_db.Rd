\name{construct_db}
\alias{construct_db}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{convert concept database to r friendly objects
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
construct_db(database.file, output.file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{database.file}{
%%     ~~Describe \code{database.file} here~~
}
  \item{output.file}{
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
function (database.file, output.file) 
{
    if (length(args) == 3) {
        old.file <- args[3]
        print(paste("Loading data from", old.file))
        load(old.file)
        rm(eg2go, n.genes, goHumanLabels)
        go2eg.all <- go2eg
    }
    else {
        print("Initializing data structures")
        go2eg.all <- list()
    }
    print(paste("Reading data from", database.file))
    X <- read.table(database.file, sep = "\t", quote = "")
    go2eg <- list()
    go2eg.database.list <- list()
    go2eg.all.database <- list()
    print("Converting the concepts to an R-friendly format")
    for (i in 1:dim(X)[1]) {
        if (as.character(X[i, 3]) == c("-")) {
            print(paste("Find illigal geneID in ", i, "th gene"))
            next
        }
        conceptName <- paste(as.character(X[i, 1]), "^", as.character(X[i, 
            2]))
        go2eg.all[[conceptName]] <- append(go2eg.all[[conceptName]], 
            as.character(X[i, 3]))
        go2eg.all.database[[conceptName]] <- as.character(X[i, 
            1])
    }
    n.genes.big <- sapply(go2eg.all, function(x) {
        length(unlist(x))
    })
    goBigLabels <- labels(go2eg.all)
    print("Filtering out concepts with too many or few genes")
    for (i in 1:length(n.genes.big)) {
        thisGO <- goBigLabels[i]
        if (n.genes.big[i] > 4) {
            go2eg[[thisGO]] <- go2eg.all[[i]]
            go2eg.database.list[[thisGO]] <- go2eg.all.database[[i]]
        }
    }
    rm(n.genes.big, go2eg.all, go2eg.all.database, goBigLabels)
    n.genes <- sapply(go2eg, function(x) {
        length(unlist(x))
    })
    eg2go <- list()
    goLabels <- labels(go2eg)
    print("Constructing the reverse map for EG->GO")
    for (i in 1:length(go2eg)) {
        thisGO <- goLabels[i]
        for (j in 1:length(go2eg[[i]])) {
            thisEG <- go2eg[[i]][[j]]
            eg2go[[as.character(thisEG)]] <- append(eg2go[[as.character(thisEG)]], 
                list(thisGO))
        }
    }
    go2eg.database <- unlist(go2eg.database.list)
    print("Writing output")
    goHumanLabels <- goLabels
    save(eg2go, go2eg, n.genes, goHumanLabels, go2eg.database, 
        file = output.file)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
