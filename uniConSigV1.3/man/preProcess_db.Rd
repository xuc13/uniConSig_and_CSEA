\name{preProcess_db}
\alias{preProcess_db}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{construct Jaccard Index matrix
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
preProcess_db(DB.file, output.file, nodes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DB.file}{
%%     ~~Describe \code{DB.file} here~~
}
  \item{output.file}{
%%     ~~Describe \code{output.file} here~~
}
  \item{nodes}{
%%     ~~Describe \code{nodes} here~~
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
function (DB.file, output.file, nodes) 
{
    load(DB.file)
    library("Matrix")
    library("doParallel")
    library("methods")
    output.intersect <- paste("intersect_", output.file, sep = "")
    n.Terms <- length(go2eg)
    if (file.exists(paste(output.intersect, ".chk.RData", sep = ""))) {
        load(paste(output.intersect, ".chk.RData", sep = ""))
        fullMatrix <- as(intMatrix, "matrix")
        startIter <- i + 1
    }
    else {
        fullMatrix <- matrix(0, nrow = n.Terms, ncol = n.Terms)
        totalTime <- 0
        startIter <- 1
    }
    hundreth <- floor(n.Terms/100)
    registerDoParallel(cores = nodes)
    for (i in startIter:n.Terms) {
        ptm <- proc.time()
        hold <- vector("numeric", n.Terms - i + 1)
        hold <- foreach(j = i:n.Terms, .combine = cbind) \%dopar\% 
            {
                hold[i] <- length(intersect(go2eg[[i]], go2eg[[j]]))
            }
        crunch.time <- as.numeric(proc.time() - ptm)[3]
        ptm <- proc.time()
        fullMatrix[i, i:n.Terms] <- hold
        copy.time <- as.numeric(proc.time() - ptm)[3]
        print(paste(i, "of", n.Terms, "terms:", n.genes[i], "genes. Took", 
            round(crunch.time, 2), "seconds to compute,", round(copy.time), 
            "seconds to copy."))
        totalTime <- totalTime + crunch.time + copy.time
        if ((i\%\%50) == 0) {
            print("Converting to sparse and saving")
            ptm <- proc.time()
            intMatrix <- as(fullMatrix, "sparseMatrix")
            print(paste("Converting took", (proc.time() - ptm)[3], 
                "seconds"))
            ptm <- proc.time()
            save(intMatrix, i, totalTime, file = paste(output.intersect, 
                ".chk.RData", sep = ""))
            print(paste("Saving took", (proc.time() - ptm)[3], 
                "seconds"))
        }
    }
    print("Done with all intersections, saving and exiting")
    intMatrix <- as(fullMatrix, "sparseMatrix")
    save(intMatrix, file = paste(output.intersect, ".RData", 
        sep = ""))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
