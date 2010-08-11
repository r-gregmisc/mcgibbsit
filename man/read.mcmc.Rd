\name{read.mcmc}
\alias{read.mcmc}
\title{Read in data from a set of MCMC runs}
\description{
  Read in data from a set of MCMC runs and create an  \code{mcmc.list} object.
}
\usage{
read.mcmc( nc, sourcepattern, ...,
               col.names, start = 1, 
               end = nrow(tmp)/numComponents*thin, 
               thin = 1, numComponents=1)
}
\arguments{
  \item{nc}{Number of MCMC sampler files to read}
  \item{sourcepattern}{MCMC data file name pattern.}
  \item{...}{Arguments to be passed to \code{read.table} when loading
    MCMC sampler data.}
  \item{col.names}{Data file column names (optional)}
  \item{start, end, thin}{See documentation for \code{mcmc}}
  \item{numComponents}{Number of component samplers.}
}
\details{
  This function reads in the states output from one or more MCMC
  samplers and creates a single \code{mcmc.list} object.
  \code{sourcepattern} will be used as a filename pattern with \code{#}
  replaced by the sampler number.  EG, \code{sourcepattern="MCMC.#.csv"}
  will be converted to "MCMC.1.csv", "MCMC.2.csv", etc.

  The function \code{read.table} is used to read in the data.  Options
  for \code{read.table} may be included as part of the call to
  \code{read.mcmc}.

  The \code{start}, \code{end}, and \code{thin} arguments can be used to
  annotate the MCMC samplers with additional information.
}
\value{
  An mcmc.list object containing \code{nc} component \code{mcmc} objects.
}
\author{Gregory R. Warnes \email{gregory\_r\_warnes@groton.pfizer.com}}
\seealso{
  \code{\link[coda]{mcmc}},
  \code{\link[coda]{mcmc.list}},
  \code{\link[utils]{read.table}}
}

\examples{

## Create set of 20 example data files
set.seed(12234)
for(i in 1:20){
  x <- matrix(rnorm(1000),ncol=4)
  x[,4] <- x[,4] + 1/3 * (x[,1] + x[,2] + x[,3])
  colnames(x) <- c("alpha","beta","gamma", "nu")
  write.csv(x, file=paste("mcmc",i,"csv",sep="."), row.names=FALSE)
}

## Read the 20 data files and convert them into a single
## MCMC object.
data <- read.mcmc(20, "mcmc.#.csv", sep=",")

## Clean up after ourselves
unlink( paste("mcmc.",1:20,".csv", sep="") )

## Now run mcgibbsit...
mcgibbsit(data)

}

\keyword{file}
