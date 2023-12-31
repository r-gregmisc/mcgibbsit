% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.mcmc.R
\name{read.mcmc}
\alias{read.mcmc}
\title{Read in data from a set of MCMC runs}
\usage{
read.mcmc(
  nc,
  sourcepattern,
  ...,
  col.names,
  start = 1,
  end = nrow(tmp)/numComponents * thin,
  thin = 1,
  numComponents = 1
)
}
\arguments{
\item{nc}{Number of MCMC sampler files to read}

\item{sourcepattern}{MCMC data file name pattern.}

\item{...}{Arguments to be passed to \code{read.table} when loading MCMC
sampler data.}

\item{col.names}{Data file column names (optional)}

\item{start, end, thin}{See documentation for \code{mcmc}}

\item{numComponents}{Number of component samplers.}
}
\value{
An mcmc.list object containing \code{nc} component \code{mcmc}
objects.
}
\description{
Read in data from a set of MCMC runs and create an \code{mcmc.list} object.
}
\details{
This function reads in the states output from one or more MCMC samplers and
creates a single \code{mcmc.list} object.  \code{sourcepattern} will be used
as a filename pattern with \code{#} replaced by the sampler number.  EG,
\code{sourcepattern="MCMC.#.csv"} will be converted to "MCMC.1.csv",
"MCMC.2.csv", etc.

The function \code{read.table} is used to read in the data.  Options for
\code{read.table} may be included as part of the call to \code{read.mcmc}.

The \code{start}, \code{end}, and \code{thin} arguments can be used to
annotate the MCMC samplers with additional information.
}
\examples{

###
# Create example data files for 20 independent chains
# with serial correlation of 0.
###

set.seed(42)
tmpdir <- tempdir()

for(i in 1:20){
  x <- matrix(rnorm(1000), ncol=4)
  
  x[,4] <- x[,4] + 1/3 * (x[,1] + x[,2] + x[,3])
  
  colnames(x) <- c("alpha","beta","gamma", "nu")
  
  write.table(
    x,
    file = file.path(
      tmpdir,
      paste("mcmc", i, "csv", sep=".")
      ),
    sep = ",",
    row.names=FALSE
  )
}

# Read them back in as an mcmc.list object
data <- read.mcmc(20, file.path(tmpdir, "mcmc.#.csv"), sep=",")

summary(data)
}
\seealso{
\code{\link[coda]{mcmc}}, \code{\link[coda]{mcmc.list}},
\code{\link[utils]{read.table}}
}
\author{
Gregory R. Warnes \email{greg@warnes.net}
}
\keyword{file}
