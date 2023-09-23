#' Read in data from a set of MCMC runs
#'
#' Read in data from a set of MCMC runs and create an \code{mcmc.list} object.
#'
#' This function reads in the states output from one or more MCMC samplers and
#' creates a single \code{mcmc.list} object.  \code{sourcepattern} will be used
#' as a filename pattern with \code{#} replaced by the sampler number.  EG,
#' \code{sourcepattern="MCMC.#.csv"} will be converted to "MCMC.1.csv",
#' "MCMC.2.csv", etc.
#'
#' The function \code{read.table} is used to read in the data.  Options for
#' \code{read.table} may be included as part of the call to \code{read.mcmc}.
#'
#' The \code{start}, \code{end}, and \code{thin} arguments can be used to
#' annotate the MCMC samplers with additional information.
#'
#' @param nc Number of MCMC sampler files to read
#' @param sourcepattern MCMC data file name pattern.
#' @param ... Arguments to be passed to \code{read.table} when loading MCMC
#'   sampler data.
#' @param col.names Data file column names (optional)
#' @param start,end,thin See documentation for \code{mcmc}
#' @param numComponents Number of component samplers.
#'
#' @return An mcmc.list object containing \code{nc} component \code{mcmc}
#'   objects.
#'
#' @author Gregory R. Warnes \email{greg@@warnes.net}
#'
#' @seealso \code{\link[coda]{mcmc}}, \code{\link[coda]{mcmc.list}},
#'   \code{\link[utils]{read.table}}
#'
#' @keywords file
#'
#' @importFrom coda mcmc
#' @importFrom coda mcmc.list
#'
#' @export
#' 
#' @examples
#'
#' ###
#' # Create example data files for 20 independent chains
#' # with serial correlation of 0.
#' ###
#' 
#' set.seed(42)
#' tmpdir <- tempdir()
#' 
#' for(i in 1:20){
#'   x <- matrix(rnorm(1000), ncol=4)
#'   
#'   x[,4] <- x[,4] + 1/3 * (x[,1] + x[,2] + x[,3])
#'   
#'   colnames(x) <- c("alpha","beta","gamma", "nu")
#'   
#'   write.table(
#'     x,
#'     file = file.path(
#'       tmpdir,
#'       paste("mcmc", i, "csv", sep=".")
#'       ),
#'     sep = ",",
#'     row.names=FALSE
#'   )
#' }
#'
#' # Read them back in as an mcmc.list object
#' data <- read.mcmc(20, file.path(tmpdir, "mcmc.#.csv"), sep=",")
#' 
#' summary(data)
read.mcmc <- function (
    nc, 
    sourcepattern, 
    ..., 
    col.names,
    start=1, 
    end=nrow(tmp)/numComponents * thin,
    thin=1,
    numComponents=1
)
{
  retval <- mcmc.list()
  for (i in 1:nc)
  {
    filename <- sub("#", as.character(i), sourcepattern)
    
    ncols <- scan(file=filename, nlines=1, what="", quiet=TRUE, ...)
    tmp   <- scan(file=filename, skip=1,   what=1,  quiet=TRUE, ...)
    
    tmp <- matrix(tmp, ncol=length(ncols), byrow=TRUE)
    
    if (missing(col.names))
      colnames(tmp) <- rep("", length(ncols))
    else
      colnames(tmp) <- col.names
    
    for (j in 1:numComponents)
    {
      indices <- seq(j, nrow(tmp), by=numComponents)
      retval[[(i - 1) * numComponents + j]] <-
        mcmc(data=tmp[indices, ], start=start, end=end, thin=thin)
    }
  }
  
  retval
}

