read.mcmc <- function (nc, sourcepattern, ..., col.names,
                         start=1, end=nrow(tmp)/numComponents * thin,
                         thin=1,
                         numComponents=1)
{
    retval <- mcmc.list()
    for (i in 1:nc)
      {
        filename <- sub("#", as.character(i), sourcepattern)
        colnames <- scan(filename, nlines=1, what="", quiet=TRUE, ...)
        tmp <- scan(file=filename, skip=1, what=1, quiet=TRUE, ...)
        tmp<-matrix(tmp, ncol=length(colnames), byrow=TRUE)
        if (missing(col.names))
          colnames(tmp) <- colnames
        else
          colnames(tmp) <- col.names
        for (j in 1:numComponents)
          {
            indices <- seq(j, nrow(tmp), by=numComponents)
            retval[[(i - 1) * numComponents + j]] <-
              mcmc(data=tmp[indices, ], start=start, end=end, thin=thin)
          }
      }
    return(retval)
}

