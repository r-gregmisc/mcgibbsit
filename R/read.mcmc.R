read.mcmc <- function( nc, sourcepattern, ..., col.names, start = 1, 
                       end = nrow(tmp)/numComponents*thin, 
                       thin = 1, numComponents=1)
  {
    retval <- mcmc.list()
    for(i in 1:nc)
      {
        filename <- sub("#",as.character(i),sourcepattern)
        tmp <-  read.table( file=filename, ...  )
        if(!missing(col.names)) colnames(tmp) <- col.names
	for (j in 1:numComponents) {
	  indices <- seq(j,nrow(tmp),by=numComponents)
          retval[[(i-1)*numComponents+j]] <- mcmc(data=as.matrix(tmp[indices,]), 
		    start=start, end=end, thin=thin)
	}
      }
    return( retval )
  }
