read.mcmc <- function( nc, sourcepattern, ..., col.names, start = 1, end = nrow(tmp), thin = 1)
  {
    retval <- mcmc.list()
    for(i in 1:nc)
      {
        
        filename <- sub("#",as.character(i),sourcepattern)
        
        tmp <-  read.table( file=filename, ...  )

        if(!missing(col.names)) colnames(tmp) <- col.names

        retval[[i]] <- mcmc(data=as.matrix(tmp),
                            start=start, end=end, thin=thin)

      }
    
    return( retval )
  }
