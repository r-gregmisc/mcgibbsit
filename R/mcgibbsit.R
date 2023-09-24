#' Warnes and Raftery's MCGibbsit MCMC diagnostic
#' 
#' mcgibbsit provides an implementation of Warnes & Raftery's MCGibbsit
#' run-length diagnostic for a set of (not-necessarily independent) MCMC
#' samplers.  It combines the estimate error-bounding approach of Raftery and
#' Lewis with the between chain variance verses within chain variance approach
#' of Gelman and Rubin.
#' 
#' 
#' \code{mcgibbsit} computes the minimum run length \eqn{N_{min}}{Nmin},
#' required burn in \eqn{M}, total run length \eqn{N}, run length inflation due
#' to \emph{auto-correlation}, \eqn{I}, and the run length inflation due to
#' \emph{between-chain} correlation, \eqn{R} for a set of exchangeable MCMC
#' simulations which need not be independent.
#' 
#' The normal usage is to perform an initial MCMC run of some pre-determined
#' length (e.g., 300 iterations) for each of a set of \eqn{k} (e.g.,
#' \eqn{k=20}) MCMC samplers.  The output from these samplers is then read in
#' to create an \code{mcmc.list} object and \code{mcgibbsit} is run specifying
#' the desired accuracy of estimation for quantiles of interest.  This will
#' return the minimum number of iterations to achieve the specified error
#' bound.  The set of MCMC samplers is now run so that the total number of
#' iterations exceeds this minimum, and \code{mcgibbsit} is again called.  This
#' should continue until the number of iterations already complete is less than
#' the minimum number computed by \code{mcgibbsit}.
#' 
#' If the initial number of iterations in \code{data} is too small to perform
#' the calculations, an error message is printed indicating the minimum pilot
#' run length.
#' 
#' The parameters \code{q}, \code{r}, \code{s}, \code{converge.eps}, and
#' \code{correct.cor} can be supplied as vectors.  This will cause
#' \code{mcgibbsit} to produce a list of results, with one element produced for
#' each set of values.  I.e., setting \code{q=(0.025,0.975), r=(0.0125,0.005)}
#' will yield a list containing two \code{mcgibbsit} objects, one computed with
#' parameters \code{q=0.025, r=0.0125}, and the other with \code{q=0.975,
#' r=0.005}.
#' 
#' @aliases mcgibbsit print.mcgibbsit
#'
#' @param data an `mcmc' object.
#' @param q quantile(s) to be estimated.
#' @param r the desired margin of error of the estimate.
#' @param s the probability of obtaining an estimate in the interval
#' @param converge.eps Precision required for estimate of time to convergence.
#' @param correct.cor should the between-chain correlation correction (R) be
#' computed and applied.  Set to false for independent MCMC chains.
#' @inheritParams base::print
#'
#' @return An \code{mcgibbsit} object with components
#' 
#' \item{call}{parameters used to call 'mcgibbsit'} \item{params}{values of r,
#' s, and q used} \item{resmatrix}{a matrix with 6 columns: \describe{
#' \item{Nmin}{The minimum required sample size for a chain with no correlation
#' between consecutive samples. Positive autocorrelation will increase the
#' required sample size above this minimum value.} \item{M}{The number of `burn
#' in' iterations to be discarded (total over all chains).} \item{N}{The number
#' of iterations after burn in required to estimate the quantile q to within an
#' accuracy of +/- r with probability p (total over all chains).}
#' \item{Total}{Overall number of iterations required (M + N).} \item{I}{An
#' estimate (the `dependence factor') of the extent to which auto-correlation
#' inflates the required sample size.  Values of `I' larger than 5 indicate
#' strong autocorrelation which may be due to a poor choice of starting value,
#' high posterior correlations, or `stickiness' of the MCMC algorithm.}
#' \item{R}{An estimate of the extent to which between-chain correlation
#' inflates the required sample size.  Large values of 'R' indicate that there
#' is significant correlation between the chains and may be indicative of a
#' lack of convergence or a poor multi-chain algorithm.} } } \item{nchains}{the
#' number of MCMC chains in the data} \item{len}{the length of each chain}
#'
#' @author Gregory R. Warnes \email{greg@@warnes.net} based on the the R
#' function \code{raftery.diag} which is part of the 'CODA' library.
#' \code{raftery.diag}, in turn, is based on the FORTRAN program `gibbsit'
#' written by Steven Lewis which is available from the Statlib archive.
#'
#' @seealso \code{\link{read.mcmc}}
#'
#' @references
#' 
#' Warnes, G.W. (2004). The Normal Kernel Coupler: An adaptive MCMC method for
#' efficiently sampling from multi-modal distributions,
#' \url{https://digital.lib.washington.edu/researchworks/handle/1773/9541}
#' 
#' Warnes, G.W. (2000).  Multi-Chain and Parallel Algorithms for Markov Chain
#' Monte Carlo. Dissertation, Department of Biostatistics, University of
#' Washington,
#' \url{https://stat.uw.edu/sites/default/files/files/reports/2001/tr395.pdf}
#' 
#' Raftery, A.E. and Lewis, S.M. (1992).  One long run with diagnostics:
#' Implementation strategies for Markov chain Monte Carlo. Statistical Science,
#' 7, 493-497.
#' 
#' Raftery, A.E. and Lewis, S.M. (1995).  The number of iterations, convergence
#' diagnostics and generic Metropolis algorithms.  In Practical Markov Chain
#' Monte Carlo (W.R. Gilks, D.J. Spiegelhalter and S. Richardson, eds.).
#' London, U.K.: Chapman and Hall.
#' 
#' @keywords models
#' 
#' @importFrom coda is.mcmc.list
#' @importFrom coda mcmc
#' @importFrom coda mcmc.list
#' @importFrom coda niter
#' @importFrom coda nvar
#' @importFrom coda thin
#' @importFrom coda varnames
#' @importFrom stats end
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @importFrom stats start
#' @importFrom stats var
#' @importFrom stats window
#'
#' @examples
#' 
#' ###
#' # Create example data files for 20 independent chains
#' # with serial correlation of 0.25
#' ###
#' 
#' set.seed(42)
#' tmpdir <- tempdir()
#' 
#' nsamples <- 1000
#' 
#' for(i in 1:20){
#'   x <- matrix(nrow = nsamples+1, ncol=4)
#'   colnames(x) <- c("alpha","beta","gamma", "nu")
#'   
#'   x[,"alpha"] <- rnorm (nsamples+1, mean=0.025, sd=0.0025)^2
#'   x[,"beta"]  <- rnorm (nsamples+1, mean=53,    sd=12)
#'   x[,"gamma"] <- rbinom(nsamples+1, 20,         p=0.25) + 1
#'   x[,"nu"]    <- rnorm (nsamples+1, mean=x[,"alpha"] * x[,"beta"], sd=1/x[,"gamma"])
#'
#'   # induce serial correlation of 0.25
#'   x <- 0.75 * x[2:(nsamples+1),] + 0.25 * x[1:nsamples,]
#'   
#'   
#'   write.table(
#'     x,
#'     file = file.path(
#'       tmpdir,
#'       paste("mcmc", i, "csv", sep=".")
#'       ),
#'     sep = ",",
#'     row.names = FALSE
#'   )
#' }
#'
#' # Read them back in as an mcmc.list object
#' data <- read.mcmc(
#'   20, 
#'   file.path(tmpdir, "mcmc.#.csv"), 
#'   sep=",",
#'   col.names=c("alpha","beta","gamma", "nu")
#'   )
#' 
#' # Summary statistics
#' summary(data)
#'
#' # Trace and Density Plots
#' plot(data)
#'
#' # And check the necessary run length 
#' mcgibbsit(data)
#'
#'
#' @export
"mcgibbsit" <- function(
    data,
    q = 0.025,
    r = 0.0125,
    s = 0.95,
    converge.eps = 0.001,
    correct.cor=TRUE
    )
{
  ## if asked for more than one q,r,s,..., combination, recurse and return
  ## a list of results
  
  parms <- cbind( q, r, s, converge.eps, correct.cor ) # so one can use
  # different r/s values
  # for each q
  
  if(dim(parms)[1] > 1 )
  {
    
    retval <- list();
    for(i in 1:length(q) )
    {
      retval[[i]] <- mcgibbsit(
        data, 
        q=parms[i,"q"], 
        r=parms[i,"r"],
        s=parms[1,"s"], 
        converge.eps=parms[i,"converge.eps"], 
        correct.cor=parms[i,"correct.cor"] 
        )
    }
    return(retval)
  }
  
  if (is.mcmc.list(data))
  {
    nchains <- length(data)
    
    ## check that all of the chains are conformant #
    for( ch in 1:nchains)
      if(
        (start(data[[1]]) != start(data[[ch]] )) ||
        (end  (data[[1]]) != end  (data[[ch]] )) ||
        (thin (data[[1]]) != thin (data[[ch]] )) ||
        (nvar (data[[1]]) != nvar (data[[ch]] )) 
      )
        stop(paste("All chains in mcmc.list must have same 'start',",
                   "'end', 'thin', and number of variables"));
    
    if(is.matrix(data[[1]]))
      combined <- mcmc(do.call("rbind",data))
    else
      combined <- mcmc(as.matrix(unlist(data)))
  }
  else
  {
    data <- mcmc.list(mcmc(as.matrix(data)))
    nchains <- 1
    combined <- mcmc(as.matrix(data))
  }
  
  
  resmatrix <- matrix(nrow = nvar(combined),
                      ncol = 6,
                      dimnames = list( varnames(data, allow.null = TRUE),
                                       c("M", "N", "Total", "Nmin",
                                         "I", "R" )) )
  
  # minimum number of iterations
  phi <- qnorm(0.5 * (1 + s))
  nmin <- as.integer(ceiling((q * (1 - q) * phi^2)/r^2))
  
  if (nmin > niter(combined))
    resmatrix <- c("Error", nmin)
  else
    for (i in 1:nvar(combined))
    {
      dichot <- list()
      
      if (is.matrix(data[[1]]))
      {
        quant <- quantile(combined[, i, drop = TRUE], probs = q)
        
        for(ch in 1:nchains)
        {
          dichot[[ch]] <- mcmc(data[[ch]][, i, drop = TRUE] <= quant,
                               start = start(data[[ch]]),
                               end   = end(data[[ch]]),
                               thin  = thin(data[[ch]]))
        }
      }
      else
      {
        quant <- quantile(combined, probs = q)
        
        for(ch in 1:nchains)
        {
          dichot[[ch]] <- mcmc(data[[ch]] <= quant,
                               start = start(data[[ch]]),
                               end = end(data[[ch]]),
                               thin = thin(data[[ch]]))
        }
      }
      kthin <- 0
      bic <- 1
      
      while (bic >= 0)
      {
        kthin <- kthin + thin(data[[1]])
        
        to.table <- function(dichot, kthin)
        {
          testres <- as.vector(window(dichot, thin = kthin))
          newdim <- length(testres)
          testtran <- table(testres[1:(newdim - 2)],
                            testres[2:(newdim - 1)],
                            testres[3:newdim])
          
          ## handle the case where one or more of the transition never
          ## happens, so that the testtran array has two few dimensions
          if( any(dim(testtran!=2)) )
          {
            tmp <- array( 0, dim=c(2,2,2),
                          dimnames=list( c("FALSE","TRUE"),
                                         c("FALSE","TRUE"),
                                         c("FALSE","TRUE") ) )
            for(t1 in dimnames(testtran)[[1]] )
              for(t2 in dimnames(testtran)[[2]])
                for(t3 in dimnames(testtran)[[3]] )
                  tmp[t1,t2,t3] <- testtran[t1,t2,t3]
            
            testtran <- tmp
          }
          
          testtran <- array(as.double(testtran), dim = dim(testtran))
          return(testtran)
        }
        
        tmp <- sapply( dichot, to.table, kthin=kthin, simplify=FALSE )
        
        ## add all of the transition matrixes together
        testtran <- tmp[[1]]
        if(nchains>1)
          for(ch in 2:nchains)
            testtran <- testtran + tmp[[ch]]
        
        ## compute the likelihoood
        g2 <- 0
        for (i1 in 1:2)
        {
          for (i2 in 1:2)
          {
            for (i3 in 1:2)
            {
              if (testtran[i1, i2, i3] != 0) {
                fitted <- (sum(testtran[i1, i2, 1:2]) *
                             sum(testtran[1:2, i2, i3])) /
                  (sum(testtran[1:2, i2, 1:2]))
                g2 <- g2 + testtran[i1, i2, i3] *
                  log(testtran[i1, i2, i3]/fitted) * 2
              }
            }
          }
        }
        
        ## compute bic
        bic <- g2 - log( sum(testtran) - 2 ) * 2
        
      }
      
      ## estimate the parameters of the first-order markov chain
      alpha <- sum(testtran[1, 2, 1:2]) / (sum(testtran[1, 1, 1:2]) +
                                             sum(testtran[1, 2, 1:2]))
      beta  <- sum(testtran[2, 1, 1:2]) / (sum(testtran[2, 1, 1:2]) +
                                             sum(testtran[2, 2, 1:2]))
      
      ## compute burn in
      tempburn <- log((converge.eps * (alpha + beta)) /
                        max(alpha, beta))/(log(abs(1 - alpha - beta)))
      
      nburn <- as.integer(ceiling(tempburn) * kthin)
      
      ## compute iterations after burn in
      tempprec <- ((2 - alpha - beta) * alpha * beta * phi^2)/
        (((alpha + beta)^3) * r^2)
      nkeep  <- ceiling(tempprec * kthin)
      
      ## compute the correlation
      if(nchains>1 && correct.cor)
      {
        dat <- do.call("cbind", dichot)
        varmat <- var(dat)
        denom <- mean(diag(varmat))  # overall variance
        diag(varmat) <- NA
        numer <- mean(c(varmat), na.rm=TRUE) # overall covariance
        rho <- numer / denom
      }
      else
        rho <- 1.0
      
      ## inflation factors
      iratio <- (nburn + nkeep)/nmin
      R      <- ( 1 + rho * (nchains - 1) )
      
      resmatrix[i, 1] <- M <- ceiling( nburn * nchains )  # M
      resmatrix[i, 2] <- N <- ceiling( nkeep * R       )  # N
      resmatrix[i, 3] <- M + N                            # Total
      resmatrix[i, 4] <- nmin                             # nmin
      resmatrix[i, 5] <- signif(iratio,  digits = 3)      # I
      if(nchains > 1 && correct.cor )
        resmatrix[i, 6] <- signif(R     ,  digits = 3)    # R
      else
        resmatrix[i, 6] <- NA                             # R
      
    }
  y <- list(params = c(r = r, s = s, q = q),
            resmatrix = resmatrix, call=match.call(),
            nchains = nchains, len=(end(data[[1]]) - start(data[[1]]) + 1) )
  class(y) <- "mcgibbsit"
  return(y)
}

#' @export
#' @rdname mcgibbsit
"print.mcgibbsit" <-
  function (x, digits = 3, ...)
  {
    cat("                  Multi-Chain Gibbsit \n")
    cat("                  ------------------- \n")
    cat("\n");
    
    cat("Call             = "); print(x$call)
    cat("\n");
    
    cat("Number of Chains =", x$nchains, "\n" )
    cat("Per-Chain Length =", x$len, "\n" )
    cat("Total Length     =", x$nchains * x$len, "\n")
    cat("\n");
    
    cat("Quantile (q)     =", x$params["q"] , "\n")
    cat("Accuracy (r)     = +/-", x$params["r"], "\n")
    cat("Probability (s)  =", x$params["s"], "\n")
    cat("\n")
    
    if (x$resmatrix[1] == "Error")
      cat("\nYou need a sample size of at least", x$resmatrix[2],
          "with these values of q, r and s\n")
    else {
      out <- x$resmatrix
      for (i in ncol(out)) out[, i] <- format(out[, i], digits = digits)
      
      maxM     <- max(x$resmatrix[,"M"])
      maxN     <- max(x$resmatrix[,"N"])
      maxTotal <- max(x$resmatrix[,"Total"])
      minBound <- max(x$resmatrix[,"Nmin"])
      
      out <- rbind(c("Burn-in ", "Estimation", "Total", "Lower bound ",
                     "Auto-Corr.", "Between-Chain"),
                   c("(M)",       "(N)",       "(M+N)", "(Nmin)",
                     "factor (I)", "Corr. factor (R)"),
                   rep('', ncol(out)),
                   out,
                   rep('-----', ncol(out) ),
                   c(maxM, maxN, maxTotal, minBound, "", "")
      )
      
      
      #    if (!is.null(rownames(x$resmatrix)) || all(rownames(x$resmatrix)==''))
      #      out <- cbind(c("", "", rownames(x$resmatrix)), out)
      
      colnames(out) <- rep("", ncol(out))
      
      
      
      print.default(out, quote = FALSE, ...)
      cat("\n")
      
      cat("NOTE: The values for M, N, and Total are combined numbers",
          "of iterations \n")
      cat("      based on using", x$nchains, "chains.\n");
      cat("\n")
      
    }
    invisible(x)
  }
