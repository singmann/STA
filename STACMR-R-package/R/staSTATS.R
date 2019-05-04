#' Statistics for state trace analysis
#' 
#' Calculates statistics for state trace analysis.
#' 
#' @param data A `list` of `lists` (returned from [gen2list]) containing
#'   `nsub` x `ncond` matrices
#' @param varnames An optional `list` of names of within-participant conditions
#' @param shrink numeric indicating amount of shrinkage to apply to the
#'   estimated covariance matrix. Generally, the covariance matrix needs to be
#'   shrunk during the bootstrap cycle to avoid ill-conditioning. If `shrink =
#'   0` then no shrinkage is applied. If `shrink = 1` then maximum shrinkage is
#'   applied. This means that the covariance matrix is diagonalized with all
#'   off-diagonal entries set to zero. If `shrink < 0` (the default) then an
#'   optimal shrinkage value is estimated for each within-participant block and
#'   applied according to an algorithm developed by Ledoit and Wolf (2004).
#' @param if warning is `TRUE` then a warning is thrown if `NA`s are detected.
#'   Default is `FALSE`.
#' 
#' @references Ledoit, O. & Wolf, M. (2004). Honey, I shrunk the sample
#'   covariance matrix. *The Journal of Portfolio Management*, 30(4), 110-119.
#' 
#' @export
staSTATS <- function(data, shrink=-1, varnames, warning=FALSE) {
  ## Calculates statistics for state trace analysis
  ##
  ## Args:
  ##   data: A list of lists containing nsub x ncond matrices
  ##   varnames: An optional list of names of within-participant conditions
  ##   shrink: 0=no shrinkage of covariance; 1=maximal shrinkage; -1=estimated optimal shrinkage
  ##   if warning is set then a warning message is printed if NAs are detected
  ##
  ## Returns a list with the following items:
  ##   means: Observed means
  ##   n: Number of subjects
  ##   cov: Observed covariance matrix
  ##   recov: Regularized (shrinked) covariance matrix 
  ##   weights: Weight matrix for monotonic regression
  ##   lm: Loftus-Masson within subjects variance
  ##   shrinkage: shrinkage parameter (estimated or returned)
  ##   nanflag: count of missing values (NAs)
  
  #if (missing(shrink)) {shrink = -1}
  #if (missing(warning)) {warning = 0}
  
  y <- data
  
  if (is(y, "data.frame")) {y = gen2list(y, varnames)} # convert from general to list format if req'd
  ngroup = length(y); nvar = length(y[[1]])
  
  output = vector("list", nvar)
  
  if ('means' %in% names(y[[1]])) {output = y # already in stats form
  } else {
    for(ivar in 1:nvar) {
      i.means = numeric()
      i.n = matrix(0,0,0)
      i.cov = matrix(0,0,0)
      i.regcov = matrix(0,0,0)
      i.shrinkage = numeric()
      i.weights = matrix(0,0,0)
      i.lm = matrix(0,0,0)
      i.bad = matrix(0,0,0)
      i.nanflag = matrix(0,0,0)
      
      for (igroup in 1:ngroup) {
        y.i <- as.matrix(y[[igroup]][[ivar]])
 #       y.i <- y.i[complete.cases(y.i), ] ## delete rows with NAs
        g.nanflag <- sum(is.na(y.i))
        g.means <- colMeans(y.i,na.rm=TRUE)
        a=y.i; a[is.na(y.i)]=0; a[!is.na(y.i)]=1; g.n=t(a)%*%a; a=g.n-1; a[a<0]=0; # no. of observations
        g.cov <- cov(y.i,use='pairwise.complete.obs')*a/g.n; g.cov[is.na(g.cov)]=0 # adjusted covariance
        
        eigcov <- eigen(g.cov)
        if ((kappa(g.cov,2) < 1e6) && (min(eigcov$values) > 0))
          {s <- shrinkDiag(y.i, shrink); g.bad <- 0}
       else
          {s <- shrinkDiag(y.i, 1); g.bad <- 1} # diagonalize ill-conditioned matrix
        g.regcov <- s$sigma; g.shrinkage = s$shrinkage
        eigRegcov <- eigen(g.regcov) # check if positive definite
        if ((kappa(g.regcov,2) > 1e6) ||(min(eigRegcov$values) <= 0))
        {g.regcov <- sum(eigRegcov$values)*diag(nrow(g.regcov))/nrow(g.regcov); g.bad <- 2}
        
        Nmin <- g.n; Nmin[which(diag(nrow(g.regcov)) > g.n)] <- 1;
        g.weights <- Nmin*solve(g.regcov)
        if (length(g.regcov)==1) {g.lm=g.regcov} else {g.lm = LoftusMasson (y.i)}
        
        # add to vectors and matrices
        i.means = c(i.means, g.means)
        i.n = magic::adiag(i.n, g.n) # requires "magic" package
        i.cov = magic::adiag(i.cov, g.cov)
        i.regcov = magic::adiag(i.regcov, g.regcov)
        i.shrinkage = c(i.shrinkage, g.shrinkage)
        i.weights = magic::adiag(i.weights, g.weights)
        i.lm = magic::adiag(i.lm, g.lm)
        i.bad = c(i.bad, g.bad)
        i.nanflag = c(i.nanflag, g.nanflag)
      }
      # add to list
      out = list(i.means, i.n, i.cov, i.regcov, i.shrinkage, i.weights, i.lm, i.nanflag, i.bad)
      output[[ivar]] = out
      names(output[[ivar]]) = c("means", "n", "cov", "regcov", "shrinkage", "weights", "lm", "nanflag", "bad")
      
      if (sum(i.nanflag > 0) && (warning))
        {warning(sum(i.nanflag), ' detected for variable',ivar)}
      if (sum(i.bad > 0) && (warning))
        {warning('Bad covariance matrix detected for variable',ivar, '. Type = ', i.bad)}
    }
  }
  return(output)
}

