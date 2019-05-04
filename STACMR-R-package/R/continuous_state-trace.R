#' CMR State-Trace Analysis
#' @rdname continuous_cmr
#' @aliases continuous_cmr
#' 
#' @description 
#' 
#' `staCMR` is the main function that conducts the CMR (state-trace) analysis
#' for continuous data. It takes a data structure or a list of structured output
#' from [staSTATS] and an optional partial order and returns the best fitting
#' values (to the data means) and the least squares fit. It fits the monotonic
#' model to the data.
#' 
#' `staMR` conducts monotonic regression on a data structure according to a given
#' partial order. We say it fits the partial order model to the data (i.e., the
#' set of dependent variables).
#' 
#' `staMRFIT` tests the fit of the partial order model.
#' 
#' `staCMRFIT`  estimates the empirical distribution (and hence *p*-value) of
#' the difference in the fit of the conjoint monotonic and the fit of the
#' partial order model.
#' 
#' 
#' @param  data Either a data structure (in list or general format, see
#'   [gen2list]) or structured output from [staSTATS].
#' @param partial is a partial order in either list or adjacency matrix format.
#' @param  shrink Shrinkage parameter (see [staSTATS]). Default calculates
#'   optimum amount of shrinkage.
#' @param approx `FALSE` (the default) uses full algorithm, `TRUE` uses an
#'   approximate algorithm.
#'   
#' @example examples/examples.delay.R
#' 
#' @export
staCMR <- function (data, partial = list(), shrink=-1, approx=FALSE) {
  # staCMR <- function (data, partial = list(), shrink=-1, approx=F)
  # wrapper function for staCMRx
  # Multidimensional CMR
  # data is cell array of data or structured output from staSTATS 
  # partial is partial order
  # shrink is parameter to control shrinkage of covariance matrix;
  # 0 = no shrinkage; 1 = diagonal matrix; -1 = calculate optimum
  # returns:
  # x = best fitting CMR values to y-means
  # fval = fit statistic
  # shrinkage = estimated shrinkage of covariance matrix
  # approx = F for full algorithm; T = approximate algorithm
  
  # *************************************************************************
  # 21 August 2017
  # modified for approx option 17 September 2018
  # *************************************************************************
  #
  # set up defaults
  tol <- 10e-5
  if (missing(partial)) {partial = list()}
  if (missing(shrink)) {shrink = -1}
  if (missing(approx)) {approx = F}
  
  output = staCMRx (data, model=NULL, E=partial, shrink=shrink, tolerance=0, proc=-1, approx=approx)
  
  return (output)
}


#' @rdname continuous_cmr 
#' @export
staMR <- function(data=list(), partial = list(), shrink=-1) {
  # function [xPrime, fit, shrinkage] = staMR (data, partial, shrink)
  # fits monotonic regression model to data according to partial order
  # data is list of lists of data or structured output from staSTATS
  # partial is partial order in list format
  # shrink is parameter to control shrinkage of covariance matrix;
  # 0 = no shrinkage; 1 = diagonal matrix; -1 = calculate optimum
  # returns:
  #   x = best fitting MR values to y-means
  #   f = total fit statistic
  #   shrinkage = shrinkage from applying staSTATS
  # *************************************************************************
  # modified from matlab 13 September 2016
  # *************************************************************************
  #

  # set up defaults
  tol <- 10e-5
  if (missing(partial)) {partial = list()}
  if (missing(shrink)) {shrink = -1}
  
  # get stats from data (depending on its form)
  if (is(data,"data.frame")) {
    y = gen2list (data) # convert from general format
    y = staSTATS (y, shrink) # get stats
  } else if (is.null(data[[1]]$means)) {y = staSTATS(data, shrink) # in list form, get stats
  } else {y = data} # already in stats form
  
  # convert partial order to list if in adjacency matrix form
  if (is(partial,"matrix")) {partial = adj2list(partial)}
  
  # extract shrinkage parameters (for information only)
  nvar = length(y)
  shrinkage = matrix(0, length(y[[1]]$shrinkage), nvar)
  for (ivar in 1:nvar) {shrinkage[,ivar] = y[[ivar]]$shrinkage}

  # do MR for each dependent variable
  xPrime = vector("list", nvar)
  fit = matrix(0, nvar, 1)
  for (ivar in 1:nvar) {
    out = jMR (y[[ivar]]$means, y[[ivar]]$weights, partial)
    xPrime[[ivar]] = out$x
    fit[ivar] = out$fval
  }
  fval = sum(fit)
  if (fval < tol) {fval = 0} # round down
  
  for (i in 1:nvar) {xPrime[[i]]=matrix(xPrime[[i]],length(xPrime[[i]]),1)}
  output = list(xPrime, fval, shrinkage)
  names(output) = c("x", "fval", "shrinkage")
  
  return(output)
}

#' @rdname continuous_cmr 
#' @param nsample no. of Monte Carlo samples (about 10000 is good)
#' @export
staMRFIT <- function (data=NULL, partial = list(), nsample=1, shrink=-1) {
# input:
  # nsample = no. of Monte Carlo samples (about 10000 is good)
  # data = data structure (cell array or general)
  # partial = optional partial order model e.g. E={[1 2] [3 4 5]} indicates that
  # condition 1 <= condition 2 and condition 3 <= condition 4 <= condition 5
  # default = none (empty)
  # shrink is parameter to control shrinkage of covariance matrix (if input is not stats form);
  # 0 = no shrinkage; 1 = diagonal matrix; -1 = calculate optimum, default = -1
# output:
  # p = empirical p-value
  # datafit = observed fit of partial order model
  # fits = nsample vector of fits of Monte Carlo samples (it is against this
  # distribution that datafit is compared to calculate p)
  # *************************************************************************
  # converted from matlab 7 February 2018
  # *************************************************************************

  if (is(data,"data.frame")) {
    y = gen2list (data) # convert from general format to list format
  } else {y = data} 
  
  tol <- 10e-6
  if (missing(partial)) {partial = list()}
  if (missing(shrink)) {shrink = -1}
  
  nvar =length(y)
  if (!is.list(partial)) {partial = adj2list(partial)} # convert from adjacency matrix to list
  
  output = jMRfits(nsample, y, partial, shrink);
  #output = jCMRfitsx(nsample, y, model, partial, shrink) # call java program
  
  output$fits[which(output$fits <= tol)] = 0;
  
  return (output)
}

#' @rdname continuous_cmr 
#' @export
staCMRFIT <- function (data=NULL, partial = list(), nsample=1, shrink=-1, approx=FALSE) {
# input:
  # nsample = no. of Monte Carlo samples (about 10000 is good)
  # data = data structure (cell array or general)
  # model is a nvar * k matrix specifying the linear model, default = ones(nvar,1))
  # partial = optional partial order model e.g. E={[1 2] [3 4 5]} indicates that
  # condition 1 <= condition 2 and condition 3 <= condition 4 <= condition 5
  # default = none (empty)
  # shrink is parameter to control shrinkage of covariance matrix (if input is not stats form);
  # 0 = no shrinkage; 1 = diagonal matrix; -1 = calculate optimum, default = -1
  # approx = approximation algorithm; F = no; T = yes
# output:
  # p = empirical p-value
  # datafit = observed fit of monotonic (1D) model
  # fits = nsample vector of fits of Monte Carlo samples (it is against this
  # distribution that datafit is compared to calculate p)
  # *************************************************************************
  # converted from matlab 18 September 2016
  # approximation option added 17 September 2018
  # *************************************************************************

  if (is(data,"data.frame")) {
    y = gen2list (data) # convert from general format to list format
  } else {y = data} 
  
  tol <- 10e-6
  if (missing(partial)) {partial = list()}
  if (missing(shrink)) {shrink = -1}
  if (missing(approx)) {approx = 0}
  proc = -1
  cheapP = F
  mrTol = 0
  seed = -1
  
  nvar =length(y[[1]])
  model = NULL
  if (missing(model) | is.null(model)) {model = matrix(1,nvar,1)} # sta default model
  
  if (!is.list(partial)) {partial = adj2list(partial)} # convert from adjacency matrix to list
  
  output = jCMRfitsx(nsample, y, model, partial, shrink, proc, cheapP, approx, mrTol, seed) # call java program
  
  output$fits[which(output$fits <= tol)] = 0;
  
  return (output)
}

