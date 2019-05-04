#' Data formats
#' @rdname data_formats
#' @aliases data_formats
#' 
#' @description \pkg{STACMR} accepts two kinds of data structure, **list
#'   format** and **general format**. `gen2list` transforms from the general
#'   format into the list format.
#' 
#' @param data `data.frame` or `matrix` in general format (see below).
#' @param varnames optional `list`  of the names of each within-participant
#'   condition.
#' 
#' @details 
#' 
#' \subsection{List Format}{In this format, the data are contained in a `b` x
#' `n` list where `b` is the number of between-participant conditions (groups)
#' and `n` is the number of dependent variables. Each component of the list is
#' itself an `N` x `w` matrix of observations where `N` is the number of
#' subjects (which may vary across groups and dependent variables) and `w` is
#' the number of within-participant conditions (fixed across groups and
#' dependent variables). The dependent variable may be either within-participant
#' or between-participant. This does not matter because the correlation between
#' dependent variables is assumed to be zero (although this might change in
#' future implementations).}
#' 
#' \subsection{General format}{This is a fixed column format organised as a
#' matrix in which each row corresponds to an observation and each column is
#' defined as follows:
#' 1. Participant number (for identification only, not used directly)
#' 2. Between-participant condition or group (if none, then set this value to 1)
#' 3.  Dependent variable (numbered 1, 2, and so on)
#' 4. column 4 to end: Values for each within-participant condition
#' }
#' 
#' @return `gen2list` returns a `ngroup` x `nvar` `list` in which each element
#'   is an `nsub` x `ncond` matrix of values

#' @rdname data_formats
#' @export
gen2list = function (data=NULL, varnames) {
# gen2cell(data)
  # R version of gen2cell.m
# converts data in "general format" to list format suitable for input to staSTATS
# general format is defined as:
  # column 1 = subject number (nsub)
  # column 2 = between-subjects condition (ngroup)
  # column 3 = dependent variable (nvar)
  # columns 4 to end = values for each within-subjects condition (ncond)
  # output is ngroup x nvar list in which each element is an nsub x ncond matrix of values
  #
  # *************************************************************************
  # written 12 September 2016
  # revised 9 March 2017 to remove missing within variables in a group
  # revised 22 August 2017 to add variable names
  # revised 28 February 2019 to repair variable names
  # *************************************************************************
  #
  if (!missing(varnames)) {colnames(data)[4:ncol(data)]=varnames}
  group = data[,2]; ugroup = sort(unique(group)); ngroup = length(ugroup)
  var = data[,3]; uvar = sort(unique(var)); nvar = length(uvar)
  within = as.matrix(data[,4:ncol(data)])
  
  y = vector("list",ngroup)
  for (igroup in 1:ngroup) {
    temp = vector("list", nvar)
    for (ivar in 1:nvar){
      k = which(group==ugroup[igroup] & var==uvar[ivar])
      a = as.matrix(within[k,])
      # delete any variables that all all missing
      n = colSums(is.na(a)); k=which(n==nrow(a)); if (length(k) > 0) {a = a[,-k]}
      # store in 2D list
      y[[igroup]][[ivar]]=a
    }
  }
  return (y)
}
  