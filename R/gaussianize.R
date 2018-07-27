#'  gaussianize null variable importances
#' @param x distr
#' @param a value to interpolate/extrapolate
#' @import stats
#' @import Hmisc


gaussianize=function(x,a){
  q=rank(x, na.last = "keep")/(length(x) + 1 - sum(is.na(x)))
  orig=x
  x2=qnorm(q)
  Hmisc::approxExtrap(orig, x2, xout = a)$y
}
