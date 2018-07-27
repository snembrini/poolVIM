#' tippett.
#' @import stats
#' @import Hmisc


#' @param p vector of pvalues
#' @param adjust if correlation has to be taken into account
#' @param R correlation matrix

tippett <- function(p, adjust, R) {

  k <- length(p)
    if (adjust == "no") {
      t <- 1 - (1 - min(p))^k
      poolp <- t

    } else if (adjust == "yes")  {
      m <- m_effective(R = R)
      t <- 1 - (1 - min(p))^m
      poolp <- t
    }

  return(poolp)

}
