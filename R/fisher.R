#' fisher
#' @import stats
#' @import Hmisc


#' @param p vector of pvalues
#' @param adjust if correlation has to be taken into account
#' @param R correlation matrix

fisher <- function(p, adjust,R) {

  k <- length(p)
  if (adjust == "no") {
    t <- -2 * sum(log(p))
    poolp <- pchisq(t, df = 2*k, lower.tail = FALSE)

  } else if (adjust == "yes")  {
    m <- m_effective(R = R)

    t <- -2 * sum(log(p))
    t <- t * (m / k)
    poolp <- pchisq(t, df = 2 * m, lower.tail = FALSE)
  }

  return(poolp)

}
