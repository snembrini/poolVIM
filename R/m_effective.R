#'  meff.
#' @param R R
#'
m_effective <- function(R) {

  evs <- eigen(R)$values
    abs.evs <- abs(evs)
  k <- length(evs)
    m <- 1 + (k - 1) * (1 - var(evs) / k)

  m <- floor(m)

  return(m)

}
