##' \code{inner_control} allows specification of control parameters for the
##' inner optimization used by \code{mpmm}
##'
##' See \code{\link{MakeADFun}} and \code{\link{newton}}
##' for details and available options. Adapted from S. Wotherspoon
##' \url{https://github.com/SWotherspoon/RWalc/blob/master/R/RWalc.R}
##' @title Control Values for \code{mpmm}
##' @param ... control parameters for the inner optimizer
##' @return Returns a list with components:
##'  \item{\code{control}}{list of control parameters for inner optimizer}
##' @seealso \code{\link{MakeADFun}}, \code{\link{newton}}
##' @examples
##' fit <- mpmm(~ ice + (ice | id),
##' data = ellie.ice.short,
##' inner.control = inner_control(tol = 1e-03))
##' @export
inner_control <- function(...) {

  dots <- list(...)
  pars <- list(maxit = 1000,
               tol = 1e-08,
               smartsearch = TRUE)

  ## Override control parameters
  pars[names(dots)] <- dots
  list(in_control = pars)
}
