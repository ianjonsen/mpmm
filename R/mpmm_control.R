##' \code{mpmm_control} selects the numerical minimizer, method, associated
##' control parameters, parameter bounds, and likelihood estimation (REML or ML)
##' used by \code{mpmm}.
##'
##' The optimizer used to minimize the objective function is
##' selected by the \code{optim} argument.  Additional control
##' parameters specific to the chosen optimizer are specified via the
##' dots argument.  See \code{\link{nlminb}} and \code{\link{optim}}
##' for available options. Adapted from S. Wotherspoon
##' \url{https://github.com/SWotherspoon/RWalc/blob/master/R/RWalc.R}
##'
##' @title Control Values for \code{mpmm}.
##' @param optim the numerical optimizer used in the fit
##' @param method optimization method to be used - one of "BFGS" or "L-BFGS-B"
##' for bounded optimization (default \code{lower} and \code{upper} bounds are
##' used if not specified
##' @param lower a list of named parameter lower bounds, if NULL then built in
##' defaults are used when \code{method = "L-BFGS-B"}, otherwise ignored
##' @param upper a list of named parameter upper bounds, if NULL then built in
##' defaults are used when \code{method = "L-BFGS-B"}, otherwise ignored
##' @param REML logical; whether to use REML (TRUE) or maximum likelihood
##' @param profile logical; option to improve speed and convergence by using
##' REML parameter estimates as initial values for ML optimization
##' @param verbose integer; report progress during minimization: 0 = silent;
##' 1 = optimizer trace; 2 = parameter trace (default))
##' @param ... control parameters for the chosen optimizer
##' @return Returns a list with components
##'   \item{\code{optim}}{the name of the numerical optimizer as a
##'   string, "nlminb" or "optim"}
##'   \item{\code{method}}{optimization method to be used}
##'   \item{\code{lower}}{named list of lower parameter bounds}
##'   \item{\code{upper}}{named list of upper parameter bounds}
##'   \item{\code{REML}}{whether REML is to be used in place of ML}
##'   \item{\code{profile}}{whether to enhance convergence robustness}
##'   \item{\code{verbose}}{level of tracing information to be reported}
##'   \item{\code{control}}{list of control parameters for the optimizer}
##' @seealso \code{\link{nlminb}}, \code{\link{optim}}.
##' @example
##' fit <- mpmm(~ ice + (ice | id),
##' data = ellie.ice.short,
##' control = mpmm_control(
##'     optim = "nlminb",
##'     REML = FALSE,
##'     eval.max = 2000)
##'     )
##' @export

mpmm_control <-
  function(optim = c("nlminb", "optim"),
           method = c("BFGS", "L-BFGS-B"),
           lower = NULL,
           upper = NULL,
           REML = FALSE,
           profile = FALSE,
           verbose = 2,
           ...) {
    optim <- match.arg(optim)
    method <- match.arg(method)

    # check for valid args
    if (!is.logical(REML))
      stop("\n'REML' argument must be either TRUE or FALSE")
    if (!is.null(lower) & !inherits(lower, "list"))
      stop("\nlower parameter bounds must be specified as a named list")
    if (!is.null(upper) & !inherits(upper, "list"))
      stop("\nupper parameter bounds must be specified as a named list")
    if (!is.logical(profile))
      stop("\n'profile' argument must be either TRUE of FALSE")

    dots <- list(...)

    ## Set default control values
    pars <- switch(
      optim,
      nlminb = {
        if(!REML) {
          list(
            eval.max = 3000,
            iter.max = 2000,
            rel.tol = 1.0e-3,
            x.tol = 1.5e-2
            )
          } else {
            list(eval.max = 1000,
                 iter.max = 1000,
                 rel.tol = 1.0e-4,
                 x.tol = 1.5e-8)
          }
        },
      optim = list(maxit = 2000, reltol = 1.0e-3)
    )
    ## Override control parameters
    pars[names(dots)] <- dots
    list(optim = optim,
         method = method,
         lower = lower,
         upper = upper,
         REML = REML,
         profile = profile,
         verbose = verbose,
         control = pars)
}
