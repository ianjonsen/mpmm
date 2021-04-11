##' @title Extract log-likelihood
##' @description extract log-likelihood from an mpmm fit object
##' @param object an mpmm model fit object
##' @param ... additional arguments to be ignored
##' @importFrom stats logLik
##' @method logLik mpmm
##' @export
logLik.mpmm <- function(object, ...) {
  if (length(list(...)) > 0) {
    warning("additional arguments ignored")
  }
  if (!is.null(object$rep)) {
    val <- if (object$rep$pdHess) {
      ifelse("objective" %in% names(object$opt), -1 * object$opt$objective, -1 * object$opt$value)
    } else {
      warning("Hessian was not positive-definite\n", immediate. = TRUE)
      NA
    }
  } else {
    ifelse("objective" %in% names(object$opt), -1 * object$opt$objective, -1 * object$opt$value)
  }

  nobs <- nrow(object$fr)
  structure(
    val,
    nobs = nobs,
    nall = nobs,
    df = length(object$tmb$par),
    class = "logLik"
  )
}

##' @title anova tables
##' @description perform likelihood ratio tests on 2 or more mpmm fit objects
##' @param object an mpmm fit object
##' @param ... additional mpmm fit objects
##' @importFrom methods is
##' @importFrom stats var getCall pchisq anova
##' @method anova mpmm
##' @export
anova.mpmm <- function (object, ...)
{
  ## borrows heavily from glmmTMB...
  mCall <- match.call(expand.dots = TRUE)
  dots <- list(...)

  ## detect multiple models, i.e. models in ...
  modp <- as.logical(vapply(dots, is, NA, "mpmm"))

  if (!any(modp))
    stop("anova() can only be used to compare multiple mpmm models")
  else {
    mods <- c(list(object), dots[modp])
    nobs <- function(object)
      nrow(object$fr)
    nobs.vec <- vapply(mods, nobs, 1L)

    if (var(nobs.vec) > 0)
      stop("models must be fit to the same data")

    mNms <- vapply(as.list(mCall), deparse, "")[-1]

    if (any(duplicated(mNms)))
      stop("model names must be unique")

    names(mods) <- sub("@env$", "", mNms)
    llks <- lapply(mods, logLik)
    ii <- order(df <- vapply(llks, attr, FUN.VALUE = numeric(1),
                             "df"))
    mods <- mods[ii]
    llks <- llks[ii]
    df <- df[ii]

    calls <- lapply(mods, getCall)
    data <- lapply(calls, `[[`, "data")

    if (!all(vapply(data, identical, NA, data[[1]])))
      stop("all models must be fit to the same data object")

    header <- paste("Data:", deparse(data[[1]]))
    subset <- lapply(calls, `[[`, "subset")

    if (!all(vapply(subset, identical, NA, subset[[1]])))
      stop("all models must use the same subset")

    if (!is.null(subset[[1]]))
      header <- c(header, paste("Subset:", deparse(subset[[1]])))

    llk <- unlist(llks)
    chisq <- 2 * pmax(0, c(NA, diff(llk)))
    dfChisq <- c(NA, diff(df))
    val <- data.frame(
      df = df,
      AIC = unlist(lapply(llks, AIC)),
      BIC = unlist(lapply(llks, BIC)),
      logLik = llk,
      deviance = -2 * llk,
      Chisq = chisq,
      `Chi_df` = dfChisq,
      `Pr(>Chisq)` = pchisq(chisq,
                            dfChisq, lower.tail = FALSE),
      row.names = names(mods),
      check.names = FALSE
    )
    class(val) <- c("anova", class(val))
    forms <- lapply(lapply(calls, `[[`, "formula"), deparse)

    structure(val, heading = c(header, "Models:",
                               paste(paste(
                                 paste(rep(names(mods), times = lengths(forms)), unlist(forms), sep = ": ")
                               ))))
  }
}

