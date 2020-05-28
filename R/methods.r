##' @importFrom stats logLik
##' @export
logLik.mpmm <- function(m, ...) {
  if (!is.null(m$rep)) {
    val <- if (m$rep$pdHess) {
      ifelse("objective" %in% names(m$opt), -1 * m$opt$objective, -1 * m$opt$value)
    } else {
      warning("Hessian was not positive-definite\n", immediate. = TRUE)
      NA
    }
  } else {
    ifelse("objective" %in% names(m$opt), -1 * m$opt$objective, -1 * m$opt$value)
  }

  nobs <- nrow(m$fr)
  structure(
    val,
    nobs = nobs,
    nall = nobs,
    df = length(m$tmb$par),
    class = "logLik"
  )
}

##' @title perform likelihood ratio tests on 2 or more mpmm model objects
##' @param m an model object with class mpmm
##' @param ... additional mpmm model objects
##' @importFrom methods is
##' @importFrom stats var getCall pchisq anova
##' @export
anova.mpmm <- function (m, ...)
{
  ## borrows heavily from glmmTMB...
  mCall <- match.call(expand.dots = TRUE)
  dots <- list(...)

  ## detect multiple models, i.e. models in ...
  modp <- as.logical(vapply(dots, is, NA, "mpmm"))

  if (!any(modp))
    stop("anova() can only be used to compare multiple mpmm models")
  else {
    mods <- c(list(m), dots[modp])
    nobs <- function(m)
      nrow(m$fr)
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
