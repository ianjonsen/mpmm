##' @importFrom stats pnorm
##' @importFrom lme4 nobars findbars
##' @importFrom lme4 nobars
##' @importFrom dplyr %>%
##' @method summary mpmm
##' @export
summary.mpmm <- function(fit, ...) {
  if (length(list(...)) > 0) {
    warning("additional arguments ignored")
  }
  nobs <-
    with(fit$data, apply(!is.na(cbind(lon, lat)), 1, min)) %>%
    sum(.)

  mkAICtab <- function(fit) {
    data.frame(
      AIC      = fit$aic,
      BIC      = fit$bic,
      logLik   = ifelse("objective" %in% names(fit$opt),  -1 * fit$opt$objective, -1 * fit$opt$value),
      deviance = ifelse("objective" %in% names(fit$opt), fit$opt$objective * 2, fit$opt$value * 2),
      df.resid = nobs - length(fit$tmb$par)
    )
  }

  mkVartab <- function(fit) {
    grpnm <- names(fit$re)[1]
    rep <- fit$tmb$env$report(fit$tmb$env$last.par.best)
    stdev <- c(rep$sd[[1]],
               log(fit$par["sigma_g", "Estimate"]))
    Vartab <-
      data.frame(
        Group = c(grpnm, rep("", length(stdev) - 1)),
        Name = c(names(fit$re)[-1], "Residual"),
        Variance = exp(2 * stdev),
        StdDev = exp(stdev)
      )

    if (length(stdev) > 2) {
      tcorr <- rep$corr[[1]]
      tcorr[upper.tri(tcorr, diag = TRUE)] <- NA
      tcorr <- round(as.vector(tcorr)[!is.na(tcorr)], 2)
      nr <- nrow(Vartab)

      excols <- nr-2
      i <- 1

      while(i <= excols) {
        if(i == 1) {
          Vartab <- data.frame(Vartab, Corr = c("", tcorr[1:(nr-2)], ""))
          tcorr <- tcorr[-c(1:(nr-2))]
        }
        else {
          if(length(tcorr) == 1) {
            Vartab <- data.frame(Vartab, c(rep("",i), tcorr, ""))
            break
          } else {
          Vartab <- data.frame(Vartab, c(rep("",i), tcorr[1:i], ""))
          tcorr <- tcorr[-c(1:i)]
          }
        }
        i <- i + 1
      }
      nc <- ncol(Vartab)
      if(nc > 5) names(Vartab)[6:nc] <- ""
    }
    Vartab
  }

  mkFixtab <- function(fit) {
    terms <- nobars(fit$formula) %>%
      terms() %>%
      attr(., "term.labels")
    terms <- c("Intercept", terms)
    val <- fit$par[rownames(fit$par) %in% terms, "Estimate"]
    stderr <- fit$par[rownames(fit$par) %in% terms, "Std. Error"]
    terms[1] <- "(Intercept)"

    Fixtab <- cbind(#terms,
                         val,
                         stderr,
                         val / stderr,
                         2 * pnorm(abs(val/stderr), lower.tail = FALSE))
    colnames(Fixtab) <- c("Value", "Std.Error", "z value", "Pr(>|z|)")
    Fixtab
  }

  nobs <- with(fit$data, apply(!is.na(cbind(lon, lat)), 1, min)) %>%
    sum(.)
  resid.df <- nobs - length(fit$opt$par)
  terms <- nobars(fit$formula) %>%
    terms() %>%
    attr(., "term.labels")
  terms <- c("Intercept", terms)
  coef <- fit$par[rownames(fit$par) %in% terms, "Estimate"]

  ranform <- lme4::findbars(fit$formula) %>%
    as.character() %>%
    paste0("(", ., ")")
  fixform <- lme4::nobars(fit$formula) %>% as.character()

  structure(
    list(
      mf = fit$mf,
      logLik = mkAICtab(fit),
      grpnm = names(fit$re)[1],
      nobs = nobs,
      coefficients = coef,
      ranform = ranform,
      fixform = fixform,
      formula = fit$formula,
      Vartab = mkVartab(fit),
      Fixtab = mkFixtab(fit)
    ),
    class = "summary.mpmm"
  )
}

##' @importFrom stats printCoefmat
##' @method print summary.mpmm
##' @export
print.summary.mpmm <- function(x, digits = 3,
                               signif.stars = getOption("show.signif.stars"),
                               ...)
{


 # print(x$logLik); cat("\n")
    cat("Formula: ~", as.character(x$formula)[2], "\n")
    cat("Data: ", x$mf$data, "\n\n")
    print(x$logLik, row.names = FALSE); cat("\n\n")

    cat("Random effects: ~", x$ranform, "\n")
    print(x$Vartab, row.names = FALSE, digits = digits)
    cat("number of obs: ", x$nobs, ", group: ", x$grpnm, "\n\n", sep = "")

    cat("fixed effects:", x$fixform, "\n")
    printCoefmat(x$Fixtab, digits = digits, signif.stars = signif.stars)

  invisible(x)
}## print.summary.mpmm
