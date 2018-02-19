##' @method print mpmm
##' @export
print.mpmm <-
  function(x, digits = max(3, getOption("digits") - 3),
           signif.stars = getOption("show.signif.stars"),
          ...)
  {
    ## Type Of Model fit --- REML? ---['class']  & Family & Call
    print(x$call); cat("\n")
    ## the 'digits' argument should have an action here
    aictab <- c(AIC = AIC(x), BIC = BIC(x), logLik = logLik(x),
                df.resid = df.residual(x))
    .prt.aictab(aictab, digits=digits+1)
    ## varcorr
    if (!all(sapply(vc <- VarCorr(x),is.null))) {
      cat("Random-effects (co)variances:\n")
      print(VarCorr(x), digits=digits, comp = ranef.comp)
    }
    ## ngroups
    gvec <- list(obs=sprintf("\nNumber of obs: %d",nobs(x)))
    ng <- ngrps.glmmTMB(x)
    for (i in seq_along(ng)) {
      if (length(ng[[i]])>0) {
        nm <- names(ng)[i]
        gvec[[nm]] <- paste0(cNames[nm],": ",
                             paste(paste(names(ng[[i]]), ng[[i]], sep=", "), collapse="; "))
      }
    }
    cat(do.call(paste,c(gvec,list(sep=" / "))),fill=TRUE)
    
    if(trivialDisp(x)) {# if trivial print here, else below(~x) or none(~0)
      printDispersion(x$modelInfo$family$family,sigma(x))  
    }
    ## Family specific parameters
    printFamily(x$modelInfo$family$family, x)
    ## Fixed effects:
    if(length(cf <- fixef(x)) > 0) {
      cat("\nFixed Effects:\n")
      print(cf, ...)
    } else
      cat("No fixed effect coefficients\n")
    invisible(x)
  }
