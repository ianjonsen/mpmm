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


##' Extract Fixed Effects
##'
##' Extract fixed effects from a fitted \code{mpmm} model.
##' @name fixef
##' @title Extract fixed-effects estimates
##' @aliases fixef fixef.mpmm
##' @docType methods
##' @param object a fitted mpmm object
##' @param \dots optional additional arguments. not implemented
##' @return an object of class \code{fixef.mpmm} comprising a list of component (\code{cond})
##' @keywords models
##' @details The print method for \code{fixef.mpmm} object \emph{only displays non-trivial components}
##' @examples
##'
##' @importFrom nlme fixef
##' @export fixef
##' @export
fixef.mpmm <- function(object, ...) {
  getXnm <- function(suffix) {
    nm <- paste0("X",suffix)
    return(colnames(getME(object, nm)))
  }
  pl <- object$tmb$env$parList(object$opt$par)
  structure(list(cond = setNames(pl$beta,   getXnm(""))),
            class = "fixef.mpmm")
}


##' Extract or Get Generalize Components from a Fitted Mixed Effects Model
##'
##' @aliases getME
##' @param object a fitted \code{mpmm} object
##' @param name of the component to be retrieved
##' @param \dots ignored, for method compatibility
##'
##' @seealso \code{\link[lme4]{getME}}
##' Get generic and re-export:
##' @importFrom lme4 getME
##' @export getME
##'
##' @method getME mpmm
##' @export
getME.mpmm <- function(object,
                       name = c("X", "Z", "theta", "beta"),
                       ...)
{
  if(missing(name)) stop("'name' must not be missing")
  ## Deal with multiple names -- "FIXME" is inefficiently redoing things
  if (length(name <- as.character(name)) > 1) {
    names(name) <- name
    return(lapply(name, getME, object = object))
  }
  if(name == "ALL") ## recursively get all provided components
    return(sapply(eval(formals()$name),
                  getME.mpmm, object=object, simplify=FALSE))

  stopifnot(inherits(object, "mpmm"))
  name <- match.arg(name)

  oo.env <- object$tmb$env
  ### Start of the switch
  allpars <- oo.env$parList(object$opt$par, object$opt$parfull)
  isSparse <- function(component) { if (is.null(om <- object$modelInfo$sparseX)) FALSE else om[[component]] }
  switch(name,
         "X"     = if (!isSparse("cond")) oo.env$data$X else oo.env$data$XS,
         "Z"     = oo.env$data$Z,
         "theta" = allpars$theta ,
         "beta"  = unlist(allpars[c("beta","betazi","betad")]),
         "..foo.." = # placeholder!
           stop(gettextf("'%s' is not implemented yet",
                         sprintf("getME(*, \"%s\")", name))),
         ## otherwise
         stop(sprintf("Mixed-Effects extraction of '%s' is not available for class \"%s\"",
                      name, class(object))))
}## {getME}


##' Calculate Variance-Covariance Matrix for a Fitted mpmm model
##'
##' @param object a \dQuote{mpmm} fit
##' @param full return a full variance-covariance matrix?
##' @param \dots ignored, for method compatibility
##' @return By default (\code{full==FALSE}).  If \code{full==TRUE}, a single square variance-covariance matrix for \emph{all} top-level model parameters (conditional and variance-covariance parameters)
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats vcov
##' @export
vcov.mpmm <- function(object, full=FALSE, ...) {
  REML <- object$method == "REML"
  if(is.null(sdr <- object$sdr)) {
    warning("Calculating sdreport. Use se=TRUE in mpmm to avoid repetitive calculation of sdreport")
    sdr <- sdreport(object$tmb, getJointPrecision=REML)
  }
  if (REML) {
    ## NOTE: This code would also work in non-REML case provided
    ## that jointPrecision is present in the object.
    Q <- sdr$jointPrecision
    whichNotRandom <- which( ! rownames(Q) %in% c("b", "lg") )
    Qm <- GMRFmarginal(Q, whichNotRandom)
    cov.all.parms <- solve(as.matrix(Qm))
  } else {
    cov.all.parms <- sdr$cov.fixed
  }
  keepTag <- if (full) { "."
  } else "beta($|[^d])"
  to_keep <- grep(keepTag,colnames(cov.all.parms)) # only keep betas
  covF <- cov.all.parms[to_keep,to_keep,drop=FALSE]

  mkNames <- function(tag) {
    X <- getME(object,paste0("X",tag))
    if (trivialFixef(nn <- colnames(X),tag)
        ## if 'full', keep disp even if trivial, if used by family
        && !(full && tag =="d" &&
             (usesDispersion(family(object)$family) && !zeroDisp(object)))) {
      return(character(0))
    }
    return(paste(tag,nn,sep="~"))
  }

  nameList <- setNames(list(colnames(getME(object,"X"))),
                       names(cNames))

  if(full) {
    ## FIXME: haven't really decided if we should drop the
    ##   trivial variance-covariance dispersion parameter ??
    ## if (trivialDisp(object))
    ##    res <- covF[-nrow(covF),-nrow(covF)]

    reNames <- function(tag) {
      re <- object$modelInfo$reStruc[[paste0(tag,"ReStruc")]]
      nn <- mapply(function(n,L) paste(n,seq(L),sep="."),
                   names(re),
                   sapply(re,"[[","blockNumTheta"))
      if (length(nn)==0) return(nn)
      return(paste("theta",gsub(" ","",nn),sep="_"))
    }
    nameList <- c(nameList,list(theta=reNames("cond")))
  }

  ## drop NA-mapped variables

  ## for matching map names vs nameList components ...
  par_components <- c("beta","theta")

  map <- object$tmb$env$map
  for (m in seq_along(map)) {
    if (length(NAmap <- which(is.na(map[[m]])))>0) {
      w <- match(names(map)[m],par_components) ##
      if (length(nameList)>=w) { ## may not exist if !full
        nameList[[w]] <- nameList[[w]][-NAmap]
      }
    }
  }

  if (full) {
    colnames(covF) <- rownames(covF) <- unlist(nameList)
    res <- covF        ## return just a matrix in this case
  } else {
    splitMat <- function(x) {
      ss <- split(seq_along(colnames(x)),
                  colnames(x))
      lapply(ss,function(z) x[z,z,drop=FALSE])
    }
    covList <- splitMat(covF)
    names(covList) <-
      names(cNames)[match(names(covList),c("beta"))]
    for (nm in names(covList)) {
      if (length(xnms <- nameList[[nm]])==0) {
        covList[[nm]] <- NULL
      }
      else dimnames(covList[[nm]]) <- list(xnms,xnms)
    }
    res <- covList
    ##  FIXME: should vcov always return a three-element list
    ## (with NULL values for trivial models)?
    class(res) <- c("vcov.mpmm","matrix")
  }
  return(res)
}

##' @method print vcov.mpmm
##' @export
print.vcov.mpmm <- function(x,...) {
  for (nm in names(x)) {
    cat(cNames[[nm]],":\n",sep="")
    print(x[[nm]])
    cat("\n")
  }
  invisible(x)
}



##' Extract Random Effects
##'
##' Extract random effects from a fitted \code{mpmm} model.
##'
##' @param object a \code{mpmm} model.
##' @param condVar whether to include conditional variances in result.
##' @param \dots not implemented
##'
##' @examples
##'
##' @aliases ranef ranef.mpmm
##' @importFrom nlme ranef
##' @export ranef
##' @export
ranef.mpmm <- function(object, condVar=TRUE, ...) {
  check_dots(...)
  ## The arrange() function converts a vector of random effects to a list of
  ## data frames, in the same way as lme4 does.
  ## FIXME: add condVar, make sure format matches lme4
  arrange <- function(x, sd, listname)
  {
    ##FIXME:: revise for mpmm object structure...
    cnms <- object$modelInfo$reTrms[[listname]]$cnms
    flist <- object$modelInfo$reTrms[[listname]]$flist
    if (!is.null(cnms)) {
      levs <- lapply(fl <- flist, levels)
      asgn <- attr(fl, "assign")
      nc <- vapply(cnms, length, 1L)     ## number of columns (terms) per RE
      nb <- nc * vapply(levs, length, 1L)[asgn] ## number of elements per RE
      nbseq <- rep.int(seq_along(nb), nb)       ## splitting vector
      ml <- split(x, nbseq)
      for (i in seq_along(ml)) {
        ml[[i]] <- matrix(ml[[i]], ncol=nc[i], byrow=TRUE,
                          dimnames=list(NULL, cnms[[i]]))
      }
      if (!is.null(sd)) {
        sd <- split(sd,nbseq)
        for (i in seq_along(sd)) {
          ii <- asgn[i]
          nr <- length(levs[[ii]])
          a <- array(NA,dim=c(nc[i],nc[i],nr))
          ## fill in diagonals: off-diagonals will stay NA (!)
          ## unless we bother to retrieve conditional covariance info
          ## from the fit
          ## when nc>1, what order is the sd vector in?
          ## guessing, level-wise
          for (j in seq(nr)) {
            a[cbind(seq(nc[i]),seq(nc[i]),j)] <-
              (sd[[i]][nc[i]*(j-1)+seq(nc[i])])^2
          }
          sd[[i]] <- a
        }
      }
      ## combine RE matrices from all terms with the same grouping factor
      x <- lapply(seq_along(fl), function(i) {
        d <- data.frame(do.call(cbind, ml[asgn==i]), row.names=levs[[i]],
                        check.names=FALSE)
        if (!is.null(sd)) {
          ## attach conditional variance info
          ## called "condVar", *not* "postVar" (contrast to lme4)
          attr(d, "condVar") <- if (length(w <- which(asgn==i))>1) {
            ## FIXME: set names?
            sd[w]  ## if more than one term, list
          } else sd[[w]]  ## else just the array
        }
        return(d)
      })
      names(x) <- names(fl)
      return(x)
    } ## if !is.null(cnms)
    else {
      list()
    }
  } ## arrange()

  pl <- getParList(object)  ## see VarCorr.R
  if (condVar && hasRandom(object))  {
    ss <- summary(object$rep,"random")
    sdl <- list(b=ss[rownames(ss)=="b","Std. Error"])
  }  else sdl <- NULL
  structure(list(cond = arrange(pl$b, sdl$b, "cond")),
            class = "ranef.mpmm")
}

##' @method print ranef.mpmm
##' @export
print.ranef.mpmm <- function(x, simplify=TRUE, ...) {
  print(if (simplify)
    unclass(x$cond) else unclass(x),
    ...)
  invisible(x)
}
