##' Fit a move persistence random walk via TMB to a pre-filtered/regularised animal
##'   track and estimate gamma as a linear function of covariates
##'
##' The input track is given as a dataframe where each row is an
##' observed location and columns
##' \describe{
##' \item{'id'}{individual animal identifier,}
##' \item{'date'}{observation time (POSIXct,GMT),}
##' \item{'lon'}{observed longitude,}
##' \item{'lat'}{observed latitude,}
##' \item{'tid'}{identifier for tracks if there are more than one track per
##' individual (optional),}
##' \item{'...'}{named covariates appended to track}
##' }
##'
##' @title Move Persistence Mixed-Effects Model
##' @param formula a right-hand-side regression formula (no response variable)
##' @param data a data frame of observations (see details)
##' @param se whether to return standard errors
##' @param control a list of control parameters for the outer optimization
##' (see \code{\link{mpmm_control}})
##' @param inner.control a list of control parameters for the inner optimization
##' (see \code{\link{TMB::MakeADFun}} and \code{\link{TMB::newton}})
##' @param model "mpmm" or "mpmm_dt"; "mpmm" (default) for data with regular
##' time intervals between locations, "mpmm_dt" is for irregular time intervals
##' @return a list with components
##' \item{\code{states}}{a dataframe of estimated states}
##' \item{\code{fitted}}{a dataframe of fitted locations}
##' \item{\code{par}}{model parameter summmary}
##' \item{\code{data}}{input dataframe}
##' \item{\code{tmb}}{the tmb object}
##' \item{\code{opt}}{the object returned by the optimizer}
##' @examples
##' data(ellie.ice)
##' fit <- mpmm(~ ice + (1 | id), data = ellie.ice)
##' summary(fit)
##'
##' @useDynLib mpmm
##' @importFrom lme4 nobars findbars subbars mkReTrms
##' @importFrom glmmTMB getReStruc splitForm
##' @importFrom Matrix t
##' @importFrom dplyr %>% arrange count mutate as_tibble tibble
##' @importFrom TMB MakeADFun sdreport newtonOption
##' @export
mpmm <- function(
                formula = NA,
                data = NULL,
                map = NULL,
                se = TRUE,
                control = mpmm_control(),
                inner.control = inner_control(),
                model = "mpmm") {
  st <- proc.time()

  call <- mf <- match.call()

  # Create a tid column if there is none specified & model == "mpmm_dt"
  if(model == "mpmm_dt" & all(colnames(data) != "tid"))
    stop("input data must have a tid column to denote sub-track id")

  # ordering the data to make sure we have continuous tracks and ids are ordered
  # NOTE: we can't do this in the fn call as model.frame (mf) pulls covars from
  #   the GlobalEnv (ie. input data _prior_ to calling mpmm())
  # solution is to check whether input data are sorted by id, date and return
  # error if not
  if(!inherits(all.equal(data, data %>% arrange(id, date)), "logical"))
    stop("input data must be sorted by id and date")
  #  data <- data %>% arrange(id, tid, date)

  # check that the formula is a formula
  is.formula <- function(x)
    tryCatch(
      inherits(x, "formula"),
      error = function(e) {
        FALSE
      }
    )
  if (!is.formula(formula))
    stop("\n'formula' must be specified as ~ x + ...")

  # check that there is no response variable in the formula
  if (attr(terms(formula), "response") != 0)
    stop("\n'formula' can not have a response variable")

  # check that formula has a random component
  if(is.null(findbars(formula)))
    stop("\n formula must include a random component; e.g., ~ (1 | id)")

  # check that covariates do not contain NA's
  ## should add proper na.action to model frame...
  ## could add prior to handle missing values...
  covars <- nobars(formula) %>% terms() %>% attr(., "term.labels")
  nas <- is.na(data[, covars]) %>% apply(., 2, sum)
  if (sum(nas) > 0)
    stop(
      paste0(
        "\n NA's detected in the following covariates: ",
        covars[which(nas > 0)],
        "\n consider imputing values"
      )
    )

  # evaluate model frame
  m <- match(c("data", "subset", "na.action"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf$formula <- formula

  f <- subbars(formula)
  environment(f) <- environment(formula)
  mf$formula <- f
  fr <- eval(mf, envir = environment(formula), enclos = parent.frame())

  # strip random part of formula
  ff <- nobars(formula)
  # num observations
  nobs <- nrow(fr)

  # build fixed effects model matrix X
  # no model matrix if fixed part of formula is empty
  if (identical(ff, ~  0) ||
      identical(ff, ~ -1)) {
    X <- NULL
  } else {
    mf$formula <- ff
    termf <- terms(eval(mf, envir = environment(ff)))
    X <- model.matrix(ff, fr)
    terms <- list(fixed = terms(termf))
  }

  # build random effects sparse matrix Z
  ref <- formula
  if (is.null(findbars(ref))) {
    Z <- matrix(0, nrow = nobs, ncol = 0)
    Z <- as(Z, "dgTMatrix")
    reTrms <- NULL
    ss <- integer(0)
  } else {
    ref <- findbars(ref)
    reTrms <- mkReTrms(ref, fr)
    mf$formula <- ref
    ss <- splitForm(formula)
    ss <- unlist(ss$reTrmClasses)
    Z <- t(reTrms$Zt)
  }

  condList  <- list(X = X, Z = Z, reTrms = reTrms, ss = ss, terms = terms)

  condReStruc <- with(condList, getReStruc(reTrms, ss))

  gnm <- names(condList$reTrms$flist)

  # Number of tracks (or individual if only one track per individual)
  switch(model,
         mpmm = {
           A <- nrow(count(data, id))
           # get index of start and end of tracks
           idx <- data$id %>%
             table() %>%
             as.numeric() %>%
             cumsum() %>%
             c(0, .)
           # create di vector
           data$di <- 1
         },
         mpmm_dt = {
           A <- nrow(count(data, id, tid))
           # get index of start and end of tracks
           data$idtid <- with(data, paste(id, tid, sep=""))
           idx <- data$idtid %>%
             table() %>%
             as.numeric() %>%
             cumsum() %>%
             c(0, .)
           # create di vector
           data$di <- c(NA, diff(data$date))
           data$di[idx[1:(length(idx)-1)] + 1] <- NA
           # Scale to median
           data$di <- data$di/median(data$di, na.rm=TRUE)
         })

  ## FIXME::the mpmm_dt code appears to be messing up the idx as at least with
  ##     some datasets (eg. ~/Dropbox/collab/vogel/data/d.all.kw.data24.7.csv)
  ##     the NA's get inserted in the wrong place - ie. in middle of an
  ##     individual, this results in di values jumping to
  ##     either -ve #'s or v big #'s

  ## build data for TMB
  data.tmb <- list(
    X     = condList$X,
    Z     = condList$Z,
    ll    = cbind(data$lon, data$lat),
    idx   = idx,
    di    = data$di,
    model = ifelse(model == "mpmm", 0, 1),
    A     = A,
    terms = condReStruc
  )

  getVal <- function(obj, component)
    vapply(obj, function(x) x[[component]], numeric(1))


  param <- with(data.tmb,
                     list(
                       lg          = rep(0, nobs),
                       beta        = rep(0, ncol(X)),
                       b           = rep(0, ncol(Z)),
                       l_sigma     = c(0,0),
                       l_rho       = 0,
                       l_sigma_g   = 0,
                       theta       = rep(0,
                                         sum(getVal(condReStruc,
                                                    "blockNumTheta")))
                     ))
  rnd <- c("lg", if(ncol(data.tmb$Z) > 0) "b")

  if(!is.null(map)) {
    names(map) <- paste0("l_", names(map))
  }

  forTMB <- list(data = data.tmb,
              param = param,
              rnd = rnd,
              gnm = gnm,
              condList = condList,
              condReStruc = condReStruc,
              allForm = list(formula),
              fr = fr,
              call = call,
              verbose = control$verbose
              )

  # integrate out the beta's from the likelihood - REML estimation;
  #  appends the beta's to the random arg.
  profl <- NULL
  if(control$REML || control$profile) profl <- "beta"

  ## TMB - create objective function
  obj <-
    with(
      forTMB,
      MakeADFun(
        data = data.tmb,
        parameters = param,
        map = map,
        random = rnd,
        profile = profl,
        DLL = "mpmm",
        hessian = TRUE,
        method = control$method,
        inner.control = inner.control,
        silent = ifelse(control$verbose == 1, FALSE, TRUE)
      )
    )

  obj$env$inner.control$trace <- ifelse(control$verbose == 1, TRUE, FALSE)
  obj$env$tracemgc <- ifelse(control$verbose == 1, TRUE, FALSE)

  ## add par values to trace if verbose = TRUE
  myfn <- function(x) {
    cat("\r", "pars: ", round(x, 3), "     ")
    flush.console()
    obj$fn(x)
  }

  if (control$method == "L-BFGS-B" & any(is.null(control$lower),
                                         is.null(control$upper))) {
    ## Set parameter bounds - most are -Inf, Inf
    L <- c(
      beta = rep(-100, ncol(X)),
      l_sigma = c(-50, -50),
      l_rho = -10,
      l_sigma_g = -50,
      theta = rep(-Inf, sum(getVal(
        condReStruc, "blockNumTheta"
      )))
    )

    U <- c(
      beta = rep(100, ncol(X)),
      l_sigma = c(100, 100),
      l_rho = 10,
      l_sigma_g = 100,
      theta = rep(Inf, sum(getVal(
        condReStruc, "blockNumTheta"
      )))
    )
  } else if (control$method == "BFGS") {
    ## Unbounded parameters - all are -Inf, Inf
    L <- c(
      beta = rep(-Inf, ncol(X)),
      l_sigma = c(-Inf, -Inf),
      l_rho = -Inf,
      l_sigma_g = -Inf,
      theta = rep(-Inf, sum(getVal(
        condReStruc, "blockNumTheta"
      )))
    )

    U <- c(
      beta = rep(Inf, ncol(X)),
      l_sigma = c(Inf, Inf),
      l_rho = Inf,
      l_sigma_g = Inf,
      theta = rep(Inf, sum(getVal(
        condReStruc, "blockNumTheta"
      )))
    )

  } else if(control$method == "L-BFGS-B" & all(!is.null(control$lower),
                                               !is.null(control$upper))) {
    L <- control$lower
    U <- control$upper
  }

  ## Minimize objective function
  cat("using", control$optim, control$method, "\n")
  opt <- suppressWarnings(switch(control$optim,
                                 nlminb = try(nlminb(
                                   start = obj$par,
                                   objective = ifelse(control$verbose == 2, myfn, obj$fn),
                                   gradient = obj$gr,
                                   control = control$control,
                                   lower = L,
                                   upper = U
    ))
    ,
    optim = try(do.call(
      optim,
      args = list(
        par = obj$par,
        fn = ifelse(control$verbose == 2, myfn, obj$fn),
        gr = obj$gr,
        method = control$method,
        control = control$control,
        lower = L,
        upper = U
      )
  ))))

  if(control$profile && !control$REML) {

    ## Sparse Schur complement (Marginal of precision matrix)
    ##' @importFrom Matrix Cholesky solve
    GMRFmarginal <- function(Q, i, ...) {
      ind <- seq_len(nrow(Q))
      i1 <- (ind)[i]
      i0 <- setdiff(ind, i1)
      if (length(i0) == 0)
        return(Q)
      Q0 <- as(Q[i0, i0, drop = FALSE], "symmetricMatrix")
      L0 <- Cholesky(Q0, ...)
      ans <- Q[i1, i1, drop = FALSE] - Q[i1, i0, drop = FALSE] %*%
        solve(Q0, as.matrix(Q[i0, i1, drop = FALSE]))
      ans
    }
    ## use glmmTMB as a guide here...

    sdr <- sdreport(obj, getJointPrecision=TRUE)
    parnames <- names(obj$env$par)
    Q <- sdr$jointPrecision; dimnames(Q) <- list(parnames, parnames)
    whichNotRandom <- which( ! parnames %in% c("b","lg"))
    Qm <- GMRFmarginal(Q, whichNotRandom)
    h <- as.matrix(Qm) ## Hessian of *all* (non-random) parameters
    parameters <- obj$env$parList(opt$par, obj$env$last.par.best)

    ## without profile, with REML estimates as inits
    obj <-
      with(
        forTMB,
        MakeADFun(
          data = data.tmb,
          parameters = parameters,
          map = map,
          random = rnd,
          profile = NULL,
          DLL = "mpmm",
          method = control$method,
          silent = ifelse(control$verbose == 1, FALSE, TRUE)
        )
      )
    ## Run up to 5 Newton iterations with fixed (off-mode) hessian
    oldpar <- par <- obj$par; iter <- 0
    ## FIXME: Make configurable ?
    max.newton.steps <- 5
    newton.tol <- 1e-10

    if (sdr$pdHess) {
      ## pdHess can be FALSE
      ##  * Happens for boundary fits (e.g. dispersion close to 0 - see 'spline' example)
      ##    * Option 1: Fall back to old method
      ##    * Option 2: Skip Newton iterations
      for (iter in seq_len(max.newton.steps)) {
        g <- as.numeric( obj$gr(par) )
        if (any(is.na(g)) || max(abs(g)) < newton.tol) break
        par <- par - solve(h, g)
      }

      if (any(is.na(g))) {
        warning("a Newton step failed in profiling")
        par <- oldpar
      }
    } else {
      warning("\n profiling not possible as Hessian was not positive-definite")
    }

    opt$par <- par
    switch(control$optim,
           nlminb = {
             opt$objective <- obj$fn(par)
           },
           optim = {
             opt$value <- obj$fn(par)
           })

    opt$newton.steps <- iter

  }

  opt$parfull <- obj$env$last.par.best[which(!names(obj$env$last.par.best) %in% "lg")]


  ## Parameters, states and the fitted values
  if(control$profile) {
    rep <- sdreport(obj, hessian.fixed = h)
  } else {
    rep <- sdreport(obj, getJointPrecision = control$REML)
  }

#  if(!rep$pdHess || !se) {
#    ## don't calculate standard errors
#    if(!rep$pdHess && method != "REML") warning("\n Hession was not positive-definite, getting fixed estimates without standard errors")
#    rep <- sdreport(obj, ignore.parm.uncertainty = TRUE)
# }
  if(!rep$pdHess && !control$REML) {
    warning("\n Hession was not positive-definite, fixed estimates do not have standard errors")
    rep <- sdreport(obj, ignore.parm.uncertainty = TRUE)
  }

  fxd <- summary(rep, "report")
  fxd_log <- summary(rep, "fixed")
  rdm <- summary(rep, "random")

  lg <- rdm[rownames(rdm) %in% "lg",]
  b <- rdm[rownames(rdm) %in% "b", ]

  ## build table of ranefs
  nms <- c(gnm, unlist(reTrms$cnms))
  ## get number of random terms
  nrt <- sapply(data.tmb$terms, function(x)
    x$blockSize) %>% sum()
  nl <- sapply(reTrms$flist, nlevels)
  ret <- matrix(b[, "Estimate"], nl, nrt, byrow = TRUE) %>%
    as.data.frame() %>%
    data.frame(unique(reTrms$flist[[1]]), .) %>%
    as_tibble()
  names(ret) <- nms

  ## build table of gamma estimates
  fitted <- tibble(
    id = data$id,
    date = data$date,
    g = plogis(lg[, 1]),
    g.se = lg[, 2]
  )

  rownames(fxd)[rownames(fxd) %in% "sigma"] <-
    c("sigma_lon", "sigma_lat")
  rownames(fxd)[rownames(fxd) %in% "beta"][1] <- "Intercept"

  ft <- attr(termf, "term.labels")
  rownames(fxd)[rownames(fxd) %in% "beta"] <- ft

  opt.time <- proc.time() - st
  cat("\n", "timing: ", opt.time, "\n")

  ## FIXME:: need to simplify and organise...
  structure(
    list(
      call = call,
      formula = formula,
      data = data,
      mf = mf,
      fr = fr,
      fitted = fitted,
      par = fxd,
      re = ret,
      tmb = obj,
      opt = opt,
      REML = control$REML,
      rep = rep,
      opt.time = opt.time
    ),
    class = "mpmm"
  )

}

