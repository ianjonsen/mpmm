##' Fit a move persistence random walk via TMB to a pre-filtered/regularised animal
##'   track and estimate gamma as a linear function of covariates
##'
##' The input track is given as a tibble grouped by id (and optionally tid) where each row is an
##' observed location and columns
##' \describe{
##' \item{'id'}{individual animal identifier,}
##' \item{'date'}{observation time (POSIXct,GMT),}
##' \item{'lon'}{observed longitude,}
##' \item{'lat'}{observed latitude,}
##' \item{'tid'}{identifier for tracks if there are more than one track per individual (optional),}
##' \item{'...'}{named covariates appended to track}
##' }
##'
##' @title Move Persistence Mixed-Effects Model
##' @param formula a right-hand-side regression formula (no response variable)
##' @param data a tibble of observations, grouped by at least 'id' (see details)
##' @param method method for maximising the log-likelihood ("ML" or "REML")
##' @param optim numerical optimizer to be used (nlminb or optim)
##' @param control a list of control parameters (currently only for nlminb)
##' @param verbose report progress during minimization
##' @param model "mpmm" or "mpmm_dt", "mpmm" is the default value and is for a model with regular time intervals between locations, "mpmm_dt" is for irregular time intervals.
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
##' @importFrom dplyr %>% arrange count mutate tbl_df data_frame is.grouped_df group_by
##' @importFrom TMB MakeADFun sdreport newtonOption
##' @export
mpmm <- function(
                formula = NA,
                data = NULL,
                method = "ML",
                optim = c("nlminb", "optim"),
                control = NULL,
                verbose = FALSE,
                model = "mpmm") {

  call <- mf <- match.call()
  optim <- match.arg(optim)

  # check that data contains an 'id' variable
  if(!"id" %in% names(data)) stop("data must contain an id variable")

  # check that data is a grouped tibble
  if(!is.grouped_df(data)) stop("data must be grouped by 'id' and (optionally) 'tid'")

  # Create a tid column if there is none specified
  if(all(colnames(data) != "tid")){
    data$tid <- NA
  }

  # ordering the data to make sure we have continuous tracks and ids are ordered
  data <- data %>% arrange(id, tid, date)

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

  # check that either ML or REML are the specified maximisation method
  if (!method %in% c("ML", "REML"))
    stop("\n'method' argument must be either ML or REML")

  # check that formula has a random component
  if(is.null(findbars(formula)))
    stop("\n formula must include a random component; e.g., ~ (1 | id)")

  # check that covariates do not contain NA's
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

  # num observations
  nobs <- nrow(fr)

  # Number of tracks (or individual if only one track per individual)
  A <- nrow(count(data, id, tid))

  # num random effects
  #nre <- sapply(reTrms$cnms, length) %>% sum()

  # get index of start and end of tracks
  data <- data %>% mutate(idtid = paste(id, tid, sep=""))
  idx <- data$idtid %>%
    table() %>%
    as.numeric() %>%
    cumsum() %>%
    c(0, .)

  # Create dt vector if model is mpmm_dt
  # dt = t_i - t_{i-1} and include in data.tmb
  if(model == "mpmm_dt"){
    data$di <- c(NA, diff(data$date))
    data$di[idx[1:(length(idx)-1)] + 1] <- NA
    # Scale to median
    data$di <- data$di/median(data$di, na.rm=TRUE)
  }else{
    data$di <- 1
  }

  ## build data for TMB
  data.tmb <- list(
    X     = condList$X,
    Z     = condList$Z,
    ll    = cbind(data$lon, data$lat),
    idx   = idx,
    di   = data$di,
    A     = A,
    terms = condReStruc
  )

  getVal <- function(obj, component)
    vapply(obj, function(x) x[[component]], numeric(1))


  param <- with(data.tmb,
                     list(
                       lg          = rep(1, nobs),
                       beta        = rep(1, ncol(X)),
                       b           = rep(1, ncol(Z)),
                       log_sigma   = c(0,0),
                       log_sigma_g = 0,
                       theta       = rep(0, sum(getVal(condReStruc, "blockNumTheta")))
                     ))
  rnd <- c("lg", if(ncol(data.tmb$Z) > 0) "b")


  forTMB <- list(data = data.tmb,
              param = param,
              rnd = rnd,
              gnm = gnm,
              condList = condList,
              condReStruc = condReStruc,
              allForm = list(formula),
              fr = fr,
              call = call,
              verbose = verbose
              )

  # integrate out the beta's from the likelihood - REML estimation; appends the beta's to the random arg.
  profile <- NULL
  if(method == "REML") profile <- "beta"

  ## TMB - create objective function
  obj <-
    with(
      forTMB,
      MakeADFun(
        data = data.tmb,
        parameters = param,
        random = rnd,
        profile = profile,
        DLL = "mpmm",
      #  hessian = TRUE,
        silent = !verbose
      )
    )

  obj$env$inner.control$trace <- verbose
  obj$env$tracemgc <- verbose

  obj$control <- list(trace = 0,
                      reltol = 1e-12,
                      maxit = 500)
  newtonOption(obj, smartsearch = TRUE)

  ## Minimize objective function
  opt.time <- system.time(opt <- suppressWarnings(switch(
    optim,
    nlminb = nlminb(
      start = obj$par,
      objective = obj$fn,
      gradient = obj$gr,
      control = control
    ),
    optim = do.call("optim", obj)
  )))

  cat(opt.time, "\n")

  ## Parameters, states and the fitted values
  rep <- sdreport(obj)
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
    tbl_df()
  names(ret) <- nms

  ## build table of gamma estimates
  fitted <- data_frame(
    id = data$id,
    date = data$date,
    g = plogis(lg[, 1]),
    g.se = lg[, 2]
  )

  if (optim == "nlminb") {
    aic <- 2 * length(opt[["par"]]) + 2 * opt[["objective"]]
    bic <-
      log(nrow(data)) * length(opt[["par"]]) + 2 * opt[["objective"]]
  }
  else {
    aic <- 2 * length(opt[["par"]]) + 2 * opt[["value"]]
    bic <- log(nrow(data)) * length(opt[["par"]]) + 2 * opt[["value"]]
  }

  rownames(fxd)[rownames(fxd) %in% "sigma"] <-
    c("sigma_lon", "sigma_lat")
  rownames(fxd)[rownames(fxd) %in% "beta"][1] <- "Intercept"

  ft <- attr(termf, "term.labels")
  rownames(fxd)[rownames(fxd) %in% "beta"] <- ft

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
      method = method,
      rep = rep,
      aic = aic,
      bic = bic,
      opt.time = opt.time
    ),
    class = "mpmm"
  )

}

