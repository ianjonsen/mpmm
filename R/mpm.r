##' Move Persistence Model
##'
##' Fit a random walk with time-varying move persistence to location data
##' with or without measurement error
##'
##' The input track is given as a dataframe where each row is an
##' observed location and columns
##' \describe{
##' \item{'id'}{individual identification,}
##' \item{'date'}{observation time (POSIXct,GMT),}
##' ##' \item{'lc'}{ARGOS location class,}
##' \item{'lon'}{observed longitude,}
##' \item{'lat'}{observed latitude.}
##' }
##'
##' @title Random Walk with autocorrelation Filter
##' @param d a data frame of observations (see details)
##' @param ts time-step in hours (used only when fitting to Argos data)
##' @param nu degrees of freedom for t-distributed measurement error (used only when fitting to Argos data)
##' @param span degree of loess smoothing for location state starting values (used only when fitting to Argos data)
##' @param optim numerical optimizer
##' @param verbose report progress during minimization
##' @return a list with components
##' \item{\code{fitted}}{a dataframe of fitted locations}
##' \item{\code{par}}{model parameter summmary}
##' \item{\code{data}}{input dataframe}
##' \item{\code{tmb}}{the tmb object}
##' \item{\code{opt}}{the object returned by the optimizer}
##' @importFrom TMB MakeADFun sdreport newtonOption
##' @importFrom dplyr %>% group_by rowwise do
##' @importFrom ssmTMB argos2tmb amfCRAWL
##' @export
mpm <- function(d,
                ts = 12,
                nu = 5,
                span = 0.4,
                optim = c("nlminb", "optim"),
                verbose = FALSE) {

  optim <- match.arg(optim)

  A <- length(unique(d$id))
  idx <- c(0, cumsum(as.numeric(table(d$id))))

  if ("lc" %in% names(d)) {
    model = "rw_err"
    compile("tmb/gamma_err.cpp")
    dyn.load(dynlib("tmb/gamma_err"))
  }
  else {
    model = "rw"
    compile("tmb/gamma.cpp")
    dyn.load(dynlib("tmb/gamma"))
  }

  ## REMOVE ALL CODE ASSOCIATED WITH GAMMA_COV...SUPERCEDED BY MPMM

  switch(model,
         rw = {
           data <- with(d,
                        list(
                          x = cbind(lon, lat),
                          A = A,
                          idx = idx
                        ))
         },
         cov = {
           data <- with(d,
                        list(
                          x = cbind(lon, lat),
                          A = A,
                          idx = idx,
                          M = M
                        ))
         },
         rw_err = {
           data.tmb <- d %>%
             group_by(id) %>%
             do(tmb = argos2tmb(., tstep = ts / 24, amf = amfCRAWL()))


           dofun <- function(dd) {
           ## Initialize parameters from loess smooths of the track
           fit.lon <- loess(
             lon ~ as.numeric(date),
             data = dd$tmb$obs,
             span = span,
             na.action = "na.exclude",
             control = loess.control(surface = "direct")
           )
           fit.lat <- loess(
             lat ~ as.numeric(date),
             data = dd$tmb$obs,
             span = span,
             na.action = "na.exclude",
             control = loess.control(surface = "direct")
           )

           ## Predict track for x starting values
           xs <-
             cbind(predict(fit.lon, newdata = data.frame(date = as.numeric(dd$tmb$ts))),
                   predict(fit.lat, newdata = data.frame(date = as.numeric(dd$tmb$ts))))

           with(dd,
                list(
                  id = id,
                  y = tmb$y,
                  K = tmb$K,
                  idx = tmb$idx,
                  ws = tmb$ws,
                  ts = tmb$ts,
                  dt = tmb$dt,
                  xs = xs
                  )
                )
           }

            prep <- data.tmb %>%
              rowwise() %>%
              do(tmp = dofun(.))

            Xidx <- c(0, cumsum(sapply(prep$tmp, function(.) length(.$ts))))
            Yidx <- c(0, cumsum(sapply(prep$tmp, function(.) nrow(.$y))))

            if(A > 1) {
              nidx <- lapply(prep$tmp, function(.) .$idx)
              idx <- c(unlist(nidx[[1]]), unlist(sapply(2:A, function(i) Xidx[i] + nidx[[i]])))
            }
            else {
              idx <- prep$tmp[[1]]$idx
            }
            xs <- do.call(rbind, lapply(prep$tmp, function(.) .$xs))

            data <- list(
              y = do.call(rbind, lapply(prep$tmp, function(.) .$y)),
              idx = idx,
              w = unlist(sapply(prep$tmp, function(.) .$ws)),
              K = do.call(rbind, lapply(prep$tmp, function(.) .$K)),
              nu = nu,
              A = A,
              Xidx = Xidx,
              Yidx = Yidx
            )
         }
      )

  switch(model,
         rw = {
           parameters <- list(
             lg = rep(1, dim(d)[1]),
             log_sigma = c(1, 1),
             log_sigma_g = 2
           )
         },
         rw_err = {
           parameters <- list(
             x = xs,
             lg = rep(1, max(Xidx)),
             log_sigma = c(1, 1),
             log_sigma_g = 1,
             log_tau = c(0, 0)
           )
         },
         cov = {
           n <- dim(M)[2]
           p <- dim(M)[1]
           parameters <- list(
             B = rep(1, n),
             lg = rep(1, p),
             log_sigma = c(1, 1),
             log_sigma_g = 2,
             log_sigma_u = 2,
#             log_sigma_b = 2,
             u = rep(0, A)
#             b = rep(0, A)
           )
         })

  ## TMB - create objective function
  switch(model,
         rw = {
           obj <-
             MakeADFun(
               data,
               parameters,
               random = c("lg"),
               DLL = "gamma",
               silent = !verbose
             )
         },
         rw_err = {
           obj <-
             MakeADFun(
               data,
               parameters,
               random = c("x", "lg"),
               DLL = "gamma_err",
               silent = !verbose
             )
         },
         cov = {
           obj <-
             MakeADFun(
               data,
               parameters,
               random = c("u", "lg"),
               DLL = "gamma_cov",
               silent = !verbose
             )
         })

  obj$env$inner.control$trace <- verbose
  obj$env$tracemgc <- verbose

  obj$control <- list(trace = 0,
                      reltol = 1e-12,
                      maxit = 500)
  obj$hessian <- TRUE
  newtonOption(obj, smartsearch = TRUE)

  ## Minimize objective function
  opt <- suppressWarnings(switch(
    optim,
    nlminb = nlminb(obj$par, obj$fn, obj$gr),
    optim = do.call("optim", obj)
  ))

  ## Parameters, states and the fitted values
  rep <- sdreport(obj)
  fxd <- summary(rep, "report")
  fxd_log <- summary(rep, "fixed")
  rdm <- summary(rep, "random")

  lg <- rdm[rownames(rdm) %in% "lg", ]

  if (model == "cov") {
    u <- rdm[rownames(rdm) %in% "u", ]
#    b <- rdm[rownames(rdm) %in% "b", ]
#    re <- data.frame(u = u[, 1], b = b[, 1]) %>%
    re <- data.frame(u = u[, 1]) %>%
      tbl_df()
  }
  if (model == "rw" || model == "cov") {
    fitted <- data_frame(
      id = d$id,
      date = d$date,
      g = plogis(lg[, 1]),
      g.se = lg[, 2]
    )
  }

  if (model == "rw_err") {
    x <- rdm[rownames(rdm) %in% "x",]

    fitted <- data_frame(
      id = rep(sapply(prep$tmp, function(.) .$id), sapply(prep$tmp, function(.) length(.$ts))),
      date = do.call(c, sapply(prep$tmp, function(.) .$ts)),
      lon = x[1:(nrow(x) / 2), 1],
      lat = x[(nrow(x) / 2 + 1):nrow(x), 1],
      lon.se = x[1:(nrow(x) / 2), 2],
      lat.se = x[(nrow(x) / 2 + 1):nrow(x), 2],
      g = plogis(lg[, 1]),
      g.se = lg[, 2]
    )
  }

  if (optim == "nlminb") {
    aic <- 2 * length(opt[["par"]]) + 2 * opt[["objective"]]
  }
  else {
    aic <- 2 * length(opt[["par"]]) + 2 * opt[["value"]]
  }

  if (model == "cov") {
    row.names(fxd)[3:4] <- c("sigma_lon", "sigma_lat")
    row.names(fxd)[5] <- "Intercept"
    row.names(fxd)[6:nrow(fxd)] <-
      attr(terms(formula), "term.labels")

    list(
      fitted = fitted,
      par = fxd,
      re = re,
      data = d,
      tmb = obj,
      opt = opt,
      rep = rep,
      aic = aic
    )

  }
  else {
    row.names(fxd)[2:3] <- c("sigma_lon", "sigma_lat")
    if(model == "rw_err") row.names(fxd)[4:5] <- c("tau_lon", "tau_lat")

    list(
      fitted = fitted,
      par = fxd,
      data = d,
      tmb = obj,
      opt = opt,
      rep = rep,
      aic = aic
    )

  }


}
