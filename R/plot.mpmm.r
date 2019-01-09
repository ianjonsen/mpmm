##' Visualise fixed and random relationships
##'
##' @title plot
##' @param m a fitted object of class mpmm
##' @param label add id labels to random effects (messy)
##' @param page 1 = plot all terms on a single page, 0 otherwise
##'
##' @importFrom lme4 nobars
##' @importFrom ggplot2 ggplot geom_line aes xlab ylab theme_bw theme ylim xlim element_text geom_text
##' @importFrom gridExtra grid.arrange
##' @importFrom stats plogis
##' @importFrom reshape2 melt
##' @method plot mpmm
##' @export
plot.mpmm <- function(m, label = FALSE, page = 1) {

terms <- attr(terms(nobars(m$formula)), "term.labels")
n <- length(terms)
rng <- sapply(1:n, function(i) range(m$fr[, terms[i]]))
xt <- sapply(1:n, function(i) seq(rng[1,i], rng[2,i], l = 200))
xt.mn <- apply(xt, 2, mean)

f_int <- m$par["Intercept","Estimate"]
betas <- sapply(1:n, function(i) m$par[terms[i],"Estimate"])

fxd <- sapply(1:n, function(i) {
  if(n > 2) {
    plogis(f_int + betas[i] * xt[, i] + betas[-i] %*% xt.mn[-i])
  } else if(n > 1){
    plogis(f_int + betas[i] * xt[, i] + betas[-i] * xt.mn[-i])
    } else {
    plogis(f_int + betas * xt)
  }
})

if(dim(m$re)[2] == 2) {
  ## intercept only random effect
  re_ints <- m$par["Intercept", "Estimate"] + m$re$`(Intercept)`
  k <- length(re_ints)

  re <- lapply(1:n, function(j) {
      if(n > 1) {
        plogis(outer(betas[j] * xt[, j] + betas[-j] * xt.mn[-j], re_ints, FUN = "+"))
      } else {
        plogis(outer(betas * xt, re_ints, FUN = "+"))
      }
})
  p <- lapply(1:n, function(j) {
    pdat <- data.frame(x = xt[, j], g = re[[j]])
    pdat <-
      melt(
        pdat,
        id.vars = "x",
        value.name = "g",
        variable.name = "Intercept"
      )
    pdat <- data.frame(id = rep(as.character(m$re$id), each = 200), pdat)
    pdat1 <- data.frame(x=xt[,j], y=fxd[,j])
    if(label) pdat.lab <- pdat[seq(1, 2001, by = 200),]

    gg <- ggplot() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 20)) +
      geom_line(aes(pdat$x, pdat$g, group = pdat$Intercept),
                size = 0.2,
                colour = "dodgerblue")
    if(label) {
      gg <- gg +
        geom_text(data = pdat.lab, aes(x, g, label = id), hjust = 0, size = 3)
    }
    gg <- gg +
      geom_line(aes(pdat1$x, pdat1$y), size = 1, colour = "firebrick") +
      xlab(terms[j]) +
      ylab(expression(gamma[t])) +
      ylim(0,1) +
      theme_bw()
  })
} else {
  ## intercept + slope(s) random effects
  re_ints <- f_int + m$re$`(Intercept)`
  k <- length(re_ints)

  rnms <- names(m$re)[!names(m$re) %in% c("id","(Intercept)")]
  # check for fixed terms not in random terms
  rmiss <- which(!terms %in% rnms)
  rpos <- which(terms %in% rnms)

  bs <- matrix(0, ncol = n, nrow = k)
  bs[, rmiss] <- 0
  bs[, rpos] <- unlist(m$re[, rnms])

  re_betas <- sapply(1:n, function(i) {
    (betas[i] + bs[, i])
  })

  re <- lapply(1:n, function(j){
    if(n > 2) {
      plogis(re_ints + re_betas[,j] %o% xt[,j] + as.vector(re_betas[,-j] %*% xt.mn[-j]))
    } else if(n > 1) {
      plogis(re_ints + re_betas[,j] %o% xt[,j] + re_betas[,-j] * xt.mn[-j])
    } else {
      plogis(re_ints + re_betas[,j] %o% xt[,j])
    }
  })
  p <- lapply(1:n, function(j) {
    pdat <- data.frame(x = xt[, j], g = t(re[[j]]))
    pdat <- melt(
      pdat,
      id.vars = "x",
      value.name = "g",
      variable.name = "re"
    )
    pdat <- data.frame(id = rep(as.character(m$re$id), each = 200), pdat)
    pdat.f <- data.frame(x=xt[,j], y=fxd[,j])
    if(label) pdat.lab <- pdat[seq(1, 2001, by = 200),]
    gg <- ggplot() +
      theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 20)) +
    geom_line(data = pdat,
                aes(x, g, group = re),
                size = 0.2,
                colour = "dodgerblue")
    if(label) {
      gg <- gg +
        geom_text(data = pdat.lab, aes(x, g, label = id), hjust = 0, size = 3)
    }
    gg <- gg +
      geom_line(data = pdat.f, aes(x, y), size = 1, colour = "firebrick") +
      xlab(terms[j]) +
      ylab(expression(gamma[t])) +
      ylim(0,1) +
      theme_bw()
    gg
})

}

if (n > 1 && page == 1) {
  grid.arrange(grobs = p, nrow = floor(sqrt(n)))
} else if (n == 1 || page == 0) {
  p
}

}



