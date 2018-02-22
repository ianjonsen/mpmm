##' Plot fixed and random components
##'
##' @title plot
##' @param x a fitted object of class mpmm
##' @param page 1 = plot all terms on a single page, 0 otherwise
##'
##'
##' @importFrom lme4 nobars
##' @importFrom ggplot2 ggplot geom_line aes xlab ylab theme_bw element_text
##' @importFrom gridExtra grid.arrange
##' @method plot mpmm
##' @export
plot.mpmm <- function(m, page = 1) {

terms <- attr(terms(nobars(m$formula)), "term.labels")

rng <- sapply(1:length(terms), function(i) range(m$fr[, terms[i]]))
xt <- sapply(1:length(terms), function(i) seq(rng[1,i], rng[2,i], l = 200))

betas <- sapply(1:dim(xt)[2], function(i) m$par[terms[i],"Estimate"])

fxd <- sapply(1:dim(xt)[2], function(i) {
  fxd_int <- m$par["Intercept","Estimate"]
  plogis(fxd_int + betas[i] * xt[, i])
})

if(dim(m$re)[2] == 2) {
  ## intercept only random effect
  re_ints <- m$par["Intercept", "Estimate"] + m$re$`(Intercept)`
  k <- length(re_ints)

  p <- lapply(1:length(terms), function(j) {
    re1 <- sapply(1:k, function(i) {
      plogis(re_ints[i] + betas[j] * xt[, j])
    })
    pdat <- data.frame(x = xt[, j], g = re1)
    pdat <-
      reshape2::melt(
        pdat,
        id.vars = "x",
        value.name = "g",
        variable.name = "Intercept"
      )

    pdat1 <- data.frame(x=xt[,j], y=fxd[,j])
    ggplot() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 20)) +
      geom_line(aes(pdat$x, pdat$g, group = pdat$Intercept),
                size = 0.2,
                colour = "dodgerblue") +
      geom_line(aes(pdat1$x, pdat1$y), size = 1, colour = "firebrick") +
      xlab(terms[j]) +
      ylab(expression(gamma[t])) +
      ylim(0,1) +

      theme_bw()
  })
} else {
  re_ints <- m$par["Intercept", "Estimate"] + m$re$`(Intercept)`
  k <- length(re_ints)
  ## below doesn't work if 2 fixed terms but only 1 random slope term...
  rnms <- names(m$re)[!names(m$re) %in% c("id","(Intercept)")]
  re_bs <- sapply(1:length(rnms), function(i) {
    m$par[rownames(m$par) %in% terms, "Estimate"][i] + m$re[, rnms]
  })

  p <- lapply(1:length(terms), function(j){
    re1 <- sapply(1:k, function(i){
      plogis(re_ints[i] + re_bs[[j]][i] * xt[, j])
    })
    pdat <- data.frame(x = xt[, j], g = re1)
    pdat <- reshape2::melt(
      pdat,
      id.vars = "x",
      value.name = "g",
      variable.name = "re"
    )
    pdat1 <- data.frame(x=xt[,j], y=fxd[,j])
    ggplot() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 20)) +
      geom_line(aes(pdat$x, pdat$g, group = pdat$re),
                size = 0.2,
                colour = "dodgerblue") +
      geom_line(aes(pdat1$x, pdat1$y), size = 1, colour = "firebrick") +
      xlab(terms[j]) +
      ylab(expression(gamma[t])) +
      ylim(0,1) +
      theme_bw()
  })

}
n <- length(p)
if (n > 1 && page == 1) {
  grid.arrange(grobs = p, nrow = floor(sqrt(n)))
} else if (n == 1 || page == 0) {
  p
}
}



