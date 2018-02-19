plot_re <- function(m) {

n <- dim(m$par)[1]

terms <- attr(terms(lme4::nobars(m$formula)), "term.labels")

rng <- sapply(1:length(terms), function(i) range(m$data[, terms[i]]))
x <- sapply(1:length(terms), function(i) seq(rng[1,i], rng[2,i], l = 200))

bs <- sapply(1:dim(x)[2], function(i) m$par[terms[i],"Estimate"])

fxd <- sapply(1:dim(x)[2], function(i) {
  fxd_int <- m$par["Intercept","Estimate"]
  plogis(fxd_int + bs[i] * x[, i])
})

if(dim(m$re)[2] == 2) {
  ## intercept only random effect
  re_ints <- m$par["Intercept", "Estimate"] + m$re$`(Intercept)`
  k <- length(re_ints)
  
  p <- lapply(1:length(terms), function(j) {
    re1 <- sapply(1:k, function(i) {
      plogis(re_ints[i] + bs[j] * x[, j])
    })
    pdat <- data.frame(x = x[, j], g = re1)
    pdat <-
      reshape2::melt(
        pdat,
        id.vars = "x",
        value.name = "g",
        variable.name = "Intercept"
      )
    
    ggplot() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 20)) +
      geom_line(aes(pdat$x, pdat$g, group = pdat$b0),
                size = 0.2,
                colour = "dodgerblue") +
      geom_line(aes(x[, j], fxd[, j]), size = 1, colour = "firebrick") +
      xlab(terms[j]) + ylab(expression(gamma[t])) +
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
      plogis(re_ints[i] + re_bs[[j]][i] * x[, j])
    })
    pdat <- data.frame(x = x[, j], g = re1)
    pdat <- reshape2::melt(
      pdat,
      id.vars = "x",
      value.name = "g",
      variable.name = "re"
    )

    ggplot() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 20)) +
      geom_line(aes(pdat$x, pdat$g, group = pdat$re),
                size = 0.2,
                colour = "dodgerblue") +
      geom_line(aes(x[, j], fxd[, j]), size = 1, colour = "firebrick") +
      xlab(terms[j]) + ylab(expression(gamma[t])) +
      theme_bw()
  })
 
}
p

}



