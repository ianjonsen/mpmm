##' @title plot
##' @description Visualise fixed and random covariate relationships from
##' an mpmm fit object
##' @param x an mpmm fit object
##' @param label add id labels to random effects
##' @param lwd a vector of regression line widths (random effect, fixed effects)
##' @param ... additional arguments to be ignored
##'
##' @importFrom lme4 nobars
##' @importFrom ggplot2 ggplot geom_line aes xlab ylab theme_bw theme ylim xlim
##' @importFrom ggplot2 element_text facet_wrap element_blank
##' @importFrom dplyr left_join mutate group_by %>% select arrange bind_cols
##' @importFrom stats plogis
##' @importFrom tidyr pivot_longer everything
##' @importFrom wesanderson wes_palette
##' @method plot mpmm
##' @export
plot.mpmm <- function(x, label = FALSE, lwd = c(0.25, 0.75), ...) {

  if (length(list(...)) > 0) {
    warning("additional arguments ignored")
  }

## set up wesanderson palette for unlabeled mixed-effects
wpal <- wes_palette("Darjeeling2", n = 5, "discrete")

nval <- 50
terms <- attr(terms(nobars(x$formula)), "term.labels")
n <- length(terms)
nid <- nrow(x$re)
fe.rng <- sapply(1:n, function(i) range(x$fr[, terms[i]]))
xt <- sapply(1:n, function(i) seq(fe.rng[1,i], fe.rng[2,i], l = nval))

## get individual ranges for random effects
fr.lst <- split(x$fr, x$fr$id)
re.rng <- lapply(fr.lst, function(x) {
  x <- x %>% select(-id)
  c(apply(x, 2, range))
}) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  data.frame(id = row.names(.), ., row.names = NULL)
names(re.rng)[-1] <- paste0(rep(terms, each=2), rep(c(".min",".max"),length(terms)))

f_int <- x$par["Intercept","Estimate"]
betas <- sapply(1:n, function(i) x$par[terms[i],"Estimate"])

fxd <- sapply(1:n, function(i) {
  if(n > 2) plogis(f_int + betas[i] * xt[, i] + sum(betas[-i] * apply(xt[, -i], 2, mean)))
  else plogis(f_int + betas[i] * xt[, i] + sum(betas[-i] * mean(xt[, -i])))
})

  re_ints <- f_int + x$re$`(Intercept)`
  k <- length(re_ints)

  rnms <- names(x$re)[!names(x$re) %in% c("id","(Intercept)")]
  # check for fixed terms not in random terms
  rmiss <- which(!terms %in% rnms)
  rpos <- which(terms %in% rnms)

  bs <- matrix(0, ncol = n, nrow = k)
  bs[, rmiss] <- 0
  bs[, rpos] <- unlist(x$re[, rnms])

   re_betas <- sapply(1:n, function(i) {
     (betas[i] + bs[, i])
   })

   xt.re <- vector(mode = "list", length = k)
   xt.re <- lapply(1:k, function(i){
     z <- list()
     q <- 0
     for(j in seq(2, 2*n, by=2)) {
       q <- q + 1
       z[[q]] <- seq(re.rng[i,j], re.rng[i,j+1], l=nval)
     }
     z <- do.call(cbind, z) %>%
       as.data.frame(.)
     #    names(z) <- terms
     z %>% mutate(id = re.rng$id[i]) %>% select(id, everything())
   })

   foo <- lapply(1:k, function(j) {
     z <- list()
     for(i in 1:n) {
       if(n > 2) z[[i]] <- plogis(re_ints[j] + xt.re[[j]][, 1+i] * re_betas[j,i] + sum(re_betas[j, -i] * apply(xt.re[[j]][, -c(1,i+1)], 2, mean)))
       else z[[i]] <- plogis(re_ints[j] + xt.re[[j]][, 1+i] * re_betas[j,i] + sum(re_betas[j, -i] * mean(xt.re[[j]][, -c(1,i+1)])))
     }
     do.call(cbind, z)
   })

   re <- lapply(1:k, function(i) {
       as.data.frame(foo[[i]]) %>%
       pivot_longer(cols = everything(), names_to = "predictor", values_to = "g") %>%
       arrange(predictor) %>%
       mutate(id = rep(re.rng$id[i])) %>%
       select(id, everything())
   }) %>%
     do.call(rbind, .)

   xt.re <- lapply(xt.re, function(x) {
     x %>%
       as.data.frame() %>%
       pivot_longer(cols = -id, names_to = "predictor", values_to = "x")
   }) %>%
     do.call(rbind, .) %>%
     arrange(id, predictor) %>%
     select(-id, -predictor)

   re.dat <- bind_cols(re, xt.re) %>%
     mutate(predictor = factor(predictor, labels = terms))

   xs <- as.data.frame(xt) %>%
     pivot_longer(cols = everything(), names_to = "predictor", values_to = "x") %>%
     mutate(predictor = factor(predictor, labels = terms)) %>%
     arrange(predictor)

   fe.dat <- as.data.frame(fxd) %>%
     pivot_longer(cols = everything(), names_to = "predictor", values_to = "g") %>%
     arrange(predictor) %>%
     select(-predictor)

   fe.dat <- data.frame(xs, fe.dat)

gg <- ggplot() + theme(axis.text = element_text(size = 14),
                          axis.title = element_text(size = 20)) +
     geom_line(data = fe.dat,
               aes(x, g),
               size = lwd[2],
               colour = wpal[2])

if(label) {
  gg <- gg +
    geom_line(data = re.dat,
               aes(x, g, group = id, colour = id), size = lwd[1], alpha = 0.7)
  } else {
    gg <- gg +
      geom_line(data = re.dat,
              aes(x, g, group = id), size = lwd[1], colour = wpal[3], alpha = 0.7)
  }

    gg <- gg +
      ylab(expression(gamma[t])) +
      facet_wrap(~ predictor, scales = "free_x") +
      ylim(0,1) +
      xlab(label = element_blank()) +
      theme_bw()

    gg
}



