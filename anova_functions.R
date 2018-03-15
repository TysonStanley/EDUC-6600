

set.seed(42)
z <- data.frame(a1 = c(rnorm(100,2), rnorm(100,1),rnorm(100,0)),
                b1 = rep(c("A", "B", "C"), each = 100),
                c1 = factor(rbinom(300, 1, .5)),
                ID = 1:300,
                a2 = c(rnorm(100,2), rnorm(100,1),rnorm(100,0)),
                a3 = c(rnorm(100,2), rnorm(100,1),rnorm(100,0)))



anoRM <- function(data, dep, rm_lab = NULL, covariates = 1, ss = "III", adjust = "none"){
  
  data <- data.frame(data)
  
  if (is.null(rm_lab)){
    repeated <- paste0("time", 1:length(dep))
  } else {
    repeated <- rm_lab
  }
  
  dep1 <- data[, dep] %>% as.matrix
  form <- as.formula(paste("dep1 ~", paste(covariates, collapse = " + ")))
  mlm1 <- lm(form, data = data)
  mlm1.aov <- car::Anova(mlm1, 
                         idata = data.frame(repeated),
                         idesign = ~repeated, 
                         type = ss)
  main <- summary(mlm1.aov, multivariate=FALSE)
  
  ps <- diff <- vector("numeric", length(dep)-1)
  nam <- vector("character", length(dep)-1)
  for (i in 1:(length(dep)-1)){
    nam[i]  <- paste(dep[i:(i+1)], collapse = " vs. ")
    diff[i] <- t.test(data[, dep[i]], 
                       data[, dep[i+1]],
                       paired = TRUE)$estimate
    ps[i]   <- t.test(data[, dep[i]], 
                      data[, dep[i+1]],
                      paired = TRUE)$p.value
  }
  nam[length(nam)+1]   <- paste(dep[c(1, length(dep))], collapse = " vs. ")
  diff[length(diff)+1] <- t.test(data[, dep[1]], 
                                 data[, dep[length(dep)]],
                                 paired = TRUE)$estimate
  ps[length(ps)+1]     <- t.test(data[, dep[1]], 
                                 data[, dep[length(dep)]],
                                 paired = TRUE)$p.value
  
  post <- switch(adjust,
                 "none" = data.frame(comp = nam, diff = diff, pvalue = ps),
                 "bonf" = data.frame(comp = nam, diff = diff, pvalue = ps/length(ps)))
  
  means <- lsmeans::lsmeans(mlm1, 
                            specs = c(paste(attr(mlm1$terms, "variables"))[3:(length(covariates)+2)]))

  list("ANOVArm" = main, "PostHoc" = post, "LSmeans" = means)
}
ano <- function(formula, data, ss = "III", posthoc = "tukey", adjust = "none", contrasts = NULL){
  call <- match.call()
  
  fit <- aov(formula, data, contrasts = contrasts)
  
  an_tab <- car::Anova(fit, type = ss)
  post <- switch(posthoc,
                 tukey = TukeyHSD(fit),
                 t.test = pairwise.t.test(fit, p.adj = adjust))
  
  out <- list("ANOVA" = an_tab, 
              "PostHoc" = post,
              "aovCall" = call,
              "residuals" = fit$residuals,
              "preds" = fit$fitted.values,
              "contrasts" = fit$contrasts,
              "levels" = fit$xlevels)
  class(out) = "ano"
  out
  
}
anoRM(z, dep = c("a1", "a2", "a3"), covariates = c("b1", "c1"))

print.ano <- function(x, ...){
  cat("-----\n",
      "Analysis of Variance\n",
      " - Formula:", deparse(x$aovCall$formula), "\n",
      " -", length(x$levels), "factors (", paste0(sapply(x$levels, length), collapse = " x "), ")\n",
      " - Contrasts:", paste0(x$contrasts, collapse = ", "), "\n",
      "\n\n")
  x$ANOVA %>%
    data.frame %>%
    round(5) %>% print
  
  cat("\n\n")
  print(x$PostHoc)
  
}

ano(a ~ b + c, x)

aov_car

library(afex)
z %>%
  aov_ez("ID", dv = "a1", between = "b1", data = .) 
z %>%
  aov_4(a1 ~ b1 + (1|ID), data = .) 


z %>%
  gather("time", "value", a1, a2, a3) %>%
  data.frame %>%
  aov_4(value ~ 1 + (time|ID),
        data = .)

z %>%
  gather("time", "value", a1, a2, a3) %>%
  lme4::lmer(value ~ time + (1|ID),
             data = .)



library(tidyverse)

glimpse(z)
