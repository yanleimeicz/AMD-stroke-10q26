# patched for loo >= 2.8 (loo_compare returns a data frame with a 'model' column)
in_sample_elpd_loo <- function(X, fits, variants, nsamps=1000){

  stopifnot(inherits(X, "cause_data"))
  if(!missing(variants)){
    if(!all(variants %in% X$snp)){
      stop("Not all `variants` are in data.", call.=FALSE)
    }
    X <- filter(X, snp %in% variants)
    X <- new_cause_data(X)
  }
  stopifnot(length(fits) >=2)

  k <- length(fits)
  for(x in fits) stopifnot(inherits(x, "cause_post"))
  mods <- expand.grid(m1 = 1:k, m2=1:k)
  mods <- mods %>% filter(m1 < m2) %>%
                   mutate(delta_elpd = NA, se_delta_elpd = NA)
  if(!is.null(names(fits))){
    mods$model1 <- names(fits)[mods$m1]
    mods$model2 <- names(fits)[mods$m2]
  }else{
    mods$model1 <- m1
    mods$model2 <- m2
  }
  #For each data point, calculate log likelihood under posterior
  lls <- lapply(fits, FUN=function(fit){
               if(is.null(fit$joint_post)){
                 llmat <- loglik_loo(0, 0, 0,
                                     fit$rho, fit$mix_grid$S1,
                                     fit$mix_grid$S2, fit$mix_grid$pi,
                                     X$beta_hat_1, X$beta_hat_2, X$seb1, X$seb2)
                 return(llmat)
               }
               samps <-samp_from_grid(fit$joint_post, c("q", "gamma", "eta"), nsamps)
               llmat <- t(loglik_loo(samps$gamma, samps$eta, samps$q, fit$rho,
                                     fit$mix_grid$S1, fit$mix_grid$S2, fit$mix_grid$pi,
                                     X$beta_hat_1, X$beta_hat_2, X$seb1, X$seb2))
               return(llmat)
              })

  loos <- lapply(lls, function(llmat){
    if(ncol(llmat) == 1) return(NULL)
    loo(llmat, r_eff = rep(1, nrow(X)))
  })

  for(j in seq(nrow(mods))){
    i1 <- mods$m1[j]
    i2 <- mods$m2[j]
    if(is.null(loos[[i1]])){
      diff <- lls[[i1]] - loos[[i2]]$pointwise[,1]
      est <- sum(diff)
      se <- sd(diff)*sqrt(nrow(X))
    }else if(is.null(loos[i2])){
      diff <- loos[[i1]]$pointwise[,1] - lls[[i2]]
      est <- sum(diff)
      se <- sd(diff)*sqrt(nrow(X))
    }else{
      comp <- as.data.frame(loo_compare(loos[[i1]], loos[[i2]]))
      rn <- if ("model" %in% names(comp)) as.character(comp$model) else rownames(comp)
      dd <- as.numeric(comp[2, "elpd_diff"]); ss <- as.numeric(comp[2, "se_diff"])
      if(rn[1]=="model2"){
        est <- dd
        se <- ss
      }else{
        est <- -1*dd
        se <- ss
      }
    }
    mods$delta_elpd[j] <- est
    mods$se_delta_elpd[j] <- se
  }
  mods <- mods %>% mutate(z = delta_elpd/se_delta_elpd) %>%
          select(model1, model2, delta_elpd, se_delta_elpd, z)
  class(mods) <- c("cause_elpd", "data.frame")
  ret <- list(mods = mods, loos = loos)
  return(ret)
}
