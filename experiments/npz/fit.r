setwd("~/dropbox/working/CONSUMER_RESOURCE/consumer_resource/experiments/npz")

mod <- stan_model("npz.stan")

data <- list(obs=obs,
             nobs=length(tobs),
             tobs=tobs,
             nobsvar=ncol(obs),
             iobsvar=iobsvar)

mcmc <- sampling(mod,data=data,open_progress=TRUE)