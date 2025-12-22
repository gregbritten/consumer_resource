setwd("~/dropbox/working/CONSUMER_RESOURCE/consumer_resource/experiments/npz")

mod <- stan_model("npz.stan")

mcmc <- sampling(mod,data=data,open_progress=TRUE)

