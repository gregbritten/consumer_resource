library(deSolve)
library(rstan)
library(viridis)
options(mc.cores = parallel::detectCores())


theta <- list(vmax       = 0.075,
              nuthalfsat = 0.3,
              graz       = 0.02,
              mort_p     = 0.02,
              mort_z     = 0.03,
              irr        = 0.8)

x <- c(N = 0.6,
       P = 0.15,
       Z = 0.23)

T  <- 2*365
dt <- 1
t  <- seq(0,T,dt)

i_n = 1
i_p = 2
i_z = 3


dxdt <- function(t,x,theta){
    with(as.list(c(x,theta)),{
        light   = 1 + 0.5*(irr*sin(pi*((t-81.25)/182.5)) - irr)
        growth  = (vmax*N/(nuthalfsat + N))*light*P
        grazing = graz*P*Z
        ploss   = mort_p*P
        zloss   = mort_z*Z*Z
        
        list(c(-growth + ploss + zloss,
               growth - grazing - ploss,
               grazing - zloss)) })}

               

x <- as.data.frame(ode(y=x, times=t, func=dxdt, parms=theta))

iobs <- sort(sample(1:length(t), 20)) 
tobs <- t[iobs]
iobsvar <- c(i_p,i_z)

sigma <- c(0.03,0.03)

obs <- cbind(x$P[iobs] + rnorm(length(iobs),sd=sigma[1]),
             x$Z[iobs] + rnorm(length(iobs),sd=sigma[2]))


light <- 1 + 0.5*(theta$irr*sin(pi*((t-81.25)/182.5)) - theta$irr)


plot(t,light/10,type='l',ylim=c(0,0.8),lty=2,col='orange',xlab='time (days)')
lines(t,x$N,col='blue')
lines(t,x$P,col='dark green')
lines(t,x$Z,col='red')

points(tobs,obs[,length(iobsvar)],col='dark green')
points(tobs,obs[,length(iobsvar)],col='red')


#####################################
## package data for Stan ############
#####################################
data <- list(obs=obs,
             nobs=length(tobs),
             tobs=tobs,
             nobsvar=ncol(obs),
             iobsvar=iobsvar)

