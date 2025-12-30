library(viridis)
#unpack from solution saved on disk
burn <- 367 
file <- "SOL_linear"
source('unpack.r')

cols  <- turbo(p$l_species)
ylims <- c(-0.04,0.04)
dts   <- c(5,10,50,100,200)
bmat  <- as.matrix(b)

#######################
## no noise 
#######################
par(mfrow=c(3,2),mar=c(2,2,2,2))
matplot(time[-1],diff(bmat),type='l',lty=2,col=cols,ylim=ylims)
matplot(time,bdot,type='l',lty=1,col=cols,add=TRUE)

for(i in 1:length(dts)){
    obs_dt    <- dts[i]
    obs_index <- seq(1,length(time),by=obs_dt)
    obs_times <- time[obs_index]

    matplot(time,bdot,type='l',lty=1,col=cols,ylim=ylims)
    #matplot(obs_times[-1],diff(as.matrix(b[obs_index,]))/obs_dt,type='l',lty=2,col=cols,add=TRUE)
    matplot(obs_times[-1] - obs_dt/2,diff(bmat[obs_index,])/obs_dt,type='l',lty=2,col=cols,add=TRUE)
}



##########################
## with noise 
##########################
noise <- matrix(rnorm(length(c(bmat)),mean=0,sd=0.05*c(bmat)),
                byrow=FALSE,
                ncol=p$l_species,
                nrow=length(time))
b_noise <- bmat + noise

#par(mfrow=c(1,1),mar=c(2,2,2,2))
#matplot(b_noise,type='l',col=cols)

ylims <- c(-0.1,0.1)

par(mfrow=c(3,2),mar=c(2,2,2,2))
matplot(time[-1],diff(b_noise),type='l',lty=2,col=cols,ylim=ylims)
matplot(time,bdot,type='l',lty=1,col=cols,add=TRUE)

for(i in 1:length(dts)){
    obs_dt    <- dts[i]
    obs_index <- seq(1,length(time),by=obs_dt)
    obs_times <- time[obs_index]

    matplot(time,bdot,type='l',lty=1,col=cols,ylim=ylims)
    #matplot(obs_times[-1],diff(as.matrix(b[obs_index,]))/obs_dt,type='l',lty=2,col=cols,add=TRUE)
    matplot(obs_times[-1] - obs_dt/2, diff(b_noise[obs_index,])/obs_dt,type='l',lty=2,col=cols,add=TRUE)
}



