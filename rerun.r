library(deSolve)
library(viridis)
library(fields)
source('src/cr_de.r')

set.seed(1)

file <- "SOL_ref"

#unpack
SOL    <- readRDS(paste0('results/',file,'.rds'))
p      <- SOL[['p']]

burn   <- 367
tf     <- length(p$times)
time   <- p$times[burn:tf]
b_ref  <- SOL[['sol']][burn:tf,2:(p$l_species+1)]
n_ref  <- t(apply(b_ref,1,function(x) p$Q*x))
n_T_ref<- rowSums(n_ref)
#r_ref <- SOL[['sol']][burn:tf,(p$l_species+2):(p$l_species+p$l_resources+1)]
#c     <- SOL[['sol']][burn:tf,(p$l_species+p$l_resources+2):(p$l_species+p$l_resources+1+p$l_resources)]
#xdot  <- SOL[['sol']][burn:tf,(p$l_species+p$l_resources+2+p$l_resources):(p$l_species+p$l_resources+1+2*p$l_resources+p$l_species)]
#bdot  <- xdot[,1:p$l_species]
#rdot  <- xdot[,(p$l_species+1):ncol(xdot)]


#######################################
#######################################

## re-run with new parameters
r0        <- rep(0.01,p$l_resources)
b0        <- rep(0.01,p$l_species)
x0        <- c(b0, r0)
names(x0) <- c(paste0("b", 1:p$l_species), paste0("r", 1:p$l_resources))

ii <- as.numeric(which(colMeans(n_ref)>1E8)[2])

## modify parameters
nfactors <- 5
factors  <- c(seq(0.75,1,length.out=3),seq(1,1.25,length.out=3))

SOL_mod=b_mod=n_mod <- list()

for(i in 1:length(factors)){
print(i)
    p_mod <- p
    #modify parameters
#    for(j in 1:p$l_species){
#        p_mod$Vmax[j,1] <- p$Vmax[j,1]*factors[i]
#        p_mod$K[j,1]    <- p$K[j,1]*factors[i]
#    }
    for(j in 1:p$l_resources){
        p_mod$Vmax[ii,j] <- p$Vmax[ii,j]*factors[i]
        #p_mod$K[ii,j]    <- p$K[ii,j]*factors[i]
    }

    SOL_mod[[i]]   <- ode(y=x0, times=p$times, func=cr_de, parms=p_mod, method="lsoda")
    b_mod[[i]]     <- SOL_mod[[i]][burn:tf,2:(p$l_species+1)]
    n_mod[[i]]     <- t(apply(b_mod[[i]],1,function(x) p$Q*x))
}

######################################################
## ANALYSIS ##########################################
######################################################

times <- 1:length(n_ref[,i])

## plot
cols <- turbo(length(factors))
par(mfrow=c(5,5),mar=c(2,2,2,2),oma=c(2,2,2,2))
cols <- turbo(length(factors))
for(i in 1:ncol(n_ref)){
    plot(times,n_ref[,i],type='l',col=cols[1])
    for(j in 2:length(factors)){
        lines(times,n_mod[[j]][,i],col=cols[j])
    }
}


#4x4 plot
pdf('~/dropbox/working/CONSUMER_RESOURCE/consumer_resource/plots/deviations.pdf',height=5,width=8)
is <- c(5,10,16,17)
par(mfrow=c(2,2),mar=c(2,2,0,1),oma=c(2,4,2,4),cex.lab=0.8,cex.axis=0.8)
cols <- turbo(length(factors))
k <- 1
ylims <- list()
ylims[[1]] <- c(1E8,6E8)
ylims[[2]] <- c(8E5,7E7)
ylims[[3]] <- c(1E2,1E9)
ylims[[4]] <- c(1E7,2E8)

for(i in is){
    plot(times/365,n_ref[,i],type='l',col='black',lwd=1,lty=2,log='y',ylim=ylims[[k]],
         ylab='',xlab='')
    for(j in 1:length(factors)){
        lines(times/365,n_mod[[j]][,i],col=cols[j])
    }
    if(k==2){
            fields::image.plot(z=factors,legend.only=TRUE)
            legend('bottomright',
                    c(expression(theta['0']),
                      expression(theta*'/'*theta['0'])),
                      lty=c(2,1),col=c('black',cols),bty='n')}
            #legend('bottomright',
             #       c(expression(theta['0']),
             #         expression(theta*'/'*theta['0']~'= 0.8'),
             #         expression(theta*'/'*theta['0']~'= 0.9'),
             #         expression(theta*'/'*theta['0']~'= 1.1')),
#                      expression(theta*'/'*theta['0']~'= 1.25')),
             #         lty=c(2,1,1,1),col=c('black',cols),bty='n')
    if(k%in%c(1,3)) mtext(side=2,'Absolute Abundance',line=2.5)
    k <- k + 1
    lines(times/365,n_mod[[20]][,i],lty=2 )
}
dev.off()






par(new=TRUE)
fields::image.plot(z=factors,legend.only=TRUE)


n10 <- n_mod[[1]] - n_ref
n20 <- n_mod[[2]] - n_ref



##################################
## generate observations 
##################################

e_T = e   <- list()
ss_T = ss <- numeric(length(factors))

for(i in 1:length(factors)){
    noise <- matrix(rnorm(length(c(n_ref)),mean=0,sd=0.0001*c(n_ref)),
                byrow=FALSE,
                ncol=ncol(n_ref),
                nrow=nrow(n_ref))
    n_noise   <- apply(n_mod[[i]] + noise, c(1,2), function(x) max(0,x))
    n_noise_T <- rowSums(n_noise) 
    e_T[[i]]  <- n_noise_T - n_T_ref
    ss_T[i]   <- sum(e_T[[i]]^2)

    e[[i]]  <- log10(n_noise+1) - log10(n_ref+1) 
    ss[i] <- sum(e[[i]]^2)
}


#sds <- c(0.1,0.2,0.5)
#fs  <- c(1,7,30)

sds <- seq(0.1,0.5,length.out=25)
fs  <- seq(1,30,length.out=26)

ll       <- array(NA,dim=c(length(fs),length(sds),length(factors)))
dll=d2ll <- matrix(NA, nrow=length(fs), ncol=length(sds)) 

for(i in 1:length(fs)){
    for(j in 1:length(sds)){
        for(k in 1:length(factors)){
            es <- e[[k]][seq(1,nrow(n_ref),fs[i]),]
            ll[i,j,k] <- sum(dnorm(es,mean=0,sd=sds[j],log=TRUE))
        }
        dll[i,j] <- abs(diff(ll[i,j,]))[19]
        d2ll[i,j] <- abs(diff(diff(ll[i,j,])))[19]
    }
}




pdf('~/dropbox/working/consumer_resource/consumer_resource/plots/expected_CI.pdf',height=4,width=5)
par(mfrow=c(1,1))
filled.contour(x=fs,y=sds*100,z=sqrt(1/d2ll)*1000)
mtext(side=1,"Sampling Period [days]",line=2.5,adj=0.1)
mtext(side=2,"Measurement noise [%]",line=2.5)
mtext(adj=0,"c) Expected estimation error")
dev.off()


pdf('~/dropbox/working/consumer_resource/consumer_resource/plots/likelihood_profile.pdf',height=3.75,width=8)
cols <- turbo(10)
par(mfrow=c(1,2),mar=c(2,2,2,2),oma=c(2,2,2,2),cex.lab=0.8,cex.axis=0.8)
plot(factors,ll[1,1,],type='l',ylim=c(-5E5,1E5),col=cols[1],xlab='',ylab='')
lines(factors,ll[7,1,],type='l',col=cols[2])
lines(factors,ll[26,1,],type='l',col=cols[3])
abline(v=1,lty=2)
mtext(side=2,expression('Log Likelihood'~italic('l('*theta*')')),line=2.5)
mtext(adj=0,'a) Sampling Frequency')
legend('bottomright', legend = c('daily', 'weekly', 'monthly'),bty='n',lty=1,col=cols[1:3])
mtext(side=1,expression(italic(theta)*'/'*theta['true']),line=2.5)

plot(factors,ll[1,1,],type='l',ylim=c(-5E5,1E5),col=cols[1],xlab='',ylab='')
lines(factors,ll[1,12,],type='l',col=cols[2])
lines(factors,ll[1,25,],type='l',col=cols[3])
abline(v=1,lty=2)
mtext(adj=0,'b) Measurement Noise')
legend('bottomright', legend = c('10%', '20%', '50%'),bty='n',lty=1,col=cols[1:3])
mtext(side=1,expression(italic(theta)*'/'*theta['true']),line=2.5)
dev.off()

