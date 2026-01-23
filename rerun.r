library(deSolve)
library(viridis)
library(fields)
source('src/cr_de.r')

burn <- 367 
file <- "SOL"
source('unpack.r')

#absolute cellular abundance
set.seed(1)
pgN_cell <- runif(p$l_species,0.01,0.2) #randomize mass per cell
Q        <- 1/1E6 * 14.01 * 1E12 * 1/pgN_cell #convert micromoles to cells [mole/umole]*[gN/moleN]*[pg/g]*[cell/pg]

#######################################
#######################################

## re-run with new parameters
r0        <- rep(0.01,p$l_resources)
b0        <- rep(0.01,p$l_species)
x0        <- c(b0, r0)
names(x0) <- c(paste0("b", 1:p$l_species), paste0("r", 1:p$l_resources))
tf        <- nrow(b)

sol_ref   <- ode(y=x0, times=p$times, func=cr_de, parms=p, method="lsoda")
b_ref     <- sol_ref[burn:tf,2:(p$l_species+1)]
n_ref     <- t(apply(b_ref,1,function(x) Q*x))
n_T_ref   <- rowSums(n_ref)

ii <- as.numeric(which(colMeans(n_ref)>1E8)[1])

## modify parameters
factors <- seq(0.1,2,length.out)

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
        p_mod$K[ii,j]    <- p$K[ii,j]*factors[i]
    }

    SOL_mod[[i]]   <- ode(y=x0, times=p$times, func=cr_de, parms=p_mod, method="lsoda")
    b_mod[[i]]     <- SOL_mod[[i]][burn:tf,2:(p$l_species+1)]
    n_mod[[i]]     <- t(apply(b_mod[[i]],1,function(x) Q*x))
}


## plot
cols <- turbo(length(factors))
par(mfrow=c(5,4),mar=c(2,2,2,2),oma=c(2,2,2,2))
cols <- turbo(length(factors))
for(i in 1:ncol(n_mod[[1]])){
    plot(time,n_mod[[1]][,i],type='l',col=cols[1])
    for(j in 2:length(factors)){
        lines(time,n_mod[[j]][,i],col=cols[j])
    }
}


##################################
## generate observations 
##################################

e_T = e   <- list()
ss_T = ss = ll1=ll2=ll3=ll4 <- numeric(length(factors))

for(i in 1:length(factors)){
    noise <- matrix(rnorm(length(c(n)),mean=0,sd=0.0001*c(n)),
                byrow=FALSE,
                ncol=ncol(n),
                nrow=nrow(n))
    n_noise   <- apply(n_mod[[i]] + noise, c(1,2), function(x) max(0,x))
    n_noise_T <- rowSums(n_noise) 
    e_T[[i]]  <- n_noise_T - n_T_ref
    ss_T[i]   <- sum(e_T[[i]]^2)

    e[[i]]  <- log10(n_noise) - log10(n_ref) 
    ss[[i]] <- sum(e[[i]]^2)
}


#sds <- c(0.1,0.2,0.5)
#fs  <- c(1,7,30)

sds <- seq(0.1,0.5,length.out=25)
fs  <- seq(1,30,length.out=20)

ll       <- array(NA,dim=c(length(fs),length(sds),length(factors)))
dll=d2ll <- matrix(NA, nrow=length(fs), ncol=length(sds)) 

for(i in 1:length(fs)){
    for(j in 1:length(sds)){
        for(k in 1:length(factors)){
            es <- e[[k]][seq(1,tf,fs[i]),]
            ll[i,j,k] <- sum(dnorm(es,mean=0,sd=sds[j],log=TRUE))
        }
        dll[i,j] <- abs(diff(ll[i,j,]))[19]
        d2ll[i,j] <- abs(diff(diff(ll[i,j,])))[19]
    }
}




pdf('~/dropbox/working/consumer_resource/consumer_resource/plots/expected_CI.pdf',height=4,width=5)
par(mfrow=c(1,1))
filled.contour(x=fs,y=sds*100,z=sqrt(1/d2ll))
mtext(side=1,"Sampling Period [days]",line=2.5,adj=0.1)
mtext(side=2,"Measurement noise [%]",line=2.5)
mtext(adj=0,"c) Expected estimation error")
dev.off()


pdf('~/dropbox/working/consumer_resource/consumer_resource/plots/likelihood_profile.pdf',height=3.75,width=8)
cols <- turbo(10)
par(mfrow=c(1,2),mar=c(2,2,2,2),oma=c(2,2,2,2),cex.lab=0.8,cex.axis=0.8)
plot(factors,ll[1,1,],type='l',ylim=c(-5E5,1E5),col=cols[1],xlab='',ylab='')
lines(factors,ll[3,1,],type='l',col=cols[2])
lines(factors,ll[10,1,],type='l',col=cols[3])
abline(v=1,lty=2)
mtext(side=2,expression('Log Likelihood'~italic('l('*theta*')')),line=2.5)
mtext(adj=0,'a) Sampling Frequency')
legend('bottomright', legend = c('daily', 'weekly', 'monthly'),bty='n',lty=1,col=cols[1:3])
mtext(side=1,expression(italic(theta)*'/'*theta['true']),line=2.5)


plot(factors,ll[1,1,],type='l',ylim=c(-5E5,1E5),col=cols[1],xlab='',ylab='')
lines(factors,ll[1,5,],type='l',col=cols[2])
lines(factors,ll[1,10,],type='l',col=cols[3])
abline(v=1,lty=2)
mtext(adj=0,'b) Measurement Noise')
legend('bottomright', legend = c('10%', '20%', '50%'),bty='n',lty=1,col=cols[1:3])
mtext(side=1,expression(italic(theta)*'/'*theta['true']),line=2.5)
dev.off()

