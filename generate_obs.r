burn <- 367 
file <- "SOL"
source('unpack.r')

#absolute cellular abundance
set.seed(1)
pgN_cell <- runif(p$l_species,0.01,0.2) #randomize mass per cell
Q        <- 1/1E6 * 14.01 * 1E12 * 1/pgN_cell #convert micromoles to cells [mole/umole]*[gN/moleN]*[pg/g]*[cell/pg]
n        <- t(apply(b,1,function(x) Q*x)) #compute relative cell abundance

#total cellular abundance
n_T      <- rowSums(n) 

#relative abundance
n_tilde  <- t(apply(n,1,function(x) x/sum(x)))

#consumption rates (c) already computed

#observed genome
Gfalse_neg <- 0.8
G_obs <- t(apply(p$G,1,function(x){
                ifelse(x==1,rbinom(n=1,size=1,p=Gfalse_neg),0)
}))

#compute likelihood of genome presence-absence
thetas <- seq(0,1,length.out=p$l_species)
plot(thetas,log10(dbinom(sum(p$G - G_obs),size=p$l_resources*p$l_species,prob=thetas)))

#### SAMPLE BGC at time points and regress against genomic context


######################################################
## extract coexistent species 
######################################################
par(mfrow=c(5,4),mar=c(2,2,2,2),oma=c(2,2,2,2))
for(i in 1:ncol(n)){
    plot(n[,i],type='l')
}

n_ce <- n[,colMeans(n)>10^4]

#plot
par(mfrow=c(5,4),mar=c(2,2,2,2),oma=c(2,2,2,2))
for(i in 1:ncol(n_ce)){
    plot(n_ce[,i],type='l')
}

#total abundance
n_T <- rowSums(n_ce)

#####################################
## add noise 
#####################################
noise <- matrix(rnorm(length(c(n_ce)),mean=0,sd=0.05*c(n_ce)),
                byrow=FALSE,
                ncol=ncol(n_ce),
                nrow=nrow(n_ce))
n_ce_noise <- n_ce + noise

#total abundance
n_T_noise <- rowSums(n_ce_noise)

#noise on the total
noise_T <- rnorm(length(n_T_obs),mean=0,sd=0.1*n_T)

n_T_noise_T <- n_T_noise + noise_T
plot(n_T_noise_T, type='l')

e <- n_T_noise_T - n_T

sum(log(dnorm(e,mean=0,sd=0.1*n_T)))



################################
## make plots 
################################
#plot absolute vs. relative abundance
par(mfrow=c(5,4),mar=c(2,2,2,2),oma=c(2,2,2,2))
for(i in 1:ncol(n_tilde)){
    plot(n_tilde[,i],type='l',col='red')
    par(new=TRUE)
    plot(b[,i],type='l',yaxt='n')
    axis(side=4,col='red')
}


