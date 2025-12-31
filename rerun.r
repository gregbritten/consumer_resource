burn <- 367 
file <- "SOL"
source('unpack.r')

#absolute cellular abundance
set.seed(1)
pgN_cell <- runif(p$l_species,0.01,0.2) #randomize mass per cell
Q        <- 1/1E6 * 14.01 * 1E12 * 1/pgN_cell #convert micromoles to cells [mole/umole]*[gN/moleN]*[pg/g]*[cell/pg]
n        <- t(apply(b,1,function(x) Q*x)) #compute relative cell abundance

which_ce <- colMeans(n)>10^4

n_ce <- n[,which_ce]

par(mfrow=c(5,4),mar=c(2,2,2,2),oma=c(2,2,2,2))
for(i in 1:ncol(n_ce)){
    plot(n_ce[,i],type='l')
}


#######################################
#######################################


## re-run with new parameters
r0        <- rep(0.01,l_resources)
b0        <- rep(0.01,l_species)
x0        <- c(b0, r0)
names(x0) <- c(paste0("b", 1:l_species), paste0("r", 1:l_resources))

p_mod1 <- p
p_mod1$Vmax[which(which_ce)[1],1] <- p$Vmax[which(which_ce)[1],1]/5
sol_mod1   <- ode(y=x0, times=times, func=cr_de, parms=p_mod1, method="lsoda")
SOL_mod1   <- list(p=p_mod1,sol=as.data.frame(sol_mod1)) #save solution

p_mod2 <- p
p_mod2$Vmax[which(which_ce)[1],1] <- p$Vmax[which(which_ce)[1],1]*3
sol_mod2   <- ode(y=x0, times=times, func=cr_de, parms=p_mod2, method="lsoda")
SOL_mod2   <- list(p=p_mod2,sol=as.data.frame(sol_mod2)) #save solution


b     <- SOL[['sol']][burn:tf,2:(p$l_species+1)]
b1    <- SOL_mod1[['sol']][burn:tf,2:(p$l_species+1)]
b2    <- SOL_mod2[['sol']][burn:tf,2:(p$l_species+1)]

Q        <- 1/1E6 * 14.01 * 1E12 * 1/pgN_cell #convert micromoles to cells [mole/umole]*[gN/moleN]*[pg/g]*[cell/pg]

n        <- t(apply(b,1,function(x) Q*x)) #compute relative cell abundance
n1       <- t(apply(b1,1,function(x) Q*x)) #compute relative cell abundance
n2       <- t(apply(b2,1,function(x) Q*x)) #compute relative cell abundance

which_ce <- colMeans(n)>10^4

n_ce  <- n[,which_ce]
n_ce1 <- n1[,which_ce]
n_ce2 <- n2[,which_ce]

par(mfrow=c(5,4),mar=c(2,2,2,2),oma=c(2,2,2,2))
for(i in 1:ncol(n_ce)){
    plot(n_ce[,i],type='l')
    lines(n_ce1[,i], col='red')
    lines(n_ce2[,i], col='red')
}

#noise
nsim  <- 2
noise <- list()

for(i in 1:nsim){
    noise[[i]] <- matrix(rnorm(length(c(n_ce)),mean=0,sd=0.05*c(n_ce)),
                    byrow=FALSE,
                    ncol=ncol(n_ce),
                    nrow=nrow(n_ce))
    n_ce_noise[[i]] <- n_ce + noise


par(mfrow=c(2,4))
for(i in 1:ncol(n_ce1)){
    hist(n_ce1[,i] - n_ce[,i],main='')
}

par(mfrow=c(2,4))
for(i in 1:ncol(n_ce1)){
    hist(n_ce2[,i] - n_ce[,i],main='')
}




