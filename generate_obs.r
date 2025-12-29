SOL <- readRDS('results/SOL.rds')
sol <- SOL[['sol']]
p   <- SOL[['p']]

## dynamics of relative abundance

B <- sol[,2:(p$n_species+1)]

#absolute cellular abundance
pgN_cell <- runif(n_species,0.01,0.2) #randomize mass per cell
Q        <- 1/1E6 * 14.01 * 1E12 * 1/pgN_cell #convert micromoles to cells [mole/umole]*[gN/moleN]*[pg/g]*[cell/pg]
C        <- t(apply(B,1,function(x) Q*x)) #compute relative cell abundance

#total cellular abundance
C_T      <- rowSums(C) 

#relative abundance
C_tilde  <- t(apply(C,1,function(x) x/sum(x)))

bgc <- sol[,(p$n_resources + p$n_species + 2):ncol(sol)]

#observed genome
Gfalse_neg <- 0.8
G_obs <- t(apply(p$G,1,function(x){
                ifelse(x==1,rbinom(n=1,size=1,p=Gfalse_neg),0)
}))

#compute likelihood of genome presence-absence
thetas <- seq(0,1,length.out=20)
plot(thetas,log10(dbinom(sum(p$G - G_obs),size=p$n_resources*p$n_species,prob=thetas)))


#### SAMPLE BGC at time points and regress against genomic context




#plot absolute vs. relative abundance
par(mfrow=c(5,4),mar=c(2,2,2,2),oma=c(2,2,2,2))
for(i in 1:ncol(C)){
    plot(C_tilde[,i],type='l',col='red')
    par(new=TRUE)
    plot(B[,i],type='l',yaxt='n')
    axis(side=4,col='red')
}


