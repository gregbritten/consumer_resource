library(viridis)

#unpack from solution saved on disk
burn <- 367 
file <- "SOL2" #choose from "SOL" (additive) or "SOL_linear" (linear)
source('unpack.r')

pdf('plots/example_simulations.pdf',height=5,width=10)
#PLOT
par(mfrow=c(2,2), mar=c(1,4,1,1), oma=c(3,3,1,1))
##########################################
#plot biomass dynamics
##########################################
cols <- turbo(p$l_species)
plot(time/365, b[,1], 
     type='l', col=cols[1], 
     ylim=c(0, max(b)*1.1),
     xlab="", ylab="", lwd=2,xlim=c(1,11))
#mtext(side=1,"Days",line=2.5)
mtext(side=2,expression("Biomass ["*mu*"M]"),line=2.5)

k <- 2
for(i in 2:p$l_species){
     lines(time/365, b[,i], col=cols[k], lwd=2)
     k <- k + 1
}
#legend("topright", 
#       legend = 1:p$l_species,
#       col=cols,lty= 1, cex = 0.6,lwd=2,bty='n')


pgN_cell <- runif(p$l_species,0.01,0.2) #randomize mass per cell
Q        <- 1/1E6 * 14.01 * 1E12 * 1/pgN_cell #convert micromoles to cells [mole/umole]*[gN/moleN]*[pg/g]*[cell/pg]
n        <- t(apply(b,1,function(x) Q*x)) #compute relative cell abundance
n_tilde  <- t(apply(n,1,function(x) x/sum(x)))

#plot bdot
cols <- turbo(p$l_species)
plot(time/365, n_tilde[,1], 
     type='l', col=cols[1], 
     ylim = c(0,0.6),
     xlab="", ylab="",lwd=2)
#mtext(side=1,"Days", line=2.5)
mtext(side=2,expression("True Relative Abundance"),line=2.5)

k <- 2
for(i in 2:(p$l_species)){
    lines(time/365, n_tilde[,i], col=cols[k], lwd=2)
    k <- k + 1
}
#legend("topright", 
#       legend=1:p$l_species,
#       col=cols, lty=1, cex=0.6, lwd=2, bty='n'



#########################################
#plot resource dynamics
#########################################
cols <- turbo(p$l_resources)
plot(time/365, r[,1], 
     type = 'l', col = cols[1], 
     ylim = c(0, max(r)*1.1),
     xlab="",ylab="",lwd=2,xlim=c(1,11))
#mtext(side=1,"Days",line=2.5)
mtext(side=2,expression("Resource ["*mu*"M]"),line=2.5)

k <- 2
for(i in 2:p$l_resources){
    lines(time/365, r[,i], col=cols[k], lwd=2)
    k <- k + 1
}
#legend("topright", 
#       legend=1:p$l_resources,
#       col=cols,lty=1, cex=0.6, lwd=2, bty='n')
mtext(side=1,"Years", line=2.5)

###########################################
#plot consumption rates
###########################################
cols <- turbo(p$l_resources)
plot(time/365, c[,1], 
     type='l', col=cols[1], 
     ylim = c(0, max(c)*1.1),
     xlab="", ylab="",lwd=2,xlim=c(1,11))
mtext(side=1,"Years", line=2.5)
mtext(side=2,expression("Consumption rate ["*mu*"M/day]"),line=2.5)

k <- 2
for(i in 2:p$l_resources){
    lines(time/365, c[,i], col=cols[k], lwd=2)
    k <- k + 1
}
#legend("topright", 
#       legend=1:p$l_resources,
#       col=cols, lty=1, cex=0.6, lwd=2, bty='n')


dev.off()

