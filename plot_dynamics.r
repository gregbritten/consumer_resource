library(viridis)

SOL <- readRDS('results/SOL.rds')
sol <- SOL[['sol']]
p   <- SOL[['p']]

#drop first year as spin-up time
sol <- sol[366:nrow(sol),]

#PLOT
par(mfrow=c(3,1), mar=c(1,3,1,1), oma=c(3,3,1,1))
##########################################
#plot biomass dynamics
##########################################
cols <- turbo(p$n_species)
plot(sol$time, sol$B1, 
     type='l', col=cols[1], 
     ylim=c(0, max(sol[,2:(p$n_species+1)])*1.1),
     xlab="", ylab="", lwd=2)
#mtext(side=1,"Days",line=2.5)
mtext(side=2,expression("Biomass ["*mu*"M]"),line=2.5)

k <- 2
for(i in 3:(1+p$n_species)){
     lines(sol$time, sol[,i], col=cols[k], lwd=2)
     k <- k + 1
}
legend("topright", 
       legend = 1:n_species,
       col=cols,lty= 1, cex = 0.6,lwd=2,bty='n')

#########################################
#plot resource dynamics
#########################################
cols <- turbo(p$n_resources)
plot(sol$time, sol$R1, 
     type = 'l', col = cols[1], 
     ylim = c(0, max(sol[, (2+p$n_species):(1+p$n_species+p$n_resources)])*1.1),
#              max(log10(sol[, (2+p$n_species):ncol(sol)]))),
     xlab="",ylab="",lwd=2)
#mtext(side=1,"Days",line=2.5)
mtext(side=2,expression("Resource ["*mu*"M]"),line=2.5)

k <- 2
for(i in (3+p$n_species):(1+p$n_species+p$n_resources)){
    lines(sol$time, sol[,i], col=cols[k], lwd=2)
    k <- k + 1
}
legend("topright", 
       legend=1:n_resources,
       col=cols,lty=1, cex=0.6, lwd=2, bty='n')

###########################################
#plot consumption rates
###########################################
cols <- turbo(p$n_resources)
plot(sol$time, sol$consumption_j1, 
     type='l', col=cols[1], 
     ylim = c(0, max(sol[,(1+p$n_species+p$n_resources):ncol(sol)])*1.1),
#              max(log10(sol[, (2+p$n_species):ncol(sol)]))),
     xlab="", ylab="",lwd=2)
mtext(side=1,"Days", line=2.5)
mtext(side=2,expression("Consumption rate ["*mu*"M/day]"),line=2.5)

k <- 2
for(i in (2+p$n_species+p$n_resources):ncol(sol)){
    lines(sol$time, sol[,i], col=cols[k], lwd=2)
    k <- k + 1
}
legend("topright", 
       legend=1:n_resources,
       col=cols, lty=1, cex=0.6, lwd=2, bty='n')




