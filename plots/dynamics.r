library(viridis)

SOL <- readRDS('results/SOL.rds')
sol <- SOL[['sol']]
p   <- SOL[['p']]

cols_R <- turbo(p$n_resources)
cols_B <- turbo(p$n_species)

#plot biomass dynamics
par(mfrow=c(2,1),mar=c(1,3,1,1),oma=c(3,1,1,1))
plot(sol$time, sol$B1, 
     type = 'l', col = 'blue', 
     ylim = c(0, max(sol[, 2:(n_species+1)])*1.1),
     xlab = "", ylab = "", 
     lwd=2)
#mtext(side=1,"Days",line=2.5)
mtext(side=2,"Biomass",line=2.5)

for(i in 3:(2+p$n_species))
lines(sol$time, sol[,i], col=i, lwd=2)

legend("topright", 
       legend = 1:n_species,
       col = 1:n_species, 
       lty = 1, cex = 0.6, lwd=2, bty='n')


#plot resource dynamics
plot(sol$time, log10(sol$R1), 
     type = 'l', col = 'blue', 
     ylim = c(min(log10(sol[, (2+p$n_species):ncol(sol)])),
              max(log10(sol[, (2+p$n_species):ncol(sol)]))),
     xlab = "", ylab = "", 
     lwd=2)
mtext(side=1,"Days",line=2.5)
mtext(side=2,"Substrate",line=2.5)

for(i in (3+p$n_species):(1+p$n_species+p$n_resources)){
    lines(sol$time, log10(sol[,i]), col=i , lwd=2)
}

legend("topright", 
       legend = 1:n_resources,
       col = 1:n_resources, 
       lty = 1, cex = 0.6, lwd=2, bty='n')
