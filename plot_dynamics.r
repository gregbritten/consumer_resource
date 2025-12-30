library(viridis)

#unpack from solution saved on disk
burn <- 367 
source('unpack.r')

#PLOT
par(mfcol=c(3,2), mar=c(1,4,1,1), oma=c(3,3,1,1))
##########################################
#plot biomass dynamics
##########################################
cols <- turbo(p$l_species)
plot(time, b[,1], 
     type='l', col=cols[1], 
     ylim=c(0, max(b)*1.1),
     xlab="", ylab="", lwd=2)
#mtext(side=1,"Days",line=2.5)
mtext(side=2,expression("Biomass ["*mu*"M]"),line=2.5)

k <- 2
for(i in 2:p$l_species){
     lines(time, b[,i], col=cols[k], lwd=2)
     k <- k + 1
}
legend("topright", 
       legend = 1:p$l_species,
       col=cols,lty= 1, cex = 0.6,lwd=2,bty='n')

#########################################
#plot resource dynamics
#########################################
cols <- turbo(p$l_resources)
plot(time, r[,1], 
     type = 'l', col = cols[1], 
     ylim = c(0, max(r)*1.1),
     xlab="",ylab="",lwd=2)
#mtext(side=1,"Days",line=2.5)
mtext(side=2,expression("Resource ["*mu*"M]"),line=2.5)

k <- 2
for(i in 2:p$l_resources){
    lines(time, r[,i], col=cols[k], lwd=2)
    k <- k + 1
}
legend("topright", 
       legend=1:p$l_resources,
       col=cols,lty=1, cex=0.6, lwd=2, bty='n')

###########################################
#plot consumption rates
###########################################
cols <- turbo(p$l_resources)
plot(time, c[,1], 
     type='l', col=cols[1], 
     ylim = c(0, max(c)*1.1),
     xlab="", ylab="",lwd=2)
mtext(side=1,"Days", line=2.5)
mtext(side=2,expression("Consumption rate ["*mu*"M/day]"),line=2.5)

k <- 2
for(i in 2:p$l_resources){
    lines(time, c[,i], col=cols[k], lwd=2)
    k <- k + 1
}
legend("topright", 
       legend=1:p$l_resources,
       col=cols, lty=1, cex=0.6, lwd=2, bty='n')


###########################################
#plot system derivatives
###########################################
#plot bdot
cols <- turbo(p$l_species)
plot(time, bdot[,1], 
     type='l', col=cols[1], 
     ylim = c(min(bdot)*1.1, max(bdot)*1.1),
     xlab="", ylab="",lwd=2)
#mtext(side=1,"Days", line=2.5)
mtext(side=2,expression("Biomass derivative ["*mu*"M/day]"),line=2.5)

k <- 2
for(i in 2:(p$l_species)){
    lines(time, bdot[,i], col=cols[k], lwd=2)
    k <- k + 1
}
legend("topright", 
       legend=1:p$l_species,
       col=cols, lty=1, cex=0.6, lwd=2, bty='n')

#plot rdot
cols <- turbo(p$l_resources)
plot(time, rdot[,1], 
     type='l', col=cols[1], 
     ylim = c(min(rdot)*1.1, max(rdot)*1.1),
     xlab="", ylab="",lwd=2)
mtext(side=1,"Days", line=2.5)
mtext(side=2,expression("Derivative ["*mu*"M/day]"),line=2.5)

k <- 2
for(i in 2:(p$l_resources)){
    lines(time, rdot[,i], col=cols[k], lwd=2)
    k <- k + 1
}
legend("topright", 
       legend=1:p$l_resources,
       col=cols, lty=1, cex=0.6, lwd=2, bty='n')

