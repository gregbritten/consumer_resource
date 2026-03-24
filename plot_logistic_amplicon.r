library(viridis)

logn_lim <- log10(1E9)
n50s  <- log10(c(1E3,1E5)) #abundance for 50% detection probability 
n     <- seq(0,logn_lim,0.1)
betas <- c(2,4,5)


pdf('plots/logistic.pdf',height=4,width=5)
plot(-999,xlim=c(0,logn_lim),ylim=c(0,1),xaxt='n',ylab='',xlab='')

cols <- turbo(5)
for(j in 1:length(n50s)){
    for(i in 1:length(betas)){
        p <- 1/(1 + exp(-betas[i]*(n-n50s[j])))
        lines(n,p,col=cols[i+1],lwd=2)
    }
}
abline(v=n50s,lty=2)

axis_labels <- parse(text = paste("10^", 0:logn_lim, sep = ""))
axis(side=1,at=seq(0,logn_lim),labels=axis_labels)
mtext(side=1,'True Cellular Abundance [cells/ml]',line=2.5)
mtext(side=2,'Probability of Detection',line=2.5)
legend('bottomright', bty='n',
       legend=c(expression('n'[50]),expression(beta~'='~2),expression(beta~'='~4),expression(beta~'='~5)),
       lty=c(2,1,1,1),
       col=c('black',cols[2:5]))
dev.off()




