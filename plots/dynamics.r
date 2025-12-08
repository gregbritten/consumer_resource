sol_df <- readRDS('results/sol_df.rds')

#plot biomass dynamics
plot(sol_df$time, sol_df$B1, 
     type = 'l', col = 'blue', 
     ylim = c(0, max(sol_df[, 2:(num_species+1)])*1.1),
     xlab = "", ylab = "", 
     lwd=2)
mtext(side=1,"Days",line=2.5)
mtext(side=2,"Biomass",line=2.5)
lines(sol_df$time, sol_df$B2, col='red', lwd=2)
lines(sol_df$time, sol_df$B3, col='green', lwd=2)

legend("bottomright", 
       legend = c("Species 1", "Species 2", "Species 3"),
       col = c("blue", "red", "green"), 
       lty = 1, cex = 0.8, lwd=2, bty='n')

