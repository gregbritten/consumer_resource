library(viridis)
library(matrixcalc)
library(glmnet)

#unpack from solution saved on disk
burn <- 367 
file <- "SOL_linear"
source('unpack.r')

#stack derivative vector
Xdot <- as.vector(vec(as.matrix(xdot)))

#combine state vector
X <- cbind(b,r)

#construct design matrix by multiplying columns and column binding
for(i in 1:p$l_species){
    for(j in 1:p$l_resources){
        XX           <- as.matrix(b[,i]*r[,j])
        colnames(XX) <- paste0("b",i,"r",j)
        X <- cbind(X,XX)
    }
}

#stack state vector that includes interactions
X <- as.matrix(X[rep(1:3650,30),])

#fit linear regression
out <- lm(Xdot ~ X)

#fit lasso regression
out_las <- glmnet(x=X, y=Xdot)




