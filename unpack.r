#unpack from solution saved on disk
SOL  <- readRDS(paste0('results/',file,'.rds'))
p    <- SOL[['p']]
tf   <- length(p$times)
time <- p$times[burn:tf]
b    <- SOL[['sol']][burn:tf,2:(p$l_species+1)]
r    <- SOL[['sol']][burn:tf,(p$l_species+2):(p$l_species+p$l_resources+1)]
c    <- SOL[['sol']][burn:tf,(p$l_species+p$l_resources+2):(p$l_species+p$l_resources+1+p$l_resources)]
xdot <- SOL[['sol']][burn:tf,(p$l_species+p$l_resources+2+p$l_resources):(p$l_species+p$l_resources+1+2*p$l_resources+p$l_species)]
bdot <- xdot[,1:p$l_species]
rdot <- xdot[,(p$l_species+1):ncol(xdot)]
