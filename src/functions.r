#generate genome matrix
genG <- function(p_fn,l_species,l_resources){
    return(
        matrix(replicate(n=l_species,  #generate 'true' genome
            exp=rbinom(n=l_resources,size=1,prob=p_fn)),
            nrow=l_species,
            ncol=l_resources,
            byrow=FALSE
        )
    )
}

#generate Vmax matrix
genVmax <- function(min,max,l_species,l_resources,G){
    Vmax <- matrix(
                runif(min,max,n=l_species*l_resources),
                nrow=l_species,
                ncol=l_resources
            )
    return(Vmax*G)
}

#generate K matrix
genK <- function(min,max,l_species,l_resources){
    return(
        matrix(
            runif(min,max,n=l_species*l_resources),
            nrow=l_species,
            ncol=l_resources)
    )
}

#generate yield matrix
genY <- function(min,max,l_species,l_resources){
    return(
        matrix(
            runif(min,max,n=l_species*l_resources),
            nrow=l_species,
            ncol=l_resources
        )
    )
}

#solve ODEs
run <- function(p,times,by,x0,filename){ 
    sol        <- ode(y=x0, times=times, func=cr_de, parms=p, method="lsoda")  #pass to ODE solver
    SOL        <- list(p=p,sol=as.data.frame(sol)) #save solution as list
 
    if(!dir.exists("results")) dir.create("results") #for ppl running for first time; results folder is not tracked in git
    saveRDS(SOL,paste0("results/",filename,".rds"))
} 
