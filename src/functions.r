genG <- function(p_fn,l_species,l_resources){
    return(
        matrix(replicate(n=l_species,  #generate 'true' genome
            exp=rbinom(n=l_resources,size=1,prob=p_fg)),
            nrow=l_species,
            ncol=l_resources,
            byrow=FALSE
        )
    )
}

genVmax <- function(min,max,l_species,l_resources,G){
    Vmax <- matrix(
        runif(min,max,n=l_species*l_resources),
        nrow=l_species,
        ncol=l_resources)
    return(Vmax*G)
}
