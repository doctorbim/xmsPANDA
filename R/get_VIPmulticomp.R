get_VIPmulticomp <-
function(p2,nvar,h) {
    
    
    b<-c(p2$loadings$Y)[1:h]
    T<-p2$variates$X[,1:h]
    SS<-b^2 * colSums(T^2)
    
    W<-p2$loadings$X[,1:h]
    # print(dim(W))
    W<-as.data.frame(W)
    nvar=nrow(W)
    h=ncol(W)
    
    Wnorm2 <- colSums(W^2)
    pls_vec<-lapply(1:nvar,function(j){
        # return(sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS)))
        
        return(sqrt(nrow(W) * sum( (W[j,]^2 / Wnorm2) * (SS/ sum(SS)))))
    })
    pls_vec<-unlist(pls_vec)
    return(pls_vec)
}
