get_plsplots <-
function(X,plsres,plsvar,samplelabels,filename=NA,ncomp=5,center=TRUE,scale=TRUE,legendcex=0.5,outloc=getwd(),col_vec=NA,sample.col.opt="default",alphacol=0.3,legendlocation="topright",class_levels=NA,pca.cex.val=3){
   
    result<-plsres
    r1<-plsvar
   
    pch.val=15
    legendlocation="bottomleft"
    
    samplelabels<-as.data.frame(samplelabels)
    samplelabels<-as.factor(samplelabels[,1])
    l2<-levels(as.factor(samplelabels))
    col_all=topo.colors(256)
    
    t1<-table(samplelabels)
    if(is.na(class_levels)==TRUE){
       
        l1<-levels(as.factor(samplelabels))
    }else{
        l1<-class_levels

        
    }

    class_labels_levels<-l1
    
    ncomp=min(dim(X)[1],dim(X)[2])
  
    if(is.na(col_vec)==TRUE){
        col_vec<-c("mediumpurple4","mediumpurple1","blueviolet","darkblue","blue","cornflowerblue","cyan4","skyblue",
        "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
        "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
        "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
        
    }
    
    if(sample.col.opt=="default"){
        
        col_vec<-c("#CC0000","#AAC000","blue","mediumpurple4","mediumpurple1","blueviolet","cornflowerblue","cyan4","skyblue",
        "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
        "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
        "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
        
    }else{
        if(sample.col.opt=="topo"){
            #col_vec<-topo.colors(256) #length(class_labels_levels))
            
            #col_vec<-col_vec[seq(1,length(col_vec),)]
            
            col_vec <- topo.colors(length(class_labels_levels), alpha=alphacol)
        }else{
            if(sample.col.opt=="heat"){
                #col_vec<-heat.colors(256) #length(class_labels_levels))
                
                col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
            }else{
                if(sample.col.opt=="rainbow"){
                    #col_vec<-heat.colors(256) #length(class_labels_levels))
                    col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
                    
                    #col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
                }else{
                    
                    if(sample.col.opt=="terrain"){
                        #col_vec<-heat.colors(256) #length(class_labels_levels))
                        #col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
                        
                        col_vec <- cm.colors(length(class_labels_levels), alpha=alphacol)
                    }
                    
                    
                }
                
            }
            
        }
    }
    #col_vec<-col_vec[sample(1:length(col_vec),length(col_vec))]
    
    l1<-gsub(x=l1,pattern="Class",replacement="")
    
    dir.create(outloc,showWarnings=FALSE)
    setwd(outloc)
  
    ## 1) raw data
    #tiff(fname,res=300, width=2000,height=2000)
    col <- rep(col_vec[1:length(t1)], t1)
    #col<-rep(col_all[1:length(l1)],t1)
    ## Choose different size of points
    cex <- rep(pca.cex.val, dim(X)[1])
    ## Choose the form of the points (square, circle, triangle and diamond-shaped
    
    
    
    pch_vec<-seq(1,50) #c(3,5,7,9,12,13,2,17,21) #seq(1,50) #
    pch <- rep(pch.val,dim(X)[1])
    cex <- rep(pca.cex.val, dim(X)[1])
    for(p1 in 1:length(l2)){
        
        pch[which(samplelabels==l2[p1])]=pch_vec[p1]
    }
 

col <- rep(col_vec[1:length(t1)], t1)
pch_vec<-seq(1,50) #c(3,5,7,9,12,13,2,17,21) #seq(1,50) #
pch <- rep(pch.val,dim(X)[1])
cex <- rep(pca.cex.val, dim(X)[1])
for(p1 in 1:length(l2)){
    
    pch[which(samplelabels==l2[p1])]=pch_vec[p1]
}


    print(plotIndiv(result, comp = c(1,2),ind.names = FALSE, group=samplelabels, cex = cex[1], pch = pch, ellipse=FALSE, ellipse.level = 0.95, X.label=paste("PLS1 (",r1[1],"% variation)",sep=""),Y.label=paste("PLS2 (",r1[2],"% variation)",sep=""),add.legend=TRUE))
    

    if(result$ncomp>2){
        
        print(plotIndiv(result, comp = c(1,3),ind.names = FALSE, group=samplelabels, cex = cex[1], pch = pch, ellipse=FALSE, ellipse.level = 0.95, X.label=paste("PLS1 (",r1[1],"% variation)",sep=""),Y.label=paste("PLS3 (",r1[3],"% variation)",sep=""),add.legend=TRUE))


        print(plotIndiv(result, comp = c(2,3),ind.names = FALSE, group=samplelabels, cex = cex[1], pch = pch, ellipse=FALSE, ellipse.level = 0.95, X.label=paste("PLS2 (",r1[2],"% variation)",sep=""),Y.label=paste("PLS3 (",r1[3],"% variation)",sep=""),add.legend=TRUE))
        
        
    }


}
