get_pca <-
function(X,samplelabels,legendlocation="topright",filename=NA,ncomp=5,center=TRUE,scale=TRUE,legendcex=0.5,outloc=getwd(),col_vec=NA,sample.col.opt="default",alphacol=0.3,class_levels=NA,pca.cex.val=3,pca.ellipse=TRUE,ellipse.conf.level=0.5,samplenames=FALSE,do_pca_anova=FALSE){
    
    X<-as.matrix(t(X))
    
    #X<-apply(X,2,scale)
    pch.val<-15
    # print("Performing PCA")
    #print(dim(X))
    #legendlocation="bottomleft"
    #  print(legendlocation)
    samplelabels<-as.data.frame(samplelabels)
 
    samplelabels<-paste("",as.factor(samplelabels[,1]),sep="")
    
    #  print(samplelabels)
    l2<-levels(as.factor(samplelabels))
    col_all=topo.colors(256)
    
    t1<-table(samplelabels)
    if(is.na(class_levels)==TRUE){
        
        l1<-levels(as.factor(samplelabels))
    }else{
        l1<-class_levels
        
        
    }
    
    # print("here")
    
    class_labels_levels<-l1
    
    ncomp=min(dim(X)[1],dim(X)[2])
    
    #   p1<-pcaMethods::pca(t(X),method="svd",center=TRUE,scale="uv",cv="q2",nPcs=10)
    metabpcaresultlog2allmetabs5pcs<-mixOmics::pca(X,ncomp=ncomp,center=TRUE,scale=TRUE)
    
    result<-metabpcaresultlog2allmetabs5pcs
    

    #save(result,file="pcares.Rda")
    
    s1<-summary(result)
    r1<-s1$importance[2,]
    r1<-round(r1,2)*100
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
    # print(paste("Generating PCA plots",sep=""))
    
    fname<-paste("PCA_eval",filename,".tiff",sep="")
    ## 1) raw data
    #tiff(fname,res=300, width=2000,height=2000)
    col <- rep(col_vec[1:length(t1)], t1)
	
    #col<-rep(col_all[1:length(l1)],t1)
    ## Choose different size of points
    cex <- rep(pca.cex.val, dim(X)[1])
    ## Choose the form of the points (square, circle, triangle and diamond-shaped
    
    pch_vec<-seq(1,50) #c(3,5,7,9,12,13,2,17,21) #seq(1,50) #
    pch <- rep(pch.val,dim(X)[1])
    
    
    
    cex <- pca.cex.val #rep(pca.cex.val, dim(X)[1])
    col_per_group<-{}
    pch_per_group<-{}
    for(p1 in 1:length(l2)){
        col[which(samplelabels==l2[p1])]=col_vec[p1]
        pch[which(samplelabels==l2[p1])]=pch_vec[p1]
        
        col_per_group<-c(col_per_group,col_vec[p1])
        pch_per_group<-c(pch_per_group,pch_vec[p1])
    }
   
   # print(table(pch))
   
   main_text=paste("Pairwise PC score plots using ",filename," features after preprocessing",sep="")
   
	   legendcex<-0.9 #0.5*pca.cex.val

   if(pca.ellipse=="car"){
      
      
      col<-col_vec[1:length(t1)]
      cex <- rep(legendcex,length(t1))
      pch <- unique(pch) #rep(15,length(t1))
      
      
      print(dataEllipse(x=result$x[,1], y=result$x[,2],groups=as.factor(samplelabels),grid=TRUE,lwd=4,levels=c(ellipse.conf.level),col=col,pch=pch,xlab=paste("PC1 (",r1[1],"% variation)",sep=""),ylab=paste("PC2 (",r1[2],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE))
       print(legend(legendlocation, l1, col = col,pch = pch, pt.cex = cex, title = "Class #", cex=legendcex))
       
       
      print(dataEllipse(x=result$x[,1], y=result$x[,3],groups=as.factor(samplelabels),grid=TRUE,lwd=4,levels=c(ellipse.conf.level),col=col,pch=pch,xlab=paste("PC1 (",r1[1],"% variation)",sep=""),ylab=paste("PC3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE))
       print(legend(legendlocation, l1, col = col,pch = pch, pt.cex = cex, title = "Class #", cex=legendcex))
       
       
       print(dataEllipse(x=result$x[,2], y=result$x[,3],groups=as.factor(samplelabels),grid=TRUE,lwd=4,levels=c(ellipse.conf.level),col=col,pch=pch,xlab=paste("PC2 (",r1[2],"% variation)",sep=""),ylab=paste("PC3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE))
       print(legend(legendlocation, l1, col = col,pch = pch, pt.cex = cex, title = "Class #", cex=legendcex))
       
   }else{
      
      # print(samplelabels)
      
      if(do_pca_anova==TRUE){
      scores_res<-result$x
      
      #  print(dim(scores_res))
      # print(length(samplelabels))
      # print(dim(samplelabels))
   
   pc1_pval<-anova(lm(cbind(scores_res[,1],scores_res[,2])~samplelabels))
   
   pc1_pval<-pc1_pval[[6]][2]
   pc2_pval<-anova(lm(cbind(scores_res[,1],scores_res[,3])~samplelabels))
   
   pc2_pval<-pc2_pval[[6]][2]
   
   pc3_pval<-anova(lm(cbind(scores_res[,2],scores_res[,3])~samplelabels))
   pc3_pval<-pc3_pval[[6]][2]
   
   pc_pval_vec<-c(pc1_pval,pc2_pval,pc3_pval)
      }

if(do_pca_anova==TRUE){
main_text=paste("Pairwise PC score plots using ",filename," features after preprocessing\np-value for overall differences between groups using PC1 and PC2 in a multivariate\n one-way ANOVA model=",round(pc1_pval,3),sep="")
}

print(plotIndiv(result, comp = c(1,2),ind.names = FALSE,col = col, cex = cex, pch = c(2), X.label=paste("PC1 (",r1[1],"% variation)",sep=""),Y.label=paste("PC2 (",r1[2],"% variation)",sep=""),title=main_text,add.legend=TRUE,plot.ellipse=pca.ellipse,ellipse.level=ellipse.conf.level,col.per.group=col_per_group,group=samplelabels))
    
    
    # scores_res<-result$x
    
    #scores_res<-as.data.frame(scores_res)
    
    #print("ggplot PCA")
    #print(head(scores_res))


#cnames_pca_res<-c("PC1","PC2")
#   colnames(scores_res)<-cnames_pca_res
    
    
    #   save(scores_res,file="scores_res.Rda")
    #res<-ggplot(scores_res, aes(PC1, PC2,colour=samplelabels)) +
    #geom_point() +
    #stat_ellipse()
    
    #print(res)

 
 if(do_pca_anova==TRUE){
 main_text=paste("Pairwise PC score plots using ",filename," features after preprocessing\np-value for overall differences between groups using PC1 and PC3 in a multivariate\n one-way ANOVA model=",round(pc2_pval,3),sep="")
 }
 
       print(plotIndiv(result, comp = c(1,3),ind.names = FALSE,col = col, cex = cex, pch = c(2), X.label=paste("PC1 (",r1[1],"% variation)",sep=""),Y.label=paste("PC3 (",r1[3],"% variation)",sep=""),title=main_text,add.legend=TRUE,plot.ellipse=pca.ellipse,ellipse.level=ellipse.conf.level,col.per.group=col_per_group,group=samplelabels))
 if(do_pca_anova==TRUE){
 
    main_text=paste("Pairwise PC score plots using ",filename," features after preprocessing\np-value for overall differences between groups using PC2 and PC3 in a multivariate\n one-way ANOVA model=",round(pc3_pval,3),sep="")
 }
    print(plotIndiv(result, comp = c(2,3),ind.names = FALSE,col = col, cex = cex, pch = c(2), X.label=paste("PC3 (",r1[3],"% variation)",sep=""),Y.label=paste("PC2 (",r1[2],"% variation)",sep=""),title=main_text,add.legend=TRUE,plot.ellipse=pca.ellipse,ellipse.level=ellipse.conf.level,col.per.group=col_per_group,group=samplelabels))
    
   
    }
    
    
    
    #print(legend(legendlocation, l1, col = col,pch = pch, pt.cex = cex, title = "Class #", cex=legendcex))
    #dataEllipse(x=pcIr@scores[,1], y=pcIr@scores[,2],groups=as.factor(iris[,5]),grid=TRUE,lwd=4,levels=c(0.95),col=c("red","blue","green"),pch=c(1,7,5))
    
    
    # print("done with PCA")
    return(result)
}
