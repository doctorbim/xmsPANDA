get_hca_child <-
function(feature_table_file,parentoutput_dir,class_labels_file,X=NA,Y=NA,heatmap.col.opt="RdBu",cor.method="spearman",is.data.znorm=FALSE,analysismode="classification",
sample.col.opt="rainbow",plots.width=8,plots.height=8,plots.res=600, plots.type="cairo", alphacol=0.3, hca_type,newdevice=FALSE,input.type="intensity",mainlab="",cexRow=1, cexCol=1)
{
    
    #print(dim(X))
    
    if(is.na(X)==TRUE){
    data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE)
    }else{
        data_matrix<-X
        rm(X)
        
    }
    
    dir.create(parentoutput_dir)
    setwd(parentoutput_dir)
    
    
    data_m<-data_matrix[,-c(1:2)]
    
    data_m<-as.matrix(data_m)
    
    col_samples<-TRUE
    
    if(is.na(Y)==TRUE){
    if(is.na(class_labels_file)==FALSE){
        classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
    }else{
        
        classlabels<-rep("classA",dim(data_m)[2])
        classlabels<-cbind(classlabels,classlabels)
        
        col_samples<-FALSE
        
    }
    }else{
        classlabels<-Y
        
    }
    
    
    #patientcolors<-rep("green",dim(data_m)[2])
    
    class_labels_levels<-levels(as.factor(classlabels[,2]))
    ordered_labels<-classlabels[,2]
    
    #class_label_alphabets<-c("A","B","C","D","E","F","G","H","I","J","K","L","M")
    class_label_alphabets<-class_labels_levels #paste("C",1:length(class_labels_levels),sep="")
    
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
                    }else{
                        col_vec <-sample.col.opt
                        
                    }
                    
                    
                }
                
            }
            
        }
    }
    col_samples=FALSE
    if(analysismode=="classification")
    {
        
        sampleclass<-{}
        patientcolors<-{}
        #print(classlabels)
        classlabels<-as.data.frame(classlabels)
        f<-factor(classlabels[,1])
        
        col_samples=TRUE
        for(c in 1:length(class_labels_levels)){
            
            num_samps_group_cur=length(which(ordered_labels==class_labels_levels[c]))
            
            #classlabels<-c(classlabels,rep(paste("Class",class_label_alphabets,sep=""),num_samps_group_cur))
            #,rep("ClassB",num_samps_group[[2]]),rep("ClassC",num_samps_group[[3]]))
            sampleclass<-c(sampleclass,rep(paste("Class",class_label_alphabets[c],sep=""),num_samps_group_cur))
            
            patientcolors <-c(patientcolors,rep(col_vec[c],num_samps_group_cur))
								}
        
        
        
    }
    
    if(heatmap.col.opt=="RdBu"){
        
        heatmap.col.opt="redblue"
    }
    
    heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
    heatmap_cols<-rev(heatmap_cols)
    
    if(heatmap.col.opt=="topo"){
        heatmap_cols<-topo.colors(256)
         heatmap_cols<-rev(heatmap_cols)
    }else{
        if(heatmap.col.opt=="heat"){
            heatmap_cols<-heat.colors(256)
             heatmap_cols<-rev(heatmap_cols)
        }else{
            
            if(heatmap.col.opt=="yellowblue"){
                
                  heatmap_cols<-colorRampPalette(c("yellow","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                #heatmap_cols<-blue2yellow(256) #colorRampPalette(c("yellow","blue"))(256)
               heatmap_cols<-rev(heatmap_cols)
            }else{
                
                if(heatmap.col.opt=="redblue"){
                    
                    heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
                    heatmap_cols<-rev(heatmap_cols)
                }else{
                    
                    #my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                    if(heatmap.col.opt=="redyellowgreen"){
                        
                        heatmap_cols <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                        heatmap_cols<-rev(heatmap_cols)
                    }else{
                        if(heatmap.col.opt=="yellowwhiteblue"){
                            
                            heatmap_cols<-colorRampPalette(c("yellow2","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                            heatmap_cols<-rev(heatmap_cols)
                        }else{
                            
                            if(heatmap.col.opt=="redwhiteblue"){
                                
                                heatmap_cols<-colorRampPalette(c("red","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                                heatmap_cols<-rev(heatmap_cols)
                            }else{
                                
                                
                                
                                    heatmap_cols <- colorRampPalette(brewer.pal(10, heatmap.col.opt))(256)
                                    heatmap_cols<-rev(heatmap_cols)
                                
                            }

                        }
                        
                    }
                    
                }
                
            }
        }
        
    }
    
    # print(classlabels)
    #print(patientcolors)
    
    if(input.type=="intensity"){
        #hc <- try(hclust(as.dist(1-WGCNA::cor(data_m,method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #samples
        hc <- try(hclust(as.dist(1-WGCNA::cor(data_m,method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #samples
        hr <- try(hclust(as.dist(1-WGCNA::cor(t(data_m),method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #metabolites
        
    }else{
        if(input.type=="correlation"){
            #hc <- try(hclust(as.dist(1-data_m)),silent=TRUE) #samples
            
             hc <- try(hclust(as.dist(1-data_m)),silent=TRUE) #samples
             hr <- try(hclust(as.dist(1-data_m)),silent=TRUE) #metabolites
             
            
        }
        
    }
    
    
   
   mainlab1<-paste("HCA using ",mainlab," significant features",sep="")
    
    if(is(hr,"try-error") || is(hc,"try-error")){
								
                                print("Hierarchical clustering can not be performed. ")
    }else{
        heatmap_file<-paste("heatmap.jpeg",sep="")
        
        if(newdevice==TRUE){
        png(heatmap_file,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        #png("test_png.png",res=600,width=8,height=8,units="in",type="cairo")
        }
        
        if(hca_type=="two-way"){
            
            if(input.type=="intensity"){
            hc <- try(hclust(as.dist(1-WGCNA::cor(data_m,method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #samples
            }else{
                if(input.type=="correlation"){
                 hc <- try(hclust(as.dist(1-data_m)),silent=TRUE) #samples
                 
                }
                
            }
            
            if(col_samples==FALSE){
                if(is.data.znorm==FALSE){
                    
                    h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1)
                    
                    
                }else{
                    h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1)
                }
                
            }else{

            
            if(is.data.znorm==FALSE){
                
                h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol, xlab="Samples",ylab="mzfeatures", main=mainlab1, ColSideColors=patientcolors)
            }else{
                h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1, ColSideColors=patientcolors)
            }
            }
   

            if(newdevice==TRUE){
            try(dev.off(),silent=TRUE)
            }
            
            
            mycl_samples <- cutree(hc, h=max(hc$height)/2)
            mycl_metabs <- cutree(hr, h=max(hr$height)/2)
            
            ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],data_matrix[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])
        }
        else{
            if(hca_type=="one-way"){
            hc<-seq(1,dim(data_m)[2])
            
            if(col_samples==FALSE){
                if(is.data.znorm==FALSE){
                    
                    h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1)
                }else{
                    h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1)
                }
                
            }else{
            if(is.data.znorm==FALSE){
                
                h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,dendrogram = c("row"),ColSideColors=patientcolors)
            }else{
                h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,dendrogram = c("row"),ColSideColors=patientcolors)
            }
            }
           
            if(newdevice==TRUE){
                try(dev.off(),silent=TRUE)
            }
            
            mycl_samples<-seq(1,dim(data_m)[2])
            #mycl_samples <- cutree(hc, h=max(hc$height)/2)
            mycl_metabs <- cutree(hr, h=max(hr$height)/2)
            
            ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],data_matrix[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])
            
            }else{
                if(hca_type=="samples"){
                hr<-seq(1,dim(data_m)[1])
                if(input.type=="intensity"){
                    hc <- try(hclust(as.dist(1-WGCNA::cor(data_m,method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #samples
                }else{
                    if(input.type=="correlation"){
                        hc <- try(hclust(as.dist(1-data_m)),silent=TRUE) #samples
                        
                    }
                    
                }
                
            
                if(col_samples==FALSE){
                    if(is.data.znorm==FALSE){
                        
                        h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1)
                    }else{
                        h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol, xlab="Samples",ylab="mzfeatures", main=mainlab1)
                    }
                    
                }else{
                    if(is.data.znorm==FALSE){
                        
                        h73<-heatmap.2(data_m, Rowv=NULL, Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,dendrogram = c("col"),ColSideColors=patientcolors)
                    }else{
                        h73<-heatmap.2(data_m, Rowv=NULL, Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,dendrogram = c("col"),ColSideColors=patientcolors)
                    }
                }
                
                # par(xpd=TRUE)
                #legend("bottomleft",legend=levels(ordered_labels),text.col=unique(patientcolors),pch=13,cex=0.4)
                #par(xpd=FALSE)
                
                if(newdevice==TRUE){
                   try(dev.off(),silent=TRUE)
                }
                
                #mycl_samples<-seq(1,dim(data_m)[2])
                mycl_samples <- cutree(hc, h=max(hc$height)/2)
                mycl_metabs <- seq(1,dim(data_m)[1]) #cutree(hr, h=max(hr$height)/2)
                
                ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],data_matrix[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])
                
                }else{
                    
                    hr<-seq(1,dim(data_m)[1])
                    if(input.type=="intensity"){
                        hc <- try(hclust(as.dist(1-WGCNA::cor(data_m,method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #samples
                    }else{
                        if(input.type=="correlation"){
                            hc <- try(hclust(as.dist(1-data_m)),silent=TRUE) #samples
                            
                        }
                        
                    }
                    
                    
                    if(col_samples==FALSE){
                        if(is.data.znorm==FALSE){
                            
                            h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,dendrogram = c("none"))
                        }else{
                            h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,dendrogram = c("none"))
                        }
                        
                    }else{
                        if(is.data.znorm==FALSE){
                            
                            #print("here")
                            h73<-heatmap.2(data_m, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,ColSideColors=patientcolors,dendrogram = c("none"))
                        }else{
                            h73<-heatmap.2(data_m, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,ColSideColors=patientcolors,dendrogram = c("none"))
                        }
                    }
                    
                    #  par(xpd=TRUE)
                    #legend("bottomleft",legend=levels(ordered_labels),text.col=unique(patientcolors),pch=13,cex=0.4)
                    #par(xpd=FALSE)
                    
                    if(newdevice==TRUE){
                        try(dev.off(),silent=TRUE)
                    }
                    
                    mycl_samples<-seq(1,dim(data_m)[2])
                    #mycl_samples <- cutree(hc, h=max(hc$height)/2)
                    mycl_metabs <- seq(1,dim(data_m)[1]) #cutree(hr, h=max(hr$height)/2)
                    
                    ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],data_matrix[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])
                    
                    
                }
            }
        }
        
        cnames1<-colnames(ord_data)
        cnames1[1]<-"mz_cluster_label"
        colnames(ord_data)<-cnames1
        fname1<-paste("Clustering_based_sorted_intensity_using_",mainlab,"features.txt",sep="")
        write.table(ord_data,file=fname1,sep="\t",row.names=FALSE)
        
        fname2<-paste("Sample_clusterlabels_using_",mainlab,"features.txt",sep="")
        
        sample_clust_num<-mycl_samples[h73$colInd]
        temp1<-classlabels[h73$colInd,1]
        temp2<-classlabels
        
        #print(head(temp2))
        #temp1<-as.data.frame(temp1)
        
        #print(dim(temp1))
        match_ind<-match(temp1,temp2[,1])
        
        temp3<-temp2[match_ind,]
        
        #print(head(temp3))
        temp4<-cbind(temp1,temp3,sample_clust_num)
       
       #write.table(temp4,file="s1.txt",sep="\t",row.names=FALSE)
    #   print(head(temp1))
        print(head(temp4))
        
        rnames1<-rownames(temp4)
        #temp4<-cbind(rnames1,temp4)
        temp4<-as.data.frame(temp4)
        temp4<-temp4[,-c(1)]
        # print(temp4[,1:4])
        
        
        
        
        if(analysismode=="regression"){
            

            
            temp3<-temp4 #[,-c(1)]
            temp3<-as.data.frame(temp3)
            temp3<-apply(temp3,2,as.numeric)
            
            
            temp_vec<-as.vector(temp3[,2])
            
            
            
            names(temp_vec)<-as.character(temp4[,1])
            
            
            # if(output.device.type!="pdf"){
                
             if(newdevice==TRUE){
                
                temp_filename_1<-"Figures/Barplot_dependent_variable_ordered_by_HCA.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
            }
            
            
            #print(temp_vec)
            #tiff("Barplot_sample_cluster_ymat.tiff", width=plots.width,height=plots.height,res=plots.res, compression="lzw")
            barplot(temp_vec,col="brown",ylab="Y",cex.axis=0.5,cex.names=0.5,main="Dependent variable levels in samples; \n ordered based on hierarchical clustering")
            #dev.off()
            
            
           
            if(newdevice==TRUE){
                try(dev.off(),silent=TRUE)
            }
            
            
            
            
            
        }
        
        write.table(temp4,file=fname2,sep="\t",row.names=FALSE)
        
        
        
        fname3<-paste("Metabolite_clusterlabels_for_",mainlab,"features.txt",sep="")
        
        mycl_metabs_ord<-mycl_metabs[rev(h73$rowInd)]
        write.table(mycl_metabs_ord,file=fname3,sep="\t",row.names=TRUE)
        
    }
    return(h73)
}
