get_pcascoredistplots_child <-
function(X,Y,feature_table_file,parentoutput_dir,class_labels_file,sample.col.opt="rainbow",plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,col_vec=NA,pairedanalysis=FALSE,pca.cex.val=4,legendlocation="topright",pca.ellipse=TRUE,ellipse.conf.level=0.5,filename="all",paireddesign=NA)
{
		
        #print("here")
			if(is.na(X)==TRUE){
			data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE,quote = "")
			
			}else{
			
				data_matrix<-X
			}
			if(is.na(Y)==TRUE){
			classlabels<-read.table(class_labels_file,sep="\t",header=TRUE,quote = "")
			
			}else{
			
				classlabels<-Y
			}
            #print("here")
dir.create(parentoutput_dir)
setwd(parentoutput_dir)

data_m<-data_matrix[,-c(1:2)]

data_m<-as.matrix(data_m)

mzvec<-data_matrix[,1]
rtvec<-data_matrix[,2]

rnames<-paste(mzvec,rtvec,sep="_")

rownames(data_m)<-as.character(rnames)


classlabels_orig<-classlabels

#print("dim classlabels_orig")
#print(dim(classlabels_orig))

#print(head(classlabels_orig))
#print("dim of X")
#print(dim(X))

#class_labels_levels<-paste("x",seq(1,length(class_labels_levels)),sep="")

#if(pairedanalysis==FALSE){

if(dim(classlabels)[2]>2){
    
    #classlabels_orig[,2]<-as.factor(paste("A",as.character(classlabels_orig[,2]),sep=""))
    #classlabels_orig[,3]<-as.factor(paste("B",as.character(classlabels_orig[,3]),sep=""))
    # print(head(classlabels_orig))
    
    classgroup<-paste(classlabels_orig[,2],":",classlabels_orig[,3],sep="") #classlabels_orig[,2]:classlabels_orig[,3]
do_pca_anova=FALSE
}else{

	classgroup<-classlabels_orig[,2]
    
    do_pca_anova=TRUE
}



class_labels_levels<-levels(as.factor(classgroup))
  ordered_labels<-classgroup
  
  class_label_alphabets<-paste("C",1:length(class_labels_levels),sep="") #c("A","B","C","D","E","F","G","H","I","J","K","L","M")
							

#if(is.na(col_vec)==TRUE)
{
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
}
			
            #	print(class_labels_levels)
			
						ordered_labels={}
							num_samps_group<-new("list")
							num_samps_group[[1]]<-0
							groupwiseindex<-new("list")
							groupwiseindex[[1]]<-0
						
                            S<-new("list")
							for(c in 1:length(class_labels_levels))
							{
								
								classlabels_index<-which(classlabels[,2]==class_labels_levels[c])
								#ordered_labels<-c(ordered_labels,as.character(classlabels[classlabels_index,2]))
								num_samps_group[[c]]<-length(classlabels_index)
								groupwiseindex[[c]]<-classlabels_index
                                
                                
                                # S[[c]]<-cov(t(data_m[,c(classlabels_index)]))
                                
							}
                            
                           
			
						sampleclass<-{}
						patientcolors<-{}
						
						classlabels<-as.data.frame(classlabels)
						
                       
						for(c in 1:length(class_labels_levels)){
                           
									num_samps_group_cur=length(which(ordered_labels==class_labels_levels[c]))
									
									
									sampleclass<-c(sampleclass,rep(paste("Class",class_label_alphabets[c],sep=""),num_samps_group_cur))

									patientcolors <-c(patientcolors,rep(col_vec[c],num_samps_group[[c]]))
								}

#print("sampleclass")
#print(sampleclass)

#print("classlabels")
#print(classlabels)


if(length(mzvec)>4){
max_per_row<-3


par_rows<-ceiling(9/max_per_row)

}else{
max_per_row<-length(mzvec)
par_rows<-1
}

file_ind<-0
			boxplots_fname<-paste("xyplots.pdf",sep="")
				#tiff(boxplots_fname, width=plots.width,height=plots.height,res=plots.res, compression="lzw")
				
		
        #print("dim")
        #print(dim(data_m))
if(dim(data_m)[1]>2){
#return(list(data_m=data_m,samplelabels=classgroup))
#res<-get_pca(X=data_m,samplelabels=classgroup,outloc=getwd(),ncomp=3,center=TRUE,scale=TRUE,col_vec=patientcolors)

#pca_res<-pca(t(data_m))

#plotIndiv(pca_res)

t1<-table(classgroup)
    l1<-levels(classgroup)

   l1<-levels(as.factor(classgroup))


patientcolors <- rep(col_vec[1:length(t1)], t1)


#print(head(classlabels_orig))


res<-get_pca(X=data_m,samplelabels=classgroup,legendlocation=legendlocation,filename=filename,ncomp=10,center=TRUE,scale=TRUE,legendcex=0.5,outloc=getwd(),col_vec=col_vec,sample.col.opt=sample.col.opt,alphacol=0.3,class_levels=NA,pca.cex.val=pca.cex.val,pca.ellipse=pca.ellipse,do_pca_anova=do_pca_anova,paireddesign=paireddesign)



#if(filename=="all")
{


fname<-paste("PCAloadings_",filename,"features.txt",sep="")

loadings_res<-res$rotation
scores_res<-res$x

loadings_res<-round(loadings_res,2)
scores_res<-round(scores_res,2)

if(dim(loadings_res)[2]>10){
    loadings_res<-loadings_res[,c(1:10)]
    scores_res<-scores_res[,c(1:10)]
}

write.table(loadings_res,file=fname,sep="\t")

fname<-paste("PCAscores_",filename,"features.txt",sep="")

scores_res1<-cbind(classlabels_orig,scores_res)


write.table(scores_res1,file=fname,sep="\t",row.names=FALSE)


}

pcnum_limit<-min(10,dim(scores_res)[2])

get_pooled_sp<-function(n1,n2,S1,S2){a<-(n1-1)*S1;b<-(n2-1)*S2;c<-(n1+n2-2);return((a+b)/c)}
#S2<-cov(scores_res[16:46,1:2])
#S1<-cov(scores_res[1:15,1:2])

 class_levels<-levels(as.factor(classlabels_orig[,2]))
 
 #if(length(class_labels_levels)==2)
 if(do_pca_anova==TRUE)
 {
     if(FALSE){
       h1<-hotelling.test(scores_res[1:num_samps_group[[1]],c(1)],scores_res[(1+num_samps_group[[1]]):(num_samps_group[[1]]+num_samps_group[[2]]),c(1)],shrinkage=FALSE)
      print("Hotelling test using PC1: ")
      print(h1)
     
       h2<-hotelling.test(scores_res[1:num_samps_group[[1]],c(2)],scores_res[(1+num_samps_group[[1]]):(num_samps_group[[1]]+num_samps_group[[2]]),c(2)],shrinkage=FALSE)
       print("Hotelling test using PC2: ")
     print(h2)
     
     h3<-hotelling.test(scores_res[1:num_samps_group[[1]],c(3)],scores_res[(1+num_samps_group[[1]]):(num_samps_group[[1]]+num_samps_group[[2]]),c(3)],shrinkage=FALSE)
     print("Hotelling test using PC3: ")
     print(h3)
     }
     
     pc_pval_vec<-{}
	
 

     for(pcnum in 1:pcnum_limit){
         
         pc1_pval<-anova(lm(cbind(scores_res[,pcnum])~classlabels_orig[,2]))
         pc1_pval<-round(pc1_pval[[5]][1],3)
         pc_pval_vec<-c(pc_pval_vec,pc1_pval)
     }
     
     
     if(FALSE){
     pc1_pval<-anova(lm(cbind(scores_res[,1])~classlabels_orig[,2]))
     
     pc1_pval<-pc1_pval[[5]][1]
     pc2_pval<-anova(lm(cbind(scores_res[,2])~classlabels_orig[,2]))
     
     pc2_pval<-pc2_pval[[5]][1]

     pc3_pval<-anova(lm(cbind(scores_res[,3])~classlabels_orig[,2]))
     pc3_pval<-pc3_pval[[5]][1]
     
     pc_pval_vec<-c(pc1_pval,pc2_pval,pc3_pval)
     }

 }
# save(num_samps_group,scores_res,file="HotellingTestInput.Rda")
#h1<-hotelling.test(scores_res[1:15,c(4:5)],scores_res[16:46,c(4:5)],shrinkage=FALSE)
#h1<-hotelling.test(scores_res[1:15,c(2:3)],scores_res[16:46,c(2:3)],shrinkage=FALSE)

 df_matrix<-{}


				if(dim(classlabels_orig)[2]>2){
                    class_levels<-levels(as.factor(classlabels_orig[,2]))
				t1<-table(classlabels_orig[,3])
                #print("Doing two factors")
                
                #print(class_levels)
               
for(pc in 1:pcnum_limit){




#print("Debug xyplot")
if(pairedanalysis==FALSE){
	 df<-cbind(res$x[,pc],classgroup,classlabels_orig[,2],classlabels_orig[,3])
     #  mzname<-paste("PC",pc," scores group-wise distribution (25th percentile, median, 75th percentile) \nusing ",filename," feats",sep="")
     mzname<-paste("PC",pc," scores distribution (25th percentile, median, 75th percentile) \n in each group using ",filename," feats vs factors",sep="")

}else{                  


	 df<-cbind(res$x[,pc],classgroup,classlabels_orig[,2],classlabels_orig[,3])
     # mzname<-paste("PC",pc," scores group-wise distribution (25th percentile, median, 75th percentile) \nusing", filename," feats",sep="")
     mzname<-paste("PC",pc," scores distribution (25th percentile, median, 75th percentile) \n in each group using ",filename," feats vs time",sep="")
     
} 

#print(class_levels)
		#pdf("pc3_scoreplots.pdf")
			fname<-paste("pc_",pc,"_scoreplot",".tiff",sep="")
			par_rows=3
			#tiff(fname)
		
						colnames(df)<-c("y","x","Category","time")
						df=as.data.frame(df)
						df$y<-as.numeric(as.character(df$y))

                        # save(df,file="df.Rda")

                        			df_fname<-paste("dftable_pc",pc,".txt",sep="")
                                    
                                    if(pc>1){
                                        if(pc==2){
                                            PC2Score=df[,1]
                                            df_matrix<-cbind(df_matrix,PC2Score)
                                        }else{
                                            
                                            
                                            PC3Score=df[,1]
                                            df_matrix<-cbind(df_matrix,PC3Score)
                                            
                                            
                                        }
                                    }else{
                                        PC1Score<-df[,1]
                                        df_matrix<-cbind(df[,-c(1)],PC1Score)
                                        
                                    }
                                    
                                    #save(list=ls(),file="pcadebug.Rda")
                                    #    write.table(df,file=df_fname,sep="\t",row.names=FALSE)
						df.summary <- df %>% group_by(x) %>%
                        
                        dplyr::summarize(ymin = quantile(y,0.25),
						      ymax = quantile(y,0.75),
						      Score = median(y),number=length(y),Category=unique(as.character(Category)),Time=unique(as.character(time)))
					

						      Category<-{}
						      for(cnum in 1:length(class_levels)){
						      
							Category<-c(Category,rep(class_levels[cnum],length(t1))) 
							
							df.summary$Category[which(df.summary$Category==cnum)]<-class_levels[cnum]
							
							
							}
					Category<-unique(Category)
       
						time.hour<- c(unique(as.character(classlabels_orig[,3])),unique(as.character(classlabels_orig[,3])))
                        #Score<-df.summary$ymean
						
						#print(df.summary)
                        #save(df.summary,file="df.summary.Rda")
                        df.summary$x<- df.summary$Time #time.hour
                        
                        #p <- ggplot(data = myData, aes(x = Class, y = Intensity, fill = Class)) + geom_bar(stat = "identity", position = dodge) + scale_fill_hue(l=40) +
                        #    geom_errorbar(limits, position = dodge, width = 0.25) + coord_cartesian(ylim=c(min_yval1,max_yval1)) +scale_y_continuous(limits=c(0,max_yval))
                        
					if(pairedanalysis==TRUE){	
						plot_res<-ggplot(df.summary, aes(x = x, y = Score,color = Category)) + geom_point(size = pca.cex.val) + geom_line(aes(group =Category))  + geom_errorbar(aes(ymin = ymin, ymax = ymax)) + scale_color_hue(l=40)
					}else{

						plot_res<-ggplot(df.summary, aes(x = x, y = Score,color = Category)) + geom_point(size = pca.cex.val) + geom_errorbar(aes(ymin = ymin, ymax = ymax)) + scale_color_hue(l=40)
					}
						
						print(plot_res + ggtitle(mzname))
				
		}
}else{

#print("Doing one factor")
#print(pairedanalysis)
				class_levels<-levels(as.factor(classgroup))
				t1<-table(classgroup)
for(pc in 1:pcnum_limit){
df<-cbind(res$x[,pc],classlabels_orig[,2],classlabels_orig[,2])

#mzname<-paste("PC",pc," scores vs factor",sep="")
#df<-cbind(res$x[,pc],classlabels[,1],classlabels[,1])

if(pairedanalysis==FALSE){
        mzname<-paste("PC",pc," scores distribution (25th percentile, median, 75th percentile) \n in each group using ",filename," feats vs factors (p=",pc_pval_vec[pc],")",sep="")
}else{ 

        mzname<-paste("PC",pc," scores distribution (25th percentile, median, 75th percentile) \n in each group using ",filename," feats vs time",sep="")

}

			fname<-paste("pc_",pc,"_scoreplot",".tiff",sep="")
			par_rows=3
		
						colnames(df)<-c("y","x","Category")
						df=as.data.frame(df)
						df$y<-as.numeric(as.character(df$y))

#print(dim(df))

					df_fname<-paste("dftable_pc",pc,".txt",sep="")
                    # write.table(df,file=df_fname,sep="\t",row.names=FALSE)
                    
                    if(pc>1){
                        if(pc==2){
                        PC2Score=df[,1]
                        df_matrix<-cbind(df_matrix,PC2Score)
                        }else{
                            
                            
                                PC3Score=df[,1]
                                df_matrix<-cbind(df_matrix,PC3Score)
                            
                            
                        }
                    }else{
                        PC1Score<-df[,1]
                        df_matrix<-cbind(df[,-c(1)],PC1Score)
                        
                    }

						df.summary <- df %>% group_by(x) %>%
					    dplyr::summarize(ymin = quantile(y,0.25),
						      ymax = quantile(y,0.75),
						      Score = median(y),number=length(y),Category=unique(as.character(Category)))
						
						      Category<-{}
						      for(cnum in 1:length(class_levels)){
						      
							Category<-c(Category,rep(class_levels[cnum],length(t1))) 
						
							if(pairedanalysis==FALSE){	
							df.summary$Category[which(df.summary$Category==cnum)]<-class_levels[cnum]
							}else{
								df.summary$Category[which(df.summary$Category==cnum)]<-class_levels[1]
							}
							
							}
						
                        # save(df.summary,file="df.summary.Rda")
                        
						time.hour<- c(unique(as.character(classlabels_orig[,2])),unique(as.character(classlabels_orig[,2])))
                        #Score<-df.summary$Score
						#df.summary$x<-time.hour
				
							
						if(pairedanalysis==FALSE){
						plot_res<-ggplot(df.summary, aes(x = x, y = Score,color = Category)) + geom_point(size = pca.cex.val)  + geom_errorbar(aes(ymin = ymin, ymax = ymax)) + scale_x_continuous(breaks=seq(1,length(class_levels))) + scale_color_hue(l=40)
						}else{
					

#plot_res<-ggplot(df.summary, aes(x = x, y = Score,color = Category)) + geom_point(size = 2) + geom_line()  #+ geom_errorbar(aes(ymin = ymin, ymax = ymax))
plot_res<-ggplot(df.summary, aes(x = x, y = Score,color = Category)) +  geom_point(size = pca.cex.val) + geom_line(aes(group =Category))  + geom_errorbar(aes(ymin = ymin, ymax = ymax)) + scale_color_hue(l=40)

#plot_res<-ggplot(df.summary, aes(x = x, y = Score,color = Category)) +  geom_point() + geom_line(aes(group =Category)) # + geom_line()

						}
						print(plot_res + ggtitle(mzname)) #
						#print(plot_res,height=2000,width=2000)
						
			#		dev.off()
		}
}
}

df_fname<-paste("PC_score_distribution_matrix_",filename,"features.txt",sep="")
#write.table(df_matrix,file=df_fname,sep="\t",row.names=TRUE)

		#dev.off()
	}
