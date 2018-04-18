get_boxplots_child <-
function(X,Y,feature_table_file,parentoutput_dir,class_labels_file,boxplot.col.opt="grey57",sample.col.opt="rainbow",alphacol=0.3,newdevice=TRUE,cex=0.8,replace.by.NA=FALSE,pairedanalysis=FALSE,filename="")
{
		

			if(is.na(X)==TRUE){
			data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE)
			
			}else{
			
				data_matrix<-X
			}
			if(is.na(Y)==TRUE){
			classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
			
			}else{
			classlabels<-Y
			
			}
dir.create(parentoutput_dir)
setwd(parentoutput_dir)

data_m<-data_matrix[,-c(1:2)]

data_m<-as.matrix(data_m)

mzvec<-data_matrix[,1]
timevec<-data_matrix[,2]
goodfeats<-data_m			
rm(data_m)



if(dim(classlabels)[2]>2){
    
    print("More than two columns found in the class labels file.")
    if(pairedanalysis==TRUE){
        
        Class<-classlabels[,3]
        
        classlabels<-classlabels[,-c(2)]
        
        
    }else{
    
        Class<-paste(classlabels[,2],":",classlabels[,3],sep="") #classlabels[,2]:classlabels[,3]
    
    }
}else{
    
    Class<-classlabels[,2]
}




#print(head(Class))
#patientcolors<-rep("green",dim(data_m)[2])

#class_labels_levels<-levels(as.factor(classlabels[,2]))
# ordered_labels<-classlabels[,2]

class_labels_levels<-levels(as.factor(Class))
ordered_labels<-Class

  class_label_alphabets<-paste("C",1:length(class_labels_levels),sep="") #c("A","B","C","D","E","F","G","H","I","J","K","L","M")
  
  if(is.na(sample.col.opt)==TRUE){
      
      col_vec<-rep(c("black"),length(class_labels_levels))
      sample.col.opt<-col_vec
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
                                        }else{
                                            
                                            col_vec <-sample.col.opt
                                        }
									
											
									}
							
						}
			
		}	
}

			
           
						ordered_labels={}
							num_samps_group<-new("list")
							num_samps_group[[1]]<-0
							groupwiseindex<-new("list")
							groupwiseindex[[1]]<-0
						
							for(c in 1:length(class_labels_levels))
							{
								
								classlabels_index<-which(Class==class_labels_levels[c])
								#ordered_labels<-c(ordered_labels,as.character(classlabels[classlabels_index,2]))
								num_samps_group[[c]]<-length(classlabels_index)
								groupwiseindex[[c]]<-classlabels_index
							}
			
			sampleclass<-{}
						patientcolors<-{}
						#print(classlabels)
                        #	classlabels<-as.data.frame(classlabels)
						
						
					
                    #sampleclass<-classlabels

			
#print(length(mzvec))

if(length(mzvec)>4){
max_per_row<-3


par_rows<-ceiling(9/max_per_row)

}else{
max_per_row<-length(mzvec)
par_rows<-1
}

#class_labels_levels<-paste("x",seq(1,length(class_labels_levels)),sep="")

file_ind<-0
			boxplots_fname<-paste("boxplots",filename,".pdf",sep="")
				#tiff(boxplots_fname, width=plots.width,height=plots.height,res=plots.res, compression="lzw")
				
				#tiff(boxplots_fname, width=2000,height=3000,res=plots.res, compression="lzw")
                
                if(newdevice==TRUE){
				pdf(boxplots_fname)
                }
				par(mfrow=c(par_rows,max_per_row))
            
		for(m in 1:dim(goodfeats)[1])
		{
            
			if(m%%9==0){
				#dev.off()
				file_ind<-file_ind+1
				boxplots_fname<-paste("boxplots_file",file_ind,".tiff",sep="")
				#tiff(boxplots_fname, width=plots.width,height=plots.height,res=plots.res, compression="lzw")
				#tiff(boxplots_fname, width=2000,height=3000,res=plots.res, compression="lzw")
				#pdf(boxplots_fname)
				par(mfrow=c(par_rows,max_per_row))
			}
			
			round_mzval<-sprintf("%.4f",mzvec[m])
            
            round_timeval<-sprintf("%.1f",timevec[m])
            
            mzname<-paste("mz_time: ",round_mzval,"_",round_timeval,sep="")

#mzname<-paste("mz ",round_mzval,sep="")
			if(length(class_labels_levels)>=2)
			{
			
					if(length(class_labels_levels)>=1)
					{
					
					t1<-table(sampleclass)
					cur_d<-{}
					for(c in 1:length(class_labels_levels))
					{
					
					
					num_samps_group[[1]]<-t1[1]
					
					cvec<-as.vector(t(goodfeats[m,c(groupwiseindex[[c]])]))
					
					cvec<-replace_outliers(cvec,replace.by.NA)
					cur_d<-cbind(cur_d,cvec)
					}
					
					cur_d<-as.data.frame(cur_d)
					
					#class_labels_boxplot<-levels(classlabelsA) #paste(seq(1,length(class_labels_levels)),sep="")
					
					#print(class_labels_boxplot)
					colnames(cur_d)<-NULL #paste(seq(1,length(class_labels_levels)),sep="") #as.character(class_labels_levels)
					cur_d<-round(cur_d,2)
					#print(dim(cur_d))
					#print(cur_d[1:4,])
					
					boxplot(cur_d,ylab="Intensity",main=mzname,xaxt="n",cex.main=cex,col=boxplot.col.opt)
					
					for(i in 1:length(class_labels_levels)){
					axis(side=1,at=c(i),labels=class_labels_levels[i], col=col_vec[i],cex.axis=cex)
					
					#text(, round(pcavar1,2), labels = round(pcavar1,2), pos = 3)
					}
					
					}
					
					
					
			}else{

			boxplot(x1,x2, ylab="Intensity",main=mzname,cex.main=cex,col=boxplot.col.opt)
		axis(side=1,at=seq(1),labels=class_labels_levels[1], col="red")
			axis(side=1,at=2,labels=class_labels_levels[2],col="green")
			}
		}
        if(newdevice==TRUE){
		dev.off()
        }
        
        par(mfrow=c(1,1))
	}
