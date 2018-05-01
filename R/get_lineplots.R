get_lineplots <-
function(X=NA,Y=NA,feature_table_file=NA,parentoutput_dir=NA,class_labels_file=NA,sample.col.opt="rainbow",alphacol=0.3,col_vec=NA,pairedanalysis=FALSE,pca.cex.val=4,legendlocation="topright",pca.ellipse=TRUE,ellipse.conf.level=0.5,filename="all",newdevice=FALSE,lineplot.col.opt=c("grey57"))
{
    
    if(is.na(parentoutput_dir)==TRUE){
        
        parentoutput_dir=getwd()
    }
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
   
    
    # print("Class group is")
    #print(classgroup)
    
    # print("HERE")
    
    if(newdevice==TRUE){
        
        fname<-paste("timeseriesplots",filename,".pdf",sep="")
        pdf(fname)
    }
    #class_labels_levels<-levels(as.factor(classgroupA))
    #ordered_labels<-classgroupA
    
    #class_label_alphabets<-paste("C",1:length(class_labels_levels),sep="") #c("A","B","C","D","E","F","G","H","I","J","K","L","M")
    
    
    #if(is.na(col_vec)==TRUE)
    if(is.na(sample.col.opt)==FALSE)
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
                        }else{
                            if(is.na(sample.col.opt)==TRUE){
                                col_vec<-c("black")
                            }else{
                                col_vec <- sample.col.opt
                            }
                        }
                        
                        
                    }
                    
                }
                
            }
        }
    }else{
        
        col_vec<-c("black")
    
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
                               
                                num_samps_group[[c]]<-length(classlabels_index)
                                groupwiseindex[[c]]<-classlabels_index
                                
                                
                               
                                
    }
    
    
    
    sampleclass<-{}
    patientcolors<-{}
    
    classlabels<-as.data.frame(classlabels)
    
    
    for(c in 1:length(class_labels_levels)){
        
        num_samps_group_cur=length(which(ordered_labels==class_labels_levels[c]))
        
        
        sampleclass<-c(sampleclass,rep(paste("Class",class_label_alphabets[c],sep=""),num_samps_group_cur))
        
        patientcolors <-c(patientcolors,rep(col_vec[c],num_samps_group[[c]]))
    }
    

    
    if(length(mzvec)>4){
        max_per_row<-3
        
        
        par_rows<-ceiling(9/max_per_row)
        
    }else{
        max_per_row<-length(mzvec)
        par_rows<-1
    }
    
    # The palette with black:
    cbPalette <- c("#E69F00","#000000","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    # To use for fills, add
    # scale_fill_manual(values=cbPalette)
    
    # To use for line and point colors, add
    #scale_colour_manual(values=cbPalette)
    
    file_ind<-0
    boxplots_fname<-paste("xyplots.pdf",sep="")
			
                if(dim(data_m)[1]>2){
                
                    
                    t1<-table(classgroup)
                    l1<-levels(classgroup)
                    
                    l1<-levels(as.factor(classgroup))
                    
                    
                    patientcolors <- rep(col_vec[1:length(t1)], t1)
                    
                    
                    
                    
                    
                    get_pooled_sp<-function(n1,n2,S1,S2){a<-(n1-1)*S1;b<-(n2-1)*S2;c<-(n1+n2-2);return((a+b)/c)}
             
                    
                    class_levels<-levels(as.factor(classlabels_orig[,2]))
                    
               
                    df_matrix<-{}
                    
                    
                    if(dim(classlabels_orig)[2]>2){
                        class_levels<-levels(as.factor(classlabels_orig[,2]))
                        t1<-table(classlabels_orig[,3])
                        #print("Doing two factors")
                        
                        #print(class_levels)
                        #length(rnames)
                        for(pc in 1:length(rnames)){
                            
                            
                            
                            
                            #print("Debug xyplot")
                            if(pairedanalysis==FALSE){
                                xvec<-(data_m[pc,])
                                
                                
                                df<-cbind(xvec,classgroup,classlabels_orig[,2])
                                mzname<-paste(rnames[pc]," distribution (25th percentile, median, 75th percentile) \n in each group using ",filename," feats vs factors",sep="")
                                
                            }else{
                                
                                
                                xvec<-(data_m[pc,])
                               
                               
                               
                                
                                mzname<-paste(rnames[pc]," distribution (25th percentile, median, 75th percentile) \n in each group using ",filename," feats vs time",sep="")
                                
                            }
                            
                            df<-cbind(xvec,classgroup,classlabels_orig[,2],classlabels_orig[,3])
                            
                            colnames(df)<-c("y","x","Category","time")
                            df=as.data.frame(df)
                            df$y<-as.numeric(as.character(df$y))
                            #dfname<-paste("df",rnames[pc],".Rda",sep="")
                            #save(df,file=dfname)
                            
                            
                            #    write.table(df,file=df_fname,sep="\t",row.names=FALSE)
                            df.summary <- df %>% group_by(x) %>%
                            
                            dplyr::summarize(ymin = quantile(y,0.25),
                            ymax = quantile(y,0.75),
                            Intensity = median(y),number=length(y),Category=unique(as.character(Category)),Time=unique(as.character(time)))
                            
                            
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
                            df.summary$x<-time.hour
                           
                            
                            if(pairedanalysis==TRUE){
                                plot_res<-ggplot(df.summary, aes(x = x, y = Intensity,color = Category)) + geom_point(size = pca.cex.val) + geom_line(aes(group =Category))  + geom_errorbar(aes(ymin = ymin, ymax = ymax)) + scale_color_hue(l=40) #+ scale_colour_manual(values=cbPalette)
                            }else{
                                
                                plot_res<-ggplot(df.summary, aes(x = x, y = Intensity,color = Category)) + geom_point(size = pca.cex.val) + geom_errorbar(aes(ymin = ymin, ymax = ymax)) + scale_color_hue(l=40) #+ scale_colour_manual(values=cbPalette)
                            }
                            
                            print(plot_res + ggtitle(mzname))
                            
                        }
                    }else{
                        
                        #print("Doing one factor")
                        #print(pairedanalysis)
                        class_levels<-levels(as.factor(classgroup))
                        t1<-table(classgroup)
                        #for(pc in 1:length(rnames)){
                        for(pc in 1:length(rnames)){
                            #df<-cbind(res$x[,pc],classlabels_orig[,2],classlabels_orig[,2])
                            
                            xvec<-(data_m[pc,])
                            df<-cbind(xvec,classlabels_orig[,2],classlabels_orig[,2])
                            
                            
                            
                            #mzname<-paste("PC",pc," scores vs factor",sep="")
                            #df<-cbind(res$x[,pc],classlabels[,1],classlabels[,1])
                            
                            if(pairedanalysis==FALSE){
                                # mzname<-paste("PC",pc," scores distribution (25th, median, 75th percentile) \n in each group using ",filename," feats vs factors (p=",pc_pval_vec[pc],")",sep="")
                                
                                mzname<-paste(rnames[pc]," distribution (25th percentile, median, 75th percentile) \n in each group using ",filename," feats vs factors",sep="")
                            }else{ 
                                
                                #mzname<-paste("PC",pc," scores distribution (25th, median, 75th percentile) \n in each group using ",filename," feats vs time",sep="")
                                
                                mzname<-paste(rnames[pc]," distribution (25th percentile, median, 75th percentile) \n in each group using ",filename," feats vs time",sep="")
                                
                            }
                            
                                                        par_rows=3
                            
                            colnames(df)<-c("y","x","Category")
                            df=as.data.frame(df)
                            df$y<-as.numeric(as.character(df$y))
                            
                            #print(dim(df))
                            
                            #df_fname<-paste("dftable_pc",pc,".txt",sep="")
                            # write.table(df,file=df_fname,sep="\t",row.names=FALSE)
                            
                            
                            df.summary <- df %>% group_by(x) %>%
                            dplyr::summarize(ymin = quantile(y,0.25),
                            ymax = quantile(y,0.75),
                            Intensity = median(y),number=length(y),Category=unique(as.character(Category)))
                            
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
                            
                            #print(Category)
                            
                            if(pairedanalysis==FALSE){
                                plot_res<-ggplot(df.summary, aes(x = x, y = Intensity,color = Category)) + geom_point(size = pca.cex.val)  + geom_errorbar(aes(ymin = ymin, ymax = ymax)) + scale_x_continuous(breaks=seq(1,length(class_levels))) + scale_color_hue(l=40) #+ scale_colour_manual(values=cbPalette)
                            }else{
                                
                                
                                #plot_res<-ggplot(df.summary, aes(x = x, y = Score,color = Category)) + geom_point(size = 2) + geom_line()  #+ geom_errorbar(aes(ymin = ymin, ymax = ymax))
                                
                                plot_res<-ggplot(df.summary, aes(x = x, y = Intensity,color = Category)) +  geom_point(size = pca.cex.val) + geom_line(aes(group =Category))  + geom_errorbar(aes(ymin = ymin, ymax = ymax)) + scale_color_hue(l=40) #+ scale_colour_manual(values=cbPalette)
                                
                                #plot_res<-ggplot(df.summary, aes(x = x, y = Intensity,color = Category)) +  geom_point(size = pca.cex.val,colour=col_vec) + geom_line(aes(group =Category),colour=col_vec)  + geom_errorbar(aes(ymin = ymin, ymax = ymax),colour=col_vec) + scale_color_hue(l=40)
                                
                                #plot_res<-ggplot(df.summary, aes(x = x, y = Score,color = Category)) +  geom_point() + geom_line(aes(group =Category)) # + geom_line()
                                
                            }
                            print(plot_res + ggtitle(mzname)) #
                            #print(plot_res,height=2000,width=2000)
                            
                            #		dev.off()
                        }
                    }
                }
                
                #df_fname<-paste("PC_score_distribution_matrix_",filename,"features.txt",sep="")
                #write.table(df_matrix,file=df_fname,sep="\t",row.names=TRUE)
                if(newdevice==TRUE){
                    
                    try(dev.off(),silent=TRUE)
                }
                
}
