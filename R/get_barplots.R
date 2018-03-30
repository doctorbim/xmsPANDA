get_barplots <-
function(feature_table_file,class_labels_file,X=NA,Y=NA,parentoutput_dir,newdevice=FALSE,ylabel="Intensity",bar.colors=NA,cex.val=0.75,barplot.col.opt=c("grey57")){
    
    if(is.na(X)==TRUE){
        data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE,quote = "")
    }else{
        data_matrix<-X
        rm(X)
        
    }
    
    dir.create(parentoutput_dir)
    setwd(parentoutput_dir)
    
    mzlabels<-data_matrix[,1]
    
    timelabels<-data_matrix[,2]

#if(is.na(bar.colors)==TRUE){
#   bar.colors<-rep(")
    
    #}

    data_m<-data_matrix[,-c(1:2)]
    
    data_m<-as.matrix(data_m)
    
    col_samples<-TRUE
    
    if(is.na(Y)==TRUE){
        if(is.na(class_labels_file)==FALSE){
            classlabels<-read.table(class_labels_file,sep="\t",header=TRUE,quote = "")
        }else{
            
            classlabels<-rep("classA",dim(data_m)[2])
            classlabels<-cbind(classlabels,classlabels)
            
            col_samples<-FALSE
            
        }
    }else{
        classlabels<-Y
        
    }
    
   
    
    classlabels<-as.data.frame(classlabels)
    
    
    
    
    if(dim(classlabels)[2]>2){
        
        Class<-paste(classlabels[,2],":",classlabels[,3],sep="") #classlabels[,2]:classlabels[,3]
    }else{
        
        Class<-classlabels[,2]
    }
    
    #      print(head(classlabels))
    #print(head(Class))
    
    #print(head(Class))
    #save(Class,file="Class.Rda")
    #save(data_m,file="mtcars.Rda")
    
    class_labels_levels<-levels(as.factor(Class))
    
    mtcars<-t(data_m)
    Class<-as.character(Class)
    
    # print("Class")
    #print(head(Class))
    
    #save(Class,file="Class.Rda")
    #save(data_m,file="mtcars.Rda")
    
    mtcars<-cbind(Class,mtcars)
    mtcars<-as.data.frame(mtcars)
    
    if(newdevice==TRUE){
    
    pdf("barplots.pdf")
    
    }
    
    par_rows=2
    max_per_row=2
    
    par(mfrow=c(par_rows,max_per_row))
    
    myData_list<-new("list")
    
    mtcars <- do.call(data.frame, mtcars)
    
    #  for(i in 2:dim(mtcars)[2]){
    lapply(2:dim(mtcars)[2],function(i){
        x<-as.numeric(as.character(mtcars[,i]))
        
        myData <- aggregate(x,
        by = list(Class = mtcars$Class),
        FUN = function(x) c(mean = mean(x), sd = sd(x),
        n = length(x)))
        
        #myData<-myData[c(2:3,1,4:5),]
        myData <- do.call(data.frame, myData)
        
        myData$se <- myData$x.sd / sqrt(myData$x.n)
        
        colnames(myData) <- c("Class", "Intensity", "sd", "n", "se")
        
        myData_list[[(i-1)]]<-myData
        
        
        # save(myData,file="myData.Rda")
        
        #mzlabel_cur<-paste("mz",mzlabels[(i-1)],sep="")
        
        round_mzval<-sprintf("%.4f",mzlabels[(i-1)])
        round_timeval<-sprintf("%.1f",timelabels[(i-1)])
        
        mzlabel_cur<-paste("mz_time: ",round_mzval,"_",round_timeval,"\n Error bar: mean +/- (1.96 x standard error)",sep="")
        
        
        max_yval<-round(max(myData$Intensity+(3*myData$se),na.rm=TRUE))
        min_yval<-max(0,round(min(myData$Intensity-(3*myData$se),na.rm=TRUE))) #round(min(myData$Intensity,na.rm=TRUE)) #
        
 
        ymax = myData$Intensity + 1.96*myData$se
 
        ymin = myData$Intensity - 1.96*myData$se
        
        below_zero_check<-which(ymin<0)
        if(length(below_zero_check)>0){
            
            ymin[below_zero_check]<-0
        }
        
        
      
        myData$Class<-as.factor(myData$Class)
        
        #95% confidence interval
        limits <- aes(ymax = ymax,
        ymin = ymin)
        
        limits1 <- c(ymax = ymax,
        ymin = ymin)
        
        
        max_yval1<-max(max_yval,max(limits1)+0.5)
        
        max_yval<-max(max_yval,max_yval1)
        
        min_yval1<-min(min_yval,(min(limits1)-0.5))
        min_yval1<-max(0,min_yval1)

        dodge <- position_dodge(width=0.9)
    
        colnames(myData) <- c("Class", "Intensity", "sd", "n", "se")
        
        #p <- ggplot(data = myData, aes(x = Class, y = Intensity, fill = Class)) + geom_bar(stat = "identity", position = dodge,width=0.8) + scale_fill_hue(l=40) +
        #  geom_errorbar(limits, position = dodge, width = 0.3,size=1) + coord_cartesian(ylim=c(min_yval1,max_yval1)) +scale_y_continuous(limits=c(0,max_yval))
        #scale_fill_manual("legend", values = c("A" = "black", "B" = "orange", "C" = "blue"))
        
        # if(length(barplot.col.opt)<length(myData$Class)){
            
            #   barplot.col.opt=rep(barplot.col.opt,length(myData$Class))
            # }
        
        if(length(barplot.col.opt)<2){
            
            barplot.col.opt1=rep(barplot.col.opt,length(myData$Class))
        }else{
            t1<-table(myData$Class)
            
            if(length(barplot.col.opt)==length(myData$Class)){
                
                barplot.col.opt1=rep(barplot.col.opt,t1)
                
                
            }else{
                
                print("Number of classes is greater than the length of the color vector. Using default colors.")
                col_clust<-topo.colors(length(t1))
                barplot.col.opt1=rep(col_clust,t1)
            }
            
            
        }
        
        p <- ggplot(data = myData, aes(x = Class, y = Intensity,fill=Class)) + geom_bar(stat = "identity", position = dodge,width=0.8)  + scale_fill_manual(values = barplot.col.opt1) +
        geom_errorbar(limits, position = dodge, width = 0.3,size=0.5) + coord_cartesian(ylim=c(min_yval1,max_yval1)) +scale_y_continuous(limits=c(0,max_yval))
        
        #+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        #axis.title.x=element_blank())
        
        print(p + ggtitle(mzlabel_cur) + theme(plot.title = element_text(face="bold", color="black", size=9)))
     
	    
    })
    if(newdevice==TRUE){
        dev.off()
    }
    
    save(myData_list,file="barplots_data.Rda")

    
}
