get_individualsampleplots <-
function(feature_table_file,class_labels_file,X=NA,Y=NA,parentoutput_dir,newdevice=FALSE,ylabel="Intensity",bar.colors=NA,cex.val=0.75,sample.col.opt=c("grey57"),plottype="barplot"){
    
    if(is.na(X)==TRUE){
        data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE,quote = "")
    }else{
        data_matrix<-X
        rm(X)
        
    }
    
    
    alphacol=0.3
    
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
    
    class_labels_levels<-levels(as.factor(Class))
    
    
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
                        
                        col_vec<-c("grey57")
                    }
                    
                    
                }
                
            }
            
        }
    }
    
    mtcars<-t(data_m)
    Class<-as.character(Class)
    
    barplot.col.opt=col_vec
    
    mtcars<-cbind(Class,mtcars)
    mtcars<-as.data.frame(mtcars)
    
    mtcars<-mtcars[order(mtcars$Class),]
    
    if(newdevice==TRUE){
        
        pdf("individual.sample.plots.pdf")
        
    }
    
    par_rows=2
    max_per_row=2
    
    par(mfrow=c(par_rows,max_per_row))
    
    myData_list<-new("list")
    
    mtcars <- do.call(data.frame, mtcars)
    
    #for(i in 2:dim(mtcars)[2]){
    lapply(2:dim(mtcars)[2],function(i){
        x<-as.numeric(as.character(mtcars[,i]))
        
        myData <- x
        
        myData<-cbind(mtcars$Class,myData)
        colnames(myData)<-c("Class","Intensity")
        #myData<-myData[c(2:3,1,4:5),]
        
        myData <- as.data.frame(myData)
        myData<-myData[order(myData$Class),]
        
        
        
        round_mzval<-sprintf("%.4f",mzlabels[(i-1)])
        round_timeval<-sprintf("%.1f",timelabels[(i-1)])
        
        mzlabel_cur<-paste("mz_time: ",round_mzval,"_",round_timeval,sep="")
        
        if(length(barplot.col.opt)<2){
            
            barplot.col.opt1=rep(barplot.col.opt,length(mtcars$Class))
        }else{
            t1<-table(mtcars$Class)
            
            if(length(barplot.col.opt)==length(levels(factor(mtcars$Class)))){
        
                    barplot.col.opt1=rep(barplot.col.opt,t1)
            
            
            }else{
                
                print("Number of classes is greater than the length of the color vector. Using default colors.")
                col_clust<-topo.colors(length(t1))
                barplot.col.opt1=rep(col_clust,t1)
            }
            
            
        }
        
        #print(barplot.col.opt1)
        #barplot.col.opt1=barplot.col.opt
        
        # print(barplot.col.opt1)
        
        myData$Class<-factor(myData$Class)
        if(plottype=="barplot"){
            #barplot(as.vector(myData[,2]),col=c(barplot.col.opt1),main=mzlabel_cur, ylab="Intensity",xlab="Sample") #,ylim=c(min(myData[,2])-1,max(myData[,2])+1),xpd=FALSE)
             plot(as.vector(myData[,2]),col=c(barplot.col.opt1),main=mzlabel_cur, ylab="Intensity",xlab="Sample",type="h",lwd=2)
        }else{
            if(plottype=="point"){
                plot(as.vector(myData[,2]),col=c(barplot.col.opt1),main=mzlabel_cur, ylab="Intensity",xlab="Sample",type="p")
            }else{
                
                if(plottype=="point_grouped"){
                    p <- ggplot(data = myData, aes(x = Class, y = Intensity,fill=Class)) + scale_fill_manual(values = c(barplot.col.opt1))
                    print(p+geom_point(size=2)+ggtitle(mzlabel_cur)) #,aes(colour=barplot.col.opt))
                }
            }
        }
        #p <- ggplot(data = myData, aes(x = Class, y = Intensity,fill=Class)) + scale_fill_manual(values = barplot.col.opt)

        #p <- ggplot(data = myData, aes(x = factor(Class), y = Intensity)) + geom_point(size=4,aes(colour=barplot.col.opt))
        #print(p + ggtitle(mzlabel_cur)) # +
        
        
    })
    if(newdevice==TRUE){
        dev.off()
    }
    
   
    
    
}
