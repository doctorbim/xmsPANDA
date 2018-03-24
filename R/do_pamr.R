do_pamr <-
function(X,Y,fdrthresh=0.1,nperms=100,pamr.threshold.select.max=FALSE){
	
	
	library(pamr)
	if(FALSE){
	x=t(leukemia$X)
    y=leukemia$Y
    x=khan.data[-c(1),-c(1:2)]
    
    x=apply(x,2,as.numeric)
    
    y=as.factor(t(khan.data[1,-c(1:2)]))
    genenames=khan.data$X1 #paste("g",as.character(1:nrow(x)),sep="")
    geneid=khan.data$X #as.character(1:nrow(x))
    
	mydata<-list(x=x,y=y,geneid=geneid,genenames=genenames)
	
    mytrain<-   pamr.train(mydata)
	pfdr<-pamr.fdr(mytrain,mydata)
	pamr.listgenes(mytrain, mydata, threshold=0)
	}
	#d1<-list(x=X,y=Y)
	
	d1 <- list(x=as.matrix(X),y=factor(Y[,1]), geneid=as.character(1:nrow(X)),
      genenames=paste("g",as.character(1:nrow(X)),sep=""))
     
	p1<-pamr.train(d1)

    set.seed(999)
	p2<-pamr.cv(data=d1,fit=p1,nfold=10)
	
    threshold_CVerror_matrix<-cbind(p2$threshold,p2$error)
    threshold_CVerror_matrix<-as.data.frame(threshold_CVerror_matrix)
    colnames(threshold_CVerror_matrix)<-c("Threshold","error")
    
    threshold_value<-max(threshold_CVerror_matrix[which(threshold_CVerror_matrix[,2]==min(threshold_CVerror_matrix[,2])),1])[1]
    
    set.seed(999)
	p3<-pamr.fdr(data=d1,p1,nperms=nperms)
	
    selected_feature_index<-{}
    max.discore<-{}
    
    if(length(which(p3$results<fdrthresh))>0){
    pamr_fdr_filt<-p3$results[which(p3$results[,5]<fdrthresh),]
    threshold_CVerror_fdrmatrix<-merge(threshold_CVerror_matrix,pamr_fdr_filt,by="Threshold")
    
    if(pamr.threshold.select.max==TRUE){
    threshold_value<-max(threshold_CVerror_fdrmatrix[which(threshold_CVerror_fdrmatrix[,2]==min(threshold_CVerror_fdrmatrix[,2])),1])[1]
    }else{
        threshold_value<-min(threshold_CVerror_fdrmatrix[which(threshold_CVerror_fdrmatrix[,2]==min(threshold_CVerror_fdrmatrix[,2])),1])[1]
        
        
    }
    
	p4<-pamr.listgenes(fit=p1,data=d1,threshold=threshold_value)
    p4<-as.data.frame(p4)
    
    selected_feature_index<-p4$id
    selected_feature_index<-as.numeric(as.character(selected_feature_index))
    
    discore_matrix<-p4[,-c(1)]
    if(nrow(discore_matrix)>1){
            discore_matrix<-apply(discore_matrix,2,as.numeric)
            abs.discore_matrix<-abs(discore_matrix)
            max.discore<-apply(abs.discore_matrix,1,max)
            
        }else{
            discore_matrix_1<-unlist(discore_matrix)
            discore_matrix<-as.numeric(as.character(discore_matrix_1))
            abs.discore_matrix<-abs(discore_matrix)
            max.discore<-max(abs.discore_matrix)
            
    }
    
    }else{
        p4<-{}
        selected_feature_index<-{}
        max.discore<-{}
    }
	
	
    #return(list("pam_model"=p1,"pam_fdr"=p3,"pam_selected_features"=selected_feature_index))
    #return(selected_feature_index)
    
    return(list("feature.list"=selected_feature_index,"max.discore"=max.discore,"pam_train_model"=p1,"pam_toplist"=p4))
}
