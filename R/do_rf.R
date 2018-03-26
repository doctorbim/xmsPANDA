do_rf <-
function(X,classlabels, ntrees=1000, analysismode="classification"){
	
	
	
	dataA<-t(X)
	dataA<-as.data.frame(dataA)
	
	
    classlabels<-as.data.frame(classlabels)
	
	dataA<-cbind(classlabels,dataA)
	dataA<-as.data.frame(dataA)
	
	colnames(dataA)<-c("Targetvar",paste("mz",seq(1,dim(X)[1]),sep=""))
	dataA<-as.data.frame(dataA)

	attach(dataA,warn.conflicts=FALSE)
	
     set.seed(290875)
     if(analysismode=="classfication"){
     rf2 <- randomForest(as.factor(Targetvar) ~ .,data=dataA, importance=TRUE,proximity=FALSE,ntree=ntrees,keep.forest=FALSE)
   	}else{
   		
   		rf2 <- randomForest(Targetvar ~ .,data=dataA, importance=TRUE,proximity=FALSE,ntree=ntrees,keep.forest=FALSE)
   		}
    
     varimp_res<-round(importance(rf2,type=2,scale=FALSE), 2)
     
     varimp_res_scaled<-round(importance(rf2,type=2,scale=TRUE), 2)
	rm(dataA)
	
	
         
    
	#return(varimp_res)
	return(list("rf_model"=rf2,"rf_varimp"=varimp_res,"rf_varimp_scaled"=varimp_res_scaled))
}
