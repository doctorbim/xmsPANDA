svmrfeFeatureRanking <-
function(x,y,svmkernel,kfold,pred.eval.method){
    n = ncol(x)
    survivingFeaturesIndexes = seq(1:n)
    featureRankedList = vector(length=n)
    rankedFeatureIndex = n
    best_acc<-0
    best_subset<-{}
    while(length(survivingFeaturesIndexes)>0){
        #train the support vector machine
        svmModel = svm(x[, survivingFeaturesIndexes], y, cachesize=500, scale=F, type="C-classification", kernel=svmkernel)
        
        subdata=x[, survivingFeaturesIndexes]
        #print(head(subdata))
        #svm_model<-try(svm_cv(v=kfold,x=x[, survivingFeaturesIndexes], y=y,kname=svmkernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
        #svm_model<-try(svm_cv(v=kfold,x=subdata,y=y,kname=svm_kernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
         svm_model<-try(svm_cv(v=kfold,x=subdata,y=y,kname=svmkernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
         
         if(is(svm_model,"try-error")){
             
             svm_acc<-0
         }else{
             
             svm_acc<-svm_model$avg_acc
         }
             
         
         kfold_acc_rand<-{}
         #for(p1 in 1:100){
         
         numcores<-round(detectCores()*0.6)
         
         set.seed(999)
         seed_vec<-round(runif(n=10,min=100,max=10000))
        kfold_acc_rand<-lapply(1:10,function(p1){
            set.seed(seed_vec[p1])
             sample_ind<-sample(1:dim(y)[1],size=dim(y)[1])
             classlabels_rand<-y[sample_ind,]
             classlabels_rand<-as.data.frame(classlabels_rand)
             svm_model_rand<-try(svm_cv(v=kfold,x=subdata,y=classlabels_rand,kname=svmkernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
             
             #svm_model_rand<-svm_cv(v=kfold,x=subdata,y=classlabels_rand,kname=svmkernel,errortype=pred.eval.method,conflevel=95)
             
             if(is(svm_model_rand,"try-error")){
                 # kfold_acc_rand<-c(kfold_acc_rand,NA)
                 return(0)
             }else{
                 #kfold_acc_rand<-c(kfold_acc_rand,svm_model$avg_acc)
                 return(svm_model_rand$avg_acc)
             }
         })
         kfold_acc_rand<-unlist(kfold_acc_rand)
         
         kfold_acc_rand<-mean(kfold_acc_rand,na.rm=TRUE)
         
         overall_acc<-0.5*svm_acc+0.5*(svm_acc-kfold_acc_rand)
         
        if(is(svm_model,"try-error")){
            
            svm_model<-0
            
        }else{
            
            if(overall_acc>best_acc){
                
                best_acc<-overall_acc
                best_subset<-survivingFeaturesIndexes
                
            }
            
            
            
        }
        
        #compute the weight vector
        w = t(svmModel$coefs)%*%svmModel$SV
        #compute ranking criteria
        rankingCriteria = w * w
        #rank the features
        ranking = sort(rankingCriteria, index.return = TRUE)$ix
        #update feature ranked list
        featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]
        rankedFeatureIndex = rankedFeatureIndex - 1
        #eliminate the feature with smallest ranking criterion
        (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
        
        
    } 
    return (list("featureRankedList"=featureRankedList,"best_subset"=best_subset))
}
