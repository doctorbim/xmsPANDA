svmrfeFeatureRankingForMulticlass <-
function(x,y,svmkernel,kfold,pred.eval.method,analysismode="classification",num_nodes=2){
    n = ncol(x)
    survivingFeaturesIndexes = seq(1:n)
    featureRankedList = vector(length=n)
    rankedFeatureIndex = n
    overall_acc_vec<-{}
    if(analysismode=="classification"){
    best_acc<-0
    }else{
        
        best_acc<-9999999999999999
    }
    best_subset<-{}
    while(length(survivingFeaturesIndexes)>0){
        #train the support vector machine
        if(analysismode=="classification"){
            svmModel = svm(x[, survivingFeaturesIndexes], y, cachesize=500, scale=F, type="C-classification", kernel=svmkernel)
            
            subdata=x[, survivingFeaturesIndexes]
            #print(head(subdata))
            #svm_model<-try(svm_cv(v=kfold,x=x[, survivingFeaturesIndexes], y=y,kname=svmkernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
            #svm_model<-try(svm_cv(v=kfold,x=subdata,y=y,kname=svm_kernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
            set.seed(999)
            svm_model<-try(svm_cv(v=kfold,x=subdata,y=y,kname=svmkernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
            
            if(is(svm_model,"try-error")){
                
                svm_acc<-0
            }else{
                
                svm_acc<-svm_model$avg_acc
            }
            
            
            
            
            overall_acc<-svm_acc
            
            if(is(svm_model,"try-error")){
                
                svm_model<-0
                
            }else{
                
                if(overall_acc>=best_acc){
                    
                    best_acc<-overall_acc
                    best_subset<-survivingFeaturesIndexes
                    
                }
                
                
                
            }
            
        }else{
            y<-as.numeric(y)
            set.seed(999)
            svmModel = svm(x[, survivingFeaturesIndexes], y, cachesize=500, scale=F, type="eps-regression", kernel=svmkernel,cross=kold)
            #save(svmModel,file="svmModel.Rda")
            overall_acc<-svmModel$tot.MSE
            
            overall_acc_vec<-c(overall_acc_vec,overall_acc)
            
            if(overall_acc<=best_acc){
                
                best_acc<-overall_acc
                best_subset<-survivingFeaturesIndexes
                
            }
        }
        
        
        #compute the weight vector
        multiclassWeights = svm.weights(svmModel)
        #compute ranking criteria
        multiclassWeights = multiclassWeights * multiclassWeights
        rankingCriteria = 0
        for(i in 1:ncol(multiclassWeights))rankingCriteria[i] = mean(multiclassWeights[,i])
        #rank the features
        (ranking = sort(rankingCriteria, index.return = TRUE)$ix)
        #update feature ranked list
        (featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]])
        rankedFeatureIndex = rankedFeatureIndex - 1
        #eliminate the feature with smallest ranking criterion
        (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
        #cat(length(survivingFeaturesIndexes),"\n")
    }
    
     save(overall_acc_vec,file="overall_acc_vec_multi.Rda")
    return (list("featureRankedList"=featureRankedList,"best_subset"=best_subset))
}
