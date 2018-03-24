svmrfeFeatureRankingForMulticlass <-
function(x,y,svmkernel,kfold,pred.eval.method){
    n = ncol(x)
    survivingFeaturesIndexes = seq(1:n)
    featureRankedList = vector(length=n)
    rankedFeatureIndex = n
    best_acc<-0
    best_subset<-{}
    while(length(survivingFeaturesIndexes)>0){
        #train the support vector machine
        svmModel = svm(x[, survivingFeaturesIndexes], y, cachesize=500,scale=F, type="C-classification", kernel=svmkernel)
        
        svm_model<-try(svm_cv(v=kfold,x=x[, survivingFeaturesIndexes], y=y,kname=svmkernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
        
        if(is(svm_model,"try-error")){
            
            svm_model<-0
            
        }else{
            
            if(svm_model$avg_acc>best_acc){
                
                best_acc<-svm_model$avg_acc
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
    return (list("featureRankedList"=featureRankedList,"best_subset"=best_subset))
}
