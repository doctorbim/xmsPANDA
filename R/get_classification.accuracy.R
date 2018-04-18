get_classification.accuracy <-
function(kfold,featuretable,classlabels,kernelname="radial",errortype="AUC",conflevel=95,classifier="svm",seednum=555,testfeaturetable=NA,testclasslabels=NA){
    
    
    v=kfold
    kname=kernelname
    
    options(warn=-1)
    mz_time<-paste(round(featuretable[,1],5),"_",round(featuretable[,2],2),sep="")
    
    x=t(featuretable[,-c(1:2)])
    
    if(is.na(testfeaturetable)==FALSE){
        
        if(nrow(featuretable)!=nrow(testfeaturetable)){
            
            stop("Number of features/variables should be same in the train and test sets.")
        }
        
        print("Note: the order of features/variables should be same in the train and test sets.")
        
        testfeaturetable<-t(testfeaturetable[,-c(1:2)])
        colnames(testfeaturetable)<-mz_time
        if(is.na(testclasslabels)==FALSE){
            
            if(is.vector(testclasslabels)==TRUE){
                testclasslabels=as.data.frame(testclasslabels)
            }else{
                testclasslabels=as.data.frame(testclasslabels)
                
                if(dim(testclasslabels)[2]>1){
                    testclasslabels<-testclasslabels[,2]
                }else{
                    testclasslabels<-testclasslabels[,1]
                }
            }
            
        }
    }
    
    colnames(x)<-mz_time
    if(is.vector(classlabels)==TRUE){
        y=as.data.frame(classlabels)
    }else{
        classlabels=as.data.frame(classlabels)
        y=classlabels
        if(dim(y)[2]>1){
            y=classlabels[,2]
        }else{
            y=classlabels[,1]
        }
    }
    y=as.data.frame(y)
    
    num_samp=dim(x)[1]
    
    if(length(num_samp)<1){
        num_samp<-length(x)
        x<-as.data.frame(x)
    }
    y<-as.data.frame(y)
    x<-as.data.frame(x)
    
    num_datasets= floor(num_samp)
    n1<-floor(num_samp/v)
    n2<-num_samp-n1*v
    n3<-v-n2
    
    ind<-rep(c(n1,n1+1),c(n3,n2))
    ind<-diffinv(ind)
    min_err=1
    best_k=1
    
    set.seed(seednum)
    group<-sample(1:num_samp,num_samp, replace=FALSE)
    
    
 
    itr=0
    #svm_acc <- matrix(0,v)  # we set K=30 before, it can be changed to any number<100.
    svm_acc<-rep(0,v)
    mod_cv_list<-new("list")
    for ( i in 1:v)
    {
        g<-group[(ind[i]+1):ind[i+1]]
        temptest<-x[g,]
        temptrain <-x[-g,]
        tempclass <-y[-g,]
        testclass<-y[g,]
        
        cv_res<-get_classification.accuracy.child(temptrain=temptrain,tempclass=tempclass,kernelname=kernelname,errortype=errortype,classifier=classifier,num_samp=num_samp,temptest=temptest,testclass=testclass,numfolds=v)

        svm_acc[i]<-cv_res$classification_acc
        
    }
    avg_acc <-mean(svm_acc,na.rm=TRUE)
    sd_acc<-sd(svm_acc,na.rm=TRUE)
    
    ##Get confidence interval
    probval<-(1-(conflevel*0.01))/2
    probval<-1-probval
    #print(probval)
    error <- qnorm(probval)*sd_acc/sqrt(length(svm_acc))
    avg_acc<-round(avg_acc,2)
    leftconfint<-avg_acc-error
    rightconfint<-avg_acc+error
    test_acc<-NA
    test_confusion_matrix<-NA
    
    leftconfint<-round(leftconfint,2)
    rightconfint<-round(rightconfint,2)
    print(paste("Training set ", kfold,"-fold CV ",errortype," classification accuracy (%):",avg_acc,sep=""))
    print(paste("Training set ", kfold,"-fold CV ",errortype," classification accuracy ", conflevel,"% confidence interval:(",leftconfint,",",rightconfint,")",sep=""))
    
    
    x<-as.data.frame(x)
    y<-y[,1]
    train_res<-get_classification.accuracy.child(temptrain=x,tempclass=y,kernelname=kernelname,errortype=errortype,classifier=classifier,num_samp=num_samp,temptest=x,testclass=y)
    print(paste("Training set ",errortype," classification accuracy (%):",train_res$classification_acc,sep=""))
    
    mod_cv<-train_res$classification_model
   
    if(is.na(testfeaturetable)==FALSE){
        
        if(is.na(testclasslabels)==FALSE){
            
            testfeaturetable<-as.data.frame(testfeaturetable)
            
              test_res<-get_classification.accuracy.child(temptrain=x,tempclass=y,kernelname=kernelname,errortype=errortype,classifier=classifier,num_samp=num_samp,temptest=testfeaturetable,testclass=testclasslabels)
              
            test_acc<-test_res$classification_acc
            test_acc<-round(test_acc,2)
            test_pred_table<-test_res$confusion_matrix
            
            print(paste("Test set ", errortype," classification accuracy (%):",test_acc,sep=""))
            print("Test set confusion matrix")
            print(test_pred_table)
            
            test_confusion_matrix<-test_pred_table
            
        }
    }
    
    options(warn=0)
    
    if(classifier=="logitreg" | classifier=="LR"){
        
            Class<-y
            
            dtemp<-cbind(Class,x)
            dtemp<-as.data.frame(dtemp)
            
            mod_cv_all<-glm(as.factor(Class)~.,data=dtemp,family=binomial(logit))
            s1<-summary(mod_cv_all)
        
         return(list(avg.train.cv.acc=avg_acc,sd.train.cv.acc=sd_acc, train.cv.acc.each.fold=svm_acc,glm_fit=s1$coefficients,train.cv.acc.confint=c(leftconfint,rightconfint),test.acc=test_acc,test.confusion.matrix=test_confusion_matrix))
    }else{
    return(list(avg.train.cv.acc=avg_acc,sd.train.cv.acc=sd_acc,train.cv.acc.each.fold=svm_acc,train.cv.acc.confint=c(leftconfint,rightconfint),test.acc=test_acc,test.confusion.matrix=test_confusion_matrix))
    }
    
    
}
