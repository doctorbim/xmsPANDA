get_classification.accuracy.child <-
function(temptrain,tempclass,kernelname="radial",errortype="AUC",classifier="svm",num_samp,temptest=NA,testclass=NA,numfolds=1)
{
    
    v=numfolds
    if(classifier=="SVM" | classifier=="svm"){
        mod_cv <- svm(x=temptrain,y=tempclass, type="C",kernel=kernelname)
       
        
        if(v==num_samp){
            
            predfit<-predict(mod_cv,(temptest))
            
        }else{
            predfit<-predict(mod_cv,(temptest))
        }
        
        
    }else{
        if(classifier=="logitreg" | classifier=="LR"){
            
            Class<-as.numeric(tempclass)-1
            
            dtemp<-cbind(Class,temptrain)
            dtemp<-as.data.frame(dtemp)
            
            save(dtemp,file="dtemp.Rda")
            save(temptest,file="temptest.Rda")
            mod_cv<-glm(as.factor(Class)~.,data=dtemp,family=binomial(logit))
            
            if(v==num_samp){
                
                predfit<-predict(mod_cv,(temptest),type="response")
            }else{
                
                
                predfit<-predict(mod_cv,temptest,type="response")
            }
            
            
            predfit <- ifelse(predfit > 0.5,1,0)
            
            testclass<-as.numeric(testclass)-1
            
        }else{
            
            
            if(classifier=="RF" | classifier=="randomforest" | classifier=="rf"){
                
                
                Class<-tempclass
                
                mod_cv<-randomForest(x=(temptrain),y=as.factor(tempclass),ntree=1000)
                
                # predfit<-predict(mod_cv,temptest)
                if(v==num_samp){
                    
                    predfit<-predict(mod_cv,(temptest))
                }else{
                    predfit<-predict(mod_cv,(temptest))
                }
                
                
            }else{
                
                if(classifier=="NaiveBayes" | classifier=="naivebayes"){
                    
                    
                    Class<-tempclass
                    
                    mod_cv<-naiveBayes(x=(temptrain),y=as.factor(tempclass))
                    
                    # predfit<-predict(mod_cv,temptest)
                    if(v==num_samp){
                        
                        predfit<-predict(mod_cv,(temptest))
                    }else{
                        predfit<-predict(mod_cv,(temptest))
                    }
                    
                    
                }else{
                    
                     if(classifier=="plsda" | classifier=="pls"){
                         set.seed(123)
                         opt_comp<-pls.lda.cv(Xtrain=temptrain, Ytrain=tempclass,  ncomp=c(1:10), nruncv=v, alpha=2/3, priors=NULL)
                         
                         predfit<-pls.lda(Xtrain=temptrain,Ytrain=tempclass,ncomp=opt_comp,nruncv=v,Xtest=temptest)
                         predfit<-predfit$predclass
                         
                        

                     }
                }
                
                
            }
            
        }
        
    }
    
    
    svm_table<-table(predfit,testclass)
    
    class_names<-rownames(svm_table)
    beracc<-{}
    auc_acc<-{}
    totacc<-length(which(predfit==testclass))/length(testclass)
    for(c in 1:dim(svm_table)[1]){
        testclass_ind<-which(testclass==class_names[c])
        beracc<-c(beracc,length(which(predfit[testclass_ind]==testclass[testclass_ind]))/length(testclass_ind))
        
    }
    if(errortype=="AUC"){
        testclass<-as.vector(testclass)
        y1<-as.vector(y[,1])
        pred_acc<-multiclass.roc(testclass,as.numeric(predfit),levels=levels(as.factor(y1)))
        pred_acc_orig<-pred_acc$auc[1]
        auc_acc<-c(auc_acc,pred_acc_orig)
    }
    
    
    
    beracc<-mean(beracc,na.rm=TRUE)
    
    if(errortype=="CV" | errortype=="total"){
        svm_acc<-(totacc*100)
    }else{
        if(errortype=="AUC"){
            svm_acc<-(auc_acc*100)
        }else{
            svm_acc<-(beracc*100)
        }
    }
    return(list(classification_acc=svm_acc,classification_model=mod_cv,confusion_matrix=svm_table))
}
