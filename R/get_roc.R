get_roc <-
function(dataA,classlabels,classifier="svm",kname="radial",rocfeatlist=seq(2,10,1),rocfeatincrement=TRUE,testset=NA,testclasslabels=NA,mainlabel=NA,col_lab=NA,legend=TRUE,newdevice=FALSE){
    
    
    d1<-dataA
    rm(dataA)
    
    
    if(newdevice==TRUE){
        
        pdf("ROC.pdf")
    }
    cnames<-colnames(d1)
    cnames[1]<-"mz"
    cnames[2]<-"time"
    colnames(d1)<-cnames
    d1<-as.data.frame(d1)
    
    #print(dim(d1))
    classlabels<-as.data.frame(classlabels)
    testclasslabels<-as.data.frame(testclasslabels)
    class_inf<-classlabels
    
    d2<-t(d1[,-c(1:2)])
    
    if(is.na(testset)==TRUE){
        testset<-d1
        testclasslabels<-classlabels
    }
    
    testset<-t(testset[,-c(1:2)])
    
    featlist<-rocfeatlist
    featincrement<-rocfeatincrement
    
    mz_names<-paste(d1$mz,d1$time,sep="_")
    
    class_vec<-as.character(class_inf[,1])
    
    class_levels<-levels(as.factor(class_vec))
    for(c in 1:length(class_levels)){
        
        classind<-(c-1)
        class_vec<-replace(class_vec,which(class_vec==class_levels[c]),classind)
    }
    
    #class_vec<-replace(class_vec,which(class_vec==class_levels[2]),1)
    
    class_vec<-as.numeric(as.character(class_vec))
    
    
    d3<-cbind(class_vec,d2)
    
    
    colnames(testset)<-as.character(mz_names)
    mz_names<-c("Class",mz_names)
    colnames(d3)<-as.character(mz_names)
    
    
    
    d3<-as.data.frame(d3)


    featlist<-unique(featlist)
   
    mod.lab<-{}
    #  col_vec<-c(heat.colors(256),topo.colors(256))
    
    #col_lab<-col_vec[sample(x=1:length(col_vec),size=length(featlist))] #c("purple","blue","darkgreen", "green","yellow","orange","red")
    
    if(is.na(col_lab)==TRUE){
    if(length(featlist)>1){
    col_lab<-palette(rainbow(length(featlist)))
    }else{
        col_lab<-c("blue")
    }
    }
    
    if(featincrement==TRUE){
        for(n in 1:length(featlist)){
            
           
            
            num_select<-featlist[n]
            
            #print(num_select)
            #print(d3[1:3,1:5])
            if(num_select>dim(d1)[1]){
                
                num_select<-dim(d1)[1]
            }
            
            if(is.na(mainlabel)==TRUE){
                
                
                if(classifier=="logitreg"){
                    mainlab<-"ROC curves using logistic regression"
                }else{
                    
                    mainlab<-"ROC curves using SVM"
                }
                
            }else{
                
                mainlab=mainlabel
            }
            
            
            if(classifier=="logitreg"){
                
                d4<-as.data.frame(d3[,c(1:num_select)])
                
                model1 <- glm(d4$Class~., data=d4, family=binomial)
                
                testset<-as.data.frame(testset)
                #model1 <- glm(d3$Class~., data=d3[,c(1:num_select)], family=binomial)
                
                
                pred<-predict(model1,testset)
                
                pred1 <- ROCR::prediction(pred, testclasslabels)
                
                
                
            }else{
                
                
                
                
                model1 <- svm(as.factor(d3$Class)~., data=d3[,c(1:num_select)], type="C",probability=TRUE,kernel=kname)
                
                
                pred<-predict(model1,testset,probability=TRUE,decision.values=TRUE,type="prob")
                pred1 <- ROCR::prediction(attributes(pred)$probabilities[,2], testclasslabels)
                
            }
            
            
           
            stats1a <- performance(pred1, 'tpr', 'fpr')
            
            roc_res<-cbind(stats1a@x.values[[1]],stats1a@y.values[[1]])
            colnames(roc_res)<-c(stats1a@x.name,stats1a@y.name)
            #  save(roc_res,file="ROCres.rda")
            x1<-seq(0,1,0.01)
            y1<-x1
            p1<-performance(pred1,"auc")
            mod.lab <-c(mod.lab,paste('using top ',(num_select-1),' m/z features: AUC ',round(p1@y.values[[1]],2),sep=""))
            if(n==1){
                plot(stats1a@x.values[[1]], stats1a@y.values[[1]], type='s', ylab=stats1a@y.name, xlab=stats1a@x.name, col=col_lab[n], lty=2, main=mainlab,cex.main=1,lwd=2)
            }else{
                lines(stats1a@x.values[[1]], stats1a@y.values[[1]], type='s', col=col_lab[n], lty=2,lwd=2)
            }
            
        }
        lines(x1, y1, col="black", type="l",lwd=2)
        if(legend==TRUE){
            
            legend('bottomright', c(mod.lab), col=col_lab, lwd=c(2,1), lty=1:3, cex=0.8, bty='n')
        }
    }else{
        
         num_select<-length(featlist)
         
        if(classifier=="logit"){
            model1 <- glm(d3$Class~., data=d3[,c(featlist)], family=binomial)
            
            
            #pred1 <- prediction(fitted(model1), d3$Class)
            pred<-predict(model1,testset)
            
            pred1 <- ROCR::prediction(pred, testclasslabels)
            
            
        }else{
            
            #model1 <- svm(d3$Class~., data=d3[,c(featlist)], type="C",kernel=kname)
            
            #pred1 <- prediction(fitted(model1), d3$Class)
            model1 <- svm(as.factor(d3$Class)~., data=d3[,c(featlist)], type="C",probability=TRUE,kernel=kname)
            
            pred<-predict(model1,testset,probability=TRUE,decision.values=TRUE,type="prob")
            pred1 <- ROCR::prediction(attributes(pred)$probabilities[,2], testclasslabels)
            
            
            
        }
        stats1a <- performance(pred1, 'tpr', 'fpr')
        roc_res<-cbind(stats1a@x.values[[1]],stats1a@y.values[[1]])
        colnames(roc_res)<-c(stats1a@x.name,stats1a@y.name)
        # save(roc_res,file="ROCres.rda")
        x1<-seq(0,1,0.01)
        y1<-x1
        p1<-performance(pred1,"auc")
        mod.lab <-c(mod.lab,paste('using selected ',(num_select-1),' m/z features: AUC ',round(p1@y.values[[1]],2),sep=""))
        if(n==1){
            plot(stats1a@x.values[[1]], stats1a@y.values[[1]], type='s', ylab=stats1a@y.name, xlab=stats1a@x.name, col=col_lab[n], lty=2, main=mainlab,cex.main=0.8,lwd=2)
        }else{
            lines(stats1a@x.values[[1]], stats1a@y.values[[1]], type='s', col=col_lab[n], lty=2,lwd=2)
        }
        
    }
    lines(x1, y1, col="black", type="l",lwd=2)
    if(legend==TRUE){
        legend('bottomright', c(mod.lab), col=col_lab, lwd=c(2,1), lty=1:3, cex=0.8, bty='n')

    }
    
    if(newdevice==TRUE){
            dev.off()
    }
}
