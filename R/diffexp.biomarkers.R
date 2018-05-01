diffexp.biomarkers <-
function(X=NA,Y=NA,feature_table_file=NA,class_labels_file=NA,feat.sel.methods=c("rfe","rf","limma","lasso","elasticnet","wilcox.test","pls","spls","ospls","opls","f.test","t.test"),num.var.sel=10,iter.quantile.thresh=0.1,split.train.test=FALSE,train.pct=0.7,outloc,kfold=10,pca.ellipse=TRUE,Xtest=NA,Ytest=NA,rsdthresh=1,pls_vip_thresh=2,seedvalue=27611,learningsetmethod="CV",confounder.matrix=NA,confounderfdrmethod="BH",confounderfdrthresh=0.05,num.methods.sel=2,globalcor=FALSE,cor.method="spearman",networktype="complete",abs.cor.thresh=0.4,cor.fdrthresh=0.05,max.cor.num=100,net_node_colors=c("green","red"), net_legend=TRUE,niter=10,output.device.type="pdf",heatmap.col.opt="RdBu",boxplot.col.opt=c("grey57"),sample.col.opt="rainbow",mz.thresh=5,time.thresh=10,svm_kernel="radial",good_feats_index=NA,pvalue.thresh=0.05,plots.width=8,plots.height=8,plots.res=600, plots.type="cairo")
{
    
    
   
    errortype="BER"
    if(is.na(X)==TRUE){
    X<-read.table(feature_table_file,sep="\t",header=TRUE)
    }
    
    if(is.na(Y)==TRUE){
    Y<-read.table(class_labels_file,sep="\t",header=TRUE)
    }
    
    
    cl<-makeSOCKcluster(2)
    
    clusterExport(cl,"do_rsd")
    
    feat_rsds<-parApply(cl,X[,-c(1:2)],1,do_rsd)
    
    stopCluster(cl)
    
    abs_feat_rsds<-abs(feat_rsds)
    
    good_metabs<-which(abs_feat_rsds>rsdthresh)
    if(length(good_metabs)>0){
        
        X<-X[good_metabs,]
    }else{
        
        stop("No features meet the defined rsd threshold.")
    }
    
    mzrt<-X[,c(1:2)]
    
    sampnames<-colnames(X[,-c(1:2)])
    
    
    Xorig<-X
    
    if(num.var.sel>dim(Xorig)[1]){
        
        num.var.sel<-dim(Xorig)[1]
    }
    
    
    #Transpose X and Xtrain; rows are samples; columns are metabs
    X<-t(X[,-c(1:2)])
    
   dir.create(outloc)
   setwd(outloc)
   
   dir.create("Figures")
   
   fname<-"InputParameters.csv"
   #sink(fname)
   #   c1<-{}
   #c1<-rbind("X","Y","feature_table_file","class_labels_file","feat.sel.methods","num.var.sel","iter.quantile.thresh","split.train.test","train.pct","outloc","kfold","pca.ellipse","Xtest","Ytest","rsdthresh","pls_vip_thresh","seedvalue","learningsetmethod","confounder.matrix","confounderfdrmethod","confounderfdrthresh","num.methods.sel","globalcor","cor.method","networktype","abs.cor.thresh","cor.fdrthresh","max.cor.num","net_node_colors","net_legend","niter","output.device.type","heatmap.col.opt","boxplot.col.opt","sample.col.opt","mz.thresh","time.thresh","svm_kernel")
   #c2<-{}
   
   #feat.sel.methods_str<-paste(feat.sel.methods,sep=";")
   
   #feat.sel.methods_str<-paste(feat.sel.methods,collapse=";")
   
   #c2<-rbind(X,Y,feature_table_file,class_labels_file,feat.sel.methods,num.var.sel,iter.quantile.thresh,split.train.test,train.pct,outloc,kfold,pca.ellipse,Xtest,Ytest,rsdthresh,pls_vip_thresh,seedvalue,learningsetmethod,confounder.matrix,confounderfdrmethod,confounderfdrthresh,num.methods.sel,globalcor,cor.method,networktype,abs.cor.thresh,cor.fdrthresh,max.cor.num,net_node_colors,net_legend,niter,output.device.type,heatmap.col.opt,boxplot.col.opt,sample.col.opt,mz.thresh,time.thresh,svm_kernel)
   # print(dim(c1))
   #print(dim(c2))
   # c3<-cbind(c1,c2)
   #colnames(c3)<-c("Parameter","Value")
   #write.csv(c3,file="InputParameters.csv",row.names=FALSE)
    
    
    nsamp<-dim(X)[1]
    
    Xtrain<-X
    Ytrain_mat<-Y
    Ytest_mat<-NA
    
    if(split.train.test==TRUE){
        
        if(nsamp<30){
            
            print("N is too small for train/test analysis. Continuing using all data as training set.")
            
            split.train.test=FALSE
            Xtrain<-X
            Y<-as.factor(Y[,2])
            
            Ytrain_mat<-Y
            Ytest_mat<-NA
            
        }else{
       
       
       
        numtrain<-round(train.pct*nsamp)
        
        set.seed(seedvalue)
        
        train_test_sets<-GenerateLearningsets(y=Y[,2],method="MCCV",ntrain=numtrain,niter=1,strat=TRUE,fold=kfold)
        
        allindex<-1:nsamp
        set.seed(seedvalue)
        # train_index<-sample(allindex,numtrain)
        
        train_index<-train_test_sets@learnmatrix[1,]
        
        check_index<-which(train_index==0)
        
        if(length(check_index)>0){
            
            train_index<-train_index[-check_index]
        }
        
        test_index<-allindex[-train_index]
        
        Ytrain_mat<-Y[train_index,]
        Ytest_mat<-Y[test_index,]
        
        Xtest<-X[test_index,]
        Ytest<-as.factor(Y[test_index,2])
        
       
        
        #   save(Xtest,file="Xtest.Rda")
        #save(Ytest,file="Ytest.Rda")
        #save(Xorig,file="Xorig.Rda")
        #save(train_index,file="train_index.Rda")
        #save(test_index,file="test_index.Rda")
        
        
        Y<-as.factor(Y[train_index,2])
        Xtrain<-X[train_index,]
        
        print("Dim of train set")
        print(dim(Xtrain))
        
        print("Dim of test set")
        print(dim(Xtest))
        
        
        sampnames<-colnames(Xorig[,c(train_index+2)])
        
        sampnames_test<-colnames(Xorig[,c(test_index+2)])
        
        if(is.na(confounder.matrix)==FALSE){
            
            confounder.matrix<-confounder.matrix[train_index,]
        }
        }
        
    }else{
        
        Xtrain<-X
        Y<-as.factor(Y[,2])
        
        if(is.na(Xtest)==FALSE){
            
            Xtestorig<-Xtest
            
            mzrt_test<-Xtest[,c(1:2)]
            
            #print(head(mzrt))
            #print(head(mzrt_test))
            
            sampnames_test<-colnames(Xtest[,-c(1:2)])
            Xtest<-t(Xtest[,-c(1:2)])
            Ytest<-as.factor(Ytest[,2])
            
            res1<-getVenn(dataA=mzrt,dataB=mzrt_test,name_a="train",name_b="test",time.thresh=time.thresh,mz.thresh=mz.thresh,xMSanalyzer.outloc=getwd())
            
            if(nrow(res1$common)<1){
                stop("No common features found.")
            }else{
                
                if(nrow(res1$common)!=nrow(mzrt)){
                    print("Not all features were common between the train and test sets. Only using the common features for further analysis.")
                    
                }
                
                print("Number of common features:")
                print(nrow(res1$common))
                
                #print(head(res1$common))
                
                Xtrain<-Xtrain[,unique(res1$common$index.A)]
                
                #matching_train_data<-matching_train_data[order(matching_train_data$mz),]
                
                Xtest<-Xtest[,unique(res1$common$index.B)]
                
            }
            
            print("Dim of train set")
            print(dim(Xtrain))
            
            print("Dim of test set")
            print(dim(Xtest))
        }
        
    }
    
    X<-(Xtrain)
    
       svm_acc<-{}
        class_levels<-levels(as.factor(Y))
        print(class_levels)
        
        if(learningsetmethod=="bootstrap"){
            
            niter=niter
        }else{
            niter=niter
        }
        
    if(is.na(Xtest)==FALSE)
    {
        
        if(ncol(Xtrain)!=ncol(Xtest)){
            
            stop("The train and test sets should have same number of variables.")
            
        }
    }
  
        
    set.seed(seedvalue)
    fiveCV10iter<-GenerateLearningsets(y=Y,method=learningsetmethod,fold=kfold,niter=niter,strat=TRUE)
    
    feat_sel_matrix<-matrix(nrow=dim(X)[2],ncol=length(feat.sel.methods),0)
     if(is.na(good_feats_index)==TRUE){
         
    if(is.na(feat.sel.methods)==FALSE){
    for(m in 1:length(feat.sel.methods)){
        
        method=feat.sel.methods[m]
    
    if(method=="pls" || method=="spls" || method=="opls" || method=="ospls"){
        
        v1<-matrix(ncol=dim(fiveCV10iter@learnmatrix)[1],nrow=dim(X)[2],0)
        
        for(i in 1:dim(fiveCV10iter@learnmatrix)[1]){
        
        temptrain<-t(X[fiveCV10iter@learnmatrix[i,],])
        tempclass<-Y[fiveCV10iter@learnmatrix[i,]]
        
        classlabels_sub<-cbind(paste("S",rep(1,length(tempclass)),sep=""),tempclass)
        
        classlabels_sub<-as.data.frame(classlabels_sub)
        
        if(method=="spls"){
            sparseselect=TRUE
        }else{
            sparseselect=FALSE
        }
    
                   numcomp<-5
                   
                   set.seed(123)
                   opt_comp<-pls.lda.cv(Xtrain=X, Ytrain=Y,  ncomp=c(1:numcomp), nruncv=kfold, alpha=2/3, priors=NULL)
                   
                   if(method=="ospls"){
                       
                       Ytemp<-as.numeric(Y)
                       leukemia.pls <- plsr(Ytemp ~ X, ncomp = opt_comp, validation = "LOO")
                       ww <- leukemia.pls$loading.weights[,1]
                       pp <- leukemia.pls$loadings[,1]
                       w.ortho <- pp - crossprod(ww, pp)/crossprod(ww) * ww
                       t.ortho <- X %*% w.ortho
                       
                       p.ortho <- crossprod(X, t.ortho) / c(crossprod(t.ortho))
                       Xcorr <- X - tcrossprod(t.ortho, p.ortho)
                       
                   
                       
                       X<-Xcorr
                       method="spls"
                   }
                  
                  if(method=="opls"){
                      
                      Ytemp<-as.numeric(Y)
                      leukemia.pls <- plsr(Ytemp ~ X, ncomp = opt_comp, validation = "LOO")
                      ww <- leukemia.pls$loading.weights[,1]
                      pp <- leukemia.pls$loadings[,1]
                      w.ortho <- pp - crossprod(ww, pp)/crossprod(ww) * ww
                      t.ortho <- X %*% w.ortho
                      
                      p.ortho <- crossprod(X, t.ortho) / c(crossprod(t.ortho))
                      Xcorr <- X - tcrossprod(t.ortho, p.ortho)
                      
                      
                      
                      X<-Xcorr
                      method="pls"
                  }


                   
                   
                   #opt_comp<-plsres1$opt_comp
                   max_comp_sel<-opt_comp
                   if(method=="spls"){
                       
                       keep_X_vec=rep(num.var.sel,opt_comp)
                       
                       linn.pls <- splsda(X, Y,ncomp=opt_comp,keepX=keep_X_vec)
                       
                       linn.vip<-linn.pls$loadings$X

                       
                       if(opt_comp>1){
                           
                           #abs
                           vip_res1<-abs(linn.vip)
                           
                           if(max_comp_sel>1){
                               vip_res1<-apply(vip_res1,1,max)
                               
                           }else{
                               
                               vip_res1<-vip_res1[,c(1)]
                           }
                       }else{
                           
                           vip_res1<-abs(linn.vip)
                       }
                       
                       pls_vip<-vip_res1 #(plsres1$vip_res)
                       
                       
                       #based on loadings for sPLS
                       #feat_sel_matrix[which(pls_vip!=0),i]<-1 #pls_vip!=0 & rand_pls_sel_fdr<fdrthresh
                       varindex<-which(pls_vip!=0)
                       if(length(varindex)>0){
                           v1[varindex,i]<-1
                       }

                       
                   }else{
                       
                       
                       
                       linn.pls <- plsda(X, Y,ncomp=opt_comp)
                       
                       linn.vip<-vip(linn.pls)
                       
                       {
                           if(opt_comp>1){
                               vip_res1<-(linn.vip)
                               if(max_comp_sel>1){
                                   vip_res1<-apply(vip_res1,1,mean)
                               }else{
                                   
                                   vip_res1<-vip_res1[,c(1)]
                               }
                           }else{
                               
                               vip_res1<-linn.vip
                           }
                           
                       }
                       
                       #vip_res1<-plsres1$vip_res
                       pls_vip<-vip_res1
                       
                       pls_vip_order<-pls_vip[order(pls_vip,decreasing=TRUE)]
                       
                       pls_vip_thresh<-min(pls_vip_order[1:num.var.sel])[1]
                       
                       #print(pls_vip_thresh)
                       #pls
                       varindex<-which(pls_vip>=pls_vip_thresh)
                       if(length(varindex)>0){
                       v1[varindex,i]<-1
                       }
                    }
        }
    }else{
    
    #set.seed(27611)
    set.seed(seedvalue)
    if(method=="rf"){
    g1<-GeneSelection(X=X,y=Y,learningsets=fiveCV10iter,method=method,trace=FALSE,seed = 100)
    }else{
        
        g1<-GeneSelection(X=X,y=Y,learningsets=fiveCV10iter,method=method,trace=FALSE)
    }
    gmatrix<-{}

    
    rank_matrix<-g1@rankings[[1]]
    
    v1<-matrix(nrow=dim(X)[2],ncol=dim(rank_matrix)[1],0)
    
 

        for(i in 1:dim(rank_matrix)[1]){
        
        varindex<-{}
        varindex1<-toplist(g1,iter=i,k=num.var.sel,show=FALSE)
        
        if(length(g1@rankings)>1){
        for(gindex in 1:length(g1@rankings)){
            
            varindex<-c(varindex,varindex1[[gindex]][,1])
        
        }
        }else{
            varindex<-varindex1[,1]
        }
        
        varindex<-unique(varindex)
        
        
        v1[varindex,i]<-1
     
        }
    }
    
    
    
    

    #hist(svm_acc,main="Inner test set accuracy distribution",col="brown")

#iter.quantile.thresh: means that value is 1 in (1-iter.quantile.thresh)% or more sets;
stability_measure<-apply(v1,1,function(x){quantile(x,iter.quantile.thresh)})
    
    stability_matrix<-stability_measure
    if(m==1){
    stability_matrix_1<-cbind(mzrt,stability_measure)
    }else{
        
        stability_matrix_1<-cbind(stability_matrix_1,stability_measure)
    }


    feat_sel_matrix[which(stability_matrix==1),m]<-1
   
    }
    
    
    
    }else{
        
        feat_sel_matrix<-matrix(nrow=dim(X)[1],ncol=length(feat.sel.methods),1)
        
    }
    
    
    if(length(feat.sel.methods)>1){
        feat_sel_matrix<-apply(feat_sel_matrix,1,sum)
        
        
    }
    
    if(num.methods.sel>length(feat.sel.methods)){
        
        num.methods.sel=length(feat.sel.methods)
    }

    
    colnames(stability_matrix_1)<-c("mz","time",feat.sel.methods)
    
    pdf("Results.pdf")
   
    #hist(stability_measure,main="Stability measure distribution",col="brown")
    write.table(stability_matrix_1,file="stability_matrix.txt",sep="\t",row.names=FALSE)
   
    good_feats_index<-which(feat_sel_matrix>=num.methods.sel)
   
   if(is.na(pvalue.thresh)==FALSE){
       
      
       
       numcores<-num_nodes #round(detectCores()*0.5)
       
       cl <- parallel::makeCluster(getOption("cl.cores", numcores))
       
       clusterExport(cl,"diffexponewayanova",envir = .GlobalEnv)
       
       clusterExport(cl,"anova",envir = .GlobalEnv)
       
       
       clusterExport(cl,"TukeyHSD",envir = .GlobalEnv)
       
       clusterExport(cl,"aov",envir = .GlobalEnv)
       
       
       #res1<-apply(data_m_fc,1,function(x){
       res1<-parApply(cl,X,2,function(x,classlabels_response_mat){
           xvec<-x
           
           
           data_mat_anova<-cbind(xvec,classlabels_response_mat)
           
           data_mat_anova<-as.data.frame(data_mat_anova)
           cnames<-colnames(data_mat_anova)
           
           cnames[1]<-"Response"
           
           colnames(data_mat_anova)<-c("Response","Factor1")
           
           #print(data_mat_anova)
           
           data_mat_anova$Factor1<-as.factor(data_mat_anova$Factor1)
           
           anova_res<-diffexponewayanova(dataA=data_mat_anova)
           
           
           
           return(anova_res)
       },Y)
       
       stopCluster(cl)
       main_pval_mat<-{}
       
       posthoc_pval_mat<-{}
       pvalues<-{}
       
       
       #print(head(res1))
       
       for(i in 1:length(res1)){
           
           
           main_pval_mat<-rbind(main_pval_mat,res1[[i]]$mainpvalues)
           pvalues<-c(pvalues,res1[[i]]$mainpvalues[1])
           posthoc_pval_mat<-rbind(posthoc_pval_mat,res1[[i]]$posthocfactor1)
           
       }
       
       pvalues<-unlist(pvalues)
       
       good_feats_index<-which(feat_sel_matrix>=num.methods.sel & pvalues<pvalue.thresh)
       if(length(good_feats_index)<1){
           stop("No features selected.")
       }
       
   }
    good_feats<-Xtrain[,good_feats_index]
    
   
   
    mzrt_sub<-mzrt[good_feats_index,]
    
    good_feats<-t(good_feats)
    
   
   
    cnames_1<-colnames(good_feats)
    
    
    colnames(good_feats)<-sampnames
    
    
    good_feats<-cbind(mzrt_sub,good_feats)
    
    if(is.na(Xtest)==FALSE){
        good_feats_test<-Xtest[,good_feats_index]
        good_feats_test<-t(good_feats_test)
        colnames(good_feats_test)<-sampnames_test
        good_feats_test<-cbind(mzrt_sub,good_feats_test)
    }else{
        good_feats_test<-NA
        Ytest_mat<-NA
    }
    
   
   
     }else{
         pdf("Results.pdf")
         good_feats<-Xtrain[,good_feats_index]
         
         
         
         mzrt_sub<-mzrt[good_feats_index,]
         
         good_feats<-t(good_feats)
         
         
         
         cnames_1<-colnames(good_feats)
         
         
         colnames(good_feats)<-sampnames
         
         
         good_feats<-cbind(mzrt_sub,good_feats)
         
         if(is.na(Xtest)==FALSE){
             good_feats_test<-Xtest[,good_feats_index]
             good_feats_test<-t(good_feats_test)
             colnames(good_feats_test)<-sampnames_test
             good_feats_test<-cbind(mzrt_sub,good_feats_test)
         }else{
             good_feats_test<-NA
             Ytest_mat<-NA
         }
         
     }
    
    print("Number of features selected")
    print(length(good_feats_index))
    
   
    if(length(good_feats_index)>1){
        
        
   
        X<-X[,good_feats_index]
        
        if(is.na(Xtest)==FALSE){
            
            Xtest<-Xtest[,good_feats_index]
        }
        
        
        
        
  
    #if(classifier=="pls")
    if(length(good_feats_index)>5)
    {
        
        #print("Sts1")
        
        #print(X[1:3,c(1:3,240:243)])
        
        # save(list=ls(),file="t1.Rda")
        tune_plslda <- CMA::tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_ldaCMA, grids = list(comp = 1:5),trace=FALSE)
        
        if(length(class_levels)==2){
            #     tune_plslr <- tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_lrCMA, grids = list(comp = 1:5))
        }
        #print("Sts2")
        tune_plsrf <- CMA::tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_rfCMA, grids = list(comp = 1:5),trace=FALSE)
    }
    #if(classifier=="scda")
    {
        
        #print("Sts3")
        tune_scda <- CMA::tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA, grids = list( ),trace=FALSE)
    }
    #if(classifier=="svm")
    {
        #print("Sts4")
        tune_svm <- CMA::tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, grids = list( ), kernel = "radial",trace=FALSE)
    }
    
    if(length(good_feats_index)>2){
    class_plslda <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_ldaCMA, tuneres = tune_plslda,trace=FALSE)
    
    class_plsrf <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_rfCMA, tuneres = tune_plsrf,trace=FALSE)
    }
    
    class_scda <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA, tuneres = tune_scda,trace=FALSE)
    
    class_svm <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, tuneres = tune_svm,kernel = "radial",trace=FALSE,probability=TRUE)
    class_rf <- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = rfCMA,trace=FALSE)
    
    class_nnet<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE)
    
    class_plr<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = plrCMA,trace=FALSE)
    
    save(class_plslda,file="class_plslda.Rda")
    save(class_plsrf,file="class_plsrf.Rda")
    save(class_scda,file="class_scda.Rda")
    save(class_svm,file="class_svm.Rda")
    save(class_rf,file="class_rf.Rda")
    save(class_nnet,file="class_nnet.Rda")
    save(class_plr,file="class_plr.Rda")
    
    
    
    if(length(class_levels)==2){
        
        # class_plslr <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_lrCMA)
        
        # class_lassoplr<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = LassoCMA)
        
        # class_elasticnetplr<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = ElasticNetCMA)
        
        #eval_plslr<-evaluation(class_plslr)
        
        # eval_lassoplr<-evaluation(class_lassoplr)
        
        #eval_elasticnetplr<-evaluation(class_elasticnetplr)
        #measures="auc"
        
        eval_measure="auc"
        #eval_measure="misclassification"
        
    }else{
        
        eval_measure="misclassification"
    }
    
    
    # print("here: evaluation")
    
     if(length(good_feats_index)>2){
    
    # if(eval_measure=="misclassification")
    
    eval_plslda<-evaluation(class_plslda,measure=eval_measure)
    eval_plsrf<-evaluation(class_plsrf,measure=eval_measure)
    eval_scda<-evaluation(class_scda,measure=eval_measure)
    eval_svm<-evaluation(class_svm,measure=eval_measure)
    eval_rf<-evaluation(class_rf,measure=eval_measure)
    eval_nnet<-evaluation(class_nnet,measure=eval_measure)
    eval_plr<-evaluation(class_plr,measure=eval_measure)

        if(eval_measure=="auc")
        {
            eval_plslda@score=1-mean(eval_plslda@score)
            eval_plsrf@score=1-mean(eval_plsrf@score)
            eval_scda@score=1-mean(eval_scda@score)
            eval_svm@score=1-mean(eval_svm@score)
            eval_rf@score=1-mean(eval_rf@score)
            eval_nnet@score=1-mean(eval_nnet@score)
            eval_plr@score=1-mean(eval_plr@score)
            
        }
    }
    
    #save(eval_svm,file="eval_svm.Rda")
    #save(class_svm,file="class_svm.Rda")
    
    
    text2<-paste(dim(v1)[2], " learning sets using training data",sep="")
    if(length(class_levels)==0){
        
         eval_mat1<-cbind(text2,100*(1-mean(eval_plslda@score)),100*(1-mean(eval_plslr@score)),100*(1-mean(eval_plsrf@score)),100*(1-mean(eval_scda@score)),100*(1-mean(eval_svm@score)),100*(1-mean(eval_rf@score)),100*(1-mean(eval_nnet@score)),100*(1-mean(eval_plr@score)),100*(1-mean(eval_lassoplr@score)),100*(1-mean(eval_elasticnetplr@score)))
        
        best_classifier<-which(eval_mat1==max(eval_mat1))
        classifier_names<-c("PLSLDA","PLSLR","PLSRF","SCDA","SVM","RF","NNet","pLR","pLRlasso","pLRelasticnet")
        
        best_classifier_name<-classifier_names[best_classifier]
        
        eval_mat1<-cbind(text2,100*(1-mean(eval_plslda@score)),100*(1-mean(eval_plslr@score)),100*(1-mean(eval_plsrf@score)),100*(1-mean(eval_scda@score)),100*(1-mean(eval_svm@score)),100*(1-mean(eval_rf@score)),100*(1-mean(eval_nnet@score)),100*(1-mean(eval_plr@score)),100*(1-mean(eval_lassoplr@score)),100*(1-mean(eval_elasticnetplr@score)))
        
        colnames(eval_mat1)<-c("Dataset","PLSLDA","PLSLR","PLSRF","SCDA","SVM","RF","NNet","pLR","pLRlasso","pLRelasticnet")
    }else{
        
        if(length(good_feats_index)>5){
        eval_mat1<-cbind(100*(1-mean(eval_plslda@score)),100*(1-mean(eval_plsrf@score)),100*(1-mean(eval_scda@score)),100*(1-mean(eval_svm@score)),100*(1-mean(eval_rf@score)),100*(1-mean(eval_nnet@score)),100*(1-mean(eval_plr@score)))
        }else{
            eval_mat1<-cbind(100*(1-1),100*(1-1),100*(1-mean(eval_scda@score)),100*(1-mean(eval_svm@score)),100*(1-mean(eval_rf@score)),100*(1-mean(eval_nnet@score)),100*(1-mean(eval_plr@score)))
            
        }
        best_classifier<-which(eval_mat1==max(eval_mat1)[1])

         classifier_names<-c("PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR")
         
         best_classifier_name<-classifier_names[best_classifier]
        
        best_innerkfold_acc<-max(eval_mat1)[1]
        
        eval_mat_1<-round(eval_mat1,2)
        eval_mat1<-cbind(text2,eval_mat_1)
        
        colnames(eval_mat1)<-c("Dataset","PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR")
        
    }
    
    if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/Barplot_classifier_comparison_CVaccuracy.png"
        
        png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
    }
    
    
    

       barplot(eval_mat_1,xaxt="n",main="Comparison of CV accuracy using different classifiers\n based on learning sets", xlab="Classifier",ylab="kfold classification accuracy(%)",col=boxplot.col.opt,type="p",ylim=c(0,100))
       axis(side=1,at=seq(1,7),labels=c("PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR"),cex.axis=0.7)
   
   
   
   if(output.device.type!="pdf"){
       
       try(dev.off(),silent=TRUE)
   }
   
   mzrt_1<-paste(round(good_feats[,1],5),round(good_feats[,2],1),sep="_")
    rownames(good_feats)<-mzrt_1
   
   
    Y1<-cbind(sampnames,as.character(Y))
    Y1<-as.data.frame(Y1)
    #print(head(Y1))
    
    if(length(good_feats_index)>3){
        
        if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/HCA_sigfeats.png"
            
            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        }
        
        
       
    try(get_hca(parentoutput_dir=outloc,X=good_feats,Y=Y1,heatmap.col.opt=heatmap.col.opt,cor.method="spearman",is.data.znorm=FALSE,analysismode="classification",
    sample.col.opt="rainbow",plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type="two-way",newdevice=FALSE),silent=TRUE)
    
   
    
    if(output.device.type!="pdf"){
        
        try(dev.off(),silent=TRUE)
    }
    
    
    if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/PCAplots_sigfeats.pdf"
        
        #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        pdf(temp_filename_1)
    }
    
    
    # save(list=ls(),file="debug1.Rda")
    try(get_pcascoredistplots(X=good_feats,Y=Y1,parentoutput_dir=outloc,sample.col.opt=sample.col.opt,alphacol=0.3,col_vec=NA,pairedanalysis=FALSE,pca.cex.val=3,legendlocation="topright",pca.ellipse=FALSE,ellipse.conf.level=0.5,filename="selected",paireddesign=paireddesign),silent=TRUE)

    }
    
  
    if(output.device.type!="pdf"){
        
        try(dev.off(),silent=TRUE)
    }
    
    
    if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/Boxplots_sigfeats.pdf"
        
        #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        pdf(temp_filename_1)
    }
    

    
    
    get_boxplots(X=good_feats,Y=Y1,parentoutput_dir=outloc,boxplot.col.opt=boxplot.col.opt,sample.col.opt=sample.col.opt,alphacol=0.3,newdevice=FALSE,cex=0.6)
    
    
    
   
    
    if(output.device.type!="pdf"){
        
        try(dev.off(),silent=TRUE)
    }
    
    num_levels<-levels(as.factor(Y))
    if(num_levels==2){
    if(split.train.test==TRUE){
        #     X<-as.matrix(X)
        #get_roc(dataA=good_feats,classlabels=Y,classifier="svm",kname=svm_kernel,rocfeatlist=seq(2,10,1),rocfeatincrement=TRUE,testset=NA,testclasslabels=NA,mainlabel="Test set")


    }else{
        #get_roc(dataA=good_feats,classlabels=Y,classifier="svm",kname=svm_kernel,rocfeatlist=seq(2,10,1),rocfeatincrement=TRUE,testset=X,testclasslabels=Y,mainlabel="Using training set based on SVM")
        
    }
    }
    
    
    #if(is.na(Xtest)==FALSE)
    {
        
       
        
        best_kfold_acc<-best_innerkfold_acc
        #outerkfold_acc<-100*(1-mean(testeval_res@score))
        
        permkfold_acc<-{}
        permkfold_acc1<-{}
        
        permkfold_acc2<-{}
        
        permkfold_acc3<-{}
        
        permkfold_acc4<-{}
        
        permkfold_acc5<-{}
        
        permkfold_acc6<-{}
        
        permkfold_acc7<-{}
        Yorig<-Y
        
         text2B<-paste("Permuted ",dim(v1)[2], " learning sets using training data",sep="")
        #if(is.na(Xtest)==FALSE)
        for(p1 in 1:100)
        {
            #set.seed(27611)
            set.seed(seedvalue)
            Y<-Yorig[sample(1:length(Yorig),size=length(Yorig))]
            #set.seed(27611)
            set.seed(seedvalue)
            fiveCV10iter<-GenerateLearningsets(y=Y,method=learningsetmethod,fold=kfold,niter=1,strat=TRUE)
            
            classifier_names<-c("PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR","PLSLR","pLRlasso","pLRelasticnet")
            
            
            if(length(good_feats_index)>5){
            
            class_plslda <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_ldaCMA,trace=FALSE)
                testeval_res1<-evaluation(class_plslda,measure=eval_measure)
            
            
             class_plsrf <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_rfCMA,trace=FALSE)
                testeval_res2<-evaluation(class_plsrf,measure=eval_measure)
                
            }
            
               class_scda <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA,trace=FALSE)
                testeval_res3<-evaluation(class_scda,measure=eval_measure)
            
            
           
              class_svm <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA,kernel = "radial",trace=FALSE,probability=TRUE)
                testeval_res4<-evaluation(class_svm,measure=eval_measure)
            
            
                class_rf <- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = rfCMA,trace=FALSE)
                testeval_res5<-evaluation(class_rf,measure=eval_measure)
           
                class_nnet<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE)
                testeval_res6<-evaluation(class_nnet,measure=eval_measure)
           
                class_plr<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = plrCMA,trace=FALSE)
                testeval_res7<-evaluation(class_plr,measure=eval_measure)
                
                
                if(eval_measure=="auc")
                {
                    eval_plslda@score=1-mean(eval_plslda@score)
                    eval_plsrf@score=1-mean(eval_plsrf@score)
                    eval_scda@score=1-mean(eval_scda@score)
                    eval_svm@score=1-mean(eval_svm@score)
                    eval_rf@score=1-mean(eval_rf@score)
                    eval_nnet@score=1-mean(eval_nnet@score)
                    eval_plr@score=1-mean(eval_plr@score)
                    
                }
                
            
            
                if(length(good_feats_index)>5){
                permkfold_acc1<-c(permkfold_acc1,100*(1-mean(testeval_res1@score)))
                
                permkfold_acc2<-c(permkfold_acc2,100*(1-mean(testeval_res2@score)))
                
                }else{
                    permkfold_acc1<-c(permkfold_acc1,0)
                    permkfold_acc2<-c(permkfold_acc2,0)
                }
         
                permkfold_acc3<-c(permkfold_acc3,100*(1-mean(testeval_res3@score)))
                
                permkfold_acc4<-c(permkfold_acc4,100*(1-mean(testeval_res4@score)))
                
                permkfold_acc5<-c(permkfold_acc5,100*(1-mean(testeval_res5@score)))
                
                permkfold_acc6<-c(permkfold_acc6,100*(1-mean(testeval_res6@score)))
                
                permkfold_acc7<-c(permkfold_acc7,100*(1-mean(testeval_res7@score)))
                
                
            }
        
        permkfold_acc1<-mean(permkfold_acc1)
        permkfold_acc2<-mean(permkfold_acc2)
        permkfold_acc3<-mean(permkfold_acc3)
        permkfold_acc4<-mean(permkfold_acc4)
        permkfold_acc5<-mean(permkfold_acc5)
        permkfold_acc6<-mean(permkfold_acc6)
        permkfold_acc7<-mean(permkfold_acc7)
        
        eval_mat_perm<-cbind(permkfold_acc1, permkfold_acc2, permkfold_acc3, permkfold_acc4, permkfold_acc5, permkfold_acc6, permkfold_acc7)
        eval_mat_perm<-round(eval_mat_perm,2)
         
        eval_mat_perm_final<-cbind(text2B,eval_mat_perm)
        
       
        
        colnames(eval_mat_perm_final)<-c("Dataset","PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR")
        
        eval_mat_actual<-as.data.frame(eval_mat_1)
        eval_mat_perm1<-as.data.frame(eval_mat_perm)
        eval_mat_perm1<-round(eval_mat_perm1,2)
        
        
        colnames(eval_mat_perm1)<-colnames(eval_mat_actual)
        emat1<-rbind(eval_mat_actual,eval_mat_perm1)
        
        emat1<-t(emat1)
        classifier_names<-c("PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR")
        
        rownames(emat1)<-classifier_names
        
        
        eval_mat_actual<-((emat1[,1]-emat1[,2])*0.5)+(0.5*(emat1[,1]))
        
        emat1<-cbind(emat1,eval_mat_actual)
        colnames(emat1)<-c("Actual.accuracy","Permuted.accuracy","Score[((Actual-Permuted)*0.5)+(0.5*Actual)]")
        print(emat1)
        
        
      
      
test_acc<-NA
test_acc_mat<-{}
 Y<-Yorig
if(is.na(Xtest)==FALSE){
    
    test_acc_mat<-{}
    class_levels_vec<-levels(as.factor(Y))
    
    
    learnmatrix <- matrix(seq(1,nrow(X)), nrow = 1)
    fiveCV10iter<-new("learningsets", learnmatrix = learnmatrix, method = "none",ntrain = ncol(learnmatrix), iter = nrow(learnmatrix))
    X<-rbind(X,Xtest)
    Y<-c(Y,Ytest)
    classifier_names<-c("PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR","PLSLR","pLRlasso","pLRelasticnet")
    
    #save(X,file="X1.Rda")
    #save(Y,file="Y1.Rda")
    #save(learnmatrix,file="learnmatrix.Rda")
    
    # save(list=ls(),file="debug3.Rda")
    
    if(length(good_feats_index)>5)
    {
        s1<-ldply(tune_plslda@tuneres,rbind)
        
        s2<-apply(s1,2,mean)
        
        
        t1<-new("list")
        confusion_matrix_list<-new("list")
        t1[[1]]<-s2
        tune_plslda1<-tune_plslda
        tune_plslda1@tuneres<-t1
        save(tune_plslda,file="tune_res.Rda")
        #
        #new("tuneres", fold = tune_plslda@fold, method = tune_plslda@method,
        #tuneres = s2, hypergrid = tune_plslda@hypergrid)
        
        class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_ldaCMA, tuneres = tune_plslda1,trace=FALSE)
        
        b1<-best(tune_plslda1)
        
        #save(b1,file="b1.Rda")
        
        learnmatrix<-as.numeric(learnmatrix)
        
        class_res2<-pls_ldaCMA(X = X, y = Y, learnind = learnmatrix, comp = median(unlist(b1)))
        
        save(class_res2,file="pls_ldaCMA.Rda")
        
        print(ftable(class_res2))
        
        save(confusion_matrix_1,file="confusion_matrix_plslda.Rda")
        confusion_matrix_list[[1]]<-table(class_res2@y,class_res2@yhat)
        
        if(length(class_levels_vec)==2){
            testeval_res_auc<-evaluation(class_res,measure = "auc")
            save(testeval_res_auc,file="testeval_res_auc.Rda")
            print(testeval_res_auc)
            test_auc<-100*(mean(testeval_res_auc@score))
            test_acc_mat<-c(test_acc_mat,test_auc)
            
        }else{
            
            test_acc<-evaluation(class_res)
            
            test_acc<-100*(1-mean(testeval_res_auc@score))
            test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
            
            
        }
        
        
        
        #if(best_classifier_name==classifier_names[2]){
        
        s1<-ldply(tune_plsrf@tuneres,rbind)
        s2<-apply(s1,2,mean)
        #tune_plsrf<-new("tuneres", fold = tune_plsrf@fold, method = tune_plsrf@method,
        #tuneres = s2, hypergrid = tune_plsrf@hypergrid)
        t1<-new("list")
        t1[[1]]<-s2
        
        tune_plsrf1<-tune_plsrf
        tune_plsrf1@tuneres<-t1
        class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_rfCMA, tuneres = tune_plsrf1,trace=FALSE)
        
        b1<-best(tune_plsrf1)
        # save(b1,file="b1_pls_rfCMA.Rda")
        
        learnmatrix<-as.numeric(learnmatrix)
        
        class_res2<-pls_rfCMA(X = X, y = Y, learnind = learnmatrix, comp = median(unlist(b1)))
        
        save(class_res2,file="pls_rfCMA.Rda")
        print(ftable(class_res2))
        
        #confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
        save(confusion_matrix_1,file="confusion_matrix_rfCMA.Rda")
        
        confusion_matrix_list[[2]]<-table(class_res2@y,class_res2@yhat)
        
        if(length(class_levels_vec)==2){
            
            testeval_res_auc<-evaluation(class_res,measure = "auc")
            save(testeval_res_auc,file="testeval_res_auc.Rda")
            
            print(testeval_res_auc)
            
            test_auc<-100*(mean(testeval_res_auc@score))
            
            test_acc_mat<-c(test_acc_mat,test_auc)
            
            
            
        }else{
            
            test_acc<-evaluation(class_res)
            
            test_acc<-100*(1-mean(test_acc@score))
            test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
            
            
        }
        
        
        #}
        
        #}
        
        #if(best_classifier_name==classifier_names[3]){
        
        s1<-ldply(tune_scda@tuneres,rbind)
        s2<-apply(s1,2,mean)
        #tune_scda<-new("tuneres", fold = tune_scda@fold, method = tune_scda@method,
        #  tuneres = s2, hypergrid = tune_scda@hypergrid)
        
        t1<-new("list")
        t1[[1]]<-s2
        
        tune_scda1<-tune_scda
        tune_scda1@tuneres<-t1
        
        class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA, tuneres = tune_scda1,trace=FALSE)
        #scdaCMA(X, y, f, learnind, delta = 0.5, models=FALSE,...)
        
        b1<-best(tune_scda1)
        
        save(b1,file="b1_scda.Rda")
        learnmatrix<-as.numeric(learnmatrix)
        class_res2<-scdaCMA(X = X, y = Y, learnind = learnmatrix, delta = median(unlist(b1)))
        save(class_res2,file="scdaCMA.Rda")
        print(ftable(class_res2))
        
        #confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
        
        save(confusion_matrix_1,file="confusion_matrix_scdaCMA.Rda")
        confusion_matrix_list[[3]]<-table(class_res2@y,class_res2@yhat)
        
        if(length(class_levels_vec)==2){
            
            testeval_res_auc<-evaluation(class_res,measure = "auc")
            save(testeval_res_auc,file="testeval_res_auc.Rda")
            
            print(testeval_res_auc)
            
            test_auc<-100*(mean(testeval_res_auc@score))
            
            test_acc_mat<-c(test_acc_mat,test_auc)
            
            
            
        }else{
            
            test_acc<-evaluation(class_res)
            
            test_acc<-100*(1-mean(test_acc@score))
            test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
            
            
        }
        
        
        #  }
        
        #if(best_classifier_name==classifier_names[4])
        {
            
            s1<-ldply(tune_svm@tuneres,rbind)
            
            s2<-apply(s1,2,mean)
            # tune_svm<-new("tuneres", fold = tune_svm@fold, method = tune_svm@method,
            #tuneres = s2, hypergrid = tune_svm@hypergrid)
            
            t1<-new("list")
            t1[[1]]<-s2
            
            tune_svm1<-tune_svm
            tune_svm1@tuneres<-t1
            
            class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, tuneres=tune_svm1, kernel ="radial",trace=FALSE,probability=TRUE)
            
            
            set.seed(999)
            class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, tuneres=tune_svm1, kernel ="radial",trace=FALSE,probability=TRUE)
            
            
            
            b1<-best(tune_svm1)
            
            save(b1,file="b1_svm.Rda")
            
            learnmatrix<-as.numeric(learnmatrix)
            
            set.seed(999)
            class_res2<-svmCMA(X = X, y = Y, learnind = learnmatrix, delta = median(unlist(b1)),probability=TRUE)
            
            save(class_res2,file="svmCMA.Rda")
            
            print("svm")
            print(ftable(class_res2))
            
            
            confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            save(confusion_matrix_1,file="confusion_matrix_svmCMA.Rda")
            confusion_matrix_list[[4]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
                
                testeval_res_auc<-evaluation(class_res,measure = "auc")
                #save(testeval_res_auc,file="testeval_res_auc.Rda")
                
                print(testeval_res_auc)
                
                test_auc<-100*(mean(testeval_res_auc@score))
                
                test_acc_mat<-c(test_acc_mat,test_auc)
                
                
                
            }else{
                
                test_acc<-evaluation(class_res)
                
                test_acc<-100*(1-mean(test_acc@score))
                test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
                
                
            }
            
        }
        
        #if(best_classifier_name==classifier_names[5])
        {
            class_res <- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = rfCMA,trace=FALSE)
            
            learnmatrix<-as.numeric(learnmatrix)
            
            # class_res2<-svmCMA(X = X, y = Y, learnind = learnmatrix, delta = median(unlist(b1)))
            
            class_res2 <- rfCMA(X =X, y = Y, learnind=learnmatrix, varimp = FALSE)
            
            save(class_res2,file="rfCMA.Rda")
            
            print(ftable(class_res2))
            
            
            confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            save(confusion_matrix_1,file="confusion_matrix_rfCMA.Rda")
            confusion_matrix_list[[5]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
                
                testeval_res_auc<-evaluation(class_res,measure = "auc")
                #save(testeval_res_auc,file="testeval_res_auc.Rda")
                
                print(testeval_res_auc)
                
                test_auc<-100*(mean(testeval_res_auc@score))
                
                test_acc_mat<-c(test_acc_mat,test_auc)
                
                
                
            }else{
                
                test_acc<-evaluation(class_res)
                
                test_acc<-100*(1-mean(test_acc@score))
                test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
                
                
            }
            
        }
        
        #if(best_classifier_name==classifier_names[6])
        {
            class_res<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE)
            
            #nnetresult <- nnetCMA(X=golubX, y=golubY, learnind=learnind, size = 3, decay = 0.01)
            class_res2 <- nnetCMA(X =X, y = Y, learnind=learnmatrix)
            
            save(class_res2,file="nnetCMA.Rda")
            
            print(ftable(class_res2))
            
            
            confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            save(confusion_matrix_1,file="confusion_matrix_nnetCMA.Rda")
            confusion_matrix_list[[6]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
                
                testeval_res_auc<-evaluation(class_res,measure = "auc")
                #save(testeval_res_auc,file="testeval_res_auc.Rda")
                
                print(testeval_res_auc)
                
                test_auc<-100*(mean(testeval_res_auc@score))
                
                test_acc_mat<-c(test_acc_mat,test_auc)
                
                
                
            }else{
                
                test_acc<-evaluation(class_res)
                
                test_acc<-100*(1-mean(test_acc@score))
                test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
                
                
            }
            
        }
        
        #if(best_classifier_name==classifier_names[7])
        {
            class_res<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = plrCMA,trace=FALSE)
            
            class_res2 <- plrCMA(X =X, y = Y, learnind=learnmatrix)
            
            save(class_res2,file="plrCMA.Rda")
            
            print(ftable(class_res2))
            
            
            confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            save(confusion_matrix_1,file="confusion_matrix_plrCMA.Rda")
            confusion_matrix_list[[7]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
                
                testeval_res_auc<-evaluation(class_res,measure = "auc")
                #save(testeval_res_auc,file="testeval_res_auc.Rda")
                
                print(testeval_res_auc)
                
                test_auc<-100*(mean(testeval_res_auc@score))
                
                test_acc_mat<-c(test_acc_mat,test_auc)
                
                
                
            }else{
                
                test_acc<-evaluation(class_res)
                
                test_acc<-100*(1-mean(test_acc@score))
                test_acc_mat<-c(test_acc_mat,test_acc) #cbind(best_kfold_acc,permkfold_acc,test_acc)
                
                
            }
            
        }
        
        
        if(length(class_levels_vec)==2){
            
            
            
            acc_mat<-cbind(emat1,test_acc_mat) #cbind(best_kfold_acc,permkfold_acc,test_acc)
            
            
            colnames(acc_mat)<-c("Training kfold CV Accuracy (AUC)","Training permuted kfold CV accuracy","Score", "Test accuracy(AUC)")
            
            print("Test set evaluation using selected features and the best classifier based on AUC measure")
            print(acc_mat)
            
            #rownames(acc_mat)<-best_classifier_name[1]
            
            write.table(acc_mat,file="Classification_evaluation_results_AUC.txt",sep="\t")
            
            
        }else{
            
            acc_mat<-cbind(emat1,test_acc_mat) #cbind(best_kfold_acc,permkfold_acc,test_acc)
            
            
            colnames(acc_mat)<-c("Training kfold CV Accuracy (Misclassification)","Training permuted kfold CV Accuracy (Misclassification)","Score", "Test accuracy (Misclassification)")
            
            print("Test set evaluation using selected features and the best classifier based on misclassification rate measure")
            print(acc_mat)
            
            #rownames(acc_mat)<-best_classifier_name[1]
            
            write.table(acc_mat,file="Classification_evaluation_results_misclassification.txt",sep="\t")
            
        }
        
        
        
        
        
        
        
        
        
        acc_mat<-acc_mat[,-c(3)]
        
        
        mainlab<-paste("Performance evaluation using classifiers and selected features",sep="")
        
        if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/Barplot_classification_accuracy.png"
            
            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        }
        
        
        barplot(acc_mat,xaxt="n",main=mainlab,col=boxplot.col.opt,ylab="Classification accuracy(%)",ylim=c(0,100))
        
        if(length(class_levels_vec)==2){
            axis(side=1,at=seq(1,4),labels=c("kfold CV","Permuted kfold CV","Test set","AUC"),cex.axis=0.9)
        }else{
            axis(side=1,at=seq(1,3),labels=c("kfold CV","Permuted kfold CV","Test set"),cex.axis=0.9)
        }
 
        
        
        if(output.device.type!="pdf"){
            
            try(dev.off(),silent=TRUE)
        }
        
        
        try(dev.off(),silent=TRUE)
     
    }
        
}else{
    
    
    acc_mat<-acc_mat[,-c(3)]
    
     colnames(acc_mat)<-c("kfold CV accuracy","Permuted kfold CV accuracy")
   mainlab<-paste("Performance evaluation using ",best_classifier_name," classifier and selected features",sep="")
   
   if(output.device.type!="pdf"){
       
       temp_filename_1<-"Figures/Barplot_classification_accuracy.png"
       
       png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
   }
   
   
  
    barplot(acc_mat,xaxt="n",main=mainlab,col=boxplot.col.opt,ylab="Classification accuracy(%)",ylim=c(0,100))
    axis(side=1,at=seq(1,2),labels=c("Best kfold CV","Permuted kfold CV"),cex.axis=0.9)
    
    
    
    if(output.device.type!="pdf"){
        
        try(dev.off(),silent=TRUE)
    }
    
    
}

    
    
 
}
    
    
    try(dev.off(),silent=TRUE)

    eval_mat1<-rbind(eval_mat1,eval_mat_perm_final)
    
   
    
    # print(eval_mat1)
    #emat1
    write.table(emat1,file="Classifier_accuracy_comparison_training_set.txt",sep="\t",row.names=TRUE)
    
    confounder_eval=NA
    if(is.na(confounder.matrix)==FALSE){
        
        response_mat<-confounder.matrix[,-c(1)]
        confounder_eval<-try(runlmreg(X=good_feats,Y=response_mat,fdrmethod=confounderfdrmethod,fdrthresh=confounderfdrthresh),silent=TRUE)
    }
    
    parentoutput_dir=outloc
    setwd(parentoutput_dir)
    
    if(globalcor==TRUE){
        
        print("##############Level 2: Correlation network analysis selected features###########")
        print(paste("Generating metabolome-wide ",cor.method," correlation network",sep=""))
        data_m_fc_withfeats<-as.data.frame(Xorig)
        
        good_feats<-as.data.frame(good_feats)
        #print(goodfeats[1:4,])
        sigfeats_index<-which(data_m_fc_withfeats$mz%in%good_feats$mz)
        sigfeats<-sigfeats_index
        
        #outloc<-paste(parentoutput_dir,"/Allcornetworksigfeats","log2fcthresh",log2.fold.change.thresh,"/",sep="")
        outloc<-paste(parentoutput_dir,"/MWASresults","/",sep="")
        
        dir.create(outloc)
        setwd(outloc)
        
        #cor.method="spearman",networktype="complete",abs.cor.thresh=0.4,cor.fdrthresh=0.05,max.cor.num=100,net_node_colors=c("green","red"), net_legend=TRUE
        
        if(networktype=="complete"){
            mwan_fdr<-do_cor(data_m_fc_withfeats,subindex=sigfeats_index,targetindex=NA,outloc,networkscope="global",cor.method,abs.cor.thresh,cor.fdrthresh,max.cor.num,net_node_colors,net_legend)
        }else{
            if(networktype=="GGM"){
                mwan_fdr<-get_partial_cornet(data_m_fc_withfeats, sigfeats.index=sigfeats_index,targeted.index=NA,networkscope="global",cor.method,abs.cor.thresh,cor.fdrthresh,outloc=outloc,net_node_colors,net_legend)
            }else{
                print("Invalid option. Please use complete or GGM.")
            }
        }
        
        print("##############Level 2: processing complete###########")
    }
    
    if(length(good_feats)>0){
        
        write.csv(good_feats,file="train_selected_data.csv",row.names=FALSE)
        try(write.csv(good_feats_test,file="test_selected_data.csv",row.names=FALSE),silent=TRUE)
        try(write.csv(confounder_eval,file="confounder_eval.csv",row.names=FALSE),silent=TRUE)
        
        write.csv(Ytrain_mat,file="train_classlabels_data.csv",row.names=FALSE)
        try(write.csv(Ytest_mat,file="test_classlabels_data.csv",row.names=FALSE),silent=TRUE)
    }
    
    
    #print(acc_mat)
    #bestcvacc=best_kfold_acc,permcvacc=permkfold_acc,testacc=test_acc,
return(list(train.selected.data=good_feats,test.selected.data=good_feats_test,classifier.comparison=emat1,feature.selection.matrix=stability_matrix_1,best.performance.measures=acc_mat,train.class=Ytrain_mat,test.class=Ytest_mat,confounder.eval.res=confounder_eval))

    }



}
