
diffRank<-function(adjMtrxSample1,adjMtrxSample2){
    
   
    #N is the number if genes in the sample.
    N=dim(adjMtrxSample1)[1];
    
    # Create a graph from the dataset; one graph for each condition.
    graph1=graph.adjacency(adjMtrxSample1, mode="undirected",weighted=TRUE)#, diag=FALSE )
    graph2=graph.adjacency(adjMtrxSample2, mode="undirected",weighted=TRUE)#, diag=FALSE )
    

    #betweennessGraph1=eigen_centrality(graph1,directed = FALSE,weights=abs(E(graph1)$weight))$vector #,normalized=TRUE);
    #betweennessGraph2=eigen_centrality(graph2,directed = FALSE,weights=abs(E(graph2)$weight))$vector #,normalized=TRUE);
    
    betweennessGraph1=betweenness(graph1,directed = FALSE,weights=abs(E(graph1)$weight),normalized=TRUE);
    betweennessGraph2=betweenness(graph2,directed = FALSE,weights=abs(E(graph2)$weight),normalized=TRUE);
    
    
    betweennessgraph1DegreeNormalized =  betweennessGraph1 #/ max(betweennessGraph1);
    betweennessgraph2DegreeNormalized =  betweennessGraph2 #/ max(betweennessGraph2);
    DBC=abs(betweennessgraph1DegreeNormalized-betweennessgraph2DegreeNormalized)
    
    
    #Create and initialize Delta_C_i
    Delta_C_i=abs(adjMtrxSample1-adjMtrxSample2);

   
    Delta_C_i<-apply(Delta_C_i,1,function(x)
    {
        
        if(sum(x)==0)
        {
            x=rep(1/N,length(x))
        }
        x=x/sum(x)
        return(x)
        #})
    })
    
    #Delta_C_i<-do.call(rbind,Delta_C_i)
    Delta_C_i<-t(Delta_C_i)
    error=100;
    count=1;
    eps=0.001
    lambda=0.5
   
    
    # Create and initialize solution array to 1/N.
    # solution contains the ranks for all genes.
    solution=array(1/N,dim=c(N,1))
    solutionEachIteration=solution;
    
    # eps (EPSILON) is a value of the difference between 2 consecutive solutions to stop the iterations.
    # Do iterate while the difference between 2 consecutive solutions is not that much big (> eps).
    while(error > eps)
    {
        count=count+1;
        
        #Keep the previous solution
        formerSoulution=solution;
        
        #Find a new solution.
        for (v in 1:N){
        
        #solution<-lapply(1:N,function(v){
            s = sum(Delta_C_i[,v] * solution); # Solution, is the (Pi ) in this case
            solution[v]=(1-lambda) * DBC[v] + lambda * s;
            #res=(1-lambda) * DBC[v] + lambda * s;
            #return(res)
            #})
       }
        #solution<-unlist(solution)
        solutionEachIteration=array(c(solutionEachIteration,solution),dim=c(N,count))
        
        #find the error to stop the iteration, the error in this case that there no significant difference between the old and the new solution
        error=sum((formerSoulution - solution)^2)
        
        #print(error)
    }
    #print(head(solution))
    
    #rnk = rank(1-solution)
    
    return(solution)
}









runlmreg<-function(X,Y,fdrmethod="BH",fdrthresh=0.05,pvalue.thresh=0.05){
    
    data_m_fc<-X[,-c(1:2)]
    
    classlabels_response_mat<-Y
    data_m_fc_withfeats<-X
    logistic_reg<-FALSE
    
    rm(X)
    fileheader="lmreg"
    
    res1<-apply(data_m_fc,1,function(x){
        xvec<-x
        
        if(dim(classlabels_response_mat)[2]>1){
            
            for(cnum in 2:dim(classlabels_response_mat)[2]){
                
                classlabels_response_mat[,cnum]<-as.numeric(classlabels_response_mat[,cnum])
                
            }
        }
        
        data_mat_anova<-cbind(xvec,classlabels_response_mat)
        
        cnames<-colnames(data_mat_anova)
        cnames[1]<-"Response"
        
        colnames(data_mat_anova)<-cnames

        
        anova_res<-diffexplmreg(dataA=data_mat_anova,logistic_reg)
        
        return(anova_res)
    })
    
    ##save(res1,file="lmregres.Rda")
    main_pval_mat<-{}
    
    posthoc_pval_mat<-{}
    pvalues<-{}
    
    all_inf_mat<-{}
    
    for(i in 1:length(res1)){
        
        main_pval_mat<-rbind(main_pval_mat,res1[[i]]$mainpvalues)
        pvalues<-c(pvalues,res1[[i]]$mainpvalues[1])
        
        cur_pvals<-t(res1[[i]]$mainpvalues)
        cur_est<-t(res1[[i]]$estimates)
        cur_stderr<-t(res1[[i]]$stderr)
        cur_tstat<-t(res1[[i]]$statistic)
        cur_res<-cbind(cur_pvals,cur_est,cur_stderr,cur_tstat)
        
        all_inf_mat<-rbind(all_inf_mat,cur_res)
        
        
        
    }
    
    cnames_1<-c(paste(colnames(cur_pvals),".pvalue",sep=""),paste(colnames(cur_pvals),".estimate",sep=""),paste(colnames(cur_pvals),".tstatistic",sep=""))
    
    
    if(fdrmethod=="BH"){
        fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
    }else{
        if(fdrmethod=="ST"){
            #fdr_adjust_pvalue<-qvalue(pvalues)
            #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
            
            fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
            
            if(is(fdr_adjust_pvalue,"try-error")){
                
                fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
            }
      
            fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
            
            
            
        }else{
            if(fdrmethod=="Strimmer"){
                pdf("fdrtool.pdf")
                
                fdr_adjust_pvalue<-fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE)
                fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
                try(dev.off(),silent=TRUE)
            }else{
                if(fdrmethod=="none"){
                    #fdr_adjust_pvalue<-pvalues
                    fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
                }else{
                    if(fdrmethod=="BY"){
                        fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                    }else{
                        if(fdrmethod=="bonferroni"){
                            fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                            
                        }
                    }
																}
            }
        }
        
        
        
    }
    
    
    if(fdrmethod=="none"){
        filename<-paste(fileheader,"_pvalall_withfeats.txt",sep="")
        
    }else{
        filename<-paste(fileheader,"_fdrall_withfeats.txt",sep="")
    }
    cnames_tab<-colnames(data_m_fc_withfeats)
    cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
    
    pvalues<-as.data.frame(pvalues)
   
    final.pvalues<-pvalues
    sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
    
    data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
    
    colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
    
    
    filename<-paste(fileheader,"_pval_coef_stderr.txt",sep="")
    
    data_allinf_withfeats<-cbind(all_inf_mat,data_m_fc_withfeats)
    
    write.table(data_allinf_withfeats, file=filename,sep="\t",row.names=FALSE)
    
    cnames_tab<-colnames(data_m_fc_withfeats)
    
    class_column_names<-colnames(classlabels_response_mat)
    
    # cnames_tab<-c(paste("P.value_var",1:dim(classlabels_response_mat)[2],sep=""),
    #paste("Estimate_var",1:dim(classlabels_response_mat)[2],sep=""), paste("StdError_var",1:dim(classlabels_response_mat)[2],sep=""),
    #paste("statistic_var",1:dim(classlabels_response_mat)[2],sep=""),cnames_tab)
    
    
    #   cnames_tab<-c(paste("P.value_",class_column_names,sep=""),
    #paste("Estimate_",class_column_names,sep=""), paste("StdError_",class_column_names,sep=""),
    #paste("statistic_",class_column_names,sep=""),cnames_tab)
    
    cnames_1<-c(paste("P.value_",colnames(cur_pvals),sep=""),paste("Estimate_",colnames(cur_pvals),sep=""),paste("t-statistic_",colnames(cur_pvals),sep=""))
    
    cnames_tab<-c(cnames_1,cnames_tab)
    colnames(data_allinf_withfeats)<-as.character(cnames_tab)
    
    
    write.table(data_allinf_withfeats, file=filename,sep="\t",row.names=FALSE)
    
    return(data_allinf_withfeats)
    
}



#c("rf","rfe","limma","lasso","elasticnet","f.test")
diffexp.biomarkers<-function(X=NA,Y=NA,feature_table_file=NA,class_labels_file=NA,feat.sel.methods=c("rfe","rf","limma","lasso","elasticnet","wilcox.test","pls","spls","ospls","opls","f.test","t.test"),num.var.sel=10,iter.quantile.thresh=0.1,split.train.test=FALSE,train.pct=0.7,outloc,kfold=10,pca.ellipse=TRUE,Xtest=NA,Ytest=NA,rsdthresh=1,pls_vip_thresh=2,seedvalue=27611,learningsetmethod="CV",confounder.matrix=NA,confounderfdrmethod="BH",confounderfdrthresh=0.05,num.methods.sel=2,globalcor=FALSE,cor.method="spearman",networktype="complete",abs.cor.thresh=0.4,cor.fdrthresh=0.05,max.cor.num=100,net_node_colors=c("green","red"), net_legend=TRUE,niter=10,output.device.type="pdf",heatmap.col.opt="RdBu",boxplot.col.opt=c("white"),barplot.col.opt=c("grey57"),sample.col.opt="rainbow",mz.thresh=5,time.thresh=10,svm_kernel="radial",good_feats_index=NA,pvalue.thresh=0.05,plots.width=8,plots.height=8,plots.res=600, plots.type="cairo",num_nodes=2,ylabel="Intensity",cex.plots=0.7,size=(5),decay=5e-4,maxit=200,rang=0.01)
{
    
    
   
   #library(randomForest)
   #library(CMA)
   # importance=randomForest::importance
   
   unlockBinding("importance", as.environment("package:ranger"))
   #assign("importance", importance, "package:randomForest")
    assign("importance", randomForest::importance, "package:ranger")
    errortype="BER"
    if(is.na(X)==TRUE){
    X<-read.table(feature_table_file,sep="\t",header=TRUE)
    }
    
    if(is.na(Y)==TRUE){
    Y<-read.table(class_labels_file,sep="\t",header=TRUE)
    }
    
    
    cl<-makeSOCKcluster(num_nodes)
    
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
        
       
        
        #   #save(Xtest,file="Xtest.Rda")
        ##save(Ytest,file="Ytest.Rda")
        ##save(Xorig,file="Xorig.Rda")
        ##save(train_index,file="train_index.Rda")
        ##save(test_index,file="test_index.Rda")
        
        
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
  
        importance<-randomForest::importance
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
        set.seed(seedvalue)
        tune_plslda <- CMA::tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_ldaCMA, grids = list(comp = 1:5),trace=FALSE)
        
        if(length(class_levels)==2){
            #     tune_plslr <- tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_lrCMA, grids = list(comp = 1:5))
        }
        #print("Sts2")
        set.seed(seedvalue)
        tune_plsrf <- CMA::tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_rfCMA, grids = list(comp = 1:5),trace=FALSE)
    }
    #if(classifier=="scda")
    {
        
       
        save(X,Y,fiveCV10iter,file="scda.Rda")
        set.seed(seedvalue)
        tune_scda <-try(CMA::tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA, grids = list( ),trace=FALSE),silent=TRUE)
        
        
        if(is(tune_scda,"try-error")){
            
            #t1<-new("list")
            #t1<-as.list(rep(0.5,nrow(fiveCV10iter@learnmatrix)))
            #  tune_scda<-new("tuningresult",tuneres=t1,method="scDA")
            
            
        }
    }
    #if(classifier=="svm")
    {
       
        set.seed(seedvalue)
        tune_svm <- CMA::tune(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, grids = list( ), kernel = "radial",trace=FALSE)
    }
    
    if(length(good_feats_index)>2){
        set.seed(seedvalue)
        class_plslda <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_ldaCMA, tuneres = tune_plslda,trace=FALSE)
        set.seed(seedvalue)
        class_plsrf <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_rfCMA, tuneres = tune_plsrf,trace=FALSE)
    }
    
    set.seed(seedvalue)
    
    if(is(tune_scda,"try-error")){
    
        class_scda <- NA
    }else{
        class_scda <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA, tuneres = tune_scda,trace=FALSE)
    }
    
   
    set.seed(seedvalue)
    class_svm <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, tuneres = tune_svm,kernel = "radial",trace=FALSE,probability=TRUE)
    set.seed(seedvalue)
    class_rf <- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = rfCMA,trace=FALSE)
    set.seed(seedvalue)
    class_nnet<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE)
    set.seed(seedvalue)
    class_plr<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = plrCMA,trace=FALSE)
    
    #save(class_plslda,file="class_plslda.Rda")
    #save(class_plsrf,file="class_plsrf.Rda")
    #save(class_scda,file="class_scda.Rda")
    #save(class_svm,file="class_svm.Rda")
    #save(class_rf,file="class_rf.Rda")
    #save(class_nnet,file="class_nnet.Rda")
    #save(class_plr,file="class_plr.Rda")
    
    
    
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
        
        if(is.na(class_scda)==FALSE){
            eval_scda<-evaluation(class_scda,measure=eval_measure)
        }else{
            eval_scda<-NA
            
        }
        eval_svm<-evaluation(class_svm,measure=eval_measure)
        eval_rf<-evaluation(class_rf,measure=eval_measure)
        eval_nnet<-evaluation(class_nnet,measure=eval_measure)
        eval_plr<-evaluation(class_plr,measure=eval_measure)

        if(eval_measure=="auc")
        {
            eval_plslda@score=1-mean(eval_plslda@score)
            eval_plsrf@score=1-mean(eval_plsrf@score)
            
            if(is.na(class_scda)==FALSE){
                eval_scda@score=1-mean(eval_scda@score)
            }else{
                
                eval_scda=eval_plslda
                eval_scda@score=NA
            }
            eval_svm@score=1-mean(eval_svm@score)
            eval_rf@score=1-mean(eval_rf@score)
            eval_nnet@score=1-mean(eval_nnet@score)
            eval_plr@score=1-mean(eval_plr@score)
            
        }
    }
    
    ##save(eval_svm,file="eval_svm.Rda")
    ##save(class_svm,file="class_svm.Rda")
    
    
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
    
    
    

       barplot(eval_mat_1,xaxt="n",main="Comparison of CV accuracy using different classifiers\n based on learning sets", xlab="Classifier",ylab="kfold classification accuracy(%)",col=barplot.col.opt,type="p",ylim=c(0,100))
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
            
            temp_filename_1<-"Figures/HCA_selectedfeats.png"
            
            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        }
        
        
       
    try(get_hca(parentoutput_dir=outloc,X=good_feats,Y=Y1,heatmap.col.opt=heatmap.col.opt,cor.method="spearman",is.data.znorm=FALSE,analysismode="classification",
    sample.col.opt="rainbow",plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type=hca_type,newdevice=FALSE),silent=TRUE)
    
   
    
    if(output.device.type!="pdf"){
        
        try(dev.off(),silent=TRUE)
    }
    
    
    if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/PCAplots_selectedfeats.pdf"
        
        #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        pdf(temp_filename_1)
    }
    
    
    # #save(list=ls(),file="debug1.Rda")
    try(get_pcascoredistplots(X=good_feats,Y=Y1,parentoutput_dir=outloc,sample.col.opt=sample.col.opt,alphacol=0.3,col_vec=NA,pairedanalysis=FALSE,pca.cex.val=3,legendlocation="topright",pca.ellipse=FALSE,ellipse.conf.level=ellipse.conf.level,filename="selected",paireddesign=paireddesign,lineplot.col.opt=lineplot.col.opt,lineplot.lty.option=lineplot.lty.option),silent=TRUE)

    }
    
  
    if(output.device.type!="pdf"){
        
        try(dev.off(),silent=TRUE)
    }
    
    
    if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/Boxplots_selectedfeats.pdf"
        
        #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        pdf(temp_filename_1)
    }
    

    #par(mfrow=c(2,2))
    par(mfrow=c(1,1),family="sans",cex=cex.plots)
    get_boxplots(X=good_feats,Y=Y1,parentoutput_dir=outloc,boxplot.col.opt=boxplot.col.opt,sample.col.opt=sample.col.opt,alphacol=0.3,newdevice=FALSE,cex=0.6,ylabel=ylabel)
    
    
    
   
    
    if(output.device.type!="pdf"){
        
        try(dev.off(),silent=TRUE)
    }
    
    num_levels<-levels(as.factor(Y))
  
    
    
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
        
        nperm<-5
        seedvalue_rand_list<-runif(nperm,1,10000)
        #for(p1 in 1:nperm)
        eval_mat_perm<-lapply(1:nperm,function(p1)
        {
            seedvalue_cur=round(seedvalue_rand_list[p1],0)
            #set.seed(27611)
            set.seed(seedvalue_cur)
            Y<-Yorig[sample(1:length(Yorig),size=length(Yorig))]
            #set.seed(27611)
            set.seed(seedvalue_cur)
            fiveCV10iter<-GenerateLearningsets(y=Y,method=learningsetmethod,fold=kfold,niter=1,strat=TRUE)
            
            classifier_names<-c("PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR") #,"PLSLR","pLRlasso","pLRelasticnet")
            
            
            if(length(good_feats_index)>5){
            
            set.seed(seedvalue)
            class_plslda <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_ldaCMA,trace=FALSE)
                testeval_res1<-evaluation(class_plslda,measure=eval_measure)
            
            set.seed(seedvalue)
             class_plsrf <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_rfCMA,trace=FALSE)
                testeval_res2<-evaluation(class_plsrf,measure=eval_measure)
                
            }
            set.seed(seedvalue)
            
            if(is.na(class_scda)==FALSE){
               class_scda <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA,trace=FALSE)
                testeval_res3<-evaluation(class_scda,measure=eval_measure)
            }else{
                testeval_res3<-NA
            }
            
           set.seed(seedvalue)
              class_svm <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA,kernel = "radial",trace=FALSE,probability=TRUE)
                testeval_res4<-evaluation(class_svm,measure=eval_measure)
            
            set.seed(seedvalue)
                class_rf <- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = rfCMA,trace=FALSE)
                testeval_res5<-evaluation(class_rf,measure=eval_measure)
           set.seed(seedvalue)
                class_nnet<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE)
                testeval_res6<-evaluation(class_nnet,measure=eval_measure)
           set.seed(seedvalue)
                class_plr<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = plrCMA,trace=FALSE)
                testeval_res7<-evaluation(class_plr,measure=eval_measure)
                
                
                if(eval_measure=="auc")
                {
                    eval_plslda@score=1-mean(eval_plslda@score)
                    eval_plsrf@score=1-mean(eval_plsrf@score)
                    
                    if(is.na(class_scda)==FALSE){
                        eval_scda@score=1-mean(eval_scda@score)
                    }
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
         
                if(is.na(class_scda)==FALSE){
                    permkfold_acc3<-c(permkfold_acc3,100*(1-mean(testeval_res3@score)))
                }else{
                    
                    permkfold_acc3<-NA
                }
                
                permkfold_acc4<-c(permkfold_acc4,100*(1-mean(testeval_res4@score)))
                
                permkfold_acc5<-c(permkfold_acc5,100*(1-mean(testeval_res5@score)))
                
                permkfold_acc6<-c(permkfold_acc6,100*(1-mean(testeval_res6@score)))
                
                permkfold_acc7<-c(permkfold_acc7,100*(1-mean(testeval_res7@score)))
                
                temp_res<-c(100*(1-mean(testeval_res1@score)),100*(1-mean(testeval_res2@score)),100*(1-mean(testeval_res3@score)),
                100*(1-mean(testeval_res4@score)),100*(1-mean(testeval_res5@score)),100*(1-mean(testeval_res6@score)),
                    100*(1-mean(testeval_res7@score)))
                    
                return(temp_res)
            })
        
        if(FALSE){
            permkfold_acc1<-mean(permkfold_acc1)
            permkfold_acc2<-mean(permkfold_acc2)
            permkfold_acc3<-mean(permkfold_acc3)
            permkfold_acc4<-mean(permkfold_acc4)
            permkfold_acc5<-mean(permkfold_acc5)
            permkfold_acc6<-mean(permkfold_acc6)
            permkfold_acc7<-mean(permkfold_acc7)
            
            eval_mat_perm<-cbind(permkfold_acc1, permkfold_acc2, permkfold_acc3, permkfold_acc4, permkfold_acc5, permkfold_acc6, permkfold_acc7)
        }
        
        eval_mat_perm<-do.call(rbind,eval_mat_perm)
        eval_mat_perm<-apply(eval_mat_perm,2,mean)
        
        eval_mat_perm<-round(eval_mat_perm,2)
         
        eval_mat_perm_final<-cbind(text2B,eval_mat_perm)
        
       
        save(eval_mat_perm_final,file="eval_mat_perm_final.Rda")
        print(dim(eval_mat_perm_final))
        print("here")
        colnames(eval_mat_perm_final)<-c("Dataset","PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR")
        print("here2")
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
    
    
    if(num_levels==2){
        if(split.train.test==TRUE){
            #get_roc(dataA=good_feats,classlabels=Y,classifier="svm",kname=svm_kernel,rocfeatlist=c(dim(good_feats)[2]),rocfeatincrement=FALSE,testset=good_feats_test,testclasslabels=Ytest,mainlabel="Test set")
           
           
           
           res<-get_classification.accuracy(kfold=kfold,featuretable=good_feats,classlabels=Y,classifier="logit",testfeaturetable=good_feats_test,testclasslabels=Ytest,errortype="BAR",kernelname=svm_kernel)

# res<-get_classification.accuracy(kfold=kfold,featuretable=good_feats,classlabels=Y,classifier="svm",testfeaturetable=good_feats_test,testclasslabels=Ytest,errortype="BAR",kernelname=svm_kernel)
#get_classification.accuracy<(kfold,featuretable,classlabels,kernelname="radial",errortype="AUC",conflevel=95,classifier="svm",seednum=555,testfeaturetable=NA,testclasslabels=NA)
            
        }else{
            #get_roc(dataA=good_feats,classlabels=Y,classifier="svm",kname=svm_kernel,rocfeatlist=seq(2,10,1),rocfeatincrement=TRUE,testset=X,testclasslabels=Y,mainlabel="Using training set based on SVM")
            
            #get_roc(dataA=good_feats,classlabels=Y,classifier="svm",kname=svm_kernel,rocfeatlist=c(dim(good_feats)[2]),rocfeatincrement=FALSE,testset=good_feats_test,testclasslabels=Ytest,mainlabel="Test set")
            
              res<-get_classification.accuracy(kfold=kfold,featuretable=good_feats,classlabels=Y,classifier="logit",testfeaturetable=good_feats_test,testclasslabels=Ytest,errortype="BAR",kernelname=svm_kernel)
              
              # res<-get_classification.accuracy(kfold=kfold,featuretable=good_feats,classlabels=Y,classifier="svm",testfeaturetable=good_feats_test,testclasslabels=Ytest,errortype="BAR",kernelname=svm_kernel)
        }
    }
    
    test_acc_mat<-{}
    class_levels_vec<-levels(as.factor(Y))
    
    
    learnmatrix <- matrix(seq(1,nrow(X)), nrow = 1)
    fiveCV10iter<-new("learningsets", learnmatrix = learnmatrix, method = "none",ntrain = ncol(learnmatrix), iter = nrow(learnmatrix))
    X<-rbind(X,Xtest)
    Y<-c(Y,Ytest)
    classifier_names<-c("PLSLDA","PLSRF","SCDA","SVM","RF","NNet","pLR","PLSLR","pLRlasso","pLRelasticnet")
    
    
    
    save(X,file="X1.Rda")
    save(Y,file="Y1.Rda")
    save(learnmatrix,file="learnmatrix.Rda")
    
    save(list=ls(),file="debug3.Rda")
    
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
        
        set.seed(seedvalue)
        class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_ldaCMA, tuneres = tune_plslda1,trace=FALSE)
        
        b1<-best(tune_plslda1)
        
        save(b1,file="b1.Rda")
        
        learnmatrix<-as.numeric(learnmatrix)
        
       set.seed(seedvalue)
        class_res2<-pls_ldaCMA(X = X, y = Y, learnind = learnmatrix, comp = median(unlist(b1)))
        
        #save(class_res2,file="pls_ldaCMA.Rda")
        
        print(ftable(class_res2))
        
        #save(confusion_matrix_1,file="confusion_matrix_plslda.Rda")
        confusion_matrix_list[[1]]<-table(class_res2@y,class_res2@yhat)
        
        if(length(class_levels_vec)==2){
            testeval_res_auc<-evaluation(class_res,measure = "auc")
            #save(testeval_res_auc,file="testeval_res_auc.Rda")
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
        set.seed(seedvalue)
        class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = pls_rfCMA, tuneres = tune_plsrf1,trace=FALSE)
        
        b1<-best(tune_plsrf1)
        # #save(b1,file="b1_pls_rfCMA.Rda")
        
        learnmatrix<-as.numeric(learnmatrix)
        
        set.seed(seedvalue)
        class_res2<-pls_rfCMA(X = X, y = Y, learnind = learnmatrix, comp = median(unlist(b1)))
        
        #save(class_res2,file="pls_rfCMA.Rda")
        print(ftable(class_res2))
        
        #confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
        #save(confusion_matrix_1,file="confusion_matrix_rfCMA.Rda")
        
        confusion_matrix_list[[2]]<-table(class_res2@y,class_res2@yhat)
        
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
        
        
       
        s1<-ldply(tune_scda@tuneres,rbind)
        s2<-apply(s1,2,mean)
       
        t1<-new("list")
        t1[[1]]<-s2
        
        tune_scda1<-tune_scda
        tune_scda1@tuneres<-t1
        
        #tune_scda<-new("tuningresult",tuneres=t1)
        
        set.seed(seedvalue)
        class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = scdaCMA, tuneres = tune_scda1,trace=FALSE)
        #scdaCMA(X, y, f, learnind, delta = 0.5, models=FALSE,...)
        
        b1<-best(tune_scda1)
        
        save(b1,file="b1_scda.Rda")
        learnmatrix<-as.numeric(learnmatrix)
        set.seed(seedvalue)
        class_res2<-scdaCMA(X = X, y = Y, learnind = learnmatrix, delta = median(unlist(b1)))
        #save(class_res2,file="scdaCMA.Rda")
        print(ftable(class_res2))
        
        #confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
        
        #save(confusion_matrix_1,file="confusion_matrix_scdaCMA.Rda")
        confusion_matrix_list[[3]]<-table(class_res2@y,class_res2@yhat)
        
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
            
            #class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, tuneres=tune_svm1, kernel ="radial",trace=FALSE,probability=TRUE)
            
            print("doing SVM")
            
            save(X,Y,fiveCV10iter,tune_svm1,file="debug_svm.Rda")
            set.seed(seedvalue)
            class_res <- classification(X = X, y = Y, learningsets = fiveCV10iter, classifier = svmCMA, tuneres=tune_svm1, kernel ="radial",trace=FALSE,probability=TRUE)
            
            
            print("Here1")
            b1<-best(tune_svm1)
            
            save(b1,file="b1_svm.Rda")
            
            learnmatrix<-as.numeric(learnmatrix)
            
          set.seed(seedvalue)
            class_res2<-svmCMA(X = X, y = Y, learnind = learnmatrix, delta = median(unlist(b1)),probability=TRUE)
            
            print("Here2")
            #save(class_res2,file="svmCMA.Rda")
            
            print("svm done")
            print(ftable(class_res2))
            
            
            confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            #save(confusion_matrix_1,file="confusion_matrix_svmCMA.Rda")
            confusion_matrix_list[[4]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
                
                testeval_res_auc<-evaluation(class_res,measure = "auc")
                ##save(testeval_res_auc,file="testeval_res_auc.Rda")
                
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
            set.seed(seedvalue)
            class_res <- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = rfCMA,trace=FALSE)
            
            learnmatrix<-as.numeric(learnmatrix)
            
          
            set.seed(seedvalue)
            class_res2 <- rfCMA(X =X, y = Y, learnind=learnmatrix, varimp = FALSE)
            
            #save(class_res2,file="rfCMA.Rda")
            
            print(ftable(class_res2))
            
            
            confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            #save(confusion_matrix_1,file="confusion_matrix_rfCMA.Rda")
            confusion_matrix_list[[5]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
                
                testeval_res_auc<-evaluation(class_res,measure = "auc")
                ##save(testeval_res_auc,file="testeval_res_auc.Rda")
                
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
            set.seed(seedvalue)
            class_res<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = nnetCMA,trace=FALSE)
            testeval_res_auc<-evaluation(class_res,measure = "auc")
            print(testeval_res_auc)
            
            set.seed(seedvalue)
            #nnetresult <- nnetCMA(X=golubX, y=golubY, learnind=learnind, size = 3, decay = 0.01)
            class_res2 <- nnetCMA(X =X, y = Y, learnind=as.numeric(learnmatrix))
            
            #save(class_res2,file="nnetCMA.Rda")
            
            print(ftable(class_res2))
            
            
            confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            #save(confusion_matrix_1,file="confusion_matrix_nnetCMA.Rda")
            confusion_matrix_list[[6]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
                
                testeval_res_auc<-evaluation(class_res,measure = "auc")
                ##save(testeval_res_auc,file="testeval_res_auc.Rda")
                
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
            set.seed(seedvalue)
            class_res<- classification(X =X, y = Y, learningsets = fiveCV10iter, classifier = plrCMA,trace=FALSE)
            
            set.seed(seedvalue)
            class_res2 <- plrCMA(X =X, y = Y, learnind=learnmatrix)
            
            #save(class_res2,file="plrCMA.Rda")
            
            print(ftable(class_res2))
            
            
            confusion_matrix_1<-table(class_res2@y,class_res2@yhat)
            
            #save(confusion_matrix_1,file="confusion_matrix_plrCMA.Rda")
            confusion_matrix_list[[7]]<-table(class_res2@y,class_res2@yhat)
            
            if(length(class_levels_vec)==2){
                
                testeval_res_auc<-evaluation(class_res,measure = "auc")
                ##save(testeval_res_auc,file="testeval_res_auc.Rda")
                
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
        
        
        barplot(acc_mat,xaxt="n",main=mainlab,col=barplot.col.opt,ylab="Classification accuracy(%)",ylim=c(0,100))
        
        if(length(class_levels_vec)==2){
            axis(side=1,at=seq(1,4),labels=c("kfold CV","Permuted kfold CV","Test set","AUC"),cex.axis=cex.plots)
        }else{
            axis(side=1,at=seq(1,3),labels=c("kfold CV","Permuted kfold CV","Test set"),cex.axis=cex.plots)
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
   
   
  
    barplot(acc_mat,xaxt="n",main=mainlab,col=barplot.col.opt,ylab="Classification accuracy(%)",ylim=c(0,100))
    axis(side=1,at=seq(1,2),labels=c("Best kfold CV","Permuted kfold CV"),cex.axis=cex.plots)
    
    
    
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

plotEigengeneNetworks_custom<-function (multiME, setLabels, letterSubPlots = FALSE, Letters = NULL,
    excludeGrey = TRUE, greyLabel = "grey", plotDendrograms = TRUE,
    plotHeatmaps = TRUE, setMargins = TRUE, marDendro = NULL,
    marHeatmap = NULL, colorLabels = TRUE, signed = TRUE, heatmapColors = NULL,
    plotAdjacency = TRUE, printAdjacency = FALSE, cex.adjacency = 0.9,
    coloredBarplot = TRUE, barplotMeans = TRUE, barplotErrors = FALSE,
    plotPreservation = "standard", zlimPreservation = c(0, 1),
    printPreservation = FALSE, cex.preservation = 0.9)
{
    size = checkSets(multiME, checkStructure = TRUE)
    if (!size$structureOK) {
        multiME = fixDataStructure(multiME)
    }
    if (is.null(Letters))
        Letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if (is.null(heatmapColors))
        if (signed) {
            heatmapColors = greenWhiteRed(50)
        }
        else {
            heatmapColors = topo.colors(30) #heat.colors(30)
        }
    heatmapColors = topo.colors(30)
    
    nSets = length(multiME)
    cex = par("cex")
    mar = par("mar")
    nPlotCols = nSets
    nPlotRows = as.numeric(plotDendrograms) + nSets * as.numeric(plotHeatmaps)
    if (nPlotRows == 0)
        stop("Nothing to plot: neither dendrograms not heatmaps requested.")
    par(mfrow = c(nPlotRows, nPlotCols))
    par(cex = cex)
    if (excludeGrey)
        for (set in 1:nSets) multiME[[set]]$data = multiME[[set]]$data[,
            substring(names(multiME[[set]]$data), 3) != greyLabel]
    plotPresTypes = c("standard", "hyperbolic", "both","change","differences")
    ipp = pmatch(plotPreservation, plotPresTypes)
    if (is.na(ipp))
        stop(paste("Invalid 'plotPreservation'. Available choices are",
            paste(plotPresTypes, sep = ", ")))
    letter.ind = 1
    if (plotDendrograms)
        for (set in 1:nSets) {
            par(mar = marDendro)
            labels = names(multiME[[set]]$data)
            uselabels = labels[substring(labels, 3) != greyLabel]
            corME = cor(multiME[[set]]$data[substring(labels,
                3) != greyLabel, substring(labels, 3) != greyLabel],
                use = "p")
            disME = as.dist(1 - corME)
            clust = flashClust(disME, method = "average")
            if (letterSubPlots) {
                main = paste(substring(Letters, letter.ind, letter.ind),
                  ". ", setLabels[set], sep = "")
            }
            else {
                main = setLabels[set]
            }
            plotLabels = uselabels
            plot(clust, main = main, sub = "", xlab = "", labels = plotLabels,
                ylab = "", ylim = c(0, 1))
            letter.ind = letter.ind + 1
        }
    if (plotHeatmaps)
        for (i.row in (1:nSets)) for (i.col in (1:nSets)) {
            letter.ind = i.row * nSets + i.col
            if (letterSubPlots) {
                letter = paste(substring(Letters, first = letter.ind,
                  last = letter.ind), ".  ", sep = "")
            }
            else {
                letter = NULL
            }
            par(cex = cex)
            if (setMargins) {
                if (is.null(marHeatmap)) {
                  if (colorLabels) {
                    par(mar = c(1, 2, 3, 4) + 0.2)
                  }
                  else {
                    par(mar = c(6, 7, 3, 5) + 0.2)
                  }
                }
                else {
                  par(mar = marHeatmap)
                }
            }
            nModules = dim(multiME[[i.col]]$data)[2]
            textMat = NULL
            if (i.row == i.col) {
                corME = cor(multiME[[i.col]]$data, use = "p")
                pME = corPvalueFisher(corME, nrow(multiME[[i.col]]$data))
                if (printAdjacency) {
                  textMat = paste(signif(corME, 2), "\n", signif(pME,
                    1))
                  dim(textMat) = dim(corME)
                }
                if (signed) {
                  if (plotAdjacency) {
                    if (printAdjacency) {
                      textMat = paste(signif((1 + corME)/2, 2),
                        "\n", signif(pME, 1))
                      dim(textMat) = dim(corME)
                    }
                    labeledHeatmap((1 + corME)/2, names(multiME[[i.col]]$data),
                      names(multiME[[i.col]]$data), main = paste(letter,
                        setLabels[[i.col]]), invertColors = FALSE,
                      zlim = c(0, 1), colorLabels = colorLabels,
                      colors = heatmapColors, setStdMargins = FALSE,
                      textMatrix = textMat, cex.text = cex.adjacency)
                  }
                  else {
                    labeledHeatmap(corME, names(multiME[[i.col]]$data),
                      names(multiME[[i.col]]$data), main = paste(letter,
                        setLabels[[i.col]]), invertColors = FALSE,
                      zlim = c(-1, 1), colorLabels = colorLabels,
                      colors = heatmapColors, setStdMargins = FALSE,
                      textMatrix = textMat, cex.text = cex.adjacency)
                  }
                }
                else {
                  labeledHeatmap(abs(corME), names(multiME[[i.col]]$data),
                    names(multiME[[i.col]]$data), main = paste(letter,
                      setLabels[[i.col]]), invertColors = FALSE,
                    zlim = c(0, 1), colorLabels = colorLabels,
                    colors = heatmapColors, setStdMargins = FALSE,
                    textMatrix = textMat, cex.text = cex.adjacency)
                }
            }
            else {
                corME1 = cor(multiME[[i.col]]$data, use = "p")
                corME2 = cor(multiME[[i.row]]$data, use = "p")
                cor.dif = (corME1 - corME2)/2
                d = tanh((corME1 - corME2)/(abs(corME1) + abs(corME2))^2)
                if (ipp == 1 | ipp == 3) {
                  dispd = cor.dif
                  main = paste(letter, "Preservation")
                  if (ipp == 3) {
                    dispd[upper.tri(d)] = d[upper.tri(d)]
                    main = paste(letter, "Hyperbolic preservation (UT)\nStandard preservation (LT)")
                  }
                }
                else {
			if (ipp == 2) {
                  dispd = d
                  main = paste(letter, "Hyperbolic preservation")
			}
			else{
				if (ipp == 4) {
				  dispd = 1-abs(cor.dif)
				  main = paste(letter, "Differences")
				  }
			}
		}
                if (i.row > i.col) {
                  if (signed) {
                    half = as.integer(length(heatmapColors)/2)
                    range = c(half:length(heatmapColors))
                    halfColors = heatmapColors[range]
                  }
                  else {
                    halfColors = heatmapColors
                  }
                  if (printPreservation) {
                    printMtx = matrix(paste(".", as.integer((1 -
                      abs(dispd)) * 100), sep = ""), nrow = nrow(dispd),
                      ncol = ncol(dispd))
                    printMtx[printMtx == ".100"] = "1"
                  }
                  else {
                    printMtx = NULL
                  }
                  if ((sum((1 - abs(dispd)) < zlimPreservation[1]) ||
                    ((1 - abs(dispd)) > zlimPreservation[2])) >
                    0)
                    warning("plotEigengeneNetworks: Correlation preservation data out of zlim range.")
                  labeledHeatmap(1 - abs(dispd), names(multiME[[i.col]]$data),
                    names(multiME[[i.col]]$data), main = main,
                    invertColors = FALSE, colorLabels = colorLabels,
                    zlim = zlimPreservation, colors = halfColors,
                    setStdMargins = FALSE, textMatrix = printMtx,
                    cex.text = cex.preservation)
                }
                else {
                  if (ipp == 2) {
                    dp = 1 - abs(d)
                    method = "Hyperbolic:"
                  }
                  else {
			  if (ipp == 4) {
				  dp = abs(cor.dif)
				  method = "Differences:"
			  }
			  else{
				    dp = 1 - abs(cor.dif)
				    method = "Preservation:"
			  }
		  }
                  diag(dp) = 0
                  #write.table(dp,file="preservation_matrix.txt",sep="\t")
                  if (barplotMeans) {
                    sum_dp = mean(dp[upper.tri(dp)])
                    means = apply(dp, 2, sum)/(ncol(dp) - 1)
                    if (barplotErrors) {
                      errors = sqrt((apply(dp^2, 2, sum)/(ncol(dp) -1) - means^2)/(ncol(dp) - 2))
                    }
                    else {
                      errors = NULL
                    }
		    Dmatrix=cbind(names(multiME[[i.col]]$data),means)
		    colnames(Dmatrix)<-c("Module","meanPreservationScore")
            
            fname<-paste("mean_preservation_matrix_",i.col,".txt",sep="")
		    write.table(Dmatrix,file=fname,sep="\t",row.names=FALSE)
                    labeledBarplot(means, names(multiME[[i.col]]$data),
                      main = paste(letter, "Pres.score vs ", i.col, signif(sum_dp,
                        2)), ylim = c(0, 1), colorLabels = colorLabels,
                      colored = coloredBarplot, setStdMargins = FALSE,
                      stdErrors = errors)
                  }
                  else {
                    sum_dp = sum(dp[upper.tri(dp)])
                    labeledBarplot(dp, names(multiME[[i.col]]$data),
                      main = paste(letter, method, "sum = ",
                        signif(sum_dp, 3)), ylim = c(0, dim(dp)[[1]]),
                      colorLabels = colorLabels, colored = coloredBarplot,
                      setStdMargins = FALSE)
                  }
                }
            }
        }
}


degree_eval<-function(feature_table_file=NA,class_labels_file=NA,X=NA,Y=NA,sigfeats=NA,sigfeatsind=NA){
   

#print("degree eval")
    if(is.na(X)==TRUE){
        data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE)
        
    }else{
        
        data_matrix<-X
    }
    if(is.na(Y)==TRUE){
        classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
        
    }else{
        
        classlabels<-Y
    }
    
    
    classlabels<-as.data.frame(classlabels)
    #print(dim(classlabels))
    # print(length(classlabels))
   
 
    if(dim(classlabels)[2]>2){
        classgroup<-paste(classlabels[,2],":",classlabels[,3],sep="") #classlabels[,2]:classlabels[,3]
    }else{
        
        classgroup<-classlabels[,2]
    }
    classlabels<-as.data.frame(classlabels)
    
    class_labels_levels<-levels(as.factor(classgroup))
    
    rnames<-paste(sprintf("%.4f",data_matrix$mz),data_matrix$time,sep="_")
    
    data_matrix_orig<-data_matrix
    data_matrix<-data_matrix[,-c(1:2)]
    #data_matrix<-na.omit(data_matrix)
    
    rnamesAB<-gsub(pattern="NA_NA",replacement=NA,x=rnames)
    rnamesAB<-na.omit(rnamesAB)
    
    nSets = length(class_labels_levels);
    multiExpr = vector(mode = "list", length = nSets)
    data_matrix_list<-new("list")
    num_samps_groups<-new("list")
    degree_list<-new("list")
    
    data_matrix_all<-t(data_matrix)
    
    
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    sft = pickSoftThreshold(data=data_matrix_all, dataIsExpr=TRUE,powerVector = powers, verbose = 0)
    power_val=sft$powerEstimate
    
    if(is.na(power_val)==TRUE){
        power_val=6
    }
    
    
    degree_overall<-softConnectivity(datExpr=data_matrix_all,power=power_val,minNSamples=10)
    degree_overall<-replace(degree_overall,which(degree_overall<0),0)
    
    
    sig_status<-rep(0,dim(data_matrix_orig)[1])
    
    colors_met<-cbind(data_matrix_orig[,c(1:2)],degree_overall,sig_status)
    
    colors_met<-as.data.frame(colors_met)
    
    colors_met$sig_status[sigfeatsind]=1
    
    colnames(colors_met)<-c("mz","time","degree_overall","sig_status")
    
    colnames_vec<-c("mz","time","degree_overall","sig_status")
    
    # print(sigfeats$mz)
    
    
    colors_met_all<-colors_met
    
    colors_met<-colors_met[order(colors_met$degree_overall,decreasing=TRUE),]
    
    sig_ind<-which(colors_met$sig_status==1)
    
    # print(summary(colors_met$degree_overall))
    
    #pdf("DICE_plots.pdf")
    
    names=paste(round(colors_met$mz,3),round(colors_met$time,0),sep="_")
    
    names=paste(round(colors_met$mz,4),sep="_")
    
    # print(names[sig_ind])
    
    if(FALSE)
    {
         plot(colors_met$degree_overall,cex=0.2,main="Overall degree distribution",ylab="Degree",xlab="Feature Index",col="orange",type="b",lwd=0.5)
       
       #for(i in 1:length(sig_ind)){
       lapply(1:length(sig_ind),function(i){
            points(y=colors_met$degree_overall[sig_ind[i]],x=sig_ind[i],col="darkgreen",cex=0.8,lwd=2)
            if(i%%2>0){
                text(y=colors_met$degree_overall[sig_ind[i]],x=sig_ind[i],names[sig_ind[i]],cex=0.31,adj=c(1,2))
            }else{
                text(y=colors_met$degree_overall[sig_ind[i]],x=sig_ind[i],names[sig_ind[i]],cex=0.31,adj=c(0,-1))
            }
            
         })
    }
    
    
    for(i in 1:length(class_labels_levels)){
        
        data_matrix_list[[i]]<-t(data_matrix[,which(classgroup==class_labels_levels[i])])
        num_samps_groups[[i]]<-dim(data_matrix_list[[i]])[1]
        #print(dim(data_matrix_list[[i]]))
        multiExpr[[i]]<-list(data = as.data.frame(data_matrix_list[[i]]));
        rownames(multiExpr[[i]]$data)=c(paste(rep(class_labels_levels[i],num_samps_groups[[i]]),seq(1,num_samps_groups[[i]]),sep=""))
        
        degree_list[[i]]<-softConnectivity(datExpr=data_matrix_list[[i]],power=power_val,minNSamples=2)
        #
        colnames_vec<-c(colnames_vec,paste("DegreeClass",i,sep=""))
        
        degree_list[[i]]<-replace(degree_list[[i]],which(degree_list[[i]]<0),0)
        
        #degree_list[[i]][which(is.na(degree_list[[i]])==TRUE)]=1
        
        colors_met<-cbind(data_matrix_orig[,c(1:2)],degree_list[[i]],sig_status)
        
        colors_met<-as.data.frame(colors_met)
        
        colors_met_all<-cbind(colors_met_all,degree_list[[i]])
        
        colors_met$sig_status[sigfeatsind]=1
        
        colnames(colors_met)<-c("mz","time","DegreeClass","sig_status")
        
        #colnames_vec<-c("mz","time","DegreeClass","sig_status")
        colors_met<-as.data.frame(colors_met)
        
        colors_met<-colors_met[order(colors_met$DegreeClass,decreasing=TRUE),]
        
        sig_ind<-which(colors_met$sig_status==1)
        
        names=paste(round(colors_met$mz,4),sep="_")
        
        # print(names[sig_ind])
        mainlab=paste("Class ",i," degree distribution",sep="")
        
        if(FALSE)
        {
            plot(colors_met$DegreeClass,cex=0.2,main=mainlab,ylab="Degree",xlab="Feature Index",col="orange",type="b",lwd=0.5)
            #for(i in 1:length(sig_ind))
            lapply(1:length(sig_ind),function(i)
            {
                points(y=colors_met$DegreeClass[sig_ind[i]],x=sig_ind[i],col="darkgreen",cex=0.8,lwd=2)
                
                if(i%%2>0){
                    text(y=colors_met$DegreeClass[sig_ind[i]],x=sig_ind[i],names[sig_ind[i]],cex=0.31,adj=c(1,2))
                }else{
                    text(y=colors_met$DegreeClass[sig_ind[i]],x=sig_ind[i],names[sig_ind[i]],cex=0.31,adj=c(0,-1))
                }
            })
            
        }
        
    }
    
    # dev.off()
    
    colors_met_all<-as.data.frame(colors_met_all)
    colnames(colors_met_all)<-colnames_vec
    
    
    #print(table(MET[[1]]$validColors))
    #print(table(classAmoduleColors))
    write.table(colors_met_all,file="Tables/Degree_eval_allfeats.txt",sep="\t",row.names=FALSE)
    
    sub_colors_met<-{}
    if(is.na(sigfeats)==FALSE){
        sub_colors_met<-colors_met_all[sigfeatsind,]
        write.table(sub_colors_met,file="Tables/Degree_eval_selectfeats.txt",sep="\t",row.names=FALSE)
    }
    
    ##save(MET,file="MET.Rda")
    
    
    
    return(list(all=colors_met_all,sigfeats=sub_colors_met))
}

#metab_data: Statistic, mz, time
#metab_annot: KEGGID, mz, time
get_fcs<-function(kegg_species_code="hsa",database="pathway",metab_data,metab_annot,reference_set=NA,type.statistic="pvalue"){
    
    if(is.na(reference_set)==TRUE){
        
        #homo sapiens: humans
        if(kegg_species_code=="hsa"){
            data(kegg_hsa)
        }else{
            
            #Mus musculus: mouse
            if(kegg_species_code=="mmu"){
                data(kegg_mmu)
            }else{
                
                #Pan troglodytes: Chimpanzee
                if(kegg_species_code=="ptr"){
                    data(kegg_ptr)
                }else{
                    
                    #Macaca mulatta: Rhesus monkey
                    if(kegg_species_code=="mcc"){
                        data(kegg_mcc)
                    }else{
                        #Bos taurus: cow
                        if(kegg_species_code=="bta"){
                            data(kegg_bta)
                        }else{
                            #Rattus norvegicus: rat
                            if(kegg_species_code=="rno"){
                                data(kegg_rno)
                            }else{
                                
                                #Danio rerio: Zebrafish
                                if(kegg_species_code=="dre"){
                                    data(kegg_dre)
                                }else{
                                    
                                    #C. elegans: nematode
                                    if(kegg_species_code=="cel"){
                                        data(kegg_cel)
                                    }else{
                                        
                                        #Drosophila melanogaster: fruit fly
                                        if(kegg_species_code=="dme"){
                                            data(kegg_dme)
                                        }
                                        
                                    }
                                    
                                }
                                
                            }
                            
                        }
                        
                    }
                }
                
            }
            
        }
            g1<-get_kegg_compounds(kegg_species_code=kegg_species_code,database=database)
    }else{
        g1<-reference_set
    }
    g1<-as.data.frame(g1)
    #colnames(hsa_module_comp)<-c("MetabSetID","KEGGID")
    colnames(metab_data)<-c("Statistic","mz","time")
    
    if(type.statistic=="pvalue"){
        
        metab_data$Statistic=(-1)*log10(metab_data$Statistic)
    }
    
    colnames(metab_annot)<-c("KEGGID","mz","time")
    
    metab_data_1<-merge(metab_data,metab_annot,by="mz")
    metab_data<-metab_data_1[,c("Statistic","KEGGID")]
    res<-get_fcs_child(metab_data=metab_data,reference_sets=g1)
    
    path_names_ids<-g1[-which(duplicated(g1$ID)==TRUE),1:2]
    
    res<-merge(res,path_names_ids,by.x="MetaboliteSet",by.y="ID")
    
    res<-res[order(as.numeric(as.character(res$pvalue)),decreasing=FALSE),]
    
    
    
    return(res)
    
}
#metab_data: KEGGID, Statistic
#reference sets: Pathway/Module ID, KEGGID
get_fcs_child<-function(metab_data,reference_sets,itrs=1000){
#g1<-get_kegg_compounds(kegg_species_code="hsa",database="module")
#hsa_module_comp<-as.data.frame(g1)
#colnames(hsa_module_comp)<-c("Module","KEGGID")
#cid_list<-names(t1);annot_res_2<-lapply(1:length(t1),function(c){sub_data<-annot_res[which(annot_res$chemical_ID==
#metab_data<-read.table("Stage2/limmasignalthresh0.8RSD1/Tables/limmaresults_allfeatures.txt",sep="\t",header=TRUE)
#metab_data<-metab_data[,c(4,6:7)]
#metab_data$P.value<-(-1)*log10(metab_data$P.value)
#colnames(metab_data)<-c("Statistic","mz","time")
#metab_data_1<-merge(metab_data,metab_annot,by="mz")
#metab_data<-metab_data_1[,c("Statistic","KEGGID")]

    metab_data_sets<-merge(metab_data,reference_sets,by="KEGGID")
    
    sid_list<-metab_data_sets[,3];
    sid_list<-unique(sid_list)
    
    metabset_res<-lapply(1:length(sid_list),function(c){
        
        sub_data<-metab_data_sets[which(as.character(metab_data_sets[,3])==sid_list[c]),];
        
        set_statistic<-sum(sub_data$Statistic)/nrow(sub_data)
        set_statistic<-round(set_statistic,3)
        res<-cbind(as.character(sid_list[c]),set_statistic)
        res<-as.data.frame(res)
        
        
        
        set.seed(500)
        
        avgrandstat <- replicate(itrs, {
            all <- sample(metab_data_sets[,2])
            randmetab_data_sets<-metab_data_sets
            randmetab_data_sets[,2]<-all
            randsub_data<-randmetab_data_sets[which(as.character(randmetab_data_sets[,3])==sid_list[c]),];
            
            randset_statistic<-sum(randsub_data$Statistic)/nrow(randsub_data)
            return(randset_statistic)
        })
        avgrandstat<-unlist(avgrandstat)
        
        pval_stat<-length(which(avgrandstat>set_statistic))/length(avgrandstat)
        
        res<-cbind(res,pval_stat)
        colnames(res)<-c("MetaboliteSet","Agg.Statistic","pvalue")
        
        return(res)
    })
    
    metabset_res_mat<-do.call(rbind,metabset_res)
    
    metabset_res_mat<-metabset_res_mat[order(as.numeric(as.character(metabset_res_mat$Agg.Statistic)),decreasing=TRUE),]
    
    
    if(FALSE){
            cid_list<-names(t1);
            filter.by=c("M+H")
            annot_res_2<-lapply(1:length(t1),function(c){
                sub_data<-annot_res[which(as.character(annot_res$chemical_ID)==cid_list[c]),];
                
                conf_level<-unique(sub_data[,2])
                
                if(conf_level==3){
                    
                    if(nrow(sub_data)<2){
                        
                        adduct_names=unique(as.character(sub_data[,13]))
                        if(length(which(adduct_names%in%filter.by))>0){
                            sub_data[,2]<-2
                            
                        }else{
                            sub_data[,2]<-0
                            
                        }
                        
                        
                    }
                }
                return(sub_data)
                
               })
            annot_res_3<-do.call(rbind,annot_res_2)
    }
    
    return(metabset_res_mat)

}

#get KEGG IDs for all pathways/modules for a given species id
get_kegg_compounds<-function(kegg_species_code="hsa",database="pathway"){
    
    
    path_list<-keggList(kegg_species_code,database=database)
    
    kegg_pathway_ids<-names(path_list)
    
    if(database=="pathway"){
        kegg_pathway_ids<-gsub(kegg_pathway_ids,pattern="path:",replacement="")
    }
    kegg_comp_list<-{}
    map_res<-{}
    kegg_module_list<-{}
    
    path_name_id_mapping<-cbind(kegg_pathway_ids,path_list)
    path_name_id_mapping<-as.data.frame(path_name_id_mapping)
    colnames(path_name_id_mapping)<-c("ID","Name")
    path_comp_mat<-{}
    path_comp_mat<-lapply(1:length(kegg_pathway_ids),function(p)
    {
        kegg_pathway_id=kegg_pathway_ids[p]
        Sys.sleep(0.01)
        k1<-keggGet(dbentries=kegg_pathway_id)
        
        kegg_comp_list<-c(kegg_comp_list,k1[[1]]$COMPOUND)
        
        if(length(kegg_comp_list)<1){
            #print(kegg_pathway_id)
            kegg_module_list<-c(kegg_module_list,k1[[1]]$MODULE)
            
        }else{
            path_comp_mat_temp<-cbind(kegg_pathway_id,names(k1[[1]]$COMPOUND))
            
            if(length(path_comp_mat_temp)>0){
                
                if(ncol(path_comp_mat_temp)==2){
                    
                    return(path_comp_mat_temp)
                }
            }
        }
    })
    path_comp_mat<-do.call(rbind,path_comp_mat)
    
    colnames(path_comp_mat)<-c("ID","KEGGID")
    
    path_comp_mat<-merge(path_name_id_mapping,path_comp_mat,by.x="ID",by.y="ID")
    return(path_comp_mat)
}



diffrank_eval<-function(feature_table_file=NA,class_labels_file=NA,X=NA,Y=NA,sigfeats=NA,sigfeatsind=NA,cor.method="pearson",num_nodes=2,abs.cor.thresh=0.4,pvalue.thresh=0.005,cor.fdrthresh=0.2,fdrmethod="Strimmer"){
    
    
    #print("degree eval")
    
    if(is.na(X)==TRUE){
        data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE)
        
    }else{
        
        data_matrix<-X
    }
    if(is.na(Y)==TRUE){
        classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
        
    }else{
        
        classlabels<-Y
    }
    
    
    classlabels<-as.data.frame(classlabels)
    #print(dim(classlabels))
    # print(length(classlabels))
    
    
    if(dim(classlabels)[2]>2){
        classgroup<-paste(classlabels[,2],":",classlabels[,3],sep="") #classlabels[,2]:classlabels[,3]
    }else{
        
        classgroup<-classlabels[,2]
    }
    classlabels<-as.data.frame(classlabels)
    
    class_labels_levels<-levels(as.factor(classgroup))
    
    rnames<-paste(sprintf("%.4f",data_matrix$mz),data_matrix$time,sep="_")
    
    data_matrix_orig<-data_matrix
    data_matrix<-data_matrix[,-c(1:2)]
    #data_matrix<-na.omit(data_matrix)
    
    rnamesAB<-gsub(pattern="NA_NA",replacement=NA,x=rnames)
    rnamesAB<-na.omit(rnamesAB)
    
    nSets = length(class_labels_levels);
    multiExpr = vector(mode = "list", length = nSets)
    data_matrix_list<-new("list")
    num_samps_groups<-new("list")
    degree_list<-new("list")
    
    data_matrix_all<-t(data_matrix)
    
    
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    sft = pickSoftThreshold(data=data_matrix_all, dataIsExpr=TRUE,powerVector = powers, verbose = 0)
    power_val=sft$powerEstimate
    
    if(is.na(power_val)==TRUE){
        power_val=6
    }
    
 
    
    #degree_overall<-softConnectivity(datExpr=data_matrix_all,power=power_val,minNSamples=10)
    #degree_overall<-replace(degree_overall,which(degree_overall<0),0)
    
    
    sig_status<-rep(0,dim(data_matrix_orig)[1])
    
    colors_met<-cbind(data_matrix_orig[,c(1:2)],sig_status)
    
    colors_met<-as.data.frame(colors_met)
    
    colors_met$sig_status[sigfeatsind]=1
    
    colnames(colors_met)<-c("mz","time","sig_status")
    
    colnames_vec<-c("mz","time","sig_status")
    
    sig_ind<-which(colors_met$sig_status==1)
    
    
    names=paste(round(colors_met$mz,4),round(colors_met$time,0),sep="_")
    
    #names=paste(round(colors_met$mz,4),sep="_")
    

    
    adj_mat_list<-lapply(1:length(class_labels_levels),function(i){
        
        data_matrix_list[[i]]<-t(data_matrix[,which(classgroup==class_labels_levels[i])])
        num_samps_groups[[i]]<-dim(data_matrix_list[[i]])[1]
        #print(dim(data_matrix_list[[i]]))
        multiExpr[[i]]<-list(data = as.data.frame(data_matrix_list[[i]]));
        rownames(multiExpr[[i]]$data)=c(paste(rep(class_labels_levels[i],num_samps_groups[[i]]),seq(1,num_samps_groups[[i]]),sep=""))
        
        cor.method="pearson"
        simmat<-WGCNA::cor(data_matrix_list[[i]],nThreads=num_nodes,method=cor.method,use = 'p')
        
        simmat<-round(simmat,5)
        pearson_resvec<-as.vector(simmat)
        
        #complete_pearsonqvalue_mat<-unlist(fdr_adjust_pvalue)
        
        cor_vec=seq(0,1,0.1)
        pvalues_vec<-corPvalueStudent(cor=cor_vec,nSamples=num_samps_groups[[i]])
        abs.cor.thresh=max(abs.cor.thresh,min(cor_vec[which(pvalues_vec<pvalue.thresh)],na.rm=TRUE),na.rm=TRUE)
        
        #dim(complete_pearsonqvalue_mat)<-dim(simmat)
        
        simmat[which(abs(simmat)<abs.cor.thresh)]<-0
        #simmat[which(complete_pearsonqvalue_mat>fdrthresh)]<-0
        
        print(dim(simmat))
        print(power_val)
        ADJdataOne<-adjacency.fromSimilarity(similarity=simmat,power=power_val)
        
        #degree_list[[i]]<-softConnectivity(datExpr=data_matrix_list[[i]],power=power_val,minNSamples=2)
        #
    
        return(ADJdataOne)
    })

    diffrank_list<-lapply(2:length(class_labels_levels),function(i){
        
        res<-diffRank(adj_mat_list[[1]],adj_mat_list[[i]])
        return(res)
    })
    
    save(diffrank_list,file="diffrank_list.Rda")
    
    diffrank_mat<-do.call("cbind",diffrank_list)
    
    
    
    if(length(class_labels_levels)>2){
        
        diffrank_res<-apply(diffrank_mat,1,function(x){max(x,na.rm=TRUE)})
        
        diffrank_res_rank<-rank(-1*diffrank_res)
        colnames_vec<-c("mz","time","sigstatus","CentralityScore","DiffRank",paste("DiffRank.",class_labels_levels[1],"vs",class_labels_levels,sep=""))
        colors_met_all<-cbind(colors_met,diffrank_res,diffrank_res_rank,diffrank_mat)
        
        
    }else{
        diffrank_res<-diffrank_mat
        
        diffrank_res_rank<-rank(-1*diffrank_res)
        
        colnames_vec<-c("mz","time","sigstatus","CentralityScore","DiffRank")
        colors_met_all<-cbind(colors_met,diffrank_res,diffrank_res_rank)
        
        
    }
    
   
   
    colors_met_all<-as.data.frame(colors_met_all)
    colnames(colors_met_all)<-colnames_vec
    
    print(head(colors_met_all))
    #print(table(MET[[1]]$validColors))
    #print(table(classAmoduleColors))
    write.table(colors_met_all,file="Tables/DiNA_eval_allfeats.txt",sep="\t",row.names=FALSE)
    
    sub_colors_met<-{}
    if(is.na(sigfeats)==FALSE){
        sub_colors_met<-colors_met_all[sigfeatsind,]
        write.table(sub_colors_met,file="Tables/DiNA_eval_selectfeats.txt",sep="\t",row.names=FALSE)
    }
    
    ##save(MET,file="MET.Rda")
    
    
    
    return(list(all=colors_met_all,sigfeats=sub_colors_met))
}



do_cv<-function(v,x,y,kname="radial",errortype="CV",conflevel=95,classifier="SVM",seednum=555){

	    
    num_samp=dim(x)[1]

    if(length(num_samp)<1){
		num_samp<-length(x)
		x<-as.data.frame(x)
    }    
    y<-as.data.frame(y)
    
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
    for ( i in 1:v)
    {
        g<-group[(ind[i]+1):ind[i+1]]
        temptest<-x[g,]
        temptrain <-x[-g,]
        tempclass <-y[-g,]
        testclass<-y[g,]
        
        
        if(classifier=="SVM"){
        mod_cv <- svm(x=temptrain,y=tempclass, type="C",kernel=kname)
        # predfit<-predict(mod_cv,temptest)
         
         if(v==num_samp){
             
             predfit<-predict(mod_cv,(temptest))
             #  predfit<-predict(mod_cv,t(temptest))
         }else{
             predfit<-predict(mod_cv,(temptest))
         }
         

        }else{
            if(classifier=="LR"){
                
                Class<-tempclass
                
                dtemp<-cbind(Class,temptrain)
                mod_cv<-glm(as.factor(Class)~.,data=dtemp,family=binomial(logit))
                
                if(v==num_samp){
                    
                    #predfit<-predict(mod_cv,t(temptest))
                    
                    predfit<-predict(mod_cv,t(temptest),type="response")
                }else{
                    # predfit<-predict(mod_cv,(temptest))
                    
                    predfit<-predict(mod_cv,temptest,type="response")
                }
                
                
                #mod_cv<-glm.fit(x=temptrain,y=tempclass,family="binomial")
                #predfit<-predict(mod_cv,temptest,type="response")
                 
                 #print(predfit)
                 
                 predfit <- ifelse(predfit > 0.5,1,0)
                

            }else{
                
               
                    if(classifier=="RF"){
                        
                        
                        Class<-tempclass
                        
                        d1<-cbind(Class,temptrain)
                        mod_cv<-randomForest(as.factor(Class)~.,data=d1)
                        #mod_cv<-randomForest(x=temptrain,y=tempclass)
                        
                        # predfit<-predict(mod_cv,temptest)
                        if(v==num_samp){
                            
                            predfit<-predict(mod_cv,t(temptest))
                        }else{
                            predfit<-predict(mod_cv,(temptest))
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
        
        if(errortype=="CV"){
            svm_acc[i]<-(totacc*100)
        }else{
            if(errortype=="AUC"){
                svm_acc[i]<-(auc_acc*100)
            }else{
                svm_acc[i]<-(beracc*100)
            }
        }
        
    }
    avg_acc <-mean(svm_acc,na.rm=TRUE)
    sd_acc<-sd(svm_acc,na.rm=TRUE)
    
    #limit<-avg_acc-(sd.error*(avg_acc) # 1 sd criterion
    #print(avg_acc)
    #print(sd_acc)
    
    #return(list(error=avg_acc,sderror=sd.error))
    probval<-(1-(conflevel*0.01))/2
    probval<-1-probval
    #print(probval)
    error <- qnorm(probval)*sd_acc/sqrt(length(svm_acc))
    
    leftconfint<-avg_acc-error
    rightconfint<-avg_acc+error
    
    #print("done")
    return(list(avg_acc=avg_acc,sd_acc=sd_acc, acc_each_fold=svm_acc,confint=c(leftconfint,rightconfint)))
    #return(list(num=best_k,error=min_err, avg=avg_acc))
}

get_classification.accuracy.child<-function(temptrain,tempclass,kernelname="radial",errortype="AUC",classifier="svm",num_samp,temptest=NA,testclass=NA,numfolds=1,plotroc=FALSE)
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
            
            #save(dtemp,file="dtemp.Rda")
            #save(temptest,file="temptest.Rda")
            mod_cv<-glm(as.factor(Class)~.,data=dtemp,family=binomial(logit))
            
            if(v==num_samp){
                
                predfit<-predict(mod_cv,(temptest),type="response")
            }else{
                
                
                predfit<-predict(mod_cv,temptest,type="response")
            }
            
            
            predfit <- ifelse(predfit > 0.5,1,0)
            
            testclass<-as.numeric(testclass)-1
            # print(predfit)
            #print(testclass)
            
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
                         
                         mod_cv<-pls.lda(Xtrain=temptrain,Ytrain=tempclass,ncomp=opt_comp,nruncv=v,Xtest=temptest)
                         predfit<-mod_cv$predclass
                         
                        

                     }
                }
                
                
            }
            
        }
        
    }
    
    #save(list=ls(),file="t1.Rda")
    
   
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
    
    roc_svm={}
    roc_lr={}
    if(dim(svm_table)[1]==2){
        
        if(plotroc==TRUE)
        {
            
            #save(list=ls(),file="t.Rda")
        temptrain<-t(temptrain)
        
        temptest<-t(temptest)
        
       
            roc_svm<-get_roc(dataA=temptrain,classlabels=tempclass,classifier="svm",kname="radial",rocfeatlist=seq(1,dim(temptrain)[2]),rocfeatincrement=FALSE,testset=temptest,testclasslabels=testclass,mainlabel="Test set (ROC);\n using SVM",mz_names=rownames(temptrain))
    
    roc_lr<-get_roc(dataA=temptrain,classlabels=tempclass,classifier="logit",kname="radial",rocfeatlist=seq(1,dim(temptrain)[2]),rocfeatincrement=FALSE,testset=temptest,testclasslabels=testclass,mainlabel="Test set (ROC);\n using LR",mz_names=rownames(temptrain))
        }
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
    return(list(classification_acc=svm_acc,classification_model=mod_cv,confusion_matrix=svm_table,roc_svm=roc_svm,roc_lr=roc_lr))
}

##This function calculates k-fold cross-validation classification accuracy on training set and overall classification accuracy for test set using balanced error rate,
##total misclassification rate, and AUC criteria.
get_classification.accuracy<-function(kfold,featuretable,classlabels,kernelname="radial",errortype="AUC",conflevel=95,classifier="svm",seednum=555,testfeaturetable=NA,testclasslabels=NA){
    
    
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
    
    error <- qnorm(probval)*sd_acc/sqrt(length(svm_acc))
    avg_acc<-round(avg_acc,2)
    leftconfint<-avg_acc-error
    rightconfint<-avg_acc+error
    test_acc<-NA
    test_confusion_matrix<-NA
    
    leftconfint<-round(leftconfint,2)
    rightconfint<-round(rightconfint,2)
    print(paste("Training set ", kfold,"-fold CV ",errortype," ",classifier," classification accuracy (%):",avg_acc,sep=""))
    print(paste("Training set ", kfold,"-fold CV ",errortype," ",classifier," classification accuracy ", conflevel,"% confidence interval:(",leftconfint,",",rightconfint,")",sep=""))
    
    
    x<-as.data.frame(x)
    y<-y[,1]
    train_res<-get_classification.accuracy.child(temptrain=x,tempclass=y,kernelname=kernelname,errortype=errortype,classifier=classifier,num_samp=num_samp,temptest=x,testclass=y)
    print(paste("Training set ",errortype," ",classifier," classification accuracy (%):",round(train_res$classification_acc,2),sep=""))
    
    mod_cv<-train_res$classification_model
   
   test_res={}
   #evaluate test set accuracy
    if(is.na(testfeaturetable)==FALSE){
        
        if(is.na(testclasslabels)==FALSE){
            
            testfeaturetable<-as.data.frame(testfeaturetable)
            
            #save(list=ls(),file="t2.Rda")
            print(dim(testfeaturetable))
            print(length(testclasslabels))
            
            
            
              test_res<-get_classification.accuracy.child(temptrain=x,tempclass=y,kernelname=kernelname,errortype=errortype,classifier=classifier,num_samp=num_samp,temptest=testfeaturetable,testclass=testclasslabels,plotroc=TRUE)
              
            test_acc<-test_res$classification_acc
            test_acc<-round(test_acc,2)
            test_pred_table<-test_res$confusion_matrix
            
            print(paste("Test set ", errortype," ",classifier," classification accuracy (%):",test_acc,sep=""))
            print(paste("Test set confusion matrix using ",classifier,sep=""))
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
    return(list(avg.train.cv.acc=avg_acc,sd.train.cv.acc=sd_acc,train.cv.acc.each.fold=svm_acc,train.cv.acc.confint=c(leftconfint,rightconfint),test.acc=test_acc,test.confusion.matrix=test_confusion_matrix,test_res=test_res))
    }
    
    
}


do_wgcna<-function(feature_table_file=NA,class_labels_file=NA,X=NA,Y=NA,sigfeats=NA){

if(is.na(X)==TRUE){
                        data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE)

                        }else{

                                data_matrix<-X
                        }
                        if(is.na(Y)==TRUE){
                        classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)

                        }else{

                                classlabels<-Y
                        }

if(dim(classlabels)[2]>2){
    classgroup<-paste(classlabels[,2],":",classlabels[,3],sep="") #classlabels[,2]:classlabels[,3]
}else{

        classgroup<-classlabels[,2]
}
classlabels<-as.data.frame(classlabels)

class_labels_levels<-levels(as.factor(classgroup))

rnames<-paste(sprintf("%.4f",data_matrix$mz),data_matrix$time,sep="_")

data_matrix_orig<-data_matrix
data_matrix<-data_matrix[,-c(1:2)]
data_matrix<-na.omit(data_matrix)

rnamesAB<-gsub(pattern="NA_NA",replacement=NA,x=rnames)
rnamesAB<-na.omit(rnamesAB)

nSets = length(class_labels_levels);
multiExpr = vector(mode = "list", length = nSets)
data_matrix_list<-new("list")
num_samps_groups<-new("list")
for(i in 1:length(class_labels_levels)){

	data_matrix_list[[i]]<-t(data_matrix[,which(classgroup==class_labels_levels[i])])
	num_samps_groups[[i]]<-dim(data_matrix_list[[i]])[1]
    #print(dim(data_matrix_list[[i]]))
	multiExpr[[i]]<-list(data = as.data.frame(data_matrix_list[[i]]));
	rownames(multiExpr[[i]]$data)=c(paste(rep(class_labels_levels[i],num_samps_groups[[i]]),seq(1,num_samps_groups[[i]]),sep=""))
}

data_matrix_all<-t(data_matrix)


# We work with two sets:
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = as.character(class_labels_levels) #c("Slow", "Rapid")
shortLabels = as.character(class_labels_levels) #c("Slow", "Rapid")

##save(class_labels_levels,file="class_labels_levels.Rda")
##save(multiExpr,file="multiExpr.Rda")

exprSize = checkSets(multiExpr)


sampleTrees = list()
for (set in 1:nSets)
{
	sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}


#pdf(file = "SampleClustering.pdf", width = 12, height = 12);


for (set in 1:nSets)
{
# Find clusters cut by the line
#labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
# Keep the largest one (labeled by the number 1)
#keep = (labels==1)
multiExpr[[set]]$data = multiExpr[[set]]$data
}
collectGarbage();


# Check the size of the leftover data
exprSize = checkSets(multiExpr)
exprSize;

traitData<-as.data.frame(classlabels)

dim(traitData)
names(traitData)

# See how big the traits are and what are the trait and sample names
# Form a multi-set structure that will hold the clinical traits.
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets)
{

	Traits[[set]] = list(data = data_matrix_list[[i]] );
	rownames(Traits[[1]]$data) = rownames(data_matrix_list[[i]]);
	#Traits[[2]] = list(data = dataexprB );
	#rownames(Traits[[2]]$data) = rownames(dataexprB);
}
collectGarbage();
# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;

##save(multiExpr, Traits, nGenes, nSamples, file = "Consensus-dataInput.RData");






powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(data=data_matrix_all, dataIsExpr=TRUE,powerVector = powers, verbose = 0)
power_val=sft$powerEstimate

if(is.na(power_val)==TRUE){
power_val=6
}

netclassA = blockwiseConsensusModules(multiExpr, power = power_val, minModuleSize = 10, deepSplit = 2,
pamRespectsDendro = FALSE,
mergeCutHeight = 0.2, numericLabels = TRUE,saveTOMs = FALSE, verbose = 0,saveIndividualTOMs=FALSE,useDiskCache=FALSE)


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(netclassA$colors)
if(FALSE){
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(netclassA$dendrograms[[1]], mergedColors[netclassA$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
}

classAmoduleLabels = netclassA$colors
classAmoduleColors = labels2colors(netclassA$colors)
classAMEs = netclassA$MEs;
classAgeneTree = netclassA$dendrograms[[1]];




# Recalculate consMEs to give them color names
consMEsC = multiSetMEs(multiExpr, universalColors = classAmoduleColors);

MET = consensusOrderMEs(consMEsC);

#save(MET,file="MET.Rda")
#save(setLabels,file="setLabels.Rda")
pdf("Module_preservation_analysis.pdf",width=10,height=8)
#tiff("module_preservation.tiff",width=2000,height=2000)
#sizeGrWindow(8,10);
#par(cex = 1)
#plotEigengeneNetworks_custom(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),zlimPreservation = c(0.5, 1), plotPreservation = "standard",plotDendrograms=FALSE)
#dev.off()

#sizeGrWindow(8,10);
#tiff("module_preservation.tiff",width=2000,height=2000)
par(cex = 0.9)
plotEigengeneNetworks_custom(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),zlimPreservation = c(0.5, 1), plotPreservation = "standard",plotDendrograms=FALSE)
try(dev.off(),silent=TRUE)
graphics.off()
graphics.off()


colors_met<-cbind(data_matrix_orig[,c(1:2)],MET[[1]]$validColors)

colors_met<-as.data.frame(colors_met)

colnames(colors_met)<-c("mz","time","Module")

#print(table(MET[[1]]$validColors))
#print(table(classAmoduleColors))
write.table(colors_met,file="Allmetabolite_modules.txt",sep="\t",row.names=FALSE)

if(is.na(sigfeats)==FALSE){
sub_colors_met<-colors_met[which(data_matrix_orig$mz%in%sigfeats$mz),]
write.table(sub_colors_met,file="Sigmetabolite_modules.txt",sep="\t",row.names=FALSE)
}

##save(MET,file="MET.Rda")


return(list(MET=MET,setLabels=setLabels))
}


get_volcanoplots<-function(xvec,yvec,up_or_down,maintext="",ythresh=0.05,y2thresh=NA,ylab,xlab,colorvec=c("darkblue","red3"),col_seq=c("brown","chocolate3","orange3","coral","pink","skyblue","blue","darkblue","purple","violet"),xincrement=1,yincrement=1,xthresh=1,pchvec=c(21,21),background.points.col="gray50",bad.feature.index=NA){
    
    # #save(list=ls(),file="volcano.Rda")
    d4<-xvec
    min_val<-round(min(d4,na.rm=TRUE)+0.5)
    max_val<-round(max(d4,na.rm=TRUE)+0.5)
    
    windowsize=xincrement
    
    d4<-as.vector(d4)
    
    logp<-as.vector(yvec)
    
    if(is.na(up_or_down)==TRUE){
        up_or_down<-rep(1,length(yvec))
    }
    
    plot(d4,logp,xaxt="n",ylab=ylab,xlab=xlab,xaxt="n",yaxt="n",cex=0.4,cex.main=0.8,main=maintext)
    axis(1, at=seq(min_val , max_val, by=xincrement) , las=2)
    axis(2, at=seq(0 , (max(logp)+2), by=yincrement) , las=2)


    points(d4,logp,col=background.points.col,cex=0.4,bg=background.points.col,pch=21)
    points(d4,logp,col=background.points.col,cex=0.4,bg=background.points.col,pch=21)
    
    
    goodip<-which(yvec>ythresh & abs(xvec)>xthresh)
   
   if(length(bad.feature.index)>0){
           
       if(is.na(bad.feature.index)==FALSE){
           if(length(goodip)>0){
               
               
               check_bad_feat_index<-which(goodip%in%bad.feature.index)
               
               if(length(check_bad_feat_index)>0){
                   goodip<-goodip[-check_bad_feat_index]
               }
               #goodip<-goodip[-which(goodip%in%bad.feature.index)] #goodip[-bad.feature.index]
           }
       }
   }
    
    for(i in goodip){
        if(up_or_down[i]>0){
            points(d4[i],logp[i],col=colorvec[1],cex=0.8,pch=pchvec[1],bg=colorvec[1]); points(d4[i],logp[i],col=colorvec[1],cex=0.4,bg=colorvec[1])
        }else{
            
            points(d4[i],logp[i],col=colorvec[2],cex=0.8,pch=pchvec[2],bg=colorvec[2]); points(d4[i],logp[i],col=colorvec[2],cex=0.4,bg=colorvec[2])
        }
    }
    if(length(bad.feature.index)>0){
        
        if(is.na(bad.feature.index)==FALSE){
            for(i in bad.feature.index){
                points(d4[i],logp[i],col=background.points.col,cex=0.4,bg=background.points.col,pch=21)
                
            }
        }
    }
    
    if(length(goodip)>0){
        
        
        abline(v=(-1)*xthresh,col="gray8",lty=2,lwd=0.8)
        abline(v=xthresh,col="gray8",lty=2,lwd=0.8)
        
        abline(h=ythresh,col="gray8",lty=2,lwd=0.8)
        
        if(is.na(y2thresh)==FALSE){
            abline(h=y2thresh,col="gray8",lty=2,lwd=0.8)
            
        }
        
    }
 
 
}


get_manhattanplots<-function(xvec,yvec,up_or_down,maintext="",ythresh=0.05,y2thresh=NA,ylab,xlab,colorvec=c("darkblue","red3"),col_seq=c("brown","chocolate3","orange3","coral","pink","skyblue","blue","darkblue","purple","violet"),xincrement=100,yincrement=1,pchvec=c(21,21),background.points.col="black",bad.feature.index=NA){
    
    d4<-xvec
    min_val<-min(c(0,d4),na.rm=TRUE)[1]
    max_val<-max(d4,na.rm=TRUE)[1]
    
    windowsize=xincrement
    
    d4<-as.vector(d4)
    pvalues<-as.vector(yvec)
    
    logp<-as.vector(yvec)
    
    if(is.na(up_or_down)==TRUE){
        up_or_down<-rep(1,length(yvec))
    }
    
    max_yval<-max(yvec,na.rm=TRUE)[1]+1.96*(sd(yvec,na.rm=TRUE)/(sqrt(length(yvec))))
    
    #,ylim=range(pretty(c(0,max(logp)+2))),xlim=range(pretty(c(min_val,max_val)))
    
    plot(d4,logp,xaxt="n",ylab=ylab,xlab=xlab,xaxt="n",yaxt="n",cex=0.4,cex.main=0.8,main=maintext,ylim=range(pretty(c(0,max(logp)))))
    
    axis(1, at=seq(min_val , max_val, by=xincrement) , las=2)
    axis(2, at=seq(0 , (max(logp)+2), by=yincrement) , las=2)
    #col_seq<-c("brown","red","orange3","coral","pink","skyblue","blue","darkblue","purple","violet")
    
    if(length(col_seq)>1){
        
            s1<-seq(windowsize,max_val,windowsize)
            points(d4[which(d4>=0 & d4<=windowsize)],logp[which(d4>=0 & d4<=windowsize)],col=col_seq[1],cex=0.4,pch=21,bg=background.points.col)
            for(i in 1:(length(s1)-1))
            {
                points(d4[which(d4>s1[i] & d4<=s1[i+1])],logp[which(d4>s1[i] & d4<=s1[i+1])],col=col_seq[i+1],cex=0.4,pch=21,bg=background.points.col)
            }
    }else{
        
        #points(d4[which(d4>=0 & d4<=windowsize)],logp[which(d4>=0 & d4<=windowsize)],col="black",bg=background.points.col,cex=0.4,pch=21)
        
        points(d4,logp,col=background.points.col,bg=background.points.col,cex=0.4,pch=21)
    }
    
    if(is.na(y2thresh)==TRUE){
        
        goodip<-which(yvec>ythresh)
    }else{
        
        goodip<-which(yvec>y2thresh)
        
    }
    
    if(length(bad.feature.index)>0){
        if(is.na(bad.feature.index)==FALSE){
            if(length(goodip)>0){
                
                
               
                check_bad_feat_index<-which(goodip%in%bad.feature.index)
                
                if(length(check_bad_feat_index)>0){
                    goodip<-goodip[-check_bad_feat_index]
                }
                
            }
        }
    }
    
    for(i in goodip){
        if(up_or_down[i]>0){
            points(d4[i],logp[i],col=colorvec[1],cex=0.8,pch=pchvec[1],bg=colorvec[1]); points(d4[i],logp[i],col=colorvec[1],cex=0.2,bg=colorvec[1])
        }else{
            
            points(d4[i],logp[i],col=colorvec[2],cex=0.8,pch=pchvec[2],bg=colorvec[2]); points(d4[i],logp[i],col=colorvec[2],cex=0.2,bg=colorvec[2])
        }
    }
    
    if(length(bad.feature.index)>0){
        if(is.na(bad.feature.index)==FALSE){
            for(i in bad.feature.index){
                
                points(d4[i],logp[i],col=background.points.col,cex=0.4,pch=pchvec[1],bg=background.points.col); #points(d4[i],logp[i],col=colorvec[1],cex=0.2,bg=)
                }
        }
    }
        
    if(length(goodip)>0){
                                #hfdrfdrthresh<-logp[which(logp==min(logp[which(yvec>ythresh)],na.rm=TRUE))]
                                #abline(h=hfdrfdrthresh,col="gray8",lty=2,lwd=2)
                                
                                abline(h=ythresh,col="gray8",lty=2,lwd=0.8)
                                if(is.na(y2thresh)==FALSE){
                                    abline(h=y2thresh,col="gray8",lty=2,lwd=0.8)
                                }
    }
   
}


#function to generate ROC curves using SVM or logistic regression classifiers
get_roc<-function(dataA,classlabels,classifier="svm",kname="radial",rocfeatlist=seq(2,10,1),rocfeatincrement=TRUE,testset=NA,testclasslabels=NA,mainlabel=NA,col_lab=NA,legend=TRUE,newdevice=FALSE,mz_names=NA){
    
    
    
    d1<-dataA
    rm(dataA)
    
    
    if(newdevice==TRUE){
        
        pdf("ROC.pdf")
    }
    cnames<-colnames(d1)
    
    
    classlabels<-as.data.frame(classlabels)
    testclasslabels<-as.data.frame(testclasslabels)
    class_inf<-classlabels
    
    if(is.na(mz_names)==TRUE){
        cnames[1]<-"mz"
        cnames[2]<-"time"
        colnames(d1)<-cnames
        d1<-as.data.frame(d1)
        
        #d1<-unique(d1)
        
        d2<-t(d1[,-c(1:2)])
        
        cnames<-colnames(d1)
    
        if(is.na(testset)==TRUE){
            testset<-d1
            testclasslabels<-classlabels
        }
        #   testset<-unique(testset)
        testset<-t(testset[,-c(1:2)])
        mz_names<-paste(d1$mz,d1$time,sep="_")
    }else{
        
       
        d1<-as.data.frame(d1)
        #d1<-unique(d1)
        
        d2<-t(d1)
        cnames<-colnames(d1)
        
        mz_names<-unique(mz_names)
        if(is.na(testset)==TRUE){
            testset<-d1
            testclasslabels<-classlabels
        }
        # testset<-unique(testset)
        testset<-t(testset)
        
    }
    
   
    featlist<-rocfeatlist
    featincrement<-rocfeatincrement
    
   
    
    class_vec<-as.character(class_inf[,1])
    
    class_levels<-levels(as.factor(class_vec))
    for(c in 1:length(class_levels)){
        
        classind<-(c-1)
        class_vec<-replace(class_vec,which(class_vec==class_levels[c]),classind)
    }
    
    #class_vec<-replace(class_vec,which(class_vec==class_levels[2]),1)
    
    class_vec<-as.numeric(as.character(class_vec))
    
    
    d3<-cbind(class_vec,d2)
    
    
    mz_names<-colnames(d2)
    # colnames(testset)<-as.character(mz_names)
    mz_names<-c("Class",mz_names)
    colnames(d3)<-as.character(mz_names)
    
    
    
    d3<-as.data.frame(d3)


    #featlist<-unique(featlist)
   
    mod.lab<-{}

    if(is.na(col_lab)==TRUE){
    if(length(featlist)>1){
    col_lab<-palette(rainbow(length(featlist)))
    }else{
        col_lab<-c("blue")
    }
    }
    
    extra_index<-which(featlist>(dim(d1)[1]+1))
    
    if(length(extra_index)>0){
        featlist<-featlist[-extra_index]
    }
    #featlist<-unique(featlist)
   
   if(is.na(mainlabel)==TRUE){
       
       
       if(classifier=="logitreg"){
           mainlab<-"ROC curves using logistic regression"
       }else{
           
           mainlab<-"ROC curves using SVM"
       }
       
   }else{
       
       mainlab=mainlabel
   }
   

    if(featincrement==TRUE){
        for(n in 1:length(featlist)){
            
           
            num_select<-featlist[n]
            
            
            
            if(classifier=="logitreg"){
                
                d4<-as.data.frame(d3[,c(1:num_select)])
                
                model1 <- glm(d4$Class~., data=d4, family=binomial)
                
                testset<-as.data.frame(testset)
               
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
            #  #save(roc_res,file="ROCres.rda")
            x1<-seq(0,1,0.01)
            y1<-x1
            p1<-performance(pred1,"auc")
            mod.lab <-c(mod.lab,paste('using top ',(num_select-1),' m/z features: AUC ',round(p1@y.values[[1]],2),sep=""))
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
    }else{
        
         num_select<-length(featlist)
         
        if(classifier=="logit"){
            
            
            d4<-as.data.frame(d3[,c(1:num_select)])
            testset<-as.data.frame(testset)
            
            model1 <- glm(d4$Class~., data=d4, family=binomial)
            
            #pred1 <- prediction(fitted(model1), d3$Class)
            pred<-predict(model1,testset)
            
            pred1 <- ROCR::prediction(pred, testclasslabels)
            
            
        }else{
            
           
          
            model1 <- svm(as.factor(d3$Class)~., data=d3[,c(featlist)], type="C",probability=TRUE,kernel=kname)
            
            pred<-predict(model1,testset,probability=TRUE,decision.values=TRUE,type="prob")
            pred1 <- ROCR::prediction(attributes(pred)$probabilities[,2], testclasslabels)
            
            
            
        }
        
        n=1
       
        stats1a <- performance(pred1, 'tpr', 'fpr')
        roc_res<-cbind(stats1a@x.values[[1]],stats1a@y.values[[1]])
        colnames(roc_res)<-c(stats1a@x.name,stats1a@y.name)
        # #save(roc_res,file="ROCres.rda")
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
            try(dev.off(),silent=TRUE)
    }
    
    return(p1@y.values[[1]])
}




svm_cv<-function(v,x,y,kname="radial",errortype="CV",conflevel=95,seednum=555){

num_samp=dim(x)[1]

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
for ( i in 1:v)
{
g<-group[(ind[i]+1):ind[i+1]]
temptest<-x[g,]
temptrain <-x[-g,]
tempclass <-y[-g,]
testclass<-y[g,]


mod_cv <- svm(x=temptrain,y=tempclass, type="C",kernel=kname)

if(v==num_samp){
    
    predfit<-predict(mod_cv,t(temptest))
}else{
    predfit<-predict(mod_cv,(temptest))
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

if(errortype=="CV"){
	svm_acc[i]<-(totacc*100)	
}else{
if(errortype=="AUC"){
	svm_acc[i]<-(auc_acc*100)
}else{
svm_acc[i]<-(beracc*100)
}
}

}
avg_acc <-mean(svm_acc,na.rm=TRUE)
sd_acc<-sd(svm_acc,na.rm=TRUE)

#limit<-avg_acc-(sd.error*(avg_acc) # 1 sd criterion
#print(avg_acc)
#print(sd_acc)

#return(list(error=avg_acc,sderror=sd.error))
probval<-(1-(conflevel*0.01))/2
probval<-1-probval
#print(probval)
error <- qnorm(probval)*sd_acc/sqrt(length(svm_acc))

leftconfint<-avg_acc-error
rightconfint<-avg_acc+error

#print("done")
return(list(avg_acc=avg_acc,sd_acc=sd_acc, acc_each_fold=svm_acc,confint=c(leftconfint,rightconfint)))
#return(list(num=best_k,error=min_err, avg=avg_acc))
}

plsda_cv<-function(v,x,y,ncomp,errortype="total",conflevel=99){

num_samp=dim(x)[1]

num_datasets= floor(num_samp)
n1<-floor(num_samp/v)
n2<-num_samp-n1*v
n3<-v-n2

ind<-rep(c(n1,n1+1),c(n3,n2))
ind<-diffinv(ind)
min_err=1
best_k=1

set.seed(555)
group<-sample(1:num_samp,num_samp, replace=FALSE)


itr=0
#plsda_error <- matrix(0,v)  # we set K=30 before, it can be changed to any number<100.
plsda_error<-rep(0,v)
for ( i in 1:v)
{
g<-group[(ind[i]+1):ind[i+1]]
temptest<-x[g,]
temptrain <-x[-g,]
tempclass <-y[-g]
testclass<-y[g]

#temptest<-as.data.frame(temptest)
#temptrain<-as.matrix(temptrain)



#print(dim(temptrain))
#print(dim(temptest))

#plsda_cv <- plsda(x=temptrain,y=tempclass, type="C",kernel=kname) 

#opt_comp<-pls.lda.cv(Xtrain=temptrain, Ytrain=tempclass,  ncomp=c(1:10), nruncv=10, alpha=2/3, priors=NULL)
predfit<-pls.lda(Xtrain=temptrain,Ytrain=tempclass,ncomp=ncomp,nruncv=v,Xtest=temptest)
#predfit<-predict(plsda_pred,temptest)

#print(length(which(plsda_pred$predclass==testclass)))
svm_table<-table(predfit$predclass,testclass)

class_names<-rownames(svm_table)
#print(testclass)
#print(predfit$predclass)
predfit<-predfit$predclass
totacc<-length(which(predfit==testclass))/length(testclass)

beracc<-{}
auc_acc<-{}
for(c in 1:dim(svm_table)[1]){
testclass_ind<-which(testclass==class_names[c])
beracc<-c(beracc,length(which(predfit[testclass_ind]==testclass[testclass_ind]))/length(testclass_ind))


	if(errortype=="AUC"){
	pred_acc<-multiclass.roc(testclass,as.numeric(predfit),levels=levels(as.factor(y)))
	pred_acc_orig<-pred_acc$auc[1]
	auc_acc<-c(auc_acc,pred_acc_orig)
	}

}

beracc<-mean(beracc,na.rm=TRUE)

if(errortype=="CV"){
	plsda_error[i]<-(totacc*100)	
}else{
if(errortype=="AUC"){
	plsda_error[i]<-(auc_acc*100)
}else{
plsda_error[i]<-(beracc*100)
}
}



}
avgacc <-mean(plsda_error)
sdacc<-sd(plsda_error)

probval<-(1-(conflevel*0.01))/2
probval<-1-probval

error <- qnorm(probval)*sdacc/sqrt(length(y))

leftconfint<-avgacc-error
rightconfint<-avgacc+error



return(list(mean_acc=avgacc,sd_acc=sdacc, acc_each_fold=plsda_error,confint=c(leftconfint,rightconfint)))


}




do_pamr<-function(X,Y,fdrthresh=0.1,nperms=100,pamr.threshold.select.max=FALSE,kfold=10){
	
	
	library(pamr)
	
	#d1<-list(x=X,y=Y)
	
	d1 <- list(x=as.matrix(X),y=factor(Y[,1]), geneid=as.character(1:nrow(X)),
      genenames=paste("g",as.character(1:nrow(X)),sep=""))
     #save(d1,file="d1.Rda")
	p1<-pamr.train(d1)

    set.seed(999)
	p2<-pamr.cv(data=d1,fit=p1,nfold=kfold)
	
    threshold_CVerror_matrix<-cbind(p2$threshold,p2$error)
    threshold_CVerror_matrix<-as.data.frame(threshold_CVerror_matrix)
    colnames(threshold_CVerror_matrix)<-c("Threshold","error")
    
    threshold_value<-max(threshold_CVerror_matrix[which(threshold_CVerror_matrix[,2]==min(threshold_CVerror_matrix[,2])),1])[1]
    
    
    set.seed(999)
	p3<-pamr.fdr(data=d1,p1,nperms=nperms)
	
    selected_feature_index<-{}
    max.discore<-{}
    
    #use median FDR
    if(length(which(p3$results[,4]<fdrthresh))>0){
    pamr_fdr_filt<-p3$results[which(p3$results[,4]<fdrthresh),]
    
    if(length(nrow(pamr_fdr_filt))>0){
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
    
    }
    
    }else{
        p4<-{}
        selected_feature_index<-{}
        max.discore<-{}
    }
    
    
    pall<-pamr.listgenes(fit=p1,data=d1,threshold=0)
    
    ##save(pall,file="pall.Rda")
    
    pall<-as.data.frame(pall)
    


    discore_matrix_all<-pall #[,-c(1)]
    discore_matrix_all<-apply(discore_matrix_all,2,as.numeric)
    discore_matrix_all<-as.data.frame(discore_matrix_all)
    discore_matrix_all<-discore_matrix_all[order(discore_matrix_all$id),]

    
    discore_matrix_all<-discore_matrix_all[,-c(1)]
    if(nrow(discore_matrix_all)>1){
        discore_matrix_all<-apply(discore_matrix_all,2,as.numeric)
        abs.discore_matrix_all<-abs(discore_matrix_all)
        max.discore.all<-apply(abs.discore_matrix_all,1,max)
        
        max.discore.all.thresh<-min(max.discore.all[selected_feature_index],na.rm=TRUE)
        
        
    }else{
        discore_matrix_all_1<-unlist(discore_matrix_all)
        discore_matrix_all<-as.numeric(as.character(discore_matrix_all_1))
        abs.discore_matrix_all<-abs(discore_matrix_all)
        max.discore.all<-max(abs.discore_matrix_all)
        
        max.discore.all.thresh<-min(max.discore.all[selected_feature_index],na.rm=TRUE)
        
    }
    
    
    # #save(list=ls(),file="debug.Rda")
   return(list("feature.list"=selected_feature_index,"max.discore.sigfeats"=max.discore,"pam_train_model"=p1,"pam_toplist"=p4,"max.discore.allfeats"=max.discore.all,"threshold_value"=threshold_value,"max.discore.all.thresh"=max.discore.all.thresh))
}

do_plsda<-function(X,Y,oscmode="pls",numcomp=3,kfold=10,evalmethod="CV",keepX=15,sparseselect=FALSE,analysismode="classification",vip.thresh=1,sample.col.opt="default",sample.col.vec=c("red","green","blue","purple"),scoreplot_legend=TRUE,feat_names=NA,pairedanalysis=FALSE,optselect=FALSE,class_labels_levels_main=NA,legendlocation="bottomleft",plotindiv=TRUE,pls.vip.selection="max",output.device.type="pdf",plots.res=600,plots.width=8,plots.height=8,plots.type="cairo")
{
    repeatmeasures=pairedanalysis
    
    
    if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/PLS_performance_plots.pdf"
        pdf(temp_filename_1)
        
        
    }
    
    
    
    num_var<-dim(X)[1]

    if(keepX>num_var){
		keepX=num_var
    }
    
   
    
    X<-t(X)
    
 
    Y<-as.data.frame(Y)

	classlabels<-Y    

    


    if(pairedanalysis==FALSE){
            Yclass<-Y[,1]
            Y<-as.numeric(Y[,1])
        
    }else{
        
            if(dim(Y)[2]>2){
                if(analysismode=="classification"){
                    
                    Yclass<-as.factor(Y[,2]):as.factor(Y[,3])
                }else{
                Yclass<-Y[,2] #:Y[,3]
                }
                Y<-as.numeric(Yclass)
            }else{
                Yclass<-Y[,2]
                Y<-as.numeric(Y[,2])
            }

	}
    
    class_labels_levels<-levels(as.factor(Yclass))
    alphacol=0.3
    if(sample.col.opt=="default"){
        
        col_vec<-c("#CC0000","#AAC000","blue","mediumpurple4","mediumpurple1","blueviolet","cornflowerblue","cyan4","skyblue",
        "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
        "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
        "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
        
    }else{
        if(sample.col.opt=="topo"){
            
            
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
                        
                        if(sample.col.opt=="colorblind"){
                            #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                            # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                            
                            if(length(class_labels_levels)<9){
                                
                                col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7", "grey57")
                                
                            }else{
                                
                                #col_vec<-colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels))
                                col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                                
                            }
                            
                            
                        }else{
                            
                            check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                            
                            if(length(check_brewer)>0){
                                sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                                col_vec <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(class_labels_levels))
                                
                            }else{
                                
                                if(sample.col.opt=="journal"){
                                    
                                    col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                    "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                    "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                    
                                    "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                    "#E64B3519","#4DBBD519","#631879E5","grey75")
                                    
                                }else{
                                    #colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                    if(length(sample.col.opt)==1){
                                        col_vec <-rep(sample.col.opt,length(class_labels_levels))
                                    }else{
                                        
                                        colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                        
                                    }
                                }
                                
                            }
                            
                        }
                    }
                    
                    
                }
                
            }
            
        }
    }
    


    
    #print("starting")
    if(dim(X)[2]>1){
        if(optselect==TRUE){
	if(analysismode=="classification")
        {
            set.seed(123)
            opt_comp<-pls.lda.cv(Xtrain=X, Ytrain=Yclass,  ncomp=c(1:numcomp), nruncv=kfold, alpha=2/3, priors=NULL)
            
          
        }else{
            if(analysismode=="regression")
            {
                
                set.seed(123)
                opt_comp<-pls.regression.cv(Xtrain=X, Ytrain=Y,  ncomp=c(1:numcomp), nruncv=kfold, alpha=2/3)
                
                
            }
            
        }
	}else{

		opt_comp<-numcomp
		keep_x_vec<-rep(keepX,opt_comp)
	}

    }
   
      suppressWarnings(dir.create("Tables"))
	if(opt_comp<2){
		opt_comp<-2
	} 
    if(oscmode=="o1pls"){
        leukemia.pls <- plsr(Y ~ X, ncomp = opt_comp, validation = "LOO")
        ww <- leukemia.pls$loading.weights[,1]
        pp <- leukemia.pls$loadings[,1]
        w.ortho <- pp - crossprod(ww, pp)/crossprod(ww) * ww
        t.ortho <- X %*% w.ortho
        
        p.ortho <- crossprod(X, t.ortho) / c(crossprod(t.ortho))
        Xcorr <- X - tcrossprod(t.ortho, p.ortho)
        
        if(analysismode=="classification")
        {
            cv_res<-plsda_cv(v=kfold,x=Xcorr,y=Yclass,ncomp=opt_comp,errortype=evalmethod)
            print(paste(kfold," CV evaluation using o1plsda",sep=""))
            print(cv_res)
        }
        
        X<-Xcorr
    }
    
    if(oscmode=="o2pls"){
        leukemia.pls <- plsr(Y ~ X, ncomp = opt_comp, validation = "LOO")
        ww <- leukemia.pls$loading.weights[,1]
        pp <- leukemia.pls$loadings[,1]
        w.ortho <- pp - crossprod(ww, pp)/crossprod(ww) * ww
        t.ortho <- X %*% w.ortho
        
        p.ortho <- crossprod(X, t.ortho) / c(crossprod(t.ortho))
        Xcorr <- X - tcrossprod(t.ortho, p.ortho)
        
        
        if(analysismode=="classification")
        {
            cv_res<-plsda_cv(v=kfold,x=Xcorr,y=Yclass,ncomp=opt_comp,errortype=evalmethod)
            
            
            print(paste(kfold," CV evaluation using o2plsda",sep=""))
            print(cv_res)
        }
        
        X<-Xcorr
    }
    
    
    
    bad_variables<-{}
    
    if(sparseselect==TRUE)
    {
        
        if(analysismode=="classification")
        {
                    if(optselect==TRUE){
                        keepx_seq<-seq(5,keepX,5)
                        best_cv_res<-c(0)
                        best_kvec<-c(5)
                        for(kvec in keepx_seq){
                                        keep_x_vec<-rep(kvec,opt_comp)

                                        if(repeatmeasures==TRUE){
                                            
                                            
                                            
                                                linn.pls <- try(multilevel(X=X, design=classlabels,ncomp = opt_comp,
                                                    keepX = keep_x_vec, method = 'splsda'),silent=TRUE)
                                                    
                                                    if(is(linn.pls,"try-error")){
                                                        
                                                       
                                                       
                                                       
                                                        linn.pls <- mixOmics::splsda(X=X,Y=classlabels[,-c(1)],ncomp=opt_comp,keepX=keep_x_vec,multilevel=classlabels[,1])
                                                        
                                                    }

                                        }else{
                                            
                                                linn.pls <- splsda(X, Yclass,ncomp=opt_comp,keepX=keep_x_vec)
                                        }


                                            #    #save(linn.pls,file="linnpls.Rda")
                                                linn.vip<-linn.pls$loadings$X
                                                #linn.vip<-linn.vip[,1]
                                                
                                                bad_variables<-linn.pls$nzv$Position
                                                
                                                good_feats<-{}
                                                for(c1 in 1:opt_comp){
                                                    good_feats<-c(good_feats,which(linn.vip[,c1]!=0))
                                                }
                                                
                                                good_feats<-unique(good_feats)
                                                
                                                
                                                
                                                if(length(good_feats)>1){
                                                    
                                                    cv_res<-try(plsda_cv(v=kfold,x=X[,good_feats],y=Yclass,ncomp=opt_comp,errortype=evalmethod),silent=TRUE)
                                                    
                                                    #cv_res<-pls.lda.cv(Xtrain=X, Ytrain=Y,  ncomp=c(1:numcomp), nruncv=kfold, alpha=2/3, priors=NULL)
                                                    
                                                    #cv_res<-cv_res$cv
                                                    print(paste(kfold," CV evaluation using spls top ",kvec," features per component",sep=""))
                                                    print(cv_res$mean_acc)
                                                    if(cv_res$mean_acc>=best_cv_res){
                                                        
                                                        best_cv_res<-cv_res$mean_acc
                                                        best_kvec<-kvec
                                                        
                                                        #best_ncomp<-cv_res$ncomp
                                                    }
                                                }else{
                                                    print("Too few variables to perfrom CV.")
                                                }
                        }
                        
                        keep_x_vec<-rep(best_kvec,opt_comp)
                    }else{

                        keep_x_vec<-rep(keepX,opt_comp)
                    }
        }
        else{
            
            keep_x_vec<-rep(keepX,opt_comp)
            
        }
        
       
        if(analysismode=="regression"){

	    if(repeatmeasures==TRUE){

                      #  #save(list=ls(),file="debugspls.Rda")

                        linn.pls <- try(multilevel(X=X, design=classlabels,ncomp = opt_comp,
                        keepX = keep_x_vec, method = 'spls'),silent=TRUE)


                        if(is(linn.pls,"try-error")){
                            
                            
                            
                            
                            linn.pls <- mixOmics::spls(X=X,Y=Y,ncomp=opt_comp,keepX=keep_x_vec,multilevel=classlabels)
                            
                        }




                }else{
                        linn.pls <- mixOmics::spls(X, Y,ncomp=opt_comp,keepX=keep_x_vec,mode="regression")
                }
	 }else{
         #analysismode classification
		if(repeatmeasures==TRUE){
            

                    linn.pls <- try(multilevel(X=X, design=classlabels,ncomp = opt_comp,
                    keepX = keep_x_vec, method = 'splsda'),silent=TRUE)

                    if(is(linn.pls,"try-error")){
                        
                        
                        linn.pls <- mixOmics::splsda(X=X,Y=classlabels[,-c(1)],ncomp=opt_comp,keepX=keep_x_vec,multilevel=classlabels[,1])
                        
                    }

                  #  #save(linn.pls,file="linnpls.Rda")



                }else{
                        linn.pls <- mixOmics::splsda(X, Y,ncomp=opt_comp,keepX=keep_x_vec)
                }
        }
        
        linn.vip<-linn.pls$loadings$X
        
        
        bad_variables<-linn.pls$nzv$Position
        
        #VIP based feature selection
        good_feats<-{}
         for(c1 in 1:opt_comp){
            good_feats<-c(good_feats,which(linn.vip[,c1]!=0))
            }
        
        good_feats<-unique(good_feats)
        
        
        
    }else{
        #PLS
        if(analysismode=="regression"){
            linn.pls <- mixOmics::pls(X, Y,ncomp=opt_comp)


        }else{
            
            if(repeatmeasures==TRUE){
                
                # print("multilevel PLS classlabels")
                #print(head(classlabels))
                
                linn.pls <- mixOmics::plsda(X, Yclass,ncomp=opt_comp,multilevel=classlabels[,1])
            }else{
                    linn.pls <- mixOmics::plsda(X, Yclass,ncomp=opt_comp)
            }
        }

        
        linn.vip<-vip(linn.pls)
   
        
        bad_variables<-linn.pls$nzv$Position
        
	
        
        good_feats<-{}
        c1<-1
        
        
        good_feats<-{}
        
        if(pls.vip.selection=="max"){
        if(opt_comp>1){
        for(c1 in 1:opt_comp){
            good_feats<-c(good_feats,which(linn.vip[,c1]>vip.thresh))
        }
	
        }else{
            good_feats<-which(linn.vip>vip.thresh)
            
        }
        
        }else{
            
            if(opt_comp>1){
                linn.vip.mean<-apply(linn.vip,1,mean)
                
                good_feats<-which(linn.vip.mean>vip.thresh)
            }else{
                good_feats<-which(linn.vip>vip.thresh)
                
            }
            
        }
        good_feats<-unique(good_feats)
        
    }
    

    v1<-{}
    cv_res<-{}
    if(length(good_feats)>1)
    {
        
        #print(paste(oscmode," PLS evaluation using selected variables",sep=""))
        
        print(paste(oscmode," PLS evaluation using all variables",sep=""))
        
        temp_d<-cbind(Y,X[,good_feats])
        temp_d<-as.data.frame(temp_d)
        
        linn.pls3 <- mixOmics::pls(X[,good_feats], Y)
        
        linn.pls2 <- mixOmics::pls(X, Y,ncomp=opt_comp)
        r2_q2_valid_res<-{}
        #print("here pls")
        ##save(linn.pls2,file="linnpls2.Rda")
        ##save(linn.pls3,file="linnpls3.Rda")
        if(length(Y)>30){
            
            if(opt_comp>1){
                
                #linn.pls2 <- pls(X, Y,ncomp=opt_comp) #pls(X, Y,ncomp=opt_comp)
                v1<-try(perf(linn.pls2,validation="loo"),silent=TRUE)
                
                if(is(v1,"try-error")){
                    
                }else{

                r2_q2_valid_res<-rbind(v1$R2,v1$Q2,v1$MSEP)
                
                if(nrow(r2_q2_valid_res)>0){
                    rownames(r2_q2_valid_res)<-c("R2","Q2","MSEP")
                }
                
                cnames_vres<-paste("PLScomp",seq(1,opt_comp),sep="")
                colnames(r2_q2_valid_res)<-cnames_vres
               
	     
		
                write.table(r2_q2_valid_res,file="Tables/pls_r2_q2_res_usingallfeatures.txt",sep="\t",row.names=TRUE)
                    if(plotindiv==TRUE){
                        w <- 0.1 #grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
                        par(omd=c(0, 1-w, 0, 1))
                barplot(r2_q2_valid_res[1:2,],beside=TRUE,main="PLS leave-one-out validation diagnostics using all features",ylab="Variation",col=c("darkgrey","lightgrey"))
                #legend("topright",c("R2","Q2"),col=c("darkgrey","lightgrey"),pch=c(20))
            
                     print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,c("R2","Q2"), col=c("darkgrey","lightgrey"),pch = c(19), pt.cex = 0.6, title = "",cex=0.8))
                     
                    }
                    
                    #linn.pls2 <- pls(X, Y,ncomp=opt_comp) #pls(X, Y,ncomp=opt_comp)
   
                }
                
            }
            
        }else{
            r2_q2_valid_res<-{}
            if(opt_comp>1){
                print("PLS loo validation diagnostics using all features")
             
                v1<-try(perf(linn.pls2,validation="loo"),silent=TRUE)
             
             if(is(v1,"try-error")){
                 
             }else{
             
                r2_q2_valid_res<-rbind(v1$R2,v1$Q2,v1$MSEP)
                
                if(nrow(r2_q2_valid_res)>0){
                    rownames(r2_q2_valid_res)<-c("R2","Q2","MSEP")
                }
                
                cnames_vres<-paste("PLScomp",seq(1,opt_comp),sep="")
                colnames(r2_q2_valid_res)<-cnames_vres
                
                if(plotindiv==TRUE){
                    
                    w <- 0.1 #grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
                    par(omd=c(0, 1-w, 0, 1))
                    barplot(r2_q2_valid_res[1:2,],beside=TRUE,main="PLS loo validation diagnostics \n using all features",ylab="Variation",col=c("darkgrey","lightgrey"))
                    #legend("topright",c("R2","Q2"),col=c("darkgrey","lightgrey"),pch=c(20))
                       print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,c("R2","Q2"), col=c("darkgrey","lightgrey"),pch = c(19), pt.cex = 0.6, title = "",cex=0.8))
                    
                }
             }
               v1<-try(perf(linn.pls3,validation="loo"),silent=TRUE)
               if(is(v1,"try-error")){
               }else{
               
               print("PLS loo validation diagnostics using selected features")

               
               r2_q2_valid_res<-rbind(v1$R2,v1$Q2,v1$MSEP)
               
               if(nrow(r2_q2_valid_res)>0){
                   rownames(r2_q2_valid_res)<-c("R2","Q2","MSEP")
               }
               
               cnames_vres<-paste("PLScomp",seq(1,dim(r2_q2_valid_res)[2]),sep="")
               colnames(r2_q2_valid_res)<-cnames_vres
		
               write.table(r2_q2_valid_res,file="Tables/pls_r2_q2_res_usingselectfeats.txt",sep="\t",row.names=TRUE)
               
               if(plotindiv==TRUE){
                   w <- 0.1 #grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
                   par(omd=c(0, 1-w, 0, 1))
                   barplot(r2_q2_valid_res[1:2,],beside=TRUE,main="PLS loo validation diagnostics \n using selected features",ylab="Variation",col=c("darkgrey","lightgrey"))
                   # legend("topright",c("R2","Q2"),col=c("darkgrey","lightgrey"),pch=c(20))
                   
                      print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,c("R2","Q2"), col=c("darkgrey","lightgrey"),pch = c(19), pt.cex = 0.6, title = "",cex=0.8))
               }
               
               }
               
            }
           
        }
        
        if(analysismode=="classification"){
           
           if(FALSE){
		if(length(Y)>10){ 
            cv_res<-try(plsda_cv(v=kfold,x=X,y=Y,ncomp=opt_comp,errortype=evalmethod),silent=TRUE)
            
            print(paste(kfold," CV ", evalmethod, " using all features: ",cv_res,sep=""))
           
            
             cv_res<-try(plsda_cv(v=kfold,x=X[,good_feats],y=Y,ncomp=opt_comp,errortype=evalmethod),silent=TRUE)
            
            print(paste(kfold," CV ", evalmethod, " using top features: ",cv_res,sep=""))
        	}
        }
	}
        
        
    }else{
        print("No variables selected.")
    }
    
    
    # print("Done with plsda")
   
   #linn.pls2 <- pls(X, Y,ncomp=opt_comp)
   
   SS<-get_plscompvar(linn.pls,nvar=dim(X)[2],opt_comp)
   
   # SS<-c(linn.pls$loadings$X)^2*colSums(linn.pls$variates$X^2)
    
    pls_var<-100*SS/(sum(SS))
    
    pls_var<-round(pls_var,2)
    
   
        
        if(output.device.type!="pdf"){
            try(dev.off(),silent=TRUE)
            
        }
        
    
    
    #barplot(pls_var,main="PLS %variation per component",cex.main=0.8)
    
    if(analysismode=="classification")
    {
        # color for plotIndiv
        col.stimu = as.numeric(Y)
    
    #print("plotting PLS")
    #print(opt_comp)
        class_labels_levels<-levels(as.factor(Yclass))
        color_vec<-col_vec #sample.col.vec #rainbow(length(class_labels_levels), start = 0, end = 0.1) #c("green","purple")
        col.stimu<-color_vec[col.stimu]
        
        class_labels_levels2<-class_labels_levels
        # pch for plots
        pch.time = rep(15, length(class_labels_levels))
        #pch.time[time == 't2'] = 4
        
        pch_vec<-seq(1,50) #c(3,5,7,9,12,13,2,17,21)
        
        Yclass2=Yclass
        
        
        samplelabels<-as.data.frame(Yclass)
        samplelabels<-as.factor(samplelabels[,1])
        l2<-levels(as.factor(samplelabels))
        col_all=topo.colors(256)
        
        t1<-table(samplelabels)
        if(is.na(class_labels_levels)==TRUE){
            
            l1<-levels(as.factor(samplelabels))
        }else{
            l1<-class_labels_levels
            
            
        }
        
        class_labels_levels<-l1
        
        col <- rep(col_vec[1:length(t1)], t1)
        #col<-rep(col_all[1:length(l1)],t1)
        ## Choose different size of points
        cex <- rep(2, length(Yclass))
        
        pch_vec<-seq(1,50) #c(3,5,7,9,12,13,2,17,21) #seq(1,50) #
        pch <- rep(15,length(Yclass))
        cex <- rep(2, length(Yclass))
        for(p1 in 1:length(l2)){
            
            pch[which(samplelabels==l2[p1])]=pch_vec[p1]
        }
        
        
        if(pairedanalysis==TRUE){
            
            if(ncol(classlabels)>2){
                
                class_labels_levels2<-levels(as.factor(classlabels[,2]):as.factor(classlabels[,3]))
                 Yclass2<-as.factor(classlabels[,2]):as.factor(classlabels[,3])
            }else{
            class_labels_levels2<-levels(as.factor(classlabels[,2]))
            Yclass2=classlabels[,2]
            }
        }
        
        #  col.stimu = as.numeric(Yclass2)
        #color_vec<-col_vec #sample.col.vec #rainbow(length(class_labels_levels), start = 0, end = 0.1) #c("green","purple")
        col.stimu<-col #color_vec[col.stimu]
        
        #pch_vec<-seq(1,length(Yclass2))
        
        for(p1 in 1:length(class_labels_levels2)){
            
                pch.time[which(Yclass2==class_labels_levels2[p1])]=pch_vec[p1]
        }
        pch.time=pch
       
        # #save(list=ls(),file="debug.Rda")
        if(plotindiv==TRUE){
	
    
    if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/PLS_pairwise_component_plots.pdf"
        
        pdf(temp_filename_1)
        #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        }
        get_plsplots(X,plsres=linn.pls,plsvar=pls_var,samplelabels=Yclass,filename=NA,ncomp=opt_comp,center=TRUE,scale=TRUE,legendcex=0.5,outloc=getwd(),col_vec=col.stimu,sample.col.opt=sample.col.opt,alphacol=0.3,legendlocation="topright",class_levels=class_labels_levels)
    #,silent=TRUE)
    
    if(output.device.type!="pdf"){
        
       try(dev.off(),silent=TRUE)
    }
    
    
    
	}
        #legend = c(class_labels_levels), cex = 0.55)
    }
    
    write.table(linn.pls$variates$X,file="Tables/pls_scores.txt",sep="\t")
    write.table(linn.pls$loadings$X,file="Tables/pls_loadings.txt",sep="\t")
    ##save(linn.pls,file="pls_res.Rda")
    return(list("model"=linn.pls,"vip_res"=linn.vip,"valid_res"=v1,"cv_res"=cv_res,"opt_comp"=opt_comp,"selected_variables"=good_feats,"bad_variables"=bad_variables))
    
}


do_plsda_rand<-function(X,Y,oscmode="pls",numcomp=3,kfold=10,evalmethod="CV",keepX=15,sparseselect=FALSE,analysismode="classification",vip.thresh=1,sample.col.opt="default",sample.col.vec=c("red","green","blue","purple"),scoreplot_legend=TRUE,feat_names=NA,pairedanalysis=FALSE,optselect=FALSE,class_labels_levels_main=NA,legendlocation="bottomleft",plotindiv=TRUE)
{
    repeatmeasures=pairedanalysis
   
   
   num_var<-dim(X)[1]
   
   if(keepX>num_var){
       keepX=num_var
   }
   
    X<-t(X)
    Y<-as.data.frame(Y)
   
    classlabels<-Y
    
    if(sample.col.opt=="default"){
        
        col_vec<-c("#CC0000","#AAC000","blue","mediumpurple4","mediumpurple1","blueviolet","cornflowerblue","cyan4","skyblue",
        "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
        "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
        "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
        
    }else{
        if(sample.col.opt=="topo"){
          
            col_vec <- topo.colors(length(class_labels_levels), alpha=alphacol)
        }else{
            if(sample.col.opt=="heat"){
            
                col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
            }else{
                if(sample.col.opt=="rainbow"){
                    
                    col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
                    
                   
                }else{
                    
                    if(sample.col.opt=="terrain"){
                        
                        
                        col_vec <- cm.colors(length(class_labels_levels), alpha=alphacol)
                    }
                    
                    
                }
                
            }
            
        }
    }
    
    if(pairedanalysis==FALSE){
        Yclass<-Y[,1]
        Y<-as.numeric(Y[,1])
    }else{
        
        if(dim(Y)[2]>2){
            if(analysismode=="classification"){
                
                Yclass<-as.factor(Y[,2]):as.factor(Y[,3])
            }else{
                Yclass<-Y[,2] #:Y[,3]
            }
            Y<-as.numeric(Yclass)
        }else{
            Yclass<-Y[,2]
            Y<-as.numeric(Y[,2])
        }
        
        
        
    }
    
   
    #Y<-as.numeric(as.factor(Y[,1]))
    #Y<-as.vector(Y)
    
    #  print(dim(X))
    #print(dim(Y))
    
    #print("starting")
    if(dim(X)[2]>1){
        if(optselect==TRUE){
            if(analysismode=="classification")
            {
                set.seed(123)
                opt_comp<-pls.lda.cv(Xtrain=X, Ytrain=Yclass,  ncomp=c(1:numcomp), nruncv=kfold, alpha=2/3, priors=NULL)
                
                #cv_res<-plsda_cv(v=kfold,x=X,y=Y,ncomp=opt_comp,errortype=evalmethod)
                
                #print(paste(kfold," CV evaluation using plsda",sep=""))
                #print(cv_res)
            }else{
                if(analysismode=="regression")
                {
                    
                    set.seed(123)
                    opt_comp<-pls.regression.cv(Xtrain=X, Ytrain=Y,  ncomp=c(1:numcomp), nruncv=kfold, alpha=2/3)
                    
                    
                }
                
            }
        }else{
            
            opt_comp<-numcomp
            keep_x_vec<-rep(keepX,opt_comp)
        }
        
    }
    
    if(opt_comp<2){
        opt_comp<-2
    }
    if(oscmode=="o1pls"){
        leukemia.pls <- plsr(Y ~ X, ncomp = opt_comp, validation = "LOO")
        ww <- leukemia.pls$loading.weights[,1]
        pp <- leukemia.pls$loadings[,1]
        w.ortho <- pp - crossprod(ww, pp)/crossprod(ww) * ww
        t.ortho <- X %*% w.ortho
        
        p.ortho <- crossprod(X, t.ortho) / c(crossprod(t.ortho))
        Xcorr <- X - tcrossprod(t.ortho, p.ortho)
        
        if(analysismode=="classification")
        {
            cv_res<-plsda_cv(v=kfold,x=Xcorr,y=Yclass,ncomp=opt_comp,errortype=evalmethod)
            print(paste(kfold," CV evaluation using o1plsda",sep=""))
            print(cv_res)
        }
        
        X<-Xcorr
    }
    
    if(oscmode=="o2pls"){
        leukemia.pls <- plsr(Y ~ X, ncomp = opt_comp, validation = "LOO")
        ww <- leukemia.pls$loading.weights[,1]
        pp <- leukemia.pls$loadings[,1]
        w.ortho <- pp - crossprod(ww, pp)/crossprod(ww) * ww
        t.ortho <- X %*% w.ortho
        
        p.ortho <- crossprod(X, t.ortho) / c(crossprod(t.ortho))
        Xcorr <- X - tcrossprod(t.ortho, p.ortho)
        
        
        if(analysismode=="classification")
        {
            cv_res<-plsda_cv(v=kfold,x=Xcorr,y=Yclass,ncomp=opt_comp,errortype=evalmethod)
            
            
            print(paste(kfold," CV evaluation using o2plsda",sep=""))
            print(cv_res)
        }
        
        X<-Xcorr
    }
    
    
    
    bad_variables<-{}
    
    if(sparseselect==TRUE)
    {
        
        if(analysismode=="classification")
        {
            if(optselect==TRUE){
                keepx_seq<-seq(5,keepX,5)
                best_cv_res<-c(0)
                best_kvec<-c(5)
                for(kvec in keepx_seq){
                    keep_x_vec<-rep(kvec,opt_comp)
                    
                    if(repeatmeasures==TRUE){
                        
                        #print("spls classlabels")
                        #   print(classlabels)
                        #print(dim(X))
                        
                        #linn.pls <- multilevel(X=X, design=classlabels,ncomp = opt_comp,
                        #keepX = keep_x_vec, method = 'splsda')
                        
                        linn.pls <- try(multilevel(X=X, design=classlabels,ncomp = opt_comp,
                        keepX = keep_x_vec, method = 'splsda'),silent=TRUE)
                        
                        if(is(linn.pls,"try-error")){
                            
                            
                            
                            
                            linn.pls <- mixOmics::splsda(X=X,Y=classlabels[,-c(1)],ncomp=opt_comp,keepX=keep_x_vec,multilevel=classlabels[,1])
                            
                        }
                        
                        #print(linn.pls)
                        
                    }else{
                        
                        linn.pls <- mixOmics::splsda(X, Yclass,ncomp=opt_comp,keepX=keep_x_vec)
                    }
                    
                    
                    
                    #linn.vip<-vip(linn.pls)
                    linn.vip<-linn.pls$loadings$X
                    #linn.vip<-linn.vip[,1]
                    
                    bad_variables<-linn.pls$nzv$Position
                    
                    good_feats<-{}
                    for(c1 in 1:opt_comp){
                        good_feats<-c(good_feats,which(linn.vip[,c1]!=0))
                    }
                    
                    good_feats<-unique(good_feats)
                    
                    
                    
                    if(length(good_feats)>1){
                        
                        cv_res<-plsda_cv(v=kfold,x=X[,good_feats],y=Yclass,ncomp=opt_comp,errortype=evalmethod)
                        
                        #cv_res<-pls.lda.cv(Xtrain=X, Ytrain=Y,  ncomp=c(1:numcomp), nruncv=kfold, alpha=2/3, priors=NULL)
                        
                        #cv_res<-cv_res$cv
                        print(paste(kfold," CV evaluation using spls top ",kvec," features per component",sep=""))
                        print(cv_res$mean_acc)
                        if(cv_res$mean_acc>best_cv_res){
                            
                            best_cv_res<-cv_res$mean_acc
                            best_kvec<-kvec
                            
                            #best_ncomp<-cv_res$ncomp
                        }
                    }else{
                        print("Too few variables to perfrom CV.")
                    }
                }
                
                keep_x_vec<-rep(best_kvec,opt_comp)
            }else{
                
                keep_x_vec<-rep(keepX,opt_comp)
            }
        }
        else{
            #keep_x_vec<-rep(dim(X)[2],opt_comp)
            keep_x_vec<-rep(keepX,opt_comp)
            
        }
        if(analysismode=="regression"){
            
            if(repeatmeasures==TRUE){
                
                #print("spls classlabels")
                #print(classlabels)
                #print(dim(X))
                #print(keep_x_vec)
                
                
                # linn.pls <- multilevel(X=X, design=classlabels,ncomp = opt_comp,
                #keepX = keep_x_vec, method = 'spls')
                
                linn.pls <- try(mixOmics::multilevel(X=X,Y=Y,design=classlabels,ncomp = opt_comp,
                keepX = keep_x_vec, method = 'spls'),silent=TRUE)
                
                if(is(linn.pls,"try-error")){
                    
                    
                    
                    
                    linn.pls <- mixOmics::spls(X=X,Y=Y,ncomp=opt_comp,keepX=keep_x_vec,multilevel=classlabels)
                    
                }
                
              
                
                
            }else{
                linn.pls <- mixOmics::spls(X, Y,ncomp=opt_comp,keepX=keep_x_vec,mode="regression")
            }
        }else{
            
            if(repeatmeasures==TRUE){
                #print("spls classlabels")
                #print(classlabels)
                #print(dim(X))
                
                
                #linn.pls <- multilevel(X=X, design=classlabels,ncomp = opt_comp,
                #keepX = keep_x_vec, method = 'splsda')
                
                linn.pls <- try(multilevel(X=X, design=classlabels,ncomp = opt_comp,
                keepX = keep_x_vec, method = 'splsda'),silent=TRUE)
                
                if(is(linn.pls,"try-error")){
                    
                    
                    
                    
                    linn.pls <- mixOmics::splsda(X=X,Y=classlabels[,-c(1)],ncomp=opt_comp,keepX=keep_x_vec,multilevel=classlabels[,1])
                    
                }
                
                if(FALSE){
                    stimu.time<-data.frame(cbind(as.character(classlabels[,2]),
                    as.character(classlabels[,3])))
                    repeat.stimu2<-classlabels[,1]
                    
                    res.2level <- multilevel(X, cond = stimu.time,
                    sample = repeat.stimu2, ncomp = 3,
                    keepX = keep_x_vec, tab.prob.gene = NULL, method = 'splsda')
                }
                
            }else{
                linn.pls <- mixOmics::splsda(X, Y,ncomp=opt_comp,keepX=keep_x_vec)
            }
        }
        
        linn.vip<-linn.pls$loadings$X
        
        
      
        
        
    }else{
        
        # print("opt comp")
        #print(opt_comp)
        if(analysismode=="regression"){
            linn.pls <- mixOmics::pls(X, Y,ncomp=opt_comp)
            
            
        }else{
            
            if(repeatmeasures==TRUE){
                linn.pls <- mixOmics::plsda(X, Yclass,ncomp=opt_comp,multileve=classlabels[,1])
            }else{
                    linn.pls <- mixOmics::plsda(X, Yclass,ncomp=opt_comp)
            }
        }
        
        linn.vip<-mixOmics::vip(linn.pls)
        
     
        
        
        
    }
    
    
    v1<-{}
    cv_res<-{}
 


 
    
return(list("model"=linn.pls,"vip_res"=linn.vip))
    
}

get_plscompvar<- function(p2,nvar,h) {
    
    
    b<-c(p2$loadings$Y)[1:h]
    T<-p2$variates$X[,1:h]
    SS<-b^2 * colSums(T^2)
    
    
    return(SS)
}



get_VIPmulticomp<- function(p2,nvar,h) {
    
    
    b<-c(p2$loadings$Y)[1:h]
    T<-p2$variates$X[,1:h]
    SS<-b^2 * colSums(T^2)
    
    W<-p2$loadings$X[,1:h]
    # print(dim(W))
    W<-as.data.frame(W)
    nvar=nrow(W)
    h=ncol(W)
    
    Wnorm2 <- colSums(W^2)
    pls_vec<-lapply(1:nvar,function(j){
        # return(sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS)))
        
        return(sqrt(nrow(W) * sum( (W[j,]^2 / Wnorm2) * (SS/ sum(SS)))))
    })
    pls_vec<-unlist(pls_vec)
    return(pls_vec)
}


#Function:find.Overlapping.mzs
#Description: This function matches features between two or more datasets using the
#following user defined criteria:
#1) Maximum m/z difference (+/-) ppm
#2) Maximum retention time difference in seconds
#Input:
#data_a->apLCMS output for dataset A,
#data_b->apLCMS output for dataset B,
#max.mz.diff->Maximum m/z difference (+/-) ppm
#max.rt.diff->Maximum retention time difference in seconds
#
#Output:
#Data frame that includes mz and retention time of common features
#
#Usage:
#common_features<-matchFeaturesmulti(data_a, data_b, max.mz.diff=10, max.rt.diff=300)
############################################
find.Overlapping.mzs<-function(dataA, dataB, mz.thresh=10, time.thresh=NA, alignment.tool=NA,use.best.match=FALSE)
{

        data_a<-as.data.frame(dataA)
        data_b<-as.data.frame(dataB)
	#data_a<-unique(data_a)
	rm(dataA)
	rm(dataB)
        
     #   data_b<-unique(data_b)
        
        com_mz_num=1
        unique_mz={}
        ppm_v={}
        rt_v={}

        commat={}

	col.names.dataA=colnames(data_a)
	col.names.dataB=colnames(data_b)

	if(is.na(alignment.tool)==FALSE){
	 if(alignment.tool=="apLCMS")
        {
              sample.col.start=5
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    sample.col.start=9
                    col.names.dataA[1]="mz"
                    col.names.dataA[2]="time"
		    col.names.dataB[1]="mz"
                    col.names.dataB[2]="time"
                    colnames(data_a)=col.names.dataA
                    colnames(data_b)=col.names.dataB
              }
	      
	}}else{
                    #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
		   
		    col.names.dataA[1]="mz"
		   
		    col.names.dataB[1]="mz"
		    if(is.na(time.thresh)==FALSE){
		     col.names.dataA[2]="time"
		     col.names.dataB[2]="time"
		      print("Using the 1st column as \"mz\" and 2nd column as \"retention time\"")
		     }else{
		      print("Using the 1st column as \"mz\"")
		     }
		    colnames(data_a)=col.names.dataA
                    colnames(data_b)=col.names.dataB
	}

       #data_a<-data_a[order(data_a$mz),]
       #data_b<-data_b[order(data_b$mz),]
       data_a<-as.data.frame(data_a)
	data_b<-as.data.frame(data_b)
	colnames(data_a)=col.names.dataA
	colnames(data_b)=col.names.dataB
        #create header for the matrix with common features
	if(is.na(time.thresh)==FALSE){
	mznames=c("index.A","mz.data.A", "time.data.A", "index.B","mz.data.B","time.data.B", "time.difference") 
        }else{
	mznames=c("index.A","mz.data.A", "index.B","mz.data.B") 
	}

        #Step 1 Group features by m/zdim(data_a)[1]
        mz_groups<-lapply(1:dim(data_a)[1],function(j){

                                commat={}
                                commzA=new("list")
                                commzB=new("list")
                                ppmb=(mz.thresh)*(data_a$mz[j]/1000000)

                                getbind_same<-which(abs(data_b$mz-data_a$mz[j])<=ppmb)

                                if(is.na(time.thresh)==FALSE){
                                  if(length(getbind_same)>0)
                                  {
                                          nearest_time_diff=10000
                                          bestmatch={}
					  rnames={}
					  temp={}
					  commat={}
                                          for (comindex in 1:length(getbind_same))
                                          {
											  tempA=cbind(j,data_a[j,c(1,2)])
											  tempB=cbind(getbind_same[comindex],data_b[getbind_same[comindex],c(1,2)])
					                                                  temp=cbind(tempA,tempB)
											
					                                                  timediff=abs(data_a[j,2]-data_b[getbind_same[comindex],2])
					                                                  
											        temp<-cbind(temp,timediff)
											 
					                                                  if(timediff<time.thresh && timediff<=nearest_time_diff)
					                                                  {
					                                                          bestmatch=as.data.frame(temp)
					                                                          nearest_time_diff=timediff
                                                                      }
												   
													
                                                    if(timediff<time.thresh)
                                                    {
												   temp<-as.data.frame(temp)
												   commat<-rbind(commat,temp)
												    rnamestemp<-paste("mz",j,"_",comindex,sep="")
												    rnames<-c(rnames,rnamestemp)
                                                    }
												
										}
					      
                          # if(use.best.match==TRUE){
                               
                               #print("HERE")
                                
                                #       commat=as.data.frame(bestmatch)
                                #}
                           
						if(length(commat)>=4){
							rownames(commat)=rnames
					        }


                                  }
                                }
                                else
                                {
                                    if(length(getbind_same)>0)
                                    {
                                    	temp1<-{}
                                    for (comindex in 1:length(getbind_same))
                                          {
						  tempA=cbind(j,data_a[j,c(1)])
						  tempB=cbind(getbind_same[comindex],data_b[getbind_same[comindex],c(1)])
                                                  temp=cbind(tempA,tempB)
                                                  #temp=cbind(data_a[j,c(1)],data_b[getbind_same[comindex],c(1)])
                                                 temp1<-rbind(temp1,temp)
                                          }
					  commat=as.data.frame(temp1)
					  rnames<-paste("mz",j,"_",seq(1,length(getbind_same)),sep="")
					  rownames(commat)=rnames
					  
                                    }
                                }
                                return(as.data.frame(commat))


                })
		
	#Step 2 Sub-group features from Step 1 by Retention time
        #find the features with RT values within the defined range as compared to the query feature

        uniqueinA={}
        uniqueinB={}
        commat=data.frame()

	
	if(length(mz_groups)>0){
        for(j in 1:length(mz_groups))
        {
                temp_diff={}
		
		if(is.list(mz_groups)==TRUE)
		{
			tempdata=mz_groups[[j]]
		}
		else
		{
			tempdata=mz_groups[j]
		}
		
                if(length(tempdata)>1)
                {
			
			colnames(tempdata)=mznames
			tempdata=as.data.frame(t(tempdata))
			
                        temp=tempdata
                        
                        temp=as.data.frame(temp)
			
                        if(is.null(commat)==TRUE)
                        {
                                commat=t(temp)

                        }
                        else
                        {

                                commat=rbind(commat,t(temp))

                        }
			

                }
	

        }

        if(is.null(dim(commat))==FALSE)
        {
                commat=as.data.frame(commat)
                
                if(use.best.match==TRUE){
                    
                    dup_index_A<-which(duplicated(commat$index.A)==TRUE)
             
                    if(length(dup_index_A)>0){
                        commat<-commat[-dup_index_A,]
                    }
                    
                    dup_index_B<-which(duplicated(commat$index.B)==TRUE)
                    
                    if(length(dup_index_B)>0){
                        commat<-commat[-dup_index_B,]
                    }
                    
                    
                }
        }
	}
	
	
        return(commat)
}



find.Unique.mzs<-function(dataA, dataB, mz.thresh=10, time.thresh=NA, alignment.tool=NA)
{

	
        data_a<-as.data.frame(dataA)
        data_b<-as.data.frame(dataB)
	data_a<-unique(data_a)
	rm(dataA)
	rm(dataB)
        
        data_b<-unique(data_b)
        
        com_mz_num=1
        unique_mz={}
        ppm_v={}
        rt_v={}

        commat={}

	col.names.dataA=colnames(data_a)
	col.names.dataB=colnames(data_b)

	if(is.na(alignment.tool)==FALSE){
	 if(alignment.tool=="apLCMS")
        {
              sample.col.start=5
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    sample.col.start=9
                    col.names.dataA[1]="mz"
                    col.names.dataA[2]="time"
		    col.names.dataB[1]="mz"
                    col.names.dataB[2]="time"
                    colnames(data_a)=col.names.dataA
                    colnames(data_b)=col.names.dataB
              }
	      
	}}else{
                    #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
		   
		    col.names.dataA[1]="mz"
		   
		    col.names.dataB[1]="mz"
		    if(is.na(time.thresh)==FALSE){
		     col.names.dataA[2]="time"
		     col.names.dataB[2]="time"
		      print("Using the 1st column as \"mz\" and 2nd column as \"retention time\"")
		     }else{
		      print("Using the 1st column as \"mz\"")
		     }
		    colnames(data_a)=col.names.dataA
                    colnames(data_b)=col.names.dataB
	}

       #data_a<-data_a[order(data_a$mz),]
       #data_b<-data_b[order(data_b$mz),]
       data_a<-as.data.frame(data_a)
	data_b<-as.data.frame(data_b)
	colnames(data_a)=col.names.dataA
	colnames(data_b)=col.names.dataB
	
        #create header for the matrix with common features
	if(is.na(time.thresh)==FALSE){
	mznames=c("index.A","mz.data.A", "time.data.A", "index.B","mz.data.B","time.data.B", "time.difference") 
        }else{
	mznames=c("index.A","mz.data.A", "index.B","mz.data.B") 
	}
	
	overlap_res<-find.Overlapping.mzs(dataA=data_a, dataB=data_b, mz.thresh, time.thresh, alignment.tool)
	
	

	if(length(overlap_res$index.A)>0){
	uniqueA<-data_a[-c(overlap_res$index.A),]
	}else{
		uniqueA<-data_a
	}
	
	if(length(overlap_res$index.B)>0){
	uniqueB<-data_b[-c(overlap_res$index.B),]
	}else{
		uniqueB<-data_b
	}
        return(list("uniqueA"=uniqueA,"uniqueB"=uniqueB))
}




#########################################################
#
#
#
#
#########################################################
getVenn<-function(dataA,name_a, dataB,name_b,mz.thresh=10,time.thresh=30,alignment.tool=NA, xMSanalyzer.outloc,use.unique.mz=FALSE,plotvenn=TRUE,use.best.match=FALSE)
{
	dir.create(xMSanalyzer.outloc,showWarnings=FALSE)
	
	data_a<-as.data.frame(dataA)
    data_b<-as.data.frame(dataB)
	rm(dataA)
	rm(dataB)

	############################################

	if(use.unique.mz==TRUE){
	data_a<-find.Unique.mzs.sameset(dataA=data_a,dataB=data_a,mz.thresh=mz.thresh,time.thresh=time.thresh,alignment.tool=alignment.tool)
	data_a<-data_a$uniqueA
	
	#print(dim(data_a))
	
	data_b<-find.Unique.mzs.sameset(dataA=data_b,dataB=data_b,mz.thresh=mz.thresh,time.thresh=time.thresh,alignment.tool=alignment.tool)
	data_b<-data_b$uniqueA
	}
	common<-find.Overlapping.mzs(data_a,data_b,mz.thresh,time.thresh=time.thresh,alignment.tool=alignment.tool,use.best.match=use.best.match)
	
	if(length(common)>0){
	commonA<-data_a[c(common$index.A),]
	commonA<-unique(commonA)
	
	commonB<-data_b[c(common$index.B),]
	commonB<-unique(commonB)
	
	data_a<-data_a[-c(common$index.A),]
	data_b<-data_b[-c(common$index.B),]
	
	

	rm_index<-which(data_a$mz%in%common$mz.data.A)
	
	if(length(rm_index)>0){
		uniqueA<-data_a[-rm_index,]
	}else{
		uniqueA<-data_a
	}
	
	rm_index<-which(data_b$mz%in%common$mz.data.B)
	
	if(length(rm_index)>0){
		uniqueB<-data_b[-rm_index,]
	}else{
		uniqueB<-data_b
	}
	
	#uniqueA<-data_a[-c(common$index.A),]
	#uniqueB<-data_b[-c(common$index.B),]
	num_commonA<-length(unique(common$index.A))
	num_commonB<-length(unique(common$index.B))
	
	num_common<-min(num_commonA,num_commonB)[1]
	}else{
	uniqueA<-data_a
	uniqueB<-data_b
	num_common<-0
	commonA<-{}
	commonB<-{}
	}
	num_commonA<-num_common
	num_commonB<-num_common
	num_uniqueA<-dim(uniqueA)[1]
	num_uniqueB<-dim(uniqueB)[1]
	
	#print(num_commonA)
	#print(num_uniqueA)
	g1 <-c(seq(1,(num_commonA+num_uniqueA)))
	g2<-c(seq(1,(num_commonB+num_uniqueB)))

	g1[1:num_commonA]=paste("x_",g1[1:num_commonA],sep="")
	g2[1:num_commonB]=paste("x_",g2[1:num_commonB],sep="")
	
	if(num_uniqueA>0){
	g1[(num_commonA+1):(num_commonA+num_uniqueA)]=paste("y_",g1[(num_commonA+1):(num_commonA+num_uniqueA)],sep="")
	}
	if(num_uniqueA>0){
		g2[(num_commonB+1):(num_commonB+num_uniqueB)]=paste("z_",g2[(num_commonB+1):(num_commonB+num_uniqueB)],sep="")
	}
	set1=as.character(g1)
	set2=as.character(g2)
	universe <- sort(unique(c(set1,set2)))
	
	Counts <- matrix(0, nrow=length(universe), ncol=2)
	colnames(Counts) <- c(name_a, name_b)
	for (i in 1:length(universe))
	{
		Counts[i,1] <- universe[i] %in% set1
		Counts[i,2] <- universe[i] %in% set2
	}
	fname<-paste(xMSanalyzer.outloc,"/Venn", name_a,"_",name_b,"_",mz.thresh,"ppm",time.thresh,"s.pdf",sep="")
	venn_counts<-vennCounts(Counts)
	if(plotvenn==TRUE){
		
	pdf(fname)
	vennDiagram(venn_counts)
	try(dev.off(),silent=TRUE)
	}
	
	return(list("common"=common,"commonA"=commonA,"uniqueA"=uniqueA,"commonB"=commonB,"uniqueB"=uniqueB,"vennCounts"=venn_counts))
	
}

	do_minmax<-function(data_m, newmin, newmax)
	{
		data_m=apply(data_m,2,function(x){
		minx=min(x,na.rm=TRUE)
		maxx=max(x,na.rm=TRUE)
		if(minx!=maxx)
		{
			(((x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))*(newmax-newmin))+newmin
		}else
		{
			(x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)+1)
		}
		})
		return(data_m)

	}
	


do_mean<-function(x){
	
	mean_val<-mean(x,na.rm=TRUE)
	return(mean_val)
}
	
do_rsd<-function(x){
	
	sd_val<-sd(x,na.rm=TRUE)
	mean_val<-mean(x,na.rm=TRUE)
	cv_val<-100*(sd_val/(mean_val))
	
	return(cv_val)
}


check_model<-function(model.fit,dependent.var){
    
    
    residuals <- resid(model.fit)
    plot(fitted(model.fit), residuals)
    abline(0,0)
    
    plot(fitted(model.fit), dependent.var)
    
    qqnorm(residuals)
    qqline(residuals)
}

getCorchild<-function(cur_mzdata,data_mt,cor.method){
	
	pearson_res<-lapply(1:dim(data_mt)[2],function(j){
		 return(WGCNA::cor(as.numeric(cur_mzdata),data_mt[,j],method=cor.method,use="pairwise.complete.obs"))
		
	
	})
	
	pvalues_list<-{}
	pearson_resmat<-{}
	pearson_list<-{}
	
	for(i in 1:length(pearson_res))
	{
		pearson_list<-c(pearson_list,pearson_res[[i]][1])
		nval<-length(which(is.na(data_mt[,i])==FALSE))
		
		pvalues_list<-c(pvalues_list,docortest(nval,pearson_res[[i]][[1]]))
	}
	
		
	return(list(cormat=pearson_list,complete_pearsonpvalue_mat=pvalues_list)) #,complete_pearsonqvalue_mat=qvalues_list))
}

#Function to perform 2-way ANOVA analysis and post-hoc comparisons using Tukey HSD method
diffexplmtwowayanova<-function(dataA){
    dataA<-as.data.frame(dataA)

    dataA$Factor1<-as.factor(dataA$Factor1)
    dataA$Factor2<-as.factor(dataA$Factor2)

    #Fit a 2-way ANOVA model
    res<-aov(as.numeric(Response) ~ (Factor1) + (Factor2) + (Factor1) * (Factor2),data=dataA)
    
    #get ANOVA results with p-values
    anova_res<-anova(res)

    #Perform post-hoc comparisons using the Tukey Honestly Significant Differences (HSD) test
    posthoc <- TukeyHSD(x=res, conf.level=0.95,test=univariate())

    num_rows<-dim(anova_res)[1]
    
    pvalues_factors<-data.frame(t(anova_res["Pr(>F)"][-c(num_rows),]))

    names(pvalues_factors)<-rownames(anova_res)[-c(num_rows)]

    interact_res<-t(c(posthoc$Factor1[,4],posthoc$Factor2[,4],posthoc$'Factor1:Factor2'[,4]))

    colnames(interact_res)<-c(rownames(posthoc$Factor1),rownames(posthoc$Factor2),rownames(posthoc$'Factor1:Factor2'))

    #return resutls
    return(list("mainpvalues"=pvalues_factors,"posthoc"=interact_res))


}


#Function to perform two-way ANOVA repeated measures analysis and post-hoc comparisons using Tukey HSD method
diffexplmtwowayanovarepeat<-function(dataA,subject_inf,modeltype="RI"){

    dataA<-as.data.frame(dataA)

    dataA$Factor1<-as.factor(dataA$Factor1)
    dataA$Factor2<-as.factor(dataA$Factor2)

    subject_inf<-as.vector(subject_inf)

    Subject<-subject_inf
    dataA<-cbind(dataA,subject_inf)

    #save(dataA,file="dataA.Rda")

    if(modeltype=="RI"){
        
        #anova(lme(as.numeric(Response) ~ Factor1 * Factor2, random=list(subject_inf=pdBlocked(list(~1,pdIdent(~Factor1-1),pdIdent(~Factor2-1)))),data=dataA))
        
        #call the lme function from the nlme package; random intercept only model
        res <- lme(as.numeric(Response) ~ Factor1 + Factor2 + Factor1 * Factor2, random = ~ 1 | subject_inf, data=dataA,control=lmeControl(opt="optim"))
    }else{
        if(modeltype=="RIRS"){
        
           #call the lme function from the nlme package;random intercept and random slope model
            res <- lme(as.numeric(Response) ~ Factor1 + Factor2 + Factor1 * Factor2, random=~ 1 + Factor1 | subject_inf, data=dataA,control=lmeControl(opt="optim"))
        }
        
    }

        if(is(res,"try-error")){
            
            return(list("mainpvalues"=NA,"posthoc"=NA))
        }else{
            
            anova_res<-anova(res)

            if(FALSE){
                #code for using glht for post-hoc comparisons
                dataA$SHD<-interaction(dataA$Factor1,dataA$Factor2)
                mod2<-lme(Response~-1+SHD, data=dataA, random=~1|subject_inf/Factor2,control=lmeControl(opt="optim"))
                posthoc<-summary(glht(mod2,linfct=mcp(SHD="Tukey")))
                posthoc_pvalues<-posthoc$test$pvalues
                names(posthoc_pvalues)<-names(posthoc$test$tstat)

            }
            
            #using lsmeans package for post-hoc comparisons since version v1.0.7.6
            means.factors=lsmeans(res,specs=c("Factor1","Factor2"))
            posthoc_res=pairs(means.factors,adjust="tukey")
            posthoc_res<-data.frame(posthoc_res)
            posthoc_pvalues<-posthoc_res$p.value
            names(posthoc_pvalues)<-as.character(posthoc_res$contrast)
            
            num_rows<-dim(anova_res)
            pvalues_factors<-data.frame(t(anova_res["p-value"][-c(1),]))

            names(pvalues_factors)<-rownames(anova_res)[-c(1)]

            return(list("mainpvalues"=pvalues_factors,"posthoc"=posthoc_pvalues))
        }

}


#Function to perform one-way ANOVA repeated measures analysis and post-hoc comparisons using Tukey HSD method
diffexplmonewayanovarepeat<-function(dataA,subject_inf,analysismode="classification",modeltype="RI"){
    
            dataA<-as.data.frame(dataA)
            
            if(analysismode=="classification"){
            dataA$Factor1<-as.factor(dataA$Factor1)
            }
            Subject<-subject_inf
            dataA<-cbind(dataA,subject_inf)
            
            save(dataA,file="dataA.Rda")

            if(modeltype=="RI"){
                #res <- try(lme(as.numeric(Response) ~ Factor1, random = ~ 1 | subject_inf, data=dataA,control=lmeControl(opt="optim")),silent=TRUE)
                res <- lme(as.numeric(Response) ~ Factor1, random = ~ 1 | subject_inf, data=dataA,control=lmeControl(opt="optim")) #,silent=TRUE)
            }else{
                if(modeltype=="RIRS"){
                 res <- try(lme(as.numeric(Response) ~ Factor1, random = ~ 1 + Factor1 | subject_inf, data=dataA,control=lmeControl(opt="optim")),silent=TRUE)
                }
            }
            if(is(res,"try-error")){
                return(list("mainpvalues"=NA,"posthoc"=NA))
            }else{
                anova_res<-anova(res)
                num_rows<-dim(anova_res)
                pvalues_factors<-data.frame(t(anova_res["p-value"][-c(1),]))

                        if(analysismode=="classification"){
                            
                            if(FALSE){
                                #using glht for post-hoc comparisons
                                posthoc<-summary(glht(res,linfct=mcp(Factor1="Tukey")))
                                posthoc_pvalues<-posthoc$test$pvalues
                                names(posthoc_pvalues)<-names(posthoc$test$tstat)
                                names(pvalues_factors)<-rownames(anova_res)[-c(1)]
                            }
                            
                            #using lsmeans package for post-hoc comparisons since version v1.0.7.6
                            means.factors=lsmeans(res,specs=c("Factor1"))
                            posthoc_res=pairs(means.factors,adjust="tukey")
                            posthoc_res<-data.frame(posthoc_res)
                            posthoc_pvalues<-posthoc_res$p.value
                            names(posthoc_pvalues)<-as.character(posthoc_res$contrast)
                            
                            return(list("mainpvalues"=pvalues_factors,"posthoc"=posthoc_pvalues))
                        }else{
                            return(list("mainpvalues"=pvalues_factors))
                        }
            }

}


#one-way ANOVA: nointeraction
diffexponewayanova<-function(dataA){
	
    dataA<-as.data.frame(dataA)

    #save(dataA,file="dataA.Rda")

    a1 <- aov(dataA$Response ~ .,data=dataA) # + chocolate$Factor1*chocolate$Factor2)

    posthoc <- TukeyHSD(x=a1, conf.level=0.95,test=univariate())
    anova_res<-anova(a1)

    num_rows<-dim(anova_res)
    pvalues_factors<-data.frame(t(anova_res["Pr(>F)"][-c(num_rows),]))

    return(list("mainpvalues"=pvalues_factors,"posthocfactor1"=posthoc$Factor1[,4]))

}


diffexpsvmrfe = function(x,y,svmkernel="radial"){
    
    #Checking for the variables
    stopifnot(!is.null(x) == TRUE, !is.null(y) == TRUE)
    
    n = ncol(x)
    survivingFeaturesIndexes = seq_len(n)
    featureRankedList = vector(length=n)
    rankedFeatureIndex = n
    
    while(length(survivingFeaturesIndexes)>0){
        #train the support vector machine
        svmModel = svm(x[, survivingFeaturesIndexes], y,
        scale=FALSE, type="C-classification", kernel=svmkernel)
        
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
    
    return (featureRankedList)
}


diffexpsvmrfemulticlass  = function(x,y,svmkernel="radial"){
    n = ncol(x)
    survivingFeaturesIndexes = seq(n)
    featureRankedList = vector(length=n)
    rankedFeatureIndex = n
    while(length(survivingFeaturesIndexes)>0){
        #train the support vector machine
        svmModel = svm(x[, survivingFeaturesIndexes],
        y,
        scale= FALSE,
        type="C-classification", kernel=svmkernel)
        #compute the weight vector
        multiclassWeights = svm.getweights(svmModel)
        #compute ranking criteria
        multiclassWeights = multiclassWeights * multiclassWeights
        rankingCriteria = 0
        for(i in 1:ncol(multiclassWeights))rankingCriteria[i] = mean(multiclassWeights[,i])
        #rank the features
        (ranking = sort(rankingCriteria, index.return = TRUE)$ix)
        
       
        ## New update feature ranked list
        #featureRankedList[rev((s-r+1):s)] <-survivingFeaturesIndexes[ranking[1:r]]
        featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]
        rankedFeatureIndex = rankedFeatureIndex - 1
        ## New to remove perc.rem
        #  rankedFeatureIndex <- rankedFeatureIndex - r
        
        #eliminate the feature with smallest ranking criterion
        #survivingFeaturesIndexes <-
        #   survivingFeaturesIndexes[-ranking[1:r]]
        
        (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
    }
    return(featureRankedList)
}


svm.getweights<-function(model){
    w=0
    if(model$nclasses==2){
        w=t(model$coefs)%*%model$SV
    }else{
        
        #when we deal with OVO svm classification
        ## compute start-index
        start <- c(1, cumsum(model$nSV)+1)
        start <- start[-length(start)]
        calcw <- function (i,j) {
            ## ranges for class i and j:
            ri <- start[i] : (start[i] + model$nSV[i] -1)
            rj <- start[j] : (start[j] + model$nSV[j] -1)
            ## coefs for (i,j):
            coef1 <- model$coefs[ri, j-1]
            coef2 <- model$coefs[rj, i]
            ## return w values:
            w=t(coef1)%*%model$SV[ri,]+t(coef2)%*%model$SV[rj,]
            return(w)
        }
        W=NULL
        for (i in 1 : (model$nclasses - 1)){
            for (j in (i + 1) : model$nclasses){
                wi=calcw(i,j) 
                W=rbind(W,wi) 
            } 
        } 
        w=W 
    } 
    return(w) 
}

#Nointeraction
diffexplmreg<-function(dataA,logistic_reg=FALSE,poisson_reg=FALSE){
	
dataA<-as.data.frame(dataA)

##save(dataA,file="lmreg_func.Rda")
if(logistic_reg==TRUE){
    cnames1<-colnames(dataA)
    cnames1[2]<-"Class"
    colnames(dataA)<-cnames1
	
    labels_1<-levels(as.factor(dataA$Class))

dataA$Class<-replace(dataA$Class,which(dataA$Class==labels_1[1]),0)
 dataA$Class<-replace(dataA$Class,which(dataA$Class==labels_1[2]),1)


	a1 <- glm(dataA$Class ~ .,family=binomial(logit),data=dataA)
}else{
    
    if(poisson_reg==TRUE){
        
        cnames1<-colnames(dataA)
        cnames1[2]<-"Class"
        colnames(dataA)<-cnames1
        
        ##save(dataA,file="temp1.Rda")
        labels_1<-levels(as.factor(dataA$Class))
        dataA$Class<-as.numeric(dataA$Class)
        
        a1 <- glm(dataA$Class ~ .,family=poisson(log),data=dataA)
        
    }else{

        a1 <- lm(dataA$Response ~ .,data=dataA) # aov(dataA$Response ~ .,data=dataA) # + chocolate$Factor1*chocolate$Factor2)
    }
}
s1<-summary(a1)


if(logistic_reg==FALSE){
    
    r2<-s1$adj.r.squared
}else{
    r2<-NA
}
if(poisson_reg==FALSE){
s1<-s1$coefficients
}else{
    cov.a1 <- vcovHC(a1, type="HC0")
    std.err <- sqrt(diag(cov.a1))
    s1 <- cbind(Estimate= coef(a1), "Robust SE" = std.err, "z value"=coef(a1)/std.err,
    "Pr(>|z|)" = 2 * pnorm(abs(coef(a1)/std.err), lower.tail=FALSE),
    LL = coef(a1) - 1.96 * std.err,
    UL = coef(a1) + 1.96 * std.err)
    
    
}

s1<-s1[-c(1),]

if(dim(dataA)[2]<3){ # && dim(dataA)[1]<3){
    #s1<-as.data.frame(s1)
s1<-t(s1)

}


    confint_lower<-s1[,1]-(1.96*s1[,2])
    confint_upper<-s1[,1]+(1.96*s1[,2])
    

return(list("mainpvalues"=s1[,4],"estimates"=s1[,1],"statistic"=s1[,3],"stderr"=s1[,2],"r2"=r2,"confint"=c(confint_lower,confint_upper)))


}

#Nointeraction
diffexpcoxph<-function(dataA,logistic_reg=FALSE,poisson_reg=FALSE){
    
    dataA<-as.data.frame(dataA)
    
   # #save(dataA,file="lmreg_func.Rda")
    
    
    if(logistic_reg==TRUE){
        cnames1<-colnames(dataA)
        cnames1[2]<-"Class"
        colnames(dataA)<-cnames1
        
        labels_1<-levels(as.factor(dataA$Class))
        
        dataA$Class<-replace(dataA$Class,which(dataA$Class==labels_1[1]),0)
        dataA$Class<-replace(dataA$Class,which(dataA$Class==labels_1[2]),1)
        
        
        a1 <- glm(dataA$Class ~ .,family=binomial(logit),data=dataA)
    }else{
        
        if(poisson_reg==TRUE){
            
            cnames1<-colnames(dataA)
            cnames1[2]<-"Class"
            colnames(dataA)<-cnames1
            
            ##save(dataA,file="temp1.Rda")
            labels_1<-levels(as.factor(dataA$Class))
            dataA$Class<-as.numeric(dataA$Class)
            
            a1 <- glm(dataA$Class ~ .,family=poisson(log),data=dataA)
            
        }else{
            
            a1 <- lm(dataA$Response ~ .,data=dataA) # aov(dataA$Response ~ .,data=dataA) # + chocolate$Factor1*chocolate$Factor2)
        }
    }
    s1<-summary(a1)
    
    
    if(logistic_reg==FALSE){
        
        r2<-s1$adj.r.squared
    }else{
        r2<-NA
    }
    if(poisson_reg==FALSE){
        s1<-s1$coefficients
    }else{
        cov.a1 <- vcovHC(a1, type="HC0")
        std.err <- sqrt(diag(cov.a1))
        s1 <- cbind(Estimate= coef(a1), "Robust SE" = std.err, "z value"=coef(a1)/std.err,
        "Pr(>|z|)" = 2 * pnorm(abs(coef(a1)/std.err), lower.tail=FALSE),
        LL = coef(a1) - 1.96 * std.err,
        UL = coef(a1) + 1.96 * std.err)
        
        
    }
    
    s1<-s1[-c(1),]
    
    if(dim(dataA)[2]<3){ # && dim(dataA)[1]<3){
        #s1<-as.data.frame(s1)
        s1<-t(s1)
        
    }
    
    
    confint_lower<-s1[,1]-(1.96*s1[,2])
    confint_upper<-s1[,1]+(1.96*s1[,2])
    
    
    return(list("mainpvalues"=s1[,4],"estimates"=s1[,1],"statistic"=s1[,3],"stderr"=s1[,2],"r2"=r2,"confint"=c(confint_lower,confint_upper)))
    
    
}
do_mars_lmreg<-function(dataA){


dataA<-as.data.frame(dataA)


#print(dim(dataA))
mars_res<-new("list")
		

	mars_res[[1]]<-try(earth(dataA$Response~.,data=as.data.frame(dataA),degree=1,ncross=10,nfold=10),silent=TRUE)
	mars_res[[2]]<-try(earth(dataA$Response~.,data=as.data.frame(dataA),degree=2,ncross=10,nfold=10),silent=TRUE)
	mars_res[[3]]<-try(earth(dataA$Response~.,data=as.data.frame(dataA),degree=3,ncross=10,nfold=10),silent=TRUE)
	
	gcv_list<-c(mars_res[[1]]$gcv,mars_res[[2]]$gcv,mars_res[[3]]$gcv)

	min_gcv<-which(gcv_list==min(gcv_list,na.rm=TRUE))

	marsfitdeg<-earth(formula=factor(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=min_gcv[1],ncross=kfold,nfold=10)

#marsfitdeg<-earth(mz~Strain+Treatment+Batch,data=as.data.frame(curdata),degree=3,ncross=10,nfold=10)
#marsfitdeg<-try(earth(dataA$Response~.,data=as.data.frame(dataA),degree=1,ncross=10,nfold=10),silent=TRUE)

if (is(marsfitdeg, "try-error")){
	lmtest_pval<-NA
	}else{
        #print(summary(marsfitdeg))
#evdeg1 <- evimp(marsfitdeg)
#plot(evdeg1)
qqnorm(marsfitdeg$residuals)
qqline(marsfitdeg$residuals)
shapiro.test(marsfitdeg$residuals)
bx1 <- model.matrix(marsfitdeg)
#deg.lm <- glm(as.vector(curdata[,1]) ~ bx1[,-1]) # -1 to drop intercept
deg.lm <- try(glm(as.numeric(dataA$Response) ~ bx1[,-1]),silent=TRUE)
if (is(deg.lm, "try-error")){
														lmtest_pval<-1
													}else{
														slm<-summary(deg.lm)
														lmtest_pval=slm$coefficients[,4]
														#print(slm$coefficients)
														
													}

#slm<-summary(deg.lm) # yields same coeffs as above summary
#aicdeg1[i]<-s1$aic
 #plot(effect('Dept:Gender', berk.mod2), multiline=TRUE)
 }
 return(list("marsummary"=summary(marsfitdeg), "pvalues"=lmtest_pval))
 
 }
 

diffexplogitreg<-function(dataA){
	
dataA<-as.data.frame(dataA)

labels_1<-labels(as.factor(dataA$Factor1))

print(labels_1)
dataA$Factor1<-replace(dataA$Factor1,which(dataA$Factor1==labels_1[1]),0)
 dataA$Factor1<-replace(dataA$Factor1,which(dataA$Factor1==labels_1[2]),1)

a1 <- glm(dataA$Factor1 ~ .,family=binomial(logit),data=dataA) # aov(dataA$Response ~ .,data=dataA) # + chocolate$Factor1*chocolate$Factor2) 

c1<-confint(a1,level=0.95)
s1<-summary(a1)
#print(summary(a1))

anova_res<-anova(a1)
num_rows<-dim(anova_res)
#pvalues_factors<-data.frame(t(anova_res["Pr(>F)"][-c(num_rows),]))

s1<-s1$coefficients
s1<-s1[-c(1),]

confint_lower<-s1[,1]-(1.96*s1[,2])
confint_upper<-s1[,1]+(1.96*s1[,2])

#print(anova_res)
#
return(list("mainpvalues"=s1[,4],"estimates"=s1[,1],"zstat"=s1[,3],"stderr"=s1[,2],"confint"=c(confint_lower,confint_upper)))


}



	
		
	do_cor<-function(data_m_fc_withfeats,subindex=NA,targetindex=NA,outloc,networkscope,cor.method,abs.cor.thresh,cor.fdrthresh,
	max.cor.num,net_node_colors,net_legend,netrandseed=555,num_nodes=6, plots.width=2000,plots.height=2000,plots.res=300){
	Sys.sleep(0.2)
dir.create(outloc)
setwd(outloc)
allsig_pcornetwork<-{}
cormat<-{}

		


allsig_pcornetwork_fdr0.05<-{}

#data_m_fc_withfeats<-apply(data_m_fc_withfeats,1,function(x){naind<-which(is.na(x)==TRUE);if(length(naind)>0){ x[naind]<-median(x,na.rm=TRUE)};return(x)})
#		data_m<-t(data_m_fc_withfeats)

allsig_pcornetwork_fdr0.05<-{}

if(dim(data_m_fc_withfeats)[1]>1)
{
	if(is.na(subindex[1])==FALSE){
		goodfeats<-data_m_fc_withfeats[subindex,]
		
	}else{
		goodfeats<-data_m_fc_withfeats
		
		}

if(is.na(targetindex[1])==FALSE){
	data_m_fc_withfeats<-data_m_fc_withfeats[targetindex,]
		
	}else{
		data_m_fc_withfeats<-data_m_fc_withfeats
		
		}
	
	if(is.na(subindex[1])==FALSE){
        if(length(subindex)==nrow(data_m_fc_withfeats)){
            
            data_m_fc_withfeats<-data_m_fc_withfeats
            
        }else{
                data_m_fc_withfeats<-data_m_fc_withfeats[-subindex,]
        }
    }
m1<-apply(data_m_fc_withfeats[,-c(1:2)],2,as.numeric)

print(dim(m1))
rnames<-paste("mzid_",seq(1,dim(data_m_fc_withfeats)[1]),sep="")


rownames(m1)=rnames

data_mt<-t(m1)
}else{
m1<-as.numeric(data_m_fc_withfeats[,-c(1:2)])
#print(dim(m1))
#rownames(m1)=rnames
data_mt<-(m1)

data_mt<-as.matrix(data_mt)

}

#listf<-list.files(".","*.txt",recursive=TRUE)
#for(l1 in listf){pmattemp<-read.table(l1,sep="\t",header=TRUE); cormat<-cbind(cormat,pmattemp[,3])}
goodfeats_inf<-as.data.frame(goodfeats[,c(1:2)])

#print(dim(goodfeats))
#print(dim(data_mt))

#for(l1 in listf){pmattemp<-read.table(l1,sep="\t",header=TRUE); cormat<-cbind(cormat,pmattemp[,3])}
goodfeats<-as.data.frame(goodfeats)

complete_pearsonpvalue_mat<-{}
complete_pearsonqvalue_mat<-{}
	l1<-list.files(".")
	

 print(paste("Computing ",cor.method," matrix",sep=""))
        #system.time(pearson_res<-parRapply(cl,goodfeats[,-c(1:2)],getCorchild,data_mt,cor.method))

		temp_mat<-rbind(goodfeats[,-c(1:2)],t(data_mt))
		
		pearson_res<-WGCNA::cor(t(temp_mat),nThreads=num_nodes,method=cor.method)


		
		pearson_res<-pearson_res[c(1:dim(goodfeats[,-c(1:2)])[1]),-c(1:dim(goodfeats[,-c(1:2)])[1])]
		
		num_samp<-dim(data_mt)[1]


num_samp<-dim(goodfeats[,-c(1:2)])[2]

		 pearson_resmat<-{}
		 
		# print("done")
		 
		 pearson_resvec<-as.vector(pearson_res)
		 
		 pearson_pvalue<-lapply(1:length(pearson_resvec),function(x){
		 	
		 	return(WGCNA::corPvalueStudent(pearson_resvec[x],num_samp))

		 	
		 })
		 
		 complete_pearsonpvalue_mat<-unlist(pearson_pvalue)
		 dim(complete_pearsonpvalue_mat)<-dim(pearson_res)
		 
			cormat<-pearson_res
			
					pearson_Res_all<-{}

		   
		   #print("here 3")
		  # print(head(cormat))
		   #print(dim(cormat[[1]]))
       
        #cormat<-t(cormat)
        cormat<-as.data.frame(cormat)
        cormat<-as.matrix(cormat)

		#print(dim(cormat))
               
               pearson_Res_all<-{}
        
        #complete_pearsonpvalue_mat<-t(complete_pearsonpvalue_mat)
        complete_pearsonpvalue_mat<-as.data.frame(complete_pearsonpvalue_mat)
        complete_pearsonpvalue_mat<-as.matrix(complete_pearsonpvalue_mat)
			
			
			complete_pearsonpvalue_mat<-t(complete_pearsonpvalue_mat)
			cormat<-t(cormat)
			
            #		print(dim(cormat))
            #print(dim(complete_pearsonpvalue_mat))
			

			
			pearson_Res_all<-{}
	#complete_pearsonqvalue_mat<-sapply(1:length(pearson_res),function(j){pearson_Res_all<-rbind(pearson_Res_all,as.matrix(pearson_res[[j]]$complete_pearsonqvalue_mat))})
	

	
	complete_pearsonqvalue_mat<-apply(cormat,2,function(x){
		
		#p.adjust(x,method="BH")
					#print(length(x))
					
					x<-as.numeric(x)
					pdf("fdrtoolB.pdf")
					
					fdr_adjust_pvalue<-try(fdrtool(x,statistic="correlation",verbose=FALSE),silent=TRUE)
					
					if(is(fdr_adjust_pvalue,"try-error")){
					fdr_adjust_pvalue<-NA
					}else{
					fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
					}
					try(dev.off(),silent=TRUE)
					return(fdr_adjust_pvalue)
                    }
		
		)
		
		
        complete_pearsonqvalue_mat<-as.data.frame(complete_pearsonqvalue_mat)
        complete_pearsonqvalue_mat<-as.matrix(complete_pearsonqvalue_mat)
        
setwd(outloc)

	
goodfeats<-as.data.frame(goodfeats)

	colnames(cormat)<-as.character(goodfeats$mz)
	

	data_m_fc_withfeats<-as.data.frame(data_m_fc_withfeats)
	

	rownames(cormat)<-as.character(data_m_fc_withfeats$mz)
	
	cormat<-round(cormat,2)


fname<-paste("correlation_matrix.txt",sep="")
write.table(cormat,file=fname,sep="\t",row.names=TRUE)

fname<-paste("correlation_pvalues.txt",sep="")

complete_pearsonpvalue_mat<-round(complete_pearsonpvalue_mat,2)

   write.table(complete_pearsonpvalue_mat,file=fname,sep="\t",row.names=TRUE)

fname<-paste("correlation_qvalues.txt",sep="")

complete_pearsonqvalue_mat<-round(complete_pearsonqvalue_mat,2)

   write.table(complete_pearsonqvalue_mat,file=fname,sep="\t",row.names=TRUE)




	#nonsig_vs_fdr0.05_pearson_mat<-
	nonsig_vs_fdr0.05_pearson_mat_bool<-cormat

	nonsig_vs_fdr0.05_pearson_mat_bool[complete_pearsonqvalue_mat>cor.fdrthresh]<-0

cormat_abs<-abs(cormat)

 nonsig_vs_fdr0.05_pearson_mat_bool[cormat_abs<abs.cor.thresh]<-0

		rm(cormat_abs)
	
	sum_mat<-apply(nonsig_vs_fdr0.05_pearson_mat_bool,1,sum)

	bad_index<-which(sum_mat==0)
	
	
	
	nonsig_vs_fdr0.05_pearson_mat_bool_filt<-nonsig_vs_fdr0.05_pearson_mat_bool
	
	mz_index<-which(data_m_fc_withfeats$mz%in%goodfeats$mz)
	
		if(length(bad_index)>0){
	nonsig_vs_fdr0.05_pearson_mat_bool_filt<-nonsig_vs_fdr0.05_pearson_mat_bool[-bad_index,]
	}
	
	
	if(length(mz_index)==length(data_m_fc_withfeats$mz))
	{
		nonsig_vs_fdr0.05_pearson_mat_bool_filt[lower.tri(nonsig_vs_fdr0.05_pearson_mat_bool_filt)==TRUE]<-0
		diag(nonsig_vs_fdr0.05_pearson_mat_bool_filt)<-0
	}else{
	
				mz_index<-which(data_m_fc_withfeats$mz%in%goodfeats$mz)
	
	if(length(mz_index)>0){
	nonsig_vs_fdr0.05_pearson_mat_bool_filt<-nonsig_vs_fdr0.05_pearson_mat_bool_filt[-mz_index,]
	}
		
	}
	
	nonsig_vs_fdr0.05_pearson_mat_bool_filt<-unique(nonsig_vs_fdr0.05_pearson_mat_bool_filt)
	
	#print(length(nonsig_vs_fdr0.05_pearson_mat_bool_filt))
	#print(dim(nonsig_vs_fdr0.05_pearson_mat_bool))
	
	
	#>=dim(goodfeats)[1]
	if(length(nonsig_vs_fdr0.05_pearson_mat_bool_filt)>0){
		
		
		nonsig_vs_fdr0.05_pearson_mat_bool_filt_1<-as.data.frame(nonsig_vs_fdr0.05_pearson_mat_bool_filt)
		
	if(length(nonsig_vs_fdr0.05_pearson_mat_bool_filt)>dim(nonsig_vs_fdr0.05_pearson_mat_bool_filt_1)[1] && (length(nonsig_vs_fdr0.05_pearson_mat_bool_filt)>1))
	#if((length(nonsig_vs_fdr0.05_pearson_mat_bool_filt)>1))
{
	
	

	check_cor<-apply(abs(nonsig_vs_fdr0.05_pearson_mat_bool_filt),1,function(x){max(x,na.rm=TRUE)})
	
	nonsig_vs_fdr0.05_pearson_mat_bool_filt<-nonsig_vs_fdr0.05_pearson_mat_bool_filt[order(check_cor,decreasing=TRUE),]
		
	#print(head(nonsig_vs_fdr0.05_pearson_mat_bool_filt))
	
	nonsig_vs_fdr0.05_pearson_mat_bool_filt<-na.omit(nonsig_vs_fdr0.05_pearson_mat_bool_filt)
	fname<-paste("significant_correlations_",networkscope,"_matrix_mzlabels.txt",sep="")
	
    #nonsig_vs_fdr0.05_pearson_mat_bool_filt<-round(nonsig_vs_fdr0.05_pearson_mat_bool_filt,2)
    
	write.table(nonsig_vs_fdr0.05_pearson_mat_bool_filt,file=fname,sep="\t",row.names=TRUE)
	

	fname<-paste("significant_correlations_",networkscope,"CIRCOSformat_mzlabels.txt",sep="")
	
	
	mz_rnames<-rownames(nonsig_vs_fdr0.05_pearson_mat_bool_filt)
	mz_cnames<-colnames(nonsig_vs_fdr0.05_pearson_mat_bool_filt)
	#rnames<-c("Data",rnames)
	
	circos_format<-abs(nonsig_vs_fdr0.05_pearson_mat_bool_filt)
    #circos_format<-round(circos_format,2)
    
	#circos_format<-rbind(cnames,circos_format)
	circos_format<-cbind(mz_rnames,circos_format)
	
    #
    
	write.table(circos_format,file=fname,sep="\t",row.names=FALSE)
	
	
	rnames<-mz_rnames
	cnames<-mz_cnames
	
	cnames<-seq(1,length(cnames))
	rnames<-seq(1,length(rnames))
	
    id_mapping_mat<-cbind(rnames,mz_rnames)
    
    # id_mapping_mat<-rbind(id_mapping_mat,cbind(cnames,mz_cnames))
    #colnames(id_mapping_mat)<-c("ID","Name")
    #write.csv(id_mapping_mat,file="node_id_mapping.csv",row.names=FALSE)
    
    
	#nonsig_vs_fdr0.05_pearson_mat_bool_filt<-nonsig_vs_fdr0.05_pearson_mat_bool_filt[-which(duplicated(rnames)==TRUE),]
	#rnames<-rnames[-which(duplicated(rnames)==TRUE)]
	
	
	
	#rnames<-rownames(nonsig_vs_fdr0.05_pearson_mat_bool_filt)
	#cnames<-colnames(nonsig_vs_fdr0.05_pearson_mat_bool_filt)
	#cnames<-round(as.numeric(cnames),5)
	#rnames<-round(as.numeric(rnames),5)
    ##save(nonsig_vs_fdr0.05_pearson_mat_bool_filt,file="nonsig_vs_fdr0.05_pearson_mat_bool_filt.Rda")
    
    
    if(length(nonsig_vs_fdr0.05_pearson_mat_bool_filt)>0)
    {
        
        if(nrow(nonsig_vs_fdr0.05_pearson_mat_bool_filt)>0){
    
    
	colnames(nonsig_vs_fdr0.05_pearson_mat_bool_filt)<-paste("Y",cnames,sep="")
	rownames(nonsig_vs_fdr0.05_pearson_mat_bool_filt)<-paste("X",rnames,sep="")

	fname<-paste("significant_correlations_",networkscope,"_matrix_rowcolnumlabels.txt",sep="")
	
    #nonsig_vs_fdr0.05_pearson_mat_bool_filt<-round(nonsig_vs_fdr0.05_pearson_mat_bool_filt,2)
    
	write.table(nonsig_vs_fdr0.05_pearson_mat_bool_filt,file=fname,sep="\t",row.names=TRUE)
        }else{
            stop("No correlations found.")
        }
    }else{
        
        stop("No correlations found.")
    }
	circos_format<-abs(nonsig_vs_fdr0.05_pearson_mat_bool_filt)
	

	circos_format<-cbind(rnames,circos_format)
 
	fname<-paste("significant_correlations_",networkscope,"CIRCOSformat_rowcolnumlabels.txt",sep="")
	
  
    write.table(circos_format,file=fname,sep="\t",row.names=FALSE)


	if(length(which(check_cor>=abs.cor.thresh))>0){
		
		
	if(is.na(max.cor.num)==FALSE){
	
	if(max.cor.num>dim(nonsig_vs_fdr0.05_pearson_mat_bool_filt)[1]){
	max.cor.num<-dim(nonsig_vs_fdr0.05_pearson_mat_bool_filt)[1]
	}
	nonsig_vs_fdr0.05_pearson_mat_bool_filt<-nonsig_vs_fdr0.05_pearson_mat_bool_filt[1:max.cor.num,]
	
	mz_rnames<-mz_rnames[1:max.cor.num]
	rnames<-rnames[1:max.cor.num]
	}

	
    
    pdfname<-paste(networkscope,"network_plot.pdf",sep="")
    pdf(pdfname,width=8,height=10)

    
    
	
	dup_ind<-which(duplicated(rnames)==TRUE)
	if(length(dup_ind)>0){
	nonsig_vs_fdr0.05_pearson_mat_bool_filt_1<-nonsig_vs_fdr0.05_pearson_mat_bool_filt[-which(duplicated(rnames)==TRUE),]
	rnames<-rnames[-which(duplicated(rnames)==TRUE)]
	mz_rnames<-mz_rnames[-which(duplicated(rnames)==TRUE)]
	
	}else{
		nonsig_vs_fdr0.05_pearson_mat_bool_filt_1<-nonsig_vs_fdr0.05_pearson_mat_bool_filt
	}
    
    dup_indB<-which(duplicated(cnames)==TRUE)
    if(length(dup_indB)>0){
        nonsig_vs_fdr0.05_pearson_mat_bool_filt_1<-nonsig_vs_fdr0.05_pearson_mat_bool_filt_1[,-which(duplicated(cnames)==TRUE)]
        cnames<-cnames[-which(duplicated(cnames)==TRUE)]
        mz_cnames<-mz_cnames[-which(duplicated(cnames)==TRUE)]
        
    }
    
	set.seed(netrandseed)
    
    ##save(nonsig_vs_fdr0.05_pearson_mat_bool_filt_1,file="nonsig_vs_fdr0.05_pearson_mat_bool_filt_1.Rda")
    
    
	net_result<-try(network(mat=as.matrix(nonsig_vs_fdr0.05_pearson_mat_bool_filt_1), threshold=abs.cor.thresh,color.node = net_node_colors,
    shape.node = c("rectangle", "circle"),
	color.edge = c("red", "blue"), lwd.edge = 1,
    show.edge.labels = FALSE, interactive = FALSE,cex.node.name=0.6),silent=TRUE)
    
    
    if(is(net_result,"try-error")){
        
       set.seed(netrandseed)
        net_result<-network(mat=as.matrix(nonsig_vs_fdr0.05_pearson_mat_bool_filt_1), cutoff=abs.cor.thresh,color.node = net_node_colors,
        shape.node = c("rectangle", "circle"),
        color.edge = c("red", "blue"), lwd.edge = 1,
        show.edge.labels = FALSE, interactive = FALSE,cex.node.name=0.6) #,silent=TRUE)
        
       
    }

   
	cytoscape_fname<-paste("network_Cytoscape_format_maxcor",max.cor.num,"_rowcolnumlabels.gml",sep="")
	write.graph(net_result$gR, file =cytoscape_fname, format = "gml")

    if(net_legend==TRUE){
    print(legend("bottomright",c("Row #","Column #"),pch=c(22,19),col=net_node_colors, cex=cex.plots,title="Network matrix values:"))
	}
	

	
    
    
	try(dev.off(),silent=TRUE)
	
	
	}else{
	if(networkscope=="all"){
	print(paste("Metabolome-wide correlation network can not be generated as the correlation threshold criteria is not met.",sep=""))
	}else{
	print(paste("Targeted correlation network can not be generated as the correlation threshold criteria is not met.",sep=""))
	}
	
	}
	
	}

	}else{
	if(length(nonsig_vs_fdr0.05_pearson_mat_bool_filt)==dim(goodfeats)[1]){
	
	#print("here")
	
		rnames<-rownames(nonsig_vs_fdr0.05_pearson_mat_bool)
		rnames1<-rnames[-bad_index]
		nonsig_vs_fdr0.05_pearson_mat_bool_filt<-as.matrix(nonsig_vs_fdr0.05_pearson_mat_bool_filt)
		nonsig_vs_fdr0.05_pearson_mat_bool_filt<-t(nonsig_vs_fdr0.05_pearson_mat_bool_filt)
		rownames(nonsig_vs_fdr0.05_pearson_mat_bool_filt)<-as.character(rnames1)
		check_cor<-max(abs(nonsig_vs_fdr0.05_pearson_mat_bool_filt))
		
		fname<-paste("significant_correlations_",networktype,"_mzlabels.txt",sep="")
	
    #nonsig_vs_fdr0.05_pearson_mat_bool_filt<-round(nonsig_vs_fdr0.05_pearson_mat_bool_filt,2)
        
		write.table(nonsig_vs_fdr0.05_pearson_mat_bool_filt,file=fname,sep="\t",row.names=TRUE)
		
	fname<-paste("significant_correlations_",networktype,"CIRCOSformat_mzlabels.txt",sep="")
	
	mz_rnames<-rownames(nonsig_vs_fdr0.05_pearson_mat_bool_filt)
	mz_cnames<-colnames(nonsig_vs_fdr0.05_pearson_mat_bool_filt)
	#rnames<-c("Data",rnames)
	
	circos_format<-abs(nonsig_vs_fdr0.05_pearson_mat_bool_filt)
	
    #circos_format<-round(circos_format,2)
    
	#circos_format<-rbind(cnames,circos_format)
	circos_format<-cbind(rnames,circos_format)
	
    
    
	write.table(circos_format,file=fname,sep="\t",row.names=FALSE)
	
	
#	print("here")
	rnames<-mz_rnames
	cnames<-mz_cnames
	
	#cnames<-round(as.numeric(cnames),5)
	#rnames<-round(as.numeric(rnames),5)
	
	cnames<-seq(1,length(cnames))
	rnames<-seq(1,length(rnames))
    
    id_mapping_mat<-cbind(rnames,mz_rnames)
    
    #id_mapping_mat<-rbind(id_mapping_mat,cbind(cnames,mz_cnames))
    # colnames(id_mapping_mat)<-c("ID","Name")
    #write.csv(id_mapping_mat,file="node_id_mapping.csv",row.names=FALSE)
	
	colnames(nonsig_vs_fdr0.05_pearson_mat_bool_filt)<-paste("Y",cnames,sep="")
	rownames(nonsig_vs_fdr0.05_pearson_mat_bool_filt)<-paste("X",rnames,sep="")
	
	fname<-paste("significant_correlations_",networktype,"CIRCOS_format_rowcolnumlabels.txt",sep="")
	
	circos_format<-abs(nonsig_vs_fdr0.05_pearson_mat_bool_filt)
    #circos_format<-round(circos_format,2)
    
	#circos_format<-rbind(cnames,circos_format)
	circos_format<-cbind(rnames,circos_format)
    
    
	
	write.table(circos_format,file=fname,sep="\t",row.names=FALSE)
	
		fname<-paste("significant_correlations_",networktype,"matrix_rowcolnumlabels.txt",sep="")
	
    #  nonsig_vs_fdr0.05_pearson_mat_bool_filt<-round(nonsig_vs_fdr0.05_pearson_mat_bool_filt,2)
    
		write.table(nonsig_vs_fdr0.05_pearson_mat_bool_filt,file=fname,sep="\t",row.names=TRUE)
	
		if(check_cor>=abs.cor.thresh){
	

	#pdf("network_plot.pdf",width=9,height=11)	
	
    pdfname<-paste(networkscope,"network_plot.pdf",sep="")
    pdf(pdfname,width=8,height=10)



#pdfname<-paste(networkscope,"network_plot.tiff",sep="")
#tiff(pdfname,res=300)

	#tiff("network_allrows_sigcols.tiff", width=plots.width,height=plots.height,res=plots.res, compression="lzw")
	#par_rows=1
	#par(mfrow=c(par_rows,1))
	nonsig_vs_fdr0.05_pearson_mat_bool_filt_1<-nonsig_vs_fdr0.05_pearson_mat_bool_filt[-which(duplicated(rnames)==TRUE),]
	set.seed(netrandseed)
	net_result<-try(network(mat=as.matrix(nonsig_vs_fdr0.05_pearson_mat_bool_filt_1), threshold=abs.cor.thresh,color.node = net_node_colors,shape.node = c("rectangle", "circle"),
	color.edge = c("red", "blue"), lwd.edge = 1,show.edge.labels = FALSE, interactive = FALSE,cex.node.name=0.6),silent=TRUE)
    
    if(is(net_result,"try-error")){
        
       
        net_result<-try(network(mat=as.matrix(nonsig_vs_fdr0.05_pearson_mat_bool_filt_1), cutoff=abs.cor.thresh,color.node = net_node_colors,
        shape.node = c("rectangle", "circle"),
        color.edge = c("red", "blue"), lwd.edge = 1,
        show.edge.labels = FALSE, interactive = FALSE,cex.node.name=0.6),silent=TRUE)
    }
    
    print(net_result)
    
	if(net_legend==TRUE){
	 print(legend("bottomright",c("Row #","Colum #"),pch=c(22,19),col=net_node_colors, cex=cex.plots,title="Network matrix values:"))
	}
	#write.graph(net_result$gR, file = "network_cytoscape_format.gml", format = "gml")
	
	cytoscape_fname<-paste("network_Cytoscape_format_maxcor",max.cor.num,"_rowcolnumlabels.gml",sep="")
	write.graph(net_result$gR, file =cytoscape_fname, format = "gml")
    #dev.off()
	
    

        
        
    try(dev.off(),silent=TRUE)
	
	}
		
	}else{
	print("No significant correlations found.")
	}
	}
	
    # #save(net_result,file="metabnet.Rda")
    
		return(nonsig_vs_fdr0.05_pearson_mat_bool_filt)
	
	
	}



docortest<-function(n,r){
	#n<-length(x)
    t<-r*sqrt((n-2)/(1-r^2))
    pvalue<-2*pt(-abs(t),df=n-2)
	
	return(pvalue)
	
}

	
	
	
	get_full_cormat<-function(data_m_fc_withfeats,targeted.index=NA,cor.method="spearman"){
	Sys.sleep(0.2)
	
	allsig_pcornetwork<-{}
	cormat<-{}

	allsig_pcornetwork_fdr0.05<-{}
	
	metab_names<-paste(data_m_fc_withfeats$mz,data_m_fc_withfeats$time,sep="_")
	rnames<-paste("mzid_",seq(1,dim(data_m_fc_withfeats)[1]),sep="")

	if(dim(data_m_fc_withfeats)[1]>1)
	{
	m1<-apply(data_m_fc_withfeats[,-c(1:2)],2,as.numeric)
	

	rownames(m1)=rnames

	data_mt<-t(m1)
	}else{

	m1<-as.numeric(data_m_fc_withfeats[,-c(1:2)])
	

	data_mt<-(m1)
	}
	data_mt<-as.matrix((data_mt))
	
	print(paste("Computing ",cor.method," correlation matrix",sep="")) 
	cormat<-WGCNA::cor(data_mt,use="pairwise.complete.obs",method=cor.method)
	colnames(cormat)<-as.character(metab_names)
	rownames(cormat)<-as.character(metab_names)
    
    num_samp<-dim(data_m_fc_withfeats[,-c(1:2)])[2]
    pearson_resvec<-as.vector(cormat)
    
    pearson_pvalue<-lapply(1:length(pearson_resvec),function(x){
        
        return(WGCNA::corPvalueStudent(pearson_resvec[x],num_samp))
        
        
    })
    
    complete_pearsonpvalue_mat<-unlist(pearson_pvalue)
    dim(complete_pearsonpvalue_mat)<-dim(cormat)
    
    
    cormat<-as.data.frame(cormat)
    cormat<-as.matrix(cormat)
    pearson_Res_all<-{}
    
    cormat<-round(cormat,2)
          
    complete_pearsonpvalue_mat<-as.data.frame(complete_pearsonpvalue_mat)
    
    complete_pearsonpvalue_mat<-round(complete_pearsonpvalue_mat,2)
    
    write.table(cormat,file="correlationmatrix.txt",sep="\t",row.names=TRUE)
    
    write.table(complete_pearsonpvalue_mat,file="pvalue_matrix.txt",sep="\t",row.names=TRUE)
    
    ##save(cormat,file="full_correlation.Rda")
	
        if(is.na(targeted.index)==FALSE){
        
            cormat<-cormat[targeted.index,]
            complete_pearsonpvalue_mat<-complete_pearsonpvalue_mat[targeted.index,]
        
             write.table(cormat,file="targetcorrelationmatrix.txt",sep="\t",row.names=TRUE)
             write.table(complete_pearsonpvalue_mat,file="targetpvalue_matrix.txt",sep="\t",row.names=TRUE)
        
        }
    
		return(cormat)
	}

	get_partial_cornet<-function(data_m_fc_withfeats, sigfeats.index=NA,targeted.index=NA,networkscope="global",cor.method="spearman",abs.cor.thresh=0.4,
	cor.fdrthresh=0.2,outloc,net_node_colors,net_legend=FALSE,netrandseed=555,pcor.method="cor2pcor"){
	
		setwd(outloc)
		l1<-list.files(".")
		
		data_m_fc_withfeats<-as.data.frame(data_m_fc_withfeats)
		
        metab_names<-paste(data_m_fc_withfeats$mz,data_m_fc_withfeats$time,sep="_")
        
		
	
    if(pcor.method=="cor2pcor"){
		
		cormat<-get_full_cormat(data_m_fc_withfeats,cor.method=cor.method)
		
		
		
        fname<-paste("complete_correlation_matrix.txt",sep="")
		
        cormat<-round(cormat,2)
        
        #write.table(cormat,file=fname,sep="\t",row.names=TRUE)
	
		
		
		p1<-cor2pcor(cormat)
    }else{
        
        p1<-pcor.shrink(data_mt)
    }
		colnames(p1)<-as.character(metab_names)
		rownames(p1)<-as.character(metab_names)
		
        p1<-round(p1,2)
        ##save(p1,file="partial_cor.Rda")
		write.table(p1,file="partial_cor.txt",sep="\t",row.names=TRUE)
		
		

		colnames(p1)<-as.character(metab_names)
		rownames(p1)<-as.character(metab_names)
		
		
		if(is.na(targeted.index[1])==FALSE){
			
			networkscope="targeted"
		}
		
		

		
		#only raw p-values
		if(is.na(cor.fdrthresh)==TRUE){
		
		#p1sig<-cor0.test(p1, kappa=16)
		
		fname<-paste("net3_corpval",cor.fdrthresh,".Rda",sep="")
					
						edge.list<-network.test.edges(p1,fdr=FALSE,verbose=FALSE,plot=FALSE)

		
						net<-extract.network(edge.list,cutoff.ggm=0.80,verbose=FALSE)
	
						net2<-net[order(net$node1),]

						net3<-net2[,c(2,3,1,4)]
		
	
						##save(net3,file=fname)
				
				
				
		
		
		
		
		netcormat<-matrix(data=0,nrow=dim(p1)[2],ncol=dim(p1)[2])
		colnames(netcormat)<-colnames(p1)
				rownames(netcormat)<-rownames(p1)
		
		for(i in 1:(dim(net3)[1])){
		
		r<-net3$node1[i]
		c<-net3$node2[i]
		p<-net3$pcor[i]
		
		netcormat[r,c]<-p
		netcormat[c,r]<-p
		}
		diag(netcormat)<-1
		#print(dim(netcormat))
			if(is.na(targeted.index[1])==TRUE){
			
				targeted.index=seq(1,dim(netcormat)[2])
			}
			
				if(is.na(sigfeats.index[1])==FALSE){
					netcormat<-netcormat[c(sigfeats.index),c(targeted.index)]
				
					}else{
					netcormat<-netcormat[,c(targeted.index)]
					}
		
		}else{
		
						
				if(cor.fdrthresh!=(-1) | (is.na(cor.fdrthresh)==TRUE)){
		
				#p1sig<-cor0.test(p1, kappa=16)
				fname<-paste("net3_corpval",cor.fdrthresh,".Rda",sep="")
					
						edge.list<-network.test.edges(p1,fdr=TRUE,verbose=FALSE,plot=FALSE)

				
						net<-extract.network(edge.list,method.ggm="qval",cutoff.ggm=(1-cor.fdrthresh),verbose=FALSE)
			
				net2<-net[order(net$node1),]

				net3<-net2[,c(2,3,1,4,5)]
				
				
				##save(net3,file=fname)
		
		
				
				netcormat<-matrix(data=0,nrow=dim(p1)[2],ncol=dim(p1)[2])
				colnames(netcormat)<-colnames(p1)
				rownames(netcormat)<-rownames(p1)
		
		for(i in 1:(dim(net3)[1])){
		netcormat[net3$node1[i],net3$node2[i]]=net3$pcor[i]
		netcormat[net3$node2[i],net3$node1[i]]=net3$pcor[i]
		}
		diag(netcormat)<-1
			if(is.na(targeted.index)==TRUE){
			
				targeted.index=seq(1,dim(netcormat)[2])
			}
					if(is.na(sigfeats.index)==FALSE){
						netcormat<-netcormat[c(sigfeats.index),c(targeted.index)]
					
						}else{
						netcormat<-netcormat[,c(targeted.index)]
						}
				}else{
			if(is.na(sigfeats.index)==FALSE){
			netcormat<-p1[c(sigfeats.index),c(targeted.index)]
		
			}else{
			netcormat<-p1[,c(targeted.index)]
			}
			}
		
		
		
		
		}
		
		if(length(sigfeats.index)<dim(data_m_fc_withfeats)[1]){
		fname<-paste("significant_correlations_for_network.txt",sep="")
		
        # netcormat<-round(netcormat,2)
        
	write.table(netcormat,file=fname,sep="\t",row.names=TRUE)
	
	net_fname<-paste("partial_corsig",cor.fdrthresh,"network.tiff",sep="")
	
	max_cor_check<-max(abs(netcormat),na.rm=TRUE,warnings=FALSE)
	
	#if(max_cor_check>0)
	
	if(max_cor_check>=abs.cor.thresh){
	
	netcormat<-t(netcormat)
	colnames(netcormat)<-paste("Y",seq(1,dim(netcormat)[2]),sep="")
	rownames(netcormat)<-paste("X",seq(1,dim(netcormat)[1]),sep="")
    #tiff(net_fname, width=plots.width,height=plots.height,res=plots.res, compression="lzw")
	pdf(net_fname,width=8,height=10)
    par_rows=1
	par(mfrow=c(par_rows,1))
	set.seed(netrandseed)
	net_result<-try(network(mat=as.matrix(netcormat), threshold=abs.cor.thresh,color.node = net_node_colors,shape.node = c("rectangle", "circle"),
	color.edge = c("red", "blue"), lwd.edge = 1,show.edge.labels = FALSE, interactive = FALSE,cex.node.name=0.6),silent=TRUE)
    
    if(is(net_result,"try-error")){
        
        net_result<-try(network(mat=as.matrix(nonsig_vs_fdr0.05_pearson_mat_bool_filt_1), cutoff=abs.cor.thresh,color.node = net_node_colors,
        shape.node = c("rectangle", "circle"),
        color.edge = c("red", "blue"), lwd.edge = 1,
        show.edge.labels = FALSE, interactive = FALSE,cex.node.name=0.6),silent=TRUE)
    }
    
	#legend("bottomright",c("Row #","Colum #"),pch=c(22,19),col=net_node_colors, cex=0.8)
	print(net_result)
    
	if(net_legend==TRUE){
	 print(legend("bottomright",c("Row #","Colum #"),pch=c(22,19),col=net_node_colors, cex=cex.plots,title="Network matrix values:"))
	}
    

# #save(net_result,"metabnet.Rda")
   
    
	write.graph(net_result$gR, file = "network_cytoscape_format.gml", format = "gml")
	try(dev.off(),silent=TRUE)
	}else{
	print("No significant correlations found.")
	}
	}else{
	
		print("Network plot can not be generated. Rows and columns need to be exclusive.")
	}
	net4<-netcormat
	
	if(FALSE){
	if(is.na(cor.fdrthresh)==FALSE){
	if(cor.fdrthresh>(-1)){
	
	net4<-cbind(data_m_fc_withfeats[net3$node1,c(1:2)],data_m_fc_withfeats[net3$node2,c(1:2)])
	net4<-cbind(net4,net3[,c(3:4)])
	
	allsig_pcornetwork_fdr0.05<-net4 #[which(net3[,7]<cor.fdrthresh),]
	
	
    #allsig_pcornetwork_fdr0.05<-round(allsig_pcornetwork_fdr0.05,2)
    
	fname<-paste("significant_correlations_",networkscope,"feats_corrsig",cor.fdrthresh,".txt",sep="")
	write.table(allsig_pcornetwork_fdr0.05,file=fname,sep="\t",row.names=FALSE)
	}else{
		net4<-netcormat
	}
	
	}else{
		net4<-cbind(data_m_fc_withfeats[net3$node1,c(1:2)],data_m_fc_withfeats[net3$node2,c(1:2)])
		net4<-cbind(net4,net3[,c(3:4)])

		allsig_pcornetwork_fdr0.05<-net4 #[which(net3[,7]<cor.fdrthresh),]
        
        #allsig_pcornetwork_fdr0.05<-round(allsig_pcornetwork_fdr0.05,2)


		fname<-paste("significant_correlations_",networkscope,"feats_corrsig",cor.fdrthresh,".txt",sep="")
		write.table(allsig_pcornetwork_fdr0.05,file=fname,sep="\t",row.names=FALSE)
	}
	}
	
		return(net4)
	}
	
do_rf<-function(X,classlabels, ntrees=1000, analysismode="classification"){
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
         set.seed(290875)
     rf2 <- randomForest(as.factor(Targetvar) ~ .,data=dataA, importance=TRUE,proximity=FALSE,ntree=ntrees,keep.forest=FALSE)
   	}else{
   		set.seed(290875)
   		rf2 <- randomForest(Targetvar ~ .,data=dataA, importance=TRUE,proximity=FALSE,ntree=ntrees,keep.forest=FALSE)
   		}
    
     #based on permutation importance;
     varimp_res<-round(randomForest::importance(rf2,type=1,scale=FALSE), 2)
     varimp_res_scaled<-round(randomForest::importance(rf2,type=1,scale=TRUE), 2)
	rm(dataA)

	#return(varimp_res)
	return(list("rf_model"=rf2,"rf_varimp"=varimp_res,"rf_varimp_scaled"=varimp_res_scaled))
}

do_rf_conditional<-function(X,classlabels, ntrees=1000, analysismode="classification"){
	
	#classlabels<-c(rep("A",6),rep("B",6))
	
	dataA<-t(X)
	dataA<-as.data.frame(dataA)
	dataA<-cbind(classlabels,dataA)
	dataA<-as.data.frame(dataA)
	
	colnames(dataA)<-c("Targetvar",paste("mz",seq(1,dim(X)[1]),sep=""))
	
	dataA<-as.data.frame(dataA)
	
	attach(dataA,warn.conflicts=FALSE)
	
	mtry_num<-if (!is.null(classlabels) && !is.factor(classlabels))
                  max(floor(ncol(X)/3), 1) else floor(sqrt(ncol(X)))

                  #return(dataA)
	set.seed(290875)
   
	if(analysismode=="classfication"){
        set.seed(290875)
        rf1<-cforest(factor(Targetvar)~.,data=dataA,control=cforest_control(teststat = "quad", testtype = "Univ",mincriterion = 0,savesplitstats = FALSE,ntree = ntrees, mtry = ceiling(sqrt(nvar)), replace = TRUE,fraction = 0.632, trace = FALSE))
        
    	}else{
   		set.seed(290875)
        #rf1<-cforest(Targetvar~.,data=dataA,control=cforest_control(teststat = "max",testtype = "Teststatistic",mincriterion = qnorm(0.9),savesplitstats = FALSE,ntree = ntrees, mtry = NULL, replace = TRUE,fraction = 0.632, trace = FALSE))
                     
                     #rf1<-cforest(Targetvar~.,data=dataA,control=cforest_unbiased(mtry = ceiling(sqrt(nvar)), ntree =ntrees))
                     #  ntree = ntrees, mtry = NULL, replace = TRUE,
                     #fraction = 0.632, trace = FALSE))
                     
                     rf1<-cforest(factor(Targetvar)~.,data=dataA,control=cforest_control(teststat = "quad", testtype = "Univ",mincriterion = 0,savesplitstats = FALSE,ntree = ntrees, mtry = ceiling(sqrt(nvar)), replace = TRUE,fraction = 0.632, trace = FALSE))
                     
   		}
	
    if(analysismode=="classfication"){
	set.seed(290875)
    varimp_res<-varimp(rf1,conditional=TRUE) #varimpAUC(rf1,conditional=TRUE)
    }else{
        varimp_res<-varimp(rf1,conditional=TRUE)
    }
	print(varimp_res)
    rm(dataA) 
        
	#return(varimp_res)
     return(list("rf_model"=rf1,"rf_varimp"=varimp_res))

}

varimp_parallel<-function (object, mincriterion = 0, conditional = FALSE, threshold = 0.2,
nperm = 1, OOB = TRUE, pre1.0_0 = conditional)
{
    response <- object@responses
    if (length(response@variables) == 1 && inherits(response@variables[[1]],
    "Surv"))
    return(varimpsurv(object, mincriterion, conditional,
    threshold, nperm, OOB, pre1.0_0))
    input <- object@data@get("input")
    xnames <- colnames(input)
    inp <- initVariableFrame(input, trafo = NULL)
    y <- object@responses@variables[[1]]
    if (length(response@variables) != 1)
    stop("cannot compute variable importance measure for multivariate response")
    if (conditional || pre1.0_0) {
        if (!all(complete.cases(inp@variables)))
        stop("cannot compute variable importance measure with missing values")
    }
    CLASS <- all(response@is_nominal)
    ORDERED <- all(response@is_ordinal)
    if (CLASS) {
        error <- function(x, oob) mean((levels(y)[sapply(x, which.max)] !=
        y)[oob])
    }
    else {
        if (ORDERED) {
            error <- function(x, oob) mean((sapply(x, which.max) !=
            y)[oob])
        }
        else {
            error <- function(x, oob) mean((unlist(x) - y)[oob]^2)
        }
    }
    w <- object@initweights
    if (max(abs(w - 1)) > sqrt(.Machine$double.eps))
    warning(sQuote("varimp"), " with non-unity weights might give misleading results")
    perror <- matrix(0, nrow = nperm * length(object@ensemble),
    ncol = length(xnames))
    colnames(perror) <- xnames
    
    l1<-lapply(1:length(object@ensemble),function(b){
        # for (b in 1:length(object@ensemble)) {
        tree <- object@ensemble[[b]]
        if (OOB) {
            oob <- object@weights[[b]] == 0
        }
        else {
            oob <- rep(TRUE, length(y))
        }
        p <- predict(tree, inp, mincriterion, -1L)
        eoob <- error(p, oob)
        #for (j in unique(varIDs(tree))) {
        l2<-lapply(unique(varIDs(tree)),function(j){
            for (per in 1:nperm) {
                if (conditional || pre1.0_0) {
                    tmp <- inp
                    ccl <- create_cond_list(conditional, threshold,
                    xnames[j], input)
                    if (length(ccl) < 1) {
                        perm <- sample(which(oob))
                    }
                    else {
                        perm <- conditional_perm(ccl, xnames, input,
                        tree, oob)
                    }
                    tmp@variables[[j]][which(oob)] <- tmp@variables[[j]][perm]
                    p <- predict(tree, tmp, mincriterion, -1L)
                }
                else {
                    p <- predict(tree, inp, mincriterion, as.integer(j))
                }
                #perror[(per + (b - 1) * nperm), j] <- (error(p,oob) - eoob)
                perror[(per), j] <- (error(p,oob) - eoob)
            }
            return(perror)
        })
        perror<-ldply(l2,rbind)
        return(perror)
    })
    perror<-ldply(l1,rbind)
    perror <- as.data.frame(perror)
    return(MeanDecreaseAccuracy = colMeans(perror))
}



do_mars<-function(X,classlabels, analysismode="classification",kfold=10){
	
	dataA<-t(X)
	dataA<-as.data.frame(dataA)
	
	
		if(analysismode=="classification"){
		
	dataA<-cbind(classlabels,dataA)
	dataA<-as.data.frame(dataA)
	
	cnames<-c("Targetvar",paste("mz",seq(1,dim(X)[1]),sep=""))
	
	colnames(dataA)<-as.character(cnames)
	
	if(dim(dataA)[1]<20){kfold=1}
	
	dataA<-as.data.frame(dataA)

	attach(dataA,warn.conflicts=FALSE)
	
	#print(dim(dataA))
	

		mars_res<-new("list")
		

	mars_res[[1]]<-earth(formula=factor(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=1,ncross=kfold,nfold=kfold)
	mars_res[[2]]<-earth(formula=factor(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=2,ncross=kfold,nfold=kfold)
	mars_res[[3]]<-earth(formula=factor(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=3,ncross=kfold,nfold=kfold)
	
	gcv_list<-c(mars_res[[1]]$gcv,mars_res[[2]]$gcv,mars_res[[3]]$gcv)

	min_gcv<-which(gcv_list==min(gcv_list,na.rm=TRUE))

	mars_res<-earth(formula=factor(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=min_gcv[1],ncross=kfold,nfold=kfold)

    	}else{
   		if(analysismode=="regression"){
		
		mars_res<-new("list")
   		
   		classlabels<-as.data.frame(classlabels)
   		#print(dim(classlabels))
	
   		if(dim(classlabels)[2]>1){
   		classlabels<-apply(classlabels,2,as.numeric)
		classlabels<-as.numeric(classlabels[,1])
   		}else{
   			classlabels<-as.numeric(classlabels[,1])
   			#classlabels<-apply(classlabels,2,as.numeric)

   			}
	
	
	dataA<-cbind(classlabels,dataA)
    #dataA<-as.data.frame(dataA)
	
	cnames<-c("Targetvar",paste("mz",seq(1,dim(X)[1]),sep=""))
	
    #cnames<-c(paste("mz",seq(1,dim(X)[1]),sep=""))
	
	colnames(dataA)<-as.character(cnames)
	
	if(dim(dataA)[1]<20){kfold=1}
	
	dataA<-as.data.frame(dataA)

	attach(dataA,warn.conflicts=FALSE)

#mars_res[[1]]<-earth(y=(classlabels),x=dataA,Use.beta.cache=FALSE,degree=1,ncross=kfold,nfold=kfold)
#		mars_res[[2]]<-earth(y=(classlabels),x=dataA,Use.beta.cache=FALSE,degree=2,ncross=kfold,nfold=kfold)
#		mars_res[[3]]<-earth(y=(classlabels),x=dataA,Use.beta.cache=FALSE,degree=3,ncross=kfold,nfold=kfold)

    mars_res[[1]]<-earth(formula=as.numeric(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=1,ncross=kfold,nfold=kfold)
    mars_res[[2]]<-earth(formula=as.numeric(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=2,ncross=kfold,nfold=kfold)
    mars_res[[3]]<-earth(formula=as.numeric(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=3,ncross=kfold,nfold=kfold)
		gcv_list<-c(mars_res[[1]]$gcv,mars_res[[2]]$gcv,mars_res[[3]]$gcv)
		min_gcv<-which(gcv_list==min(gcv_list,na.rm=TRUE))

		mars_res<-earth(formula=as.numeric(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=min_gcv[1],ncross=kfold,nfold=kfold)
		#mars_res<-earth(formula=(Targetvar) ~ .,data=dataA,Use.beta.cache=FALSE,degree=3,ncross=kfold,nfold=kfold)
        #mars_res<-earth(y=(classlabels),x=dataA,Use.beta.cache=FALSE,degree=min_gcv[1],ncross=kfold,nfold=kfold)
	
	
   		}else{
   			stop("Invalid analysismode entered. Please use classification or regression.")
   			}
		}
	dataA<-NULL
	detach(dataA)
	rm(dataA)

    #save(mars_res,file="mars_res.Rda")
	#varimp_res <- evimp(mars_res)

	#if(FALSE)
	{
	varimp_res <- evimp(mars_res,trim=FALSE)

	#print(varimp_res[1:10,])
    varimp_marsres1<-varimp_res #[order(varimp_res[,2],decreasing=TRUE),]
		
		mars_mznames<-rownames(varimp_marsres1)
		
		g1<-grep(pattern="NA",x=mars_mznames)
		if(length(g1)>0){
			varimp_marsres1<-varimp_marsres1[-g1,]
		}
	}

	
return(list("mars_model"=mars_res,"mars_varimp"=varimp_marsres1))
	
}

getSumreplicates<-function(curdata,alignment.tool,numreplicates,numcluster,rep.max.missing.thresh,summary.method="mean",summary.na.replacement="zeros",missing.val=0)
{
		 mean_replicate_difference<-{}
		sd_range_duplicate_pairs<-{}
		  #print(alignment.tool)
		if(alignment.tool=="apLCMS")
		{
		      col_end=2
		}
		else
		{
		      if(alignment.tool=="XCMS")
		      {
			    col_end=2
		      }
		      else
		      {
			    stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
		      }
		}
		
		curdata_mz_rt_info=curdata[,c(1:col_end)]
		curdata=curdata[,-c(1:col_end)]
		
		
		
		cl<-parallel::makeCluster(numcluster)
		numfeats=dim(curdata)[1]
		numsamp=dim(curdata)[2]	

		clusterEvalQ(cl, "getSumreplicateschild")
		sub_samp_list<-list()
	
		sampcount=1
		for(samp in seq(1,(numsamp),numreplicates))
                {
                        i=samp
                        j=i+numreplicates-1
			if(dim(curdata[,c(i:j)])[1]>0){
                        sub_samp_list[[sampcount]]=curdata[,c(i:j)]
                        }
			sampcount=sampcount+1
                }
		
		avg.res<-parSapply(cl,sub_samp_list,getSumreplicateschild,alignment.tool=alignment.tool,numreplicates=numreplicates,rep.max.missing.thresh=rep.max.missing.thresh,method=summary.method,missing.val=missing.val)
		#avg.res<-getAvgreplicateschild(sub_samp_list[[1]],alignment.tool,numreplicates)
		#print("done")
		
		
		stopCluster(cl)
		
		
			
			final_set<-as.data.frame(avg.res)
			colnames_data<-colnames(curdata)
			colnames_data<-colnames_data[seq(1,(numsamp),numreplicates)]
			colnames(final_set)<-colnames_data
			rownames(final_set)=NULL
			#final_set<-cbind(curdata_mz_rt_info,final_set)
		
		final_set<-apply(final_set,2,as.numeric)
        #	write.table(final_set,file="final_Set.txt",sep="\t",row.names=FALSE)
		
		if(summary.na.replacement=="zeros"){
			
			final_set<-replace(final_set,which(is.na(final_set)==TRUE),0)
		}else{
			if(summary.na.replacement=="halfsamplemin"){
				
			
		final_set<-apply(final_set,2,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-min(x,na.rm=TRUE)/2}; return(x)})
		}else{
			
			if(summary.na.replacement=="halfdatamin"){
				
				
			min_val<-min(final_set,na.rm=TRUE)*0.5
			final_set<-replace(final_set,which(is.na(final_set)==TRUE),min_val)

			#data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-min(x,na.rm=TRUE)/2}; return(x)})
			}else{
			if(summary.na.replacement=="halffeaturemin"){
				
				
			final_set<-apply(final_set,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-min(x,na.rm=TRUE)/2}; return(x)})
			final_set<-t(final_set)
			}
			}
		}
			
			
		}
		final_set<-as.data.frame(final_set)
		return(final_set)
}


getSumreplicateschild<-function(curdata,alignment.tool,numreplicates,rep.max.missing.thresh,method="mean",missing.val)
{
       #curdata<-t(curdata)
       #write.table(curdata,file="test.txt",sep="\t",row.names=FALSE)
        numfeats=dim(curdata)[1]
        numsamp=dim(curdata)[2]
     # if(FALSE){
        resvec_1<-lapply(1:numfeats,function(r)
        {
                newrow={}
                finalmat={}
                #for(samp in seq(1,(numsamp),numreplicates))
                {
                       # i=samp
                        #j=i+numreplicates-1

                        curdata_int=curdata[r,]
			
            #if(is.na(missing.val)==FALSE){
            #           check_zeros=which(curdata_int==missing.val)
            #           }else{
                        	check_zeros=which(is.na(curdata_int)==TRUE)
                            #           	}
                        na_thresh=round(rep.max.missing.thresh*numreplicates)
                        
                        
                        if(length(check_zeros)>na_thresh)
                        {
                                meanval<-missing.val
                        }
                        else
                        {
                                #temporarily replace the missing intensities, set to 0 in apLCMS,
                                #with mean intensity value of the corresponding replicates (with non-zero values)
                                #curdata_int[check_zeros]=mean(t(curdata_int[-c(check_zeros)]))
                                 if(length(check_zeros)>0)
                                {
                                		if(method=="mean"){
                                			 meanval<-mean(t(curdata_int[-check_zeros]),na.rm=TRUE)
                                		}else{
                                        meanval<-median(t(curdata_int[-check_zeros]),na.rm=TRUE)
                                		}
                                }
                                else
                                {
                                		if(method=="mean"){
                                				meanval<-mean(t(curdata_int),na.rm=TRUE)
                                			}else{
                                        meanval<-median(t(curdata_int),na.rm=TRUE)
                                		}
                                }

                        }
                        newrow<-cbind(newrow,meanval)
                }


                finalmat<-rbind(finalmat, newrow)
                return(finalmat)
        })
        
        #colnames(final_set)<-colnames_data
        #rownames(final_set)=NULL
        return(resvec_1)
	
	
	
}


metabnet<-function(feature_table_file,target.metab.file,sig.metab.file,class_labels_file=NA,parentoutput_dir,num_replicates=3,cor.method="spearman",abs.cor.thresh=0.4,cor.fdrthresh=0.05,
target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=100,feat.filt.thresh=NA,summarize.replicates=TRUE,summary.method="mean",all.missing.thresh=0.5, group.missing.thresh=0.7,
log2transform=TRUE,medcenter=TRUE,znormtransform=FALSE,quantile_norm=TRUE,lowess_norm=FALSE,madscaling=FALSE,missing.val=0, networktype="complete", samplermindex=NA,
rep.max.missing.thresh=0.3,summary.na.replacement="zeros",net_node_colors=c("pink","skyblue"),net_legend=FALSE,netrandseed=555){

		options(warn=-1)
		
		

		data_matrix<-data_preprocess(Xmat=NA,Ymat=NA,feature_table_file,parentoutput_dir,class_labels_file,num_replicates=num_replicates,feat.filt.thresh,summarize.replicates,summary.method, all.missing.thresh, group.missing.thresh,
log2transform,medcenter,znormtransform,quantile_norm,lowess_norm,madscaling,missing.val, samplermindex,rep.max.missing.thresh,summary.na.replacement)

	data_matrix<-data_matrix$data_matrix_afternorm_scaling
	


		
		dir.create(parentoutput_dir)
		setwd(parentoutput_dir)
		
	data_m<-data_matrix[,-c(1:2)]
	
	#print(data_matrix[1:3,1:10])
	
	
	
	if(is.na(sig.metab.file)==FALSE){
						goodfeats<-read.table(sig.metab.file,sep="\t",header=TRUE)
			}else{
		goodfeats<-data_matrix
		}
		
	
	if(is.na(target.metab.file)==FALSE){
	dataA<-read.table(target.metab.file,sep="\t",header=TRUE)
	
	dataA<-as.data.frame(dataA)
	outloc<-paste(parentoutput_dir,"/Stage4","/",sep="")
	dir.create(outloc)
	print(paste("Searching for metabolites matching target list",sep=""))
	
g1<-getVenn(dataA=dataA,name_a="TargetSet",name_b="ExperimentalSet",dataB=data_matrix[,c(1:2)],mz.thresh=target.mzmatch.diff,time.thresh=target.rtmatch.diff,
xMSanalyzer.outloc=outloc,alignment.tool=NA)
#names(g1)


	
	if(length(g1$common)>1){
	
	if(is.na(sig.metab.file)==FALSE){
	com_mzs<-find.Overlapping.mzs(dataA=data_matrix,dataB=goodfeats,mz.thresh=1,time.thresh=1,alignment.tool=NA)
	
	sigfeats.index<-com_mzs$index.A #which(data_matrix$mz%in%goodfeats$mz)
	print(paste(length(unique(sigfeats.index))," selected features",sep=""))
	

	}else{
		sigfeats.index<-NA #seq(1,dim(data_matrix)[1])
		#sigfeats.index<-seq(1,50)
		
		}
	
	print(paste(length(unique(g1$common$index.B))," metabolites matched the target list",sep=""))
	
		print(paste("Generating targeted network",sep=""))
	
	#print(data_matrix[1:3,1:5])
	#print(length(sigfeats.index))
if(networktype=="complete"){
	targetedan_fdr<-do_cor(data_matrix,subindex=sigfeats.index,targetindex=g1$common$index.B,outloc,networkscope="targeted",cor.method,
	abs.cor.thresh,cor.fdrthresh,max.cor.num,net_node_colors,net_legend,netrandseed)
    
    #  pdf("Cornetworkplot.pdf")
    #load("metabnet.Rda")
    #print(plot(net_result))
    #dev.off()
    
	}else{
	if(networktype=="GGM"){
	targetedan_fdr<-get_partial_cornet(data_matrix, sigfeats.index,targeted.index=g1$common$index.B,networkscope="targeted",cor.method,abs.cor.thresh,
	cor.fdrthresh,outloc=outloc,net_node_colors,net_legend,netrandseed)
    
   
    
	}else{
		print("Invalid option. Please use complete or GGM.")
	}
	}
	}else{
	print(paste("Targeted metabolites were not found.",sep=""))
	}
	}else{
		
		outloc<-paste(outloc,"/MWASresults","/",sep="")
        dir.create(outloc)
        setwd(outloc)
		
			if(is.na(sig.metab.file)==FALSE){
	com_mzs<-find.Overlapping.mzs(dataA=data_matrix,dataB=goodfeats,mz.thresh=1,time.thresh=1,alignment.tool=NA)
	
	sigfeats.index<-com_mzs$index.A #which(data_matrix$mz%in%goodfeats$mz)
	}else{
		sigfeats.index<-NA #seq(1,dim(data_matrix)[1])
		}
		
		if(networktype=="complete"){
		
	
	#targetedan_fdr<-do_cor(data_matrix,sigfeats.index,outloc,networkscope="global",cor.method,abs.cor.thresh,cor.fdrthresh,max.cor.num)
	targetedan_fdr<-do_cor(data_matrix,subindex=sigfeats.index,targetindex=NA,outloc,networkscope="global",cor.method,abs.cor.thresh,cor.fdrthresh,
	max.cor.num,net_node_colors,net_legend,netrandseed)

if(FALSE){
pdf("Cornetworkplot.pdf")
load("metabnet.Rda")
print(plot(net_result))
dev.off()
        }

	}else{
	if(networktype=="GGM"){
	targetedan_fdr<-get_partial_cornet(data_matrix, sigfeats.index,targeted.index=NA,networkscope="global",cor.method,abs.cor.thresh,
	cor.fdrthresh,outloc=outloc,net_node_colors,net_legend,netrandseed)
	
    if(FALSE){
    pdf("GGMnetworkplot.pdf")
    load("metabnet.Rda")
    print(plot(net_result))
    dev.off()
    }
    
	}else{
		print("Invalid option")
	}
	}
			
		}
		
		
		print("Processing complete.")
		return(targetedan_fdr)
}



get_boxplots<-function(X=NA,Y=NA,feature_table_file,parentoutput_dir,class_labels_file,boxplot.col.opt="white",sample.col.opt="rainbow",alphacol=0.3,newdevice=TRUE,cex=0.6,replace.by.NA=FALSE,pairedanalysis=FALSE,filename="",ylabel="Intensity")
{

get_boxplots_child(X=X,Y=Y,feature_table_file=feature_table_file,parentoutput_dir=parentoutput_dir,class_labels_file=class_labels_file,boxplot.col.opt,sample.col.opt=sample.col.opt,alphacol=alphacol,newdevice=newdevice,cex=cex,replace.by.NA=replace.by.NA,pairedanalysis=pairedanalysis,filename=filename,ylabel=ylabel)

}

get_boxplots_child<-function(X,Y,feature_table_file,parentoutput_dir,class_labels_file,boxplot.col.opt="journal",alphacol=0.3,newdevice=TRUE,cex=0.6,replace.by.NA=FALSE,pairedanalysis=FALSE,filename="",ylabel="Intensity",...)
{
		

			if(is.na(X)==TRUE){
			data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE)
			
			}else{
			
				data_matrix<-X
			}
			if(is.na(Y)==TRUE){
			classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
			
			}else{
			classlabels<-Y
			
			}
            rm(X)
            rm(Y)
            
dir.create(parentoutput_dir)
setwd(parentoutput_dir)

data_m<-data_matrix[,-c(1:2)]

data_m<-as.matrix(data_m)

mzvec<-data_matrix[,1]
timevec<-data_matrix[,2]
goodfeats<-data_m			
rm(data_m)



if(dim(classlabels)[2]>2){
    
    print("More than two columns found in the class labels file.")
    if(pairedanalysis==TRUE){
        
        Class<-classlabels[,3]
        
        classlabels<-classlabels[,-c(2)]
        
        
    }else{
    
        Class<-paste(classlabels[,2],":",classlabels[,3],sep="") #classlabels[,2]:classlabels[,3]
    
    }
}else{
    
    Class<-classlabels[,2]
}




class_labels_levels<-levels(as.factor(Class))
ordered_labels<-Class

  class_label_alphabets<-paste("C",1:length(class_labels_levels),sep="") #c("A","B","C","D","E","F","G","H","I","J","K","L","M")
  
  if(is.na(boxplot.col.opt)==TRUE){
      
      col_vec<-rep(c("white"),length(class_labels_levels))
      boxplot.col.opt<-col_vec
  }

    if(boxplot.col.opt=="default"){

    col_vec<-c("#CC0000","#AAC000","blue","mediumpurple4","mediumpurple1","blueviolet","cornflowerblue","cyan4","skyblue",
    "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
    "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
    "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
  
		   }else{ 
		if(boxplot.col.opt=="topo"){
			#col_vec<-topo.colors(256) #length(class_labels_levels)) 
			
			#col_vec<-col_vec[seq(1,length(col_vec),)]
			
			col_vec <- topo.colors(length(class_labels_levels), alpha=alphacol)
		}else{
				if(boxplot.col.opt=="heat"){
					#col_vec<-heat.colors(256) #length(class_labels_levels))
					
					col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
					}else{
								if(boxplot.col.opt=="rainbow"){
									#col_vec<-heat.colors(256) #length(class_labels_levels))
									col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
									
									#col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
								}else{
									
										if(boxplot.col.opt=="terrain"){
									#col_vec<-heat.colors(256) #length(class_labels_levels))
									#col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
									
									col_vec <- cm.colors(length(class_labels_levels), alpha=alphacol)
                                        }else{
                                            
                                            if(boxplot.col.opt=="colorblind"){
                                                #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                                                # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                                                
                                                if(length(class_labels_levels)<9){
                                                    
                                                    col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7", "grey57")
                                                    
                                                }else{
                                                    
                                                    #   col_vec<-colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels))
                                                    col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                                    "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                                                    
                                                }
                                                
                                                
                                            }else{
                                                
                                                check_brewer<-grep(pattern="brewer",x=boxplot.col.opt)
                                                
                                                if(length(check_brewer)>0){
                                                    
                                                     boxplot.col.opt=gsub(x=boxplot.col.opt,pattern="brewer.",replacement="")
                                                    col_vec <- colorRampPalette(brewer.pal(10, boxplot.col.opt))(length(class_labels_levels))
                                                    
                                                }else{
                                                    
                                                    if(boxplot.col.opt=="journal"){
                                                        
                                                        col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                                        "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                                        "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                                        
                                                        "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                                        "#E64B3519","#4DBBD519","#631879E5","grey75")
                                                        
                                                    }else{
                                                        col_vec <-boxplot.col.opt
                                                    }
                                                    
                                                }
                                                
                                            }
                                        }
									
											
									}
							
						}
			
		}	
}

			
           
						ordered_labels={}
							num_samps_group<-new("list")
							num_samps_group[[1]]<-0
							groupwiseindex<-new("list")
							groupwiseindex[[1]]<-0
						
							for(c in 1:length(class_labels_levels))
							{
								
								classlabels_index<-which(Class==class_labels_levels[c])
								#ordered_labels<-c(ordered_labels,as.character(classlabels[classlabels_index,2]))
								num_samps_group[[c]]<-length(classlabels_index)
								groupwiseindex[[c]]<-classlabels_index
							}
			
			sampleclass<-{}
						patientcolors<-{}
						
if(length(mzvec)>4){
max_per_row<-3


par_rows<-ceiling(9/max_per_row)

}else{
max_per_row<-length(mzvec)
par_rows<-1
}

#class_labels_levels<-paste("x",seq(1,length(class_labels_levels)),sep="")

file_ind<-0
			boxplots_fname<-paste("boxplots",filename,".pdf",sep="")
				#tiff(boxplots_fname, width=plots.width,height=plots.height,res=plots.res, compression="lzw")
				
				#tiff(boxplots_fname, width=2000,height=3000,res=plots.res, compression="lzw")
                
                if(newdevice==TRUE){
				pdf(boxplots_fname)
                }
                #par(mfrow=c(par_rows,max_per_row))
            
            par(mfrow=c(1,1),family="sans",cex=cex)
            
		for(m in 1:dim(goodfeats)[1])
		{
            
			if(m%%9==0){
				
				file_ind<-file_ind+1
				boxplots_fname<-paste("boxplots_file",file_ind,".tiff",sep="")
		
			}
			
			round_mzval<-sprintf("%.4f",mzvec[m])
            
            round_timeval<-sprintf("%.1f",timevec[m])
            
            mzname<-paste("mz_time: ",round_mzval,"_",round_timeval,sep="")


			if(length(class_labels_levels)>=2)
			{
			
					if(length(class_labels_levels)>=1)
					{
					
					t1<-table(sampleclass)
					cur_d<-{}
					for(c in 1:length(class_labels_levels))
					{
                
					num_samps_group[[1]]<-t1[1]
					
					cvec<-as.vector(t(goodfeats[m,c(groupwiseindex[[c]])]))
					
					cvec<-replace_outliers(cvec,replace.by.NA)
					cur_d<-cbind(cur_d,cvec)
					}
					
					cur_d<-as.data.frame(cur_d)
					
					
					colnames(cur_d)<-NULL
					cur_d<-round(cur_d,2)
					
                   
            
                    
                   
                    
					par(mfrow=c(1,1),family="sans",cex=cex)
                    #,ylim=c(0,max(cur_d,na.rm=TRUE))
                    max_yval=max(cur_d,na.rm=TRUE)[1]
                    
                    w <- 0.1
                    par(omd=c(0, 1-w, 0, 1))
                    
                    boxplot(cur_d,ylab=ylabel,main=mzname,xaxt="n",cex.main=0.8,col=col_vec) #,ylim=range(pretty(c(0,max_yval))))
					
					for(i in 1:length(class_labels_levels)){
                        axis(side=1,at=c(i),labels=class_labels_levels[i], col=col_vec[i],cex.axis=cex,cex.names=cex)
					
					
					}
                    
                    print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec[1:length(class_labels_levels)],pch = rep(19,length(col_vec[1:length(class_labels_levels)])), pt.cex = 0.6, title = "Class",cex=0.8))
                    
					
                }
					
					
					
			}else{
                par(mfrow=c(1,1),family="sans",cex=cex)
                #,ylim=c(0,max(c(x1,x2),na.rm=TRUE))
               max_yval=max(c(x1,x2),na.rm=TRUE)[1]
               
               w <- 0.1
               par(omd=c(0, 1-w, 0, 1))
               
               boxplot(x1,x2, ylab=ylabel,main=mzname,cex.main=0.8,col=col_vec) #,ylim=range(pretty(c(0,max_yval))))
            axis(side=1,at=seq(1),labels=class_labels_levels[1], col=col_vec[1])
			axis(side=1,at=2,labels=class_labels_levels[2],col=col_vec[2])
            
            print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec[1:length(class_labels_levels)],pch = rep(19,length(col_vec[1:length(class_labels_levels)])), pt.cex = 0.6, title = "Class",cex=0.8))
            
			}
		}
        if(newdevice==TRUE){
		try(dev.off(),silent=TRUE)
        }
        
        par(mfrow=c(1,1))
	}



get_pcascoredistplots<-function(X=NA,Y=NA,feature_table_file,parentoutput_dir,class_labels_file,sample.col.opt="rainbow",plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,col_vec,pairedanalysis=FALSE,pca.cex.val=3,legendlocation="topright",pca.ellipse=TRUE,ellipse.conf.level=0.95,filename="all",paireddesign=NA,error.bar=TRUE,lineplot.col.opt="black",lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"))
{


    res<-get_pcascoredistplots_child(X=X,Y=Y,feature_table_file=feature_table_file,parentoutput_dir=parentoutput_dir,class_labels_file=class_labels_file,sample.col.opt=sample.col.opt,plots.width=plots.width,plots.height=plots.height,plots.res=plots.res, alphacol=alphacol,col_vec=col_vec,pairedanalysis=pairedanalysis,pca.cex.val=pca.cex.val,legendlocation=legendlocation,pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,filename=filename,paireddesign=paireddesign,error.bar=error.bar,lineplot.col.opt=lineplot.col.opt,lineplot.lty.option=lineplot.lty.option)

    return(res)
}

get_pcascoredistplots_child<-function(X,Y,feature_table_file,parentoutput_dir,class_labels_file,sample.col.opt="rainbow",plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,col_vec=NA,pairedanalysis=FALSE,pca.cex.val=4,legendlocation="topright",pca.ellipse=TRUE,ellipse.conf.level=0.95,filename="all",paireddesign=NA,error.bar=TRUE,lineplot.col.opt="black",lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"))
{
    
                if(is.na(X)==TRUE){
                data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE)
                
                }else{
                
                    data_matrix<-X
                }
                if(is.na(Y)==TRUE){
                classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
                
                }else{
                
                    classlabels<-Y
                }
                #print("here")
    dir.create(parentoutput_dir)
    setwd(parentoutput_dir)

    data_m<-data_matrix[,-c(1:2)]

    data_m<-as.matrix(data_m)

    mzvec<-data_matrix[,1]
    rtvec<-data_matrix[,2]

    rnames<-paste(mzvec,rtvec,sep="_")

    rownames(data_m)<-as.character(rnames)


    classlabels_orig<-classlabels


    if(dim(classlabels)[2]>2){

        classgroup<-paste(classlabels_orig[,2],":",classlabels_orig[,3],sep="") #classlabels_orig[,2]:classlabels_orig[,3]
    do_pca_anova=FALSE
    }else{

        classgroup<-classlabels_orig[,2]
        
        do_pca_anova=TRUE
    }

 
   print(dim(classlabels_orig))
   print(head(classlabels_orig))
   
    class_labels_levels<-levels(as.factor(classgroup))
    ordered_labels<-classgroup
    
    class_label_alphabets<-paste("C",1:length(class_labels_levels),sep="") #c("A","B","C","D","E","F","G","H","I","J","K","L","M")
    
    class_col_vec=col_vec

    #if(is.na(col_vec)==TRUE)
    {
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
                                                
                                                if(sample.col.opt=="colorblind"){
                                                    #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                                                    # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                                                    
                                                    if(length(class_labels_levels)<9){
                                                        
                                                        col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7", "grey57")
                                                        
                                                    }else{
                                                        
                                                        col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                                        "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                                                        
                                                        #colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels)))
                                                    }
                                                    
                                                    
                                                }else{
                                                    
                                                    check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                                                    
                                                    if(length(check_brewer)>0){
                                                        
                                                        sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                                                        
                                                        col_vec <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(class_labels_levels))
                                                        
                                                    }else{
                                                        
                                                        if(sample.col.opt=="journal"){
                                                            
                                                            col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                                            "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                                            "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                                            
                                                            "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                                            "#E64B3519","#4DBBD519","#631879E5","grey75")
                                                            
                                                        }else{
                                                            #colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                                            if(length(sample.col.opt)==1){
                                                                col_vec <-rep(sample.col.opt,length(class_labels_levels))
                                                            }else{
                                                                
                                                                colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                                                
                                                            }
                                                        }
                                                    }
                                                    
                                                }
                                        
                                            }
                                        }
                                
                            }
                
            }
    }
    }
                    class_col_vec=rep("red",nrow(classlabels))
                #	print(class_labels_levels)
                
                            ordered_labels={}
                                num_samps_group<-new("list")
                                num_samps_group[[1]]<-0
                                groupwiseindex<-new("list")
                                groupwiseindex[[1]]<-0
                            
                                S<-new("list")
                                for(c in 1:length(class_labels_levels))
                                {
                                    
                                    classlabels_index<-which(classlabels[,2]==class_labels_levels[c])
                                    #ordered_labels<-c(ordered_labels,as.character(classlabels[classlabels_index,2]))
                                    num_samps_group[[c]]<-length(classlabels_index)
                                    groupwiseindex[[c]]<-classlabels_index
                                    
                                    class_col_vec[which(classlabels[,2]==class_labels_levels[c])]<-col_vec[c]
                                    # S[[c]]<-cov(t(data_m[,c(classlabels_index)]))
                                    
                                }
                                
                               
                
                            sampleclass<-{}
                            patientcolors<-{}
                            
                            classlabels<-as.data.frame(classlabels)
                            
                           
                            for(c in 1:length(class_labels_levels)){
                               
                                        num_samps_group_cur=length(which(ordered_labels==class_labels_levels[c]))
                                        
                                        
                                        sampleclass<-c(sampleclass,rep(paste("Class",class_label_alphabets[c],sep=""),num_samps_group_cur))

                                        patientcolors <-c(patientcolors,rep(col_vec[c],num_samps_group[[c]]))
                            }

    if(length(mzvec)>4){
    max_per_row<-3


    par_rows<-ceiling(9/max_per_row)

    }else{
    max_per_row<-length(mzvec)
    par_rows<-1
    }

    file_ind<-0
    boxplots_fname<-paste("xyplots.pdf",sep="")
    
    suppressWarnings(dir.create("Tables"))

    if(dim(data_m)[1]>2){

    t1<-table(classgroup)

    l1<-levels(as.factor(classgroup))


    patientcolors <- rep(col_vec[1:length(t1)], t1)


    res<-get_pca(X=data_m,samplelabels=classgroup,legendlocation=legendlocation,filename=filename,ncomp=10,center=TRUE,scale=TRUE,legendcex=0.5,outloc=getwd(),col_vec=col_vec,sample.col.opt=sample.col.opt,alphacol=0.3,class_levels=NA,pca.cex.val=pca.cex.val,pca.ellipse=pca.ellipse,do_pca_anova=do_pca_anova,paireddesign=paireddesign,classlabelsorig=classlabels_orig)

    pc_pval_vec<-res$pca_pval_vec
    res<-res$result
    
    fname<-paste("Tables/PCAloadings_",filename,"features.txt",sep="")

    loadings_res<-res$rotation
    scores_res<-res$x

    loadings_res<-round(loadings_res,2)
    scores_res<-round(scores_res,2)

    if(dim(loadings_res)[2]>10){
        loadings_res<-loadings_res[,c(1:10)]
        scores_res<-scores_res[,c(1:10)]
    }

    write.table(loadings_res,file=fname,sep="\t")

    fname<-paste("Tables/PCAscores_",filename,"features.txt",sep="")

    scores_res1<-cbind(classlabels_orig,scores_res)


    write.table(scores_res1,file=fname,sep="\t",row.names=FALSE)

    pcnum_limit<-min(5,dim(scores_res)[2])

    get_pooled_sp<-function(n1,n2,S1,S2){a<-(n1-1)*S1;b<-(n2-1)*S2;c<-(n1+n2-2);return((a+b)/c)}
   

     class_levels<-levels(as.factor(classlabels_orig[,2]))
     
     sample.col.opt=lineplot.col.opt
     
     {
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
                             
                             if(sample.col.opt=="colorblind"){
                                 #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                                 # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                                 
                                 if(length(class_labels_levels)<9){
                                     
                                     col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7", "grey57")
                                     
                                 }else{
                                     
                                     col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                     "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                                     
                                     #colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels)))
                                 }
                                 
                                 
                             }else{
                                 
                                 check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                                 
                                 if(length(check_brewer)>0){
                                     
                                     sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                                     
                                     col_vec <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(class_labels_levels))
                                     
                                 }else{
                                     
                                     if(sample.col.opt=="journal"){
                                         
                                         col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                         "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                         "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                         
                                         "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                         "#E64B3519","#4DBBD519","#631879E5","grey75")
                                         
                                     }else{
                                         #colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                         if(length(sample.col.opt)==1){
                                             col_vec <-rep(sample.col.opt,length(class_labels_levels))
                                         }else{
                                             
                                             colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                             
                                         }
                                     }
                                 }
                                 
                             }
                             
                         }
                     }
                     
                 }
                 
             }
         }
     }
     df_matrix<-{}

     if(pairedanalysis==TRUE)
     {
         
            if(dim(classlabels_orig)[2]>2)
            {
                    class_levels<-levels(as.factor(classlabels_orig[,2]))
                    t1<-table(classlabels_orig[,3])
                    
                    for(c in 1:length(class_levels))
                    {
                        
                        
                        class_col_vec[which(classlabels[,2]==class_levels[c])]<-col_vec[c]
                        
                        
                    }
                    #print("Doing two factors")
                   
                   # for(pc in 1:pcnum_limit)
                   lapply(1:pcnum_limit,function(pc)
                   {

                            ylabel_text=paste("PC",pc,"score",sep="")
                            
                            if(pairedanalysis==FALSE){
                                 df<-cbind(res$x[,pc],classgroup,classlabels_orig[,2],classlabels_orig[,3])
                                 
                                 mzname<-paste("PC",pc," scores distribution with 95% confidence interval \n in each group using ",filename," feats vs factors (p=",pc_pval_vec[pc],")",sep="")

                            }else{


                                 df<-cbind(res$x[,pc],classgroup,classlabels_orig[,2],classlabels_orig[,3])
                                 
                                 mzname<-paste("PC",pc," scores distribution with 95% confidence interval \n in each group using ",filename," feats vs time (p=",pc_pval_vec[pc],")",sep="")
                                 
                            }

                            fname<-paste("pc_",pc,"_scoreplot",".tiff",sep="")
                            par_rows=3
                            
                            colnames(df)<-c("y","x","Class","time")
                            df=as.data.frame(df)
                            df$y<-as.numeric(as.character(df$y))

                          #  #save(df,file="df.Rda")

                                        df_fname<-paste("dftable_pc",pc,".txt",sep="")
                                        
                                        if(pc>1){
                                            if(pc==2){
                                                PC2Score=df[,1]
                                                df_matrix<-cbind(df_matrix,PC2Score)
                                            }else{
                                                
                                                
                                                PC3Score=df[,1]
                                                df_matrix<-cbind(df_matrix,PC3Score)
                                                
                                                
                                            }
                                        }else{
                                            PC1Score<-df[,1]
                                            df_matrix<-cbind(df[,-c(1)],PC1Score)
                                            
                                        }
                                        
                                     #   #save(list=ls(),file="pcadebug.Rda")
                                        #    write.table(df,file=df_fname,sep="\t",row.names=FALSE)
                                        #df.summary <- df %>% group_by(x) %>%
                            
                            #      dplyr::summarize(ymin = quantile(y,0.25),
                            #     ymax = quantile(y,0.75),
                            #      Score = median(y),number=length(y),Class=unique(as.character(Class)),Time=unique(as.character(time)))
                        
                    #    df.summary <- aggregate(df$y,
                   #     by = list(df$Class,df$time),
                   #     FUN = function(y) c(ymin = quantile(y,0.25),
                    #    ymax = quantile(y,0.75),
                    #    Score = median(y),number=length(y)))
			
			  df.summary <- aggregate(df$y,
                        by = list(df$Class,df$time),
                        FUN = function(y) c(ysd = sd(y),
                        yse =sd(y)/sqrt(length(y)),
                        Score = mean(y),number=length(y)))
                        
                        df.summary<-do.call(data.frame, df.summary)
                        
                        df.summary<-cbind(df.summary[,c(1,3,4,5,6,1,2)])
                        
                        colnames(df.summary)<-c("x","sd","se","PCscore","number","Class","time")
			
			ymax = df.summary$PCscore + 1.96* df.summary$se
 
			ymin =  df.summary$PCscore - 1.96* df.summary$se
			
			max_yval<-ceiling(max((df.summary$PCscore + (2.5* df.summary$se)),na.rm=TRUE)) #round(max( df.summary$Intensity+(4* df.summary$se),na.rm=TRUE))
			 
			
			min_yval<-floor(min((df.summary$PCscore - (2.5*df.summary$se)),na.rm=TRUE))
					
                        class_levels_time<-levels(as.factor(classlabels_orig[,3]))
                        
                        df.summary$Class<-as.numeric(df.summary$Class)
                        
                                  Class<-{}
                                  for(cnum in 1:length(class_levels)){
                                      
                                        Class<-c(Class,rep(class_levels[cnum],length(t1)))
                                        
                                        df.summary$Class[which(df.summary$Class==cnum)]<-class_levels[cnum]
                                        
                                  }
                                  
                         df.summary$time<-as.numeric(df.summary$time)
                         
                         for(cnum in 1:length(class_levels_time)){
                             
                             
                             df.summary$time[which(df.summary$time==cnum)]<-class_levels_time[cnum]
                             
                         }
                         
                            Class<-unique(Class)
           
           
           
                            #save(df.summary,file="df.summary.Rda")
                            #df.summary$x<- df.summary$time
                            
                            #save(class_levels,file="class_levels.Rda")
                            #save(Class,file="Class.Rda")
                            
                            #save(list=ls(),file="pcadebugfactor.Rda")
                         
                        
                        if(pairedanalysis==TRUE){
                            
                            if(FALSE){
                            plot_res<-ggplot(df.summary, aes(x = x, y = PCscore,color = Class)) + geom_point(size = pca.cex.val) + geom_line(aes(group =Class))  + geom_errorbar(aes(ymin = ymin, ymax = ymax)) + xlab("") + scale_color_manual(values=unique(class_col_vec)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.margin=unit(c(10,5,5,5),"mm"),
                            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                            strip.text = element_text(face="bold"))
                            }
                            
                            colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_levels))
                            
                            
                            if(length(class_levels)>6){
                                
                                lineplot.lty.option=c(lineplot.lty.option,rep("solid",(length(class_levels)-6)))
                            }
                            
                            lapply(1:length(class_levels),function(clevelind){
                             
                                rindvec<-which(df.summary$Class==class_levels[clevelind])
                                
                                lty_val=lineplot.lty.option[clevelind]
                                if(clevelind==1){
                                    barCenters=plot(df.summary$PCscore[rindvec]~as.numeric(factor(df.summary$time[rindvec])),type="o",col=col_vec[clevelind],pch=19,,xaxt="n",lwd=2,ylab=ylabel_text,xlab="TimePoint",main=mzname,cex.main=0.8,
				    ylim=range(pretty(c(min_yval,max_yval))),lty=lineplot.lty.option[clevelind])
				    
                                    axis(1, labels=c(class_levels_time),at=1:length(class_levels_time))
                                    
                                }else{
                                    lines(df.summary$PCscore[rindvec]~as.numeric(factor(df.summary$time[rindvec])),lty=lineplot.lty.option[clevelind],col=col_vec[clevelind],pch=19,type="o",lwd=2)
                                    
                                }
                                
                                if(error.bar==TRUE){
                                    arrows(as.numeric(factor(df.summary$time[rindvec])),df.summary$PCscore[rindvec],as.numeric(factor(df.summary$time[rindvec])),ymin[rindvec],angle=90,lty=1,lwd=1.25,length=0.05,col=col_vec[clevelind],code=3)
                                    arrows(as.numeric(factor(df.summary$time[rindvec])),df.summary$PCscore[rindvec],as.numeric(factor(df.summary$time[rindvec])),ymax[rindvec],angle=90,lty=1,lwd=1.25,length=0.05,col=col_vec[clevelind],code=3)
                                }
                                
                                print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,unique(df.summary$Class), col = col_vec[1:length(t1)],pch = rep(NA,length(col_vec[1:length(t1)])), lty=lineplot.lty.option,pt.cex = 0.6, title = "Class",cex=0.8,lwd=2))
                                
                            })
                           
                           
                        }else{

                            lapply(1:length(class_levels),function(clevelind){
                                
                                rindvec<-which(df.summary$Class==class_levels[clevelind])
                                if(clevelind==1){
                                    barCenters=plot(df.summary$PCscore[rindvec]~as.numeric(factor(df.summary$time[rindvec])),type="o",col=col_vec[clevelind],pch=19,xaxt="n",lwd=2,ylab=ylabel_text,xlab="",lty=0,main=mzname,cex.main=0.8,
				    ylim=range(pretty(c(min_yval,max_yval))))
				    
                                    axis(1, labels=c(class_levels_time),at=1:length(class_levels_time))
                                    
                                }else{
                                    lines(df.summary$PCscore[rindvec]~as.numeric(factor(df.summary$time[rindvec])),col=col_vec[clevelind],pch=19,type="o",lwd=2,lty=0)
                                    
                                }
                                
                                if(error.bar==TRUE){
                                    arrows(as.numeric(factor(df.summary$time[rindvec])),df.summary$PCscore[rindvec],as.numeric(factor(df.summary$time[rindvec])),ymin[rindvec],angle=90,lty=1,lwd=1.25,length=0.05,col=col_vec[clevelind],code=3)
                                    arrows(as.numeric(factor(df.summary$time[rindvec])),df.summary$PCscore[rindvec],as.numeric(factor(df.summary$time[rindvec])),ymax[rindvec],angle=90,lty=1,lwd=1.25,length=0.05,col=col_vec[clevelind],code=3)
                                }
                                
                                print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,unique(df.summary$Class), col = col_vec[1:length(t1)],pch = rep(19,length(col_vec[1:length(t1)])),lty=0,pt.cex = 0.6, title = "Class",cex=0.8,lwd=2))
                                
                            })
                           
                        }
                        
                            #print(plot_res + ggtitle(mzname))
                    
            })

            ##############
            }else{

                            print("Doing one factor")
                            #print(pairedanalysis)
                            class_levels<-levels(as.factor(classgroup))
                            t1<-table(classgroup)
                            # for(pc in 1:pcnum_limit)
                            lapply(1:pcnum_limit,function(pc)
                            {
                                    ylabel_text<-paste("PC",pc,"score",sep="")
                                    df<-cbind(res$x[,pc],classlabels_orig[,2],classlabels_orig[,2])

                                    if(pairedanalysis==FALSE){
                                            mzname<-paste("PC",pc," scores distribution (25th percentile, median, 75th percentile) \n in each group using ",filename," feats vs factors (p=",pc_pval_vec[pc],")",sep="")
                                    }else{

                                            mzname<-paste("PC",pc," scores distribution (25th percentile, median, 75th percentile) \n in each group using ",filename," feats vs time (p=",pc_pval_vec[pc],")",sep="")

                                    }

                                    fname<-paste("pc_",pc,"_scoreplot",".tiff",sep="")
                                    par_rows=3
                                
                                    colnames(df)<-c("y","x","Class")
                                    df=as.data.frame(df)
                                    df$y<-as.numeric(as.character(df$y))



					
					
                                    df_fname<-paste("dftable_pc",pc,".txt",sep="")
                                # write.table(df,file=df_fname,sep="\t",row.names=FALSE)
                                
                                    if(pc>1){
                                        if(pc==2){
                                        PC2Score=df[,1]
                                        df_matrix<-cbind(df_matrix,PC2Score)
                                        }else{
                                            
                                            
                                                PC3Score=df[,1]
                                                df_matrix<-cbind(df_matrix,PC3Score)
                                            
                                            
                                        }
                                    }else{
                                        PC1Score<-df[,1]
                                        df_matrix<-cbind(df[,-c(1)],PC1Score)
                                        
                                    }

                                
                                if(FALSE){
                                    df.summary <- aggregate(df$y,
                                    by = list(df$Class),
                                    FUN = function(y) c(ymin = quantile(y,0.25),
                                    ymax = quantile(y,0.75),
                                    Score = median(y),number=length(y)))
				    }
				      df.summary <- aggregate(df$y,
                        by = list(df$Class),
                        FUN = function(y) c(ysd = sd(y),
                        yse =sd(y)/sqrt(length(y)),
                        Score = mean(y),number=length(y)))
                        
                        df.summary<-do.call(data.frame, df.summary)
		
                                    
			df.summary<-cbind(df.summary[,c(1,2:5,1)])
                        
                       
                        colnames(df.summary)<-c("x","sd","se","PCscore","number","Class")
			
			ymax = df.summary$PCscore + 1.96* df.summary$se
 
			ymin =  df.summary$PCscore - 1.96* df.summary$se
			
			max_yval<-ceiling(max((df.summary$PCscore + (2.5* df.summary$se)),na.rm=TRUE)) #round(max( df.summary$Intensity+(4* df.summary$se),na.rm=TRUE))
			 
			
			min_yval<-floor(min((df.summary$PCscore - (2.5*df.summary$se)),na.rm=TRUE))
					
				    
                                    
                                    
                                          Class<-{}
                                          for(cnum in 1:length(class_levels)){
                                          
                                                Class<-c(Class,rep(class_levels[cnum],length(t1)))
                                            
                                        #        if(pairedanalysis==FALSE){
                                         #       df.summary$Class[which(df.summary$Class==cnum)]<-class_levels[cnum]
                                         #       }else{
                                           #         df.summary$Class[which(df.summary$Class==cnum)]<-class_levels[1]
                                            #    }
                                        
                                          }
                                    
                                   # #save(df.summary,file="df.summary.Rda")
                                    
                                    time.hour<- c(unique(as.character(classlabels_orig[,2])),unique(as.character(classlabels_orig[,2])))

                                    ##save(list=ls(),file="pcadebugonefc.Rda")
                            
                            if(FALSE){
                                    if(pairedanalysis==FALSE){
                                        plot_res<-ggplot(df.summary, aes(x = Class, y = PCscore,color = Class)) + geom_point(size = pca.cex.val)  + geom_errorbar(aes(ymin = ymin, ymax = ymax)) + xlab("") + scale_color_manual(values=unique(class_col_vec)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.margin=unit(c(10,5,5,5),"mm"),
                                        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                        strip.text = element_text(face="bold"))
                                    }else{
                                        
                                        plot_res<-ggplot(df.summary, aes(x = Class, y = PCscore,color = Class)) +  geom_point(size = pca.cex.val) + geom_line(aes(group =Class)) + xlab("") + geom_errorbar(aes(ymin = ymin, ymax = ymax)) + scale_color_manual(values=unique(class_col_vec)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.margin=unit(c(10,5,5,5),"mm"),
                                        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                        strip.text = element_text(face="bold"))

                                    }
                            }
                            
                            
                                    
                                    if(pairedanalysis==TRUE){
                                     #   ylim=c(0,max_yval),
                                        clevelind=1
                                        barCenters=plot(df.summary$PCscore~as.numeric(df.summary$Class),type="o",col=col_vec[1],pch=19,xaxt="n",lwd=2,ylab=ylabel_text,cex.main=0.8,main=mzname,xlab="TimePoint",ylim=range(pretty(c(min_yval,max_yval))))
                                        axis(1, labels=c(class_levels),at=1:length(class_levels))
                                        
                                        
                                        
                                        if(error.bar==TRUE){
                                            arrows(df.summary$Class,df.summary$PCscore,df.summary$Class,ymin,angle=90,lty=1,lwd=1.25,length=0.05,col=col_vec[1],code=3)
                                            arrows(df.summary$Class,df.summary$PCscore,df.summary$Class,ymax,angle=90,lty=1,lwd=1.25,length=0.05,col=col_vec[clevelind],code=3)
                                        }
                                        
                                        #print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,unique(df.summary$Class), col = col_vec[1:length(t1)],pch = rep(NA,length(col_vec[1:length(t1)])), lty=rep(1,length(t1)),pt.cex = 0.6, title = "Class",cex=0.8,lwd=2))
                                        
                                        
                                    }
                                    #print(plot_res + ggtitle(mzname))
                            
                            })
                }
		
	}
    }

    df_fname<-paste("PC_score_distribution_matrix_",filename,"features.txt",sep="")

    return(res)
}




get_lineplots<-function(X=NA,Y=NA,feature_table_file=NA,parentoutput_dir=NA,class_labels_file=NA,alphacol=0.3,col_vec=NA,pairedanalysis=FALSE,pca.cex.val=4,legendlocation="topright",pca.ellipse=TRUE,ellipse.conf.level=0.95,filename="all",newdevice=FALSE,lineplot.col.opt=c("grey57"),ylabel="Intensity",error.bar=TRUE,cex.val=0.7,lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"))
{
    cex.plots=cex.val
    if(is.na(parentoutput_dir)==TRUE){
        
        parentoutput_dir=getwd()
    }
    if(is.na(X)==TRUE){
        data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE)
        
    }else{
        
        data_matrix<-X
    }
    if(is.na(Y)==TRUE){
        classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
        
    }else{
        
        classlabels<-Y
    }
    
    sample.col.opt=lineplot.col.opt
    
    par(mfrow=c(1,1),family="sans",cex=cex.val)
    
    #print("here")
    dir.create(parentoutput_dir)
    setwd(parentoutput_dir)
    
    data_m<-data_matrix[,-c(1:2)]
    
    data_m<-as.matrix(data_m)
    
    mzvec<-data_matrix[,1]
    rtvec<-data_matrix[,2]
    
    rnames<-paste(mzvec,rtvec,sep="_")
    
    rownames(data_m)<-as.character(rnames)
    
    
    classlabels_orig<-classlabels
    
    
   
   
   if(dim(classlabels)[2]>2){
       
       #classlabels_orig[,2]<-as.factor(paste("A",as.character(classlabels_orig[,2]),sep=""))
       #classlabels_orig[,3]<-as.factor(paste("B",as.character(classlabels_orig[,3]),sep=""))
       # print(head(classlabels_orig))
       
       classgroup<-paste(classlabels_orig[,2],":",classlabels_orig[,3],sep="") #classlabels_orig[,2]:classlabels_orig[,3]
       do_pca_anova=FALSE
       
       col_class_levels<-levels(as.factor(classlabels_orig[,2]))
   }else{
       
       classgroup<-classlabels_orig[,2]
       
       do_pca_anova=TRUE
   }
   
   
   
   class_labels_levels<-levels(as.factor(classgroup))
   ordered_labels<-classgroup
   
   class_label_alphabets<-paste("C",1:length(class_labels_levels),sep="") #c("A","B","C","D","E","F","G","H","I","J","K","L","M")
   
    if(newdevice==TRUE){
        
        fname<-paste("timeseriesplots",filename,".pdf",sep="")
        pdf(fname)
    }
    
    if(is.na(sample.col.opt)==FALSE)
    {
        if(sample.col.opt=="default"){
            
            col_vec<-c("#CC0000","#AAC000","blue","mediumpurple4","mediumpurple1","blueviolet","cornflowerblue","cyan4","skyblue",
            "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
            "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
            "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
            
        }else{
            if(sample.col.opt=="topo"){
                #col_vec<-topo.colors(256) #length(class_labels_levels))
                
                #col_vec<-col_vec[seq(1,length(col_vec),)]
                
                col_vec <- topo.colors(length(col_class_levels), alpha=alphacol)
            }else{
                if(sample.col.opt=="heat"){
                    #col_vec<-heat.colors(256) #length(class_labels_levels))
                    
                    col_vec <- heat.colors(length(col_class_levels), alpha=alphacol)
                }else{
                    if(sample.col.opt=="rainbow"){
                        #col_vec<-heat.colors(256) #length(class_labels_levels))
                        col_vec<-rainbow(length(col_class_levels), start = 0, end = alphacol)
                        
                        #col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
                    }else{
                        
                        if(sample.col.opt=="terrain"){
                            #col_vec<-heat.colors(256) #length(class_labels_levels))
                            #col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
                            
                            col_vec <- cm.colors(length(col_class_levels), alpha=alphacol)
                        }else{
                            if(is.na(sample.col.opt)==TRUE){
                                col_vec<-c("black")
                            }else{
                                
                                if(sample.col.opt=="colorblind"){
                                    #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                                    # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                                    
                                    if(length(col_class_levels)<9){
                                        
                                        col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7", "grey57")
                                        
                                    }else{
                                        
                                        #col_vec<-colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels))
                                        col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                        "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                                        
                                    }
                                    
                                    
                                }else{
                                    
                                    check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                                    
                                    if(length(check_brewer)>0){
                                        
                                        sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                                        
                                        col_vec <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(col_class_levels))
                                        
                                    }else{
                                        
                                        #colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                        if(sample.col.opt=="journal"){
                                            
                                            col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                            "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                            "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                            
                                            "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                            "#E64B3519","#4DBBD519","#631879E5","grey75")
                                            
                                        }else{
                                            #colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                            
                                            if(length(sample.col.opt)==1){
                                                col_vec <-rep(sample.col.opt,length(col_class_levels))
                                            }else{
                                                
                                               colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(col_class_levels))
                                               
                                               save(colfunc,file="colfunc.Rda")
                                                
                                            }
                                        }
                                        
                                    }
                                    
                                }
                            }
                        }
                        
                        
                    }
                    
                }
                
            }
        }
    }else{
        
        col_vec<-c("black")
    
    }
    
    #	print(class_labels_levels)
    
    ordered_labels={}
    num_samps_group<-new("list")
    num_samps_group[[1]]<-0
    groupwiseindex<-new("list")
    groupwiseindex[[1]]<-0
    
    S<-new("list")
    for(c in 1:length(class_labels_levels))
    {
								
                                classlabels_index<-which(classlabels[,2]==class_labels_levels[c])
                               
                                num_samps_group[[c]]<-length(classlabels_index)
                                groupwiseindex[[c]]<-classlabels_index
                                
                                
                               
                                
    }
    
    
    
    sampleclass<-{}
    patientcolors<-{}
    
    classlabels<-as.data.frame(classlabels)
    
    
    for(c in 1:length(class_labels_levels)){
        
        num_samps_group_cur=length(which(ordered_labels==class_labels_levels[c]))
        
        
        sampleclass<-c(sampleclass,rep(paste("Class",class_label_alphabets[c],sep=""),num_samps_group_cur))
        
        patientcolors <-c(patientcolors,rep(col_vec[c],num_samps_group[[c]]))
    }
    

    
    if(length(mzvec)>4){
        max_per_row<-3
        
        
        par_rows<-ceiling(9/max_per_row)
        
    }else{
        max_per_row<-length(mzvec)
        par_rows<-1
    }
    
    # The palette with black:
    #col_vec <- c("#E69F00","#000000","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    # To use for fills, add
    # scale_fill_manual(values=col_vec)
    
    # To use for line and point colors, add
    #scale_colour_manual(values=col_vec)
    
    file_ind<-0
    boxplots_fname<-paste("xyplots.pdf",sep="")
			
                if(dim(data_m)[1]>2){
                
                    
                    t1<-table(classgroup)
                    l1<-levels(classgroup)
                    
                    l1<-levels(as.factor(classgroup))
                    
                    
                    patientcolors <- rep(col_vec[1:length(t1)], t1)
                    
                    
                    
                    
                    
                    get_pooled_sp<-function(n1,n2,S1,S2){a<-(n1-1)*S1;b<-(n2-1)*S2;c<-(n1+n2-2);return((a+b)/c)}
             
                    
                    class_levels<-levels(as.factor(classlabels_orig[,2]))
                    
               
                    df_matrix<-{}
                    
                    
                    if(dim(classlabels_orig)[2]>2){
                        class_levels<-levels(as.factor(classlabels_orig[,2]))
                        t1<-table(classlabels_orig[,3])
                        
                        
                        
                         
                        myData<-cbind(classgroup,classlabels_orig[,2],classlabels_orig[,3],t(data_m))
                 
                        myData<-as.data.frame(myData)
                        
                        myData_sum <- do.call(data.frame,aggregate(list(myData[,-c(1:3)]),
                        by = list(myData[,2],myData[,3]),
                        FUN = function(x){
                         
                        x<-as.numeric(as.character(x))
                        c(mean = mean(x), sd = sd(x),n = length(x),se=sd(x)/sqrt(length(x)))
                            
                        }))
                       label_inc_list<-seq(3,dim(myData_sum)[2],4)
                       
                       #save(myData_sum,file="myDatasum.Rda")
                       
                       
                       #save(label_inc_list,file="label_inc_list.Rda")
                       #save(rnames,file="rnames.Rda")
                       #save(data_m,file="data_m.Rda")
                        lapply(seq(3,dim(myData_sum)[2],4),function(pc)
                        {
                         
                         df.summary<-myData_sum[,c(pc:(pc+3),1,2)] #cbind(xvec,classgroup,classlabels_orig[,2],classlabels_orig[,3])
                         
                         get_label_ind<-which(label_inc_list==pc)
                         if(pairedanalysis==FALSE){
                             
                             
                             if(error.bar==FALSE){
                                 mzname<-paste(rnames[get_label_ind]," distribution \n in each group using ",filename," feats vs factors",sep="")
                             }else{
                                 mzname<-paste(rnames[get_label_ind]," distribution with 95% confidence interval\n in each group using ",filename," feats vs factors",sep="")
                             }
                         }else{
                             
                             if(error.bar==FALSE){
                                 mzname<-paste(rnames[get_label_ind]," distribution \n in each group using ",filename," feats vs time",sep="")
                             }else{
                                 
                                 mzname<-paste(rnames[get_label_ind]," distribution 95% confidence interval\n in each group using ",filename," feats vs time",sep="")
                             }
                             
                         }
                         
                         colnames(df.summary)<-c("Intensity","sd","number","se","Class","time")
                         
                         df.summary$Class<-as.numeric(df.summary$Class)
                         df.summary$time<-as.numeric(df.summary$time)
                         
                         ymax = df.summary$Intensity + 1.96*df.summary$se
                         
                         ymin = df.summary$Intensity - 1.96*df.summary$se
                         
                         max_yval<-ceiling(max(df.summary$Intensity+(2.5* df.summary$se),na.rm=TRUE)) #max(ymax,na.rm=TRUE)
                         
                         min_yval<-max(0,floor(min(df.summary$Intensity-(2.5* df.summary$se),na.rm=TRUE)))
                         
                         class_levels_time<-levels(as.factor(classlabels_orig[,3]))
                            
                            Class<-{}
                            for(cnum in 1:length(class_levels)){
                                
                                Class<-c(Class,rep(class_levels[cnum],length(t1)))
                        
                                df.summary$Class[which(df.summary$Class==cnum)]<-class_levels[cnum]

                            }
                            
                        
                            Class<-unique(Class)
                            
                            t1<-table(df.summary$Class)
                            
                            
                            
                            #print(df.summary)
                            
                            df.summary$x<-df.summary$time # time.hour
                            save(df.summary,file="df.summary.Rda")
                            save(col_vec,file="col_vec.Rda")
                            save(class_levels,file="class_levels.Rda")
                            save(list=ls(),file="debuglinep.Rda")
                            
                            col_vec<-colfunc(length(t1))
                            
                            if(pairedanalysis==TRUE){
                                
                                if(length(class_levels)>6){
                                    
                                    lineplot.lty.option=c(lineplot.lty.option,rep("solid",(length(class_levels)-6)))
                                }
                                
                                lapply(1:length(class_levels),function(clevelind){
                                        rindvec<-which(df.summary$Class==class_levels[clevelind])
                                        
                                        lty_val=lineplot.lty.option[clevelind]
                                        
                                        if(clevelind==1){
                                    
                                           w <- 0.1
                                           par(omd=c(0, 1-w, 0, 1))
                                           barCenters=plot(df.summary$Intensity[rindvec]~as.numeric(factor(df.summary$time[rindvec])),type="o",col=col_vec[clevelind],pch=19,ylim=c(min_yval,max_yval),xaxt="n",lwd=2,ylab=ylabel,xlab="TimePoint",main=mzname,cex.main=0.8,lty=lty_val)
                                                axis(1, labels=c(class_levels_time),at=1:length(class_levels_time))
                                  
                                        }else{
                                                lines(df.summary$Intensity[rindvec]~as.numeric(factor(df.summary$time[rindvec])),lty=lty_val,col=col_vec[clevelind],pch=19,type="o",lwd=2)
                                            
                                        }
                                        
                                       
                                        
                                        if(error.bar==TRUE){
                                        arrows(as.numeric(factor(df.summary$time[rindvec])),df.summary$Intensity[rindvec],as.numeric(factor(df.summary$time[rindvec])),ymin[rindvec],angle=90,lty=1,lwd=1.25,length=0.05,col=col_vec[clevelind],code=3)
                                        arrows(as.numeric(factor(df.summary$time[rindvec])),df.summary$Intensity[rindvec],as.numeric(factor(df.summary$time[rindvec])),ymax[rindvec],angle=90,lty=1,lwd=1.25,length=0.05,col=col_vec[clevelind],code=3)
                                        }
                                        
                                })
                                print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,unique(df.summary$Class), col = col_vec[1:length(t1)],pch = rep(NA,length(col_vec[1:length(t1)])), lty=lineplot.lty.option,pt.cex = 0.6, title = "Class",cex=cex.plots,lwd=2))
                                
                                
                            }else{
                                
                                    lapply(1:length(class_levels),function(clevelind){
                                        
                                        rindvec<-which(df.summary$Class==class_levels[clevelind])
                                        if(clevelind==1){
                                            w <- 0.1
                                            par(omd=c(0, 1-w, 0, 1))
                                           barCenters=plot(df.summary$Intensity[rindvec]~as.numeric(factor(df.summary$time[rindvec])),type="o",col=col_vec[clevelind],pch=19,ylim=c(min_yval,max_yval),xaxt="n",lwd=2,ylab=ylabel,xlab="Factor2",lty=0)
                                            axis(1, labels=c(class_levels_time),at=1:length(class_levels_time))
                                            
                                        }else{
                                            lines(df.summary$Intensity[rindvec]~as.numeric(factor(df.summary$time[rindvec])),col=col_vec[clevelind],pch=19,type="o",lwd=2,lty=0)
                                            
                                        }
                                       
                                       if(error.bar==TRUE){
                                        arrows(as.numeric(factor(df.summary$time[rindvec])),df.summary$Intensity[rindvec],as.numeric(factor(df.summary$time[rindvec])),ymin[rindvec],angle=90,lty=1,lwd=1.25,length=0.05,col=col_vec[clevelind],code=3)
                                        arrows(as.numeric(factor(df.summary$time[rindvec])),df.summary$Intensity[rindvec],as.numeric(factor(df.summary$time[rindvec])),ymax[rindvec],angle=90,lty=1,lwd=1.25,length=0.05,col=col_vec[clevelind],code=3)
                                       }
                                       
                                        print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,unique(df.summary$Class), col = col_vec[1:length(t1)],pch = rep(19,length(col_vec[1:length(t1)])),lty=0,pt.cex = 0.6, title = "Factor1",cex=cex.plots,lwd=2))
                                        
                                    })
                                
                            }
                            
                            # print(plot_res + ggtitle(mzname))
                            
                        })
                    }else{
                        
                      
                      #one-way anova or one-way with repeated measures
                        class_levels<-levels(as.factor(classgroup))
                        t1<-table(classgroup)
                        
                        
                     
                     myData<-cbind(classgroup,classlabels_orig[,2],t(data_m))
                     
                     myData<-as.data.frame(myData)
                     
                     myData_sum <- do.call(data.frame,aggregate(list(myData[,-c(1:2)]),
                     by = list(myData[,2]),
                     FUN = function(x){
                         
                         x<-as.numeric(as.character(x))
                         c(mean = mean(x), sd = sd(x),n = length(x),se=sd(x)/sqrt(length(x)))
                         
                     }))
                     
                     #save(myData_sum,file="myDatasum.Rda")
                     label_inc_list<-seq(2,dim(myData_sum)[2],4)
                     
                      #save(label_inc_list,file="label_inc_list.Rda")
                      #save(rnames,file="rnames.Rda")
                      
                     lapply(seq(2,dim(myData_sum)[2],4),function(pc)
                     {
                         
                         df.summary<-myData_sum[,c(pc:(pc+3),1)] #cbind(xvec,classgroup,classlabels_orig[,2],classlabels_orig[,3])
                         
                         get_label_ind<-which(label_inc_list==pc)
                         if(pairedanalysis==FALSE){
                             
                             
                             if(error.bar==FALSE){
                                 mzname<-paste(rnames[get_label_ind]," distribution \n in each group using ",filename," feats vs factors",sep="")
                             }else{
                                 mzname<-paste(rnames[get_label_ind]," distribution with 95% confidence interval\n in each group using ",filename," feats vs factors",sep="")
                             }
                         }else{
                             
                             if(error.bar==FALSE){
                                 mzname<-paste(rnames[get_label_ind]," distribution \n in each group using ",filename," feats vs time",sep="")
                             }else{
                                 
                                 mzname<-paste(rnames[get_label_ind]," distribution 95% confidence interval\n in each group using ",filename," feats vs time",sep="")
                             }
                             
                         }
                         
                         colnames(df.summary)<-c("Intensity","sd","number","se","Class")
                         
                         df.summary$Class<-as.numeric(df.summary$Class)
                        
                         
                         ymax = df.summary$Intensity + 1.96*df.summary$se
                         
                         ymin = df.summary$Intensity - 1.96*df.summary$se
                         
                        # max_yval<-ceiling(max(ymax,na.rm=TRUE)) #round(max(myData$Intensity+(4*myData$se),na.rm=TRUE))
                         
                       #  min_yval<-max(0,floor(min(ymin,na.rm=TRUE)))
                         
                         
                         max_yval<-ceiling(max(df.summary$Intensity+(2.5* df.summary$se),na.rm=TRUE)) #max(ymax,na.rm=TRUE)
                         
                         min_yval<-max(0,floor(min(df.summary$Intensity-(2.5* df.summary$se),na.rm=TRUE)))
                         Class<-{}
                         for(cnum in 1:length(class_levels)){
                             
                             Class<-c(Class,rep(class_levels[cnum],length(t1)))
                             
                           #  df.summary$Class[which(df.summary$Class==cnum)]<-class_levels[cnum]
                             
                         }
                         
                       
                         Class<-unique(Class)
                         
                         t1<-table(df.summary$Class)
                         
                         
                         #print(df.summary)
                         ##save(df.summary,file="df.summarydotp1.Rda")
                        
                         
                        # #save(list=ls(),file="debuglinep1.Rda")
                         
                         
                         
                         if(pairedanalysis==TRUE){
                             
                             clevelind=1
                             barCenters=plot(df.summary$Intensity~as.numeric(df.summary$Class),type="o",col=col_vec[1],pch=19,ylim=c(min_yval,max_yval),xaxt="n",lwd=2,ylab=ylabel,xlab="TimePoint",main=mzname,main.cex=cex.plots)
                                     axis(1, labels=c(class_levels),at=1:length(class_levels))
                                     
                                     
                                 
                                 if(error.bar==TRUE){
                                     arrows(df.summary$Class,df.summary$Intensity,df.summary$Class,ymin,angle=90,lty=1,lwd=1.25,length=0.05,col=col_vec[1],code=3)
                                     arrows(df.summary$Class,df.summary$Intensity,df.summary$Class,ymax,angle=90,lty=1,lwd=1.25,length=0.05,col=col_vec[clevelind],code=3)
                                 }
                                 
                                 #print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,unique(df.summary$Class), col = col_vec[1:length(t1)],pch = rep(NA,length(col_vec[1:length(t1)])), lty=rep(1,length(t1)),pt.cex = 0.6, title = "Class",cex=0.8,lwd=2))
                           
                     
                         }else{
                             
                             clevelind=1
                             barCenters=plot(df.summary$Intensity~as.numeric(df.summary$Class),type="o",col=col_vec[1],pch=19,ylim=c(min_yval,max_yval),xaxt="n",lwd=2,ylab=ylabel,xlab="Class",lty=0)
                             axis(1, labels=c(class_levels),at=1:length(class_levels))
                             
                             
                             
                             if(error.bar==TRUE){
                                 arrows(df.summary$Class,df.summary$Intensity,df.summary$Class,ymin,angle=90,lwd=1.25,length=0.05,col=col_vec[1],code=3)
                                 arrows(df.summary$Class,df.summary$Intensity,df.summary$Class,ymax,angle=90,lwd=1.25,length=0.05,col=col_vec[1],code=3)
                             }
                             
                         }
                })
                    }
            }
                #df_fname<-paste("PC_score_distribution_matrix_",filename,"features.txt",sep="")
                #write.table(df_matrix,file=df_fname,sep="\t",row.names=TRUE)
                if(newdevice==TRUE){
                    
                    try(dev.off(),silent=TRUE)
                }
            


}


get_pca<-function(X,samplelabels,legendlocation="topright",filename=NA,ncomp=5,center=TRUE,scale=TRUE,legendcex=0.5,outloc=getwd(),col_vec=NA,sample.col.opt="default",alphacol=0.3,class_levels=NA,pca.cex.val=3,pca.ellipse=TRUE,ellipse.conf.level=0.95,samplenames=FALSE,do_pca_anova=FALSE,paireddesign=NA,classlabelsorig=NA){
    
    X<-as.matrix(t(X))
    
    par(mfrow=c(1,1),family="sans",cex=0.9)
    pch.val<-19
 
    samplelabels<-as.data.frame(samplelabels)
 
    samplelabels<-paste("",as.factor(samplelabels[,1]),sep="")
    

    l2<-levels(as.factor(samplelabels))
    col_all=topo.colors(256)
    
    t1<-table(samplelabels)
    if(is.na(class_levels)==TRUE){
        
        l1<-levels(as.factor(samplelabels))
    }else{
        l1<-class_levels
        
        
    }
    
    class_labels_levels<-l1
    
    ncomp=min(c(10,dim(X)[1],dim(X)[2]))
    
    #   p1<-pcaMethods::pca(t(X),method="svd",center=TRUE,scale="uv",cv="q2",nPcs=10)
    if(is.na(paireddesign)==TRUE){
    metabpcaresultlog2allmetabs5pcs<-mixOmics::pca(X,ncomp=ncomp,center=TRUE,scale=TRUE)
    }else{
        metabpcaresultlog2allmetabs5pcs<-mixOmics::pca(X,ncomp=ncomp,center=TRUE,scale=TRUE) #,multilevel=paireddesign)
          
    }
    result<-metabpcaresultlog2allmetabs5pcs
    

    ##save(result,file="pcares.Rda")
    
    s1<-summary(result)
    r1<-s1$importance[2,]
    r1<-round(r1,2)*100
    
      barplot(r1,beside=TRUE,main="% variation explained by each component",ylab="% variation",col=c("#0072B2"),cex.main=0.8,ylim=range(pretty(c(0,max(r1)))))
    
    
    if(is.na(col_vec)==TRUE){
        col_vec<-c("mediumpurple4","mediumpurple1","blueviolet","darkblue","blue","cornflowerblue","cyan4","skyblue",
        "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
        "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
        "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
        
    }
    
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
                        
                        if(sample.col.opt=="colorblind"){
                            #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                            # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                            
                            if(length(class_labels_levels)<9){
                                
                                col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7", "grey57")
                                
                            }else{
                                
                                #col_vec<-colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels))
                                col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                                
                            }
                            
                            
                        }else{
                            
                            check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                            
                            if(length(check_brewer)>0){
                                
                                 sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                                col_vec <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(class_labels_levels))
                                
                            }else{
                                
                                if(sample.col.opt=="journal"){
                                    
                                    col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                    "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                    "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                    
                                    "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                    "#E64B3519","#4DBBD519","#631879E5","grey75")
                                    
                                }else{
                                    #colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                    if(length(sample.col.opt)==1){
                                        col_vec <-rep(sample.col.opt,length(class_labels_levels))
                                    }else{
                                        
                                        colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                        
                                    }
                                }
                                
                            }
                            
                        }
                        
                    }
                    
                    
                }
                
            }
            
        }
    }
    #col_vec<-col_vec[sample(1:length(col_vec),length(col_vec))]
    
    l1<-gsub(x=l1,pattern="Class",replacement="")
    
    dir.create(outloc,showWarnings=FALSE)
    setwd(outloc)
    # print(paste("Generating PCA plots",sep=""))
    
    fname<-paste("PCA_eval",filename,".tiff",sep="")
    ## 1) raw data
    #tiff(fname,res=300, width=2000,height=2000)
    col <- rep(col_vec[1:length(t1)], t1)
	
    #col<-rep(col_all[1:length(l1)],t1)
    ## Choose different size of points
    cex <- rep(pca.cex.val, dim(X)[1])
    ## Choose the form of the points (square, circle, triangle and diamond-shaped
    
    #pch_vec<-seq(1,50) #c(3,5,7,9,12,13,2,17,21) #seq(1,50) #
    pch <- rep(pch.val,dim(X)[1])
    pch_vec <- rep(pch.val,dim(X)[1])
    
    
    cex <- pca.cex.val #rep(pca.cex.val, dim(X)[1])
    col_per_group<-{}
    pch_per_group<-{}
    for(p1 in 1:length(l2)){
        col[which(samplelabels==l2[p1])]=col_vec[p1]
        pch[which(samplelabels==l2[p1])]=pch_vec[p1]
        
        col_per_group<-c(col_per_group,col_vec[p1])
        pch_per_group<-c(pch_per_group,pch_vec[p1])
    }
   
   # print(table(pch))
   
   main_text=paste("Pairwise PC score plots using ",filename," features after preprocessing",sep="")
   
	   legendcex<-0.7 #0.5*pca.cex.val

       #if(pca.ellipse=="car")
   {
      
      
      pc_pval_single<-{}
      print("here")
      print(dim(classlabelsorig))
      print(paireddesign)
      
      do_pca_anova=TRUE
      
      #save(classlabelsorig,file="classlabelsorigpc.Rda")
      #save(result,file="result.Rda")
      #save(samplelabels,file="samplelabels.Rda")
      #if(do_pca_anova==TRUE)
      {
          scores_res<-result$x
          
          # #save(scores_res,file="scores_res.Rda")
          ##save(samplelabels,file="samplelabels.Rda")
          ##save(classlabelsorig,file="classlabelsorig.Rda")
          ##save(paireddesign,file="paireddesign.Rda")
          pc1_pval<-anova(lm(cbind(scores_res[,1],scores_res[,2])~samplelabels))
          
          pc1_pval<-pc1_pval[[6]][2]
          pc2_pval<-anova(lm(cbind(scores_res[,1],scores_res[,3])~samplelabels))
          
          pc2_pval<-pc2_pval[[6]][2]
          
          pc3_pval<-anova(lm(cbind(scores_res[,2],scores_res[,3])~samplelabels))
          pc3_pval<-pc3_pval[[6]][2]
          
          if(dim(classlabelsorig)[2]==2){
              dtemp<-cbind(classlabelsorig,scores_res)
              dtemp<-as.data.frame(dtemp)
              
              #save(dtemp,file="pcdtemp.Rda")
              if(is.na(paireddesign)==TRUE){
                  
                 testname="one-way ANOVA"
                  pc_pval_single<-lapply(1:5,function(pcn1){
                      pc1_only<-anova(lm(scores_res[,pcn1]~samplelabels))
                      pc1_only<-pc1_only[[5]][1]
                      return(pc1_only)
                      
                  })
                  
                  pc_pval_single<-unlist(pc_pval_single)
              }else{
                  
                   testname="one-way ANOVA with repeated measures"
                     #one-way ANOVA repeat
                      pc_pval_single<-lapply(3:ncol(dtemp),function(pcn1){
                          dataA<-cbind(dtemp[,pcn1],dtemp[,c(2)])
                          
                          colnames(dataA)<-c("Response","Factor1")
                          
                          pc1_res<-diffexplmonewayanovarepeat(dataA=dataA,subject_inf=paireddesign)
                          
                          return(pc1_res$mainpvalues)
                          
                      })
                      
                      pc_pval_single<-unlist(pc_pval_single)
                   
              }
          
          }else{
              
              dtemp<-cbind(classlabelsorig,scores_res)
              dtemp<-as.data.frame(dtemp)
              
              #save(dtemp,file="pcdtemp2.Rda")
              
              if(dim(classlabelsorig)[2]>=2){
                  
                  if(is.na(paireddesign)==FALSE){
                      
                      testname="two-way ANOVA with repeated measures in one factor"
                          #two-way ANOVA repeat
                          pc_pval_single<-lapply(4:ncol(dtemp),function(pcn1){
                              dataA<-cbind(dtemp[,pcn1],dtemp[,c(2:3)])
                              
                              colnames(dataA)<-c("Response","Factor1","Factor2")
                              
                              dataA$Response<-dtemp[,pcn1]
                              
                              #save(dataA,file="pcdataA.rda")
                              pc1_res<-diffexplmtwowayanovarepeat(dataA=dataA,subject_inf=paireddesign)

                               return(pc1_res$mainpvalues)
                              
                          })
                          pc_pval_single<-unlist(pc_pval_single)
                  }else{
                      
                      testname="two-way ANOVA"
                      #two-way ANOVA
                      pc_pval_single<-lapply(4:ncol(dtemp),function(pcn1){
                          dataA<-cbind(dtemp[,pcn1],dtemp[,c(2:3)])
                          
                          colnames(dataA)<-c("Response","Factor1","Factor2")
                          
                          dataA$Response<-dtemp[,pcn1]
                          
                          #save(dataA,file="pcdataA.rda")
                          pc1_res<-diffexplmtwowayanova(dataA=dataA)
                          
                          return(pc1_res$mainpvalues)
                          
                      })
                      pc_pval_single<-unlist(pc_pval_single)
                      
                  }
                  
              }
              
          }
          
          pc_pval_vec<-c(pc1_pval,pc2_pval,pc3_pval)
      }

      
      pc_pval_single<-round(pc_pval_single)
     
     
     if(do_pca_anova==TRUE){
         main_text=paste("Pairwise PC score plots using ",filename," features after preprocessing\np-value for overall differences between groups using PC1 and PC3 in a multivariate\n one-way ANOVA model=",round(pc1_pval,3),sep="")
     }
    
      
      # plot(c(1,1),plot=FALSE)
      #l <- legend(0, 0, bty='n', l1,plot=FALSE, pch = pch_per_group, pt.cex = 0.6)
      # calculate right margin width in ndc
      w <- 0.1 #grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
      par(omd=c(0, 1-w, 0, 1))
     
     iqr_xlim=4*sd(result$variates$X[,1],na.rm=TRUE) #4*(summary(result$variates$X[,1])[5]-summary(result$variates$X[,1])[3])
     iqr_ylim=4*sd(result$variates$X[,2],na.rm=TRUE) #4*(summary(result$variates$X[,2])[5]-summary(result$variates$X[,2])[3])
     
     
     
      print(dataEllipse(x=result$x[,1], y=result$x[,2],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=c(ellipse.conf.level),col=col_per_group,pch=pch,xlab=paste("PC1 (",r1[1],"% variation)",sep=""),ylab=paste("PC2 (",r1[2],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE,cex.main=0.8,center.pch=FALSE,ylim=range(pretty(c(floor(min(result$variates$X[,2])-iqr_ylim),ceiling(max(result$variates$X[,2])+iqr_ylim)))),xlim=range(pretty(c(floor(min(result$variates$X[,1])-iqr_xlim),ceiling(max(result$variates$X[,1])+iqr_xlim))))))
      
      
      # print(legend(legendlocation, l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.7, title = "Class", cex=legendcex))
      
      # print(dataEllipse(x=result$x[,1], y=result$x[,2],groups=as.factor(samplelabels),grid=FALSE,lwd=2,levels=c(ellipse.conf.level),col=col_per_group,pch=pch,xlab=paste("PC1 (",r1[1],"% variation)",sep=""),ylab=paste("PC2 (",r1[2],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE,cex.main=0.7,center.pch=FALSE,ylim=c(min(result$x[,2])-10,max(result$x[,2])+10),xlim=c(min(result$x[,1])-10,max(result$x[,1])+10)))
      
       print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
       
       
       #      plot(c(1,1),plot=FALSE)
       #l <- legend(0, 0, bty='n', l1,plot=FALSE, pch = pch_per_group, pt.cex = 0.6)
       # calculate right margin width in ndc
       #w <- grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
       par(omd=c(0, 1-w, 0, 1))
       
       print(dataEllipse(x=result$x[,1], y=result$x[,2],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=NULL,col=col_per_group,pch=pch,xlab=paste("PC1 (",r1[1],"% variation)",sep=""),ylab=paste("PC2 (",r1[2],"% variation)",sep=""),group.labels="",main=main_text,fill=FALSE,add=FALSE,cex.main=0.8,ylim=range(pretty(c(floor(min(result$variates$X[,2])-iqr_ylim),ceiling(max(result$variates$X[,2])+iqr_ylim)))),xlim=range(pretty(c(floor(min(result$variates$X[,1])-iqr_xlim),ceiling(max(result$variates$X[,1])+iqr_xlim))))))
       
       #print(legend(legendlocation, l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.7, title = "Class", cex=legendcex))
       
       print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
       
       
       if(ncomp>2){
           if(do_pca_anova==TRUE){
               main_text=paste("Pairwise PC score plots using ",filename," features after preprocessing\np-value for overall differences between groups using PC1 and PC3 in a multivariate\n one-way ANOVA model=",round(pc2_pval,3),sep="")
           }
           
           #     plot(c(1,1),plot=FALSE)
           #l <- legend(0, 0, bty='n', l1,plot=FALSE, pch = pch_per_group, pt.cex = 0.6)
           # calculate right margin width in ndc
           #w <- grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
            w <- 0.1
           par(omd=c(0, 1-w, 0, 1))
           
           iqr_xlim=4*sd(result$variates$X[,1],na.rm=TRUE) #(summary(result$variates$X[,1])[5]-summary(result$variates$X[,1])[3])
           iqr_ylim=4*sd(result$variates$X[,3],na.rm=TRUE) #(summary(result$variates$X[,3])[5]-summary(result$variates$X[,3])[3])
           
           
      print(dataEllipse(x=result$x[,1], y=result$x[,3],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=c(ellipse.conf.level),col=col_per_group,pch=pch,xlab=paste("PC1 (",r1[1],"% variation)",sep=""),ylab=paste("PC3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE,cex.main=0.8,center.pch=FALSE,ylim=c(min(result$variates$X[,3])-iqr_ylim,max(result$variates$X[,3])+iqr_ylim),xlim=c(min(result$variates$X[,1])-iqr_xlim,max(result$variates$X[,1])+iqr_xlim)))
      
      #print(legend(legendlocation, l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.7, title = "Class", cex=legendcex))
      print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
      
      
       w <- 0.1
      par(omd=c(0, 1-w, 0, 1))
      
      print(dataEllipse(x=result$x[,1], y=result$x[,3],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=NULL,col=col_per_group,pch=pch,xlab=paste("PC1 (",r1[1],"% variation)",sep=""),ylab=paste("PC3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=FALSE,add=FALSE,cex.main=0.8,ylim=c(min(result$variates$X[,3])-iqr_ylim,max(result$variates$X[,3])+iqr_ylim),xlim=c(min(result$variates$X[,1])-iqr_xlim,max(result$variates$X[,1])+iqr_xlim)))
      
      # print(legend(legendlocation, l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.7, title = "Class", cex=legendcex))
      print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
      
      
      w <- 0.1
      par(omd=c(0, 1-w, 0, 1))
      
      iqr_xlim=4*sd(result$variates$X[,2],na.rm=TRUE) #2.5*(summary(result$variates$X[,2])[5]-summary(result$variates$X[,2])[3])
      iqr_ylim=4*sd(result$variates$X[,3],na.rm=TRUE) #2.5*(summary(result$variates$X[,3])[5]-summary(result$variates$X[,3])[3])
      
      
       print(dataEllipse(x=result$x[,2], y=result$x[,3],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=c(ellipse.conf.level),col=col_per_group,pch=pch,xlab=paste("PC2 (",r1[2],"% variation)",sep=""),ylab=paste("PC3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE,cex.main=0.8,center.pch=FALSE,ylim=c(min(result$variates$X[,3])-iqr_ylim,max(result$variates$X[,3])+iqr_ylim),xlim=c(min(result$variates$X[,2])-iqr_xlim,max(result$variates$X[,2])+iqr_xlim)))
       
       #print(legend(legendlocation, l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.7, title = "Class", cex=legendcex))
       print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
       
       
       if(do_pca_anova==TRUE){
           
           main_text=paste("Pairwise PC score plots using ",filename," features after preprocessing\np-value for overall differences between groups using PC2 and PC3 in a multivariate\n one-way ANOVA model=",round(pc3_pval,3),sep="")
       }
       
      
       
        w <- 0.1
       par(omd=c(0, 1-w, 0, 1))
       
       print(dataEllipse(x=result$x[,2], y=result$x[,3],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=NULL,col=col_per_group,pch=pch,xlab=paste("PC2 (",r1[2],"% variation)",sep=""),ylab=paste("PC3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=FALSE,add=FALSE,cex.main=0.8,ylim=c(min(result$variates$X[,3])-iqr_ylim,max(result$variates$X[,3])+iqr_ylim),xlim=c(min(result$variates$X[,2])-iqr_xlim,max(result$variates$X[,2])+iqr_xlim)))
       #print(legend(legendlocation, l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.7, title = "Class", cex=legendcex))
       print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
       
       
       
       }
       
       
      

}
   
   
   
   Class<-samplelabels
   
   if(ncol(result$x)>5){
       
       melted <- cbind(Class, melt(result$x[,1:5]))
       
   }else{
       
       melted <- cbind(Class, melt(result$x))
   }
   
   colnames(melted)<-c("Class","Samples","Var2","PCscore")
   melted$Samples<-seq(1,nrow(result$x))
   
   # for(i in 1:5)
   
   testname=""
   lapply(1:5,function(i){
       
       pcname<-paste("PC",i,sep="")
       melted_pc1<-melted[which(melted$Var2==pcname),]
       
       if(do_pca_anova==TRUE){
           
          
       main_text=paste("PC score plots using ",filename," features after preprocessing\np-value for differences between groups using ",pcname,"\n",testname," p=",round(pc_pval_single[i],3),sep="")
       }else{
           
            main_text=paste("PC score plots using ",filename," features after preprocessing",sep="")
       }
       
   
       w <- 0.1
       par(omd=c(0, 1-w, 0, 1))
       
       # barCenters=barplot(myData$Intensity,ylab="Intensity",xlab="Class",main=mzlabel_cur,col=barplot.col.opt1,ylim=c(0,max_yval))
       #arrows(barCenters, ymin,barCenters, ymax,angle=90,code=3,lty=1,lwd=1.25,length=0.05)
       ylab1<-paste(pcname,"score",sep="")
       plot(as.vector(melted_pc1[,4]),col=c(col),main=main_text, ylab=ylab1,xlab="Sample",type="h",lwd=2,cex.main=0.8)
       
       print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,unique(melted_pc1$Class), col = col_per_group,pch = rep(19,length(col_per_group)), pt.cex = 0.6, title = "Class",cex=0.8))
       
       
       #print(barplot1)
       
       
       
   })
   
   # #save(list=ls(),file="getpcadebug.Rda")

 
    return(list(result=result,pca_pval_vec=pc_pval_single))
}

get_plsplots<-function(X,plsres,plsvar,samplelabels,filename=NA,ncomp=5,center=TRUE,scale=TRUE,legendcex=0.5,outloc=getwd(),col_vec=NA,sample.col.opt="default",alphacol=0.3,legendlocation="topright",class_levels=NA,pca.cex.val=3,ellipse.conf.level=0.95,main_text="PLS-DA score plots"){
   
    result<-plsres
    r1<-plsvar
   
    pch.val=19
    legendlocation="bottomleft"
    
    samplelabels<-as.data.frame(samplelabels)
    samplelabels<-as.factor(samplelabels[,1])
    l2<-levels(as.factor(samplelabels))
    col_all=topo.colors(256)
    
    t1<-table(samplelabels)
    if(is.na(class_levels)==TRUE){
       
        l1<-levels(as.factor(samplelabels))
    }else{
        l1<-class_levels

        
    }

    class_labels_levels<-l1
    
    ncomp=min(dim(X)[1],dim(X)[2])
  
    if(is.na(col_vec)==TRUE){
        col_vec<-c("mediumpurple4","mediumpurple1","blueviolet","darkblue","blue","cornflowerblue","cyan4","skyblue",
        "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
        "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
        "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
        
    }
    
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
                        if(sample.col.opt=="colorblind"){
                            #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                            # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                            
                            if(length(class_labels_levels)<9){
                                
                                col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7", "grey57")
                                
                            }else{
                                
                                #col_vec<-colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels))
                                col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                                
                            }
                            
                            
                        }else{
                            
                            check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                            
                            if(length(check_brewer)>0){
                                
                                sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                                col_vec <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(class_labels_levels))
                                
                            }else{
                                
                                #colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                if(sample.col.opt=="journal"){
                                    
                                    col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                    "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                    "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                    
                                    "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                    "#E64B3519","#4DBBD519","#631879E5","grey75")
                                    
                                }else{
                                    
                                    if(length(sample.col.opt)==1){
                                        col_vec <-rep(sample.col.opt,length(class_labels_levels))
                                    }else{
                                    
                                        colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                        
                                    }
                                }
                                
                            }
                            
                        }
                        
                    }
                    
                    
                }
                
            }
            
        }
    }
    #col_vec<-col_vec[sample(1:length(col_vec),length(col_vec))]
    
    l1<-gsub(x=l1,pattern="Class",replacement="")
    
    dir.create(outloc,showWarnings=FALSE)
    setwd(outloc)
  
    ## 1) raw data
    #tiff(fname,res=300, width=2000,height=2000)
    col <- rep(col_vec[1:length(t1)], t1)
    #col<-rep(col_all[1:length(l1)],t1)
    ## Choose different size of points
    cex <- rep(pca.cex.val, dim(X)[1])
    ## Choose the form of the points (square, circle, triangle and diamond-shaped
    
    
    
    #pch_vec<-seq(1,50) #c(3,5,7,9,12,13,2,17,21) #seq(1,50) #
    pch_vec <- rep(pch.val,dim(X)[1])
    pch=pch_vec
    cex <- rep(pca.cex.val, dim(X)[1])
    col_per_group<-{}
    pch_per_group<-{}
    for(p1 in 1:length(l2)){
        
        pch[which(samplelabels==l2[p1])]=pch_vec[p1]
        col[which(samplelabels==l2[p1])]=col_vec[p1]
        col_per_group<-c(col_per_group,col_vec[p1])
        pch_per_group<-c(pch_per_group,pch_vec[p1])
    }
 

 #pch_vec<-seq(1,50) #c(3,5,7,9,12,13,2,17,21) #seq(1,50) #
 #pch_vec <- rep(pch.val,dim(X)[1])
cex <- rep(pca.cex.val, dim(X)[1])


##save(list=ls(),file="plsdebug.Rda")

#print(plotIndiv(result, comp = c(1,2),ind.names = FALSE, group=samplelabels, cex = cex[1], pch = pch, ellipse=FALSE, ellipse.level = 0.95, X.label=paste("PLS1 (",r1[1],"% variation)",sep=""),Y.label=paste("PLS2 (",r1[2],"% variation)",sep=""),add.legend=TRUE))


     if(result$ncomp>1){
    
      w <- 0.1
    par(omd=c(0, 1-w, 0, 1))
    # print(dataEllipse(x=result$variates$X[,1], y=result$variates$X[,2],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=c(ellipse.conf.level),pch=pch,col=col_per_group,center.pch=NULL,xlab=paste("PLS1 (",r1[1],"% variation)",sep=""),ylab=paste("PLS2 (",r1[2],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE,cex.main=0.7))
    
    
    iqr_xlim=3*(sd(result$variates$X[,1],na.rm=TRUE)) #2.5*(summary(result$variates$X[,1])[5]-summary(result$variates$X[,1])[3])
    iqr_ylim=3*(sd(result$variates$X[,2],na.rm=TRUE)) #2.5*(summary(result$variates$X[,2])[5]-summary(result$variates$X[,2])[3])
    
    main_text1=paste(main_text,"\n with elliptical contours at 95% confidence interval",sep="")
     print(dataEllipse(x=result$variates$X[,1], y=result$variates$X[,2],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=c(ellipse.conf.level),col=col_per_group,pch=pch,
     xlab=paste("PLS1 (",r1[1],"% variation)",sep=""),ylab=paste("PLS2 (",r1[2],"% variation)",sep=""),group.labels="",main=main_text1,fill=TRUE,cex.main=0.8,center.pch=FALSE,
     ylim=c(min(result$variates$X[,2])-iqr_ylim,max(result$variates$X[,2])+iqr_ylim),xlim=c(min(result$variates$X[,1])-iqr_xlim,max(result$variates$X[,1])+iqr_xlim)))
     
     # axis(side = 1, at=seq(min(result$variates$X[,1])-iqr_xlim,max(result$variates$X[,1])+iqr_xlim,10))
     #axis(side = 2, at=seq(min(result$variates$X[,2])-iqr_ylim,max(result$variates$X[,2])+iqr_ylim,10))
     
    
     print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
     
     
     w <- 0.1
     par(omd=c(0, 1-w, 0, 1))
     #   print(dataEllipse(x=result$variates$X[,1], y=result$variates$X[,2],groups=as.factor(samplelabels),grid=FALSE,lwd=1,levels=NULL,pch=pch,col=col_per_group,center.pch=NULL,xlab=paste("PLS1 (",r1[1],"% variation)",sep=""),ylab=paste("PLS2 (",r1[2],"% variation)",sep=""),group.labels="",main=main_text,fill=FALSE,add=FALSE,cex.main=0.7))
     
 #    iqr_xlim=2.5*(summary(result$variates$X[,1])[5]-summary(result$variates$X[,1])[3])
    # iqr_ylim=2.5*(summary(result$variates$X[,2])[5]-summary(result$variates$X[,2])[3])
	
	  iqr_xlim=3*(sd(result$variates$X[,1],na.rm=TRUE)) #2.5*(summary(result$variates$X[,1])[5]-summary(result$variates$X[,1])[3])
    iqr_ylim=3*(sd(result$variates$X[,2],na.rm=TRUE)) #2.5*(summary(result$variates$X[,2])[5]-summary(result$variates$X[,2])[3])
   


     print(dataEllipse(x=result$variates$X[,1], y=result$variates$X[,2],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=NULL,col=col_per_group,center.pch=NULL,pch=pch,xlab=paste("PLS1 (",r1[1],"% variation)",sep=""),ylab=paste("PLS2 (",r1[2],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE,cex.main=0.8,add=FALSE,ylim=c(min(result$variates$X[,2])-iqr_ylim,max(result$variates$X[,2])+iqr_ylim),xlim=c(min(result$variates$X[,1])-iqr_xlim,max(result$variates$X[,1])+iqr_xlim)))
     
     print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
     
     #  print(legend(legendlocation, l1, col = col_per_group,pch = pch_per_group, pt.cex = cex[1], title = "Class", cex=cex[1]))
     #mtext(l1,side=3,line=1,cex=0.6,adj=NA,col = col_per_group,pch = pch_per_group)
     
     }
    if(result$ncomp>2){
        
        #      print(plotIndiv(result, comp = c(1,3),ind.names = FALSE, group=samplelabels, cex = cex[1], pch = pch, ellipse=FALSE, ellipse.level = 0.95, X.label=paste("PLS1 (",r1[1],"% variation)",sep=""),Y.label=paste("PLS3 (",r1[3],"% variation)",sep=""),add.legend=TRUE))

         w <- 0.1
        par(omd=c(0, 1-w, 0, 1))
        #  print(dataEllipse(x=result$variates$X[,1], y=result$variates$X[,3],groups=as.factor(samplelabels),grid=FALSE,lwd=1,levels=c(ellipse.conf.level),pch=pch,col=col_per_group,center.pch=NULL,xlab=paste("PLS1 (",r1[1],"% variation)",sep=""),ylab=paste("PLS3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE,cex.main=0.7))
       
  #     iqr_xlim=2.5*(summary(result$variates$X[,1])[5]-summary(result$variates$X[,1])[3])
     #  iqr_ylim=2.5*(summary(result$variates$X[,3])[5]-summary(result$variates$X[,3])[3])
     
       iqr_xlim=3*(sd(result$variates$X[,1],na.rm=TRUE)) #2.5*(summary(result$variates$X[,1])[5]-summary(result$variates$X[,1])[3])
      iqr_ylim=3*(sd(result$variates$X[,3],na.rm=TRUE)) #2.5*(summary(result$variates$X[,2])[5]-summary(result$variates$X[,2])[3])
   
       
       main_text1=paste(main_text,"\n with elliptical contours at 95% confidence interval",sep="")
       print(dataEllipse(x=result$variates$X[,1], y=result$variates$X[,3],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=c(ellipse.conf.level),col=col_per_group,pch=pch,xlab=paste("PLS1 (",r1[1],"% variation)",sep=""),ylab=paste("PLS3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text1,fill=TRUE,cex.main=0.8,center.pch=FALSE,ylim=c(min(result$variates$X[,3])-iqr_ylim,max(result$variates$X[,3])+iqr_ylim),xlim=c(min(result$variates$X[,1])-iqr_xlim,max(result$variates$X[,1])+iqr_xlim)))
       
       print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
        
  
  w <- 0.1
  par(omd=c(0, 1-w, 0, 1))
  #  print(dataEllipse(x=result$variates$X[,1], y=result$variates$X[,3],groups=as.factor(samplelabels),grid=FALSE,lwd=1,levels=NULL,pch=pch,col=col_per_group,center.pch=NULL,xlab=paste("PLS1 (",r1[1],"% variation)",sep=""),ylab=paste("PLS3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=FALSE,add=FALSE,cex.main=0.7))
#  iqr_xlim=2.5*(summary(result$variates$X[,1])[5]-summary(result$variates$X[,1])[3])
 # iqr_ylim=2.5*(summary(result$variates$X[,3])[5]-summary(result$variates$X[,3])[3])
 
   iqr_xlim=3*(sd(result$variates$X[,1],na.rm=TRUE)) #2.5*(summary(result$variates$X[,1])[5]-summary(result$variates$X[,1])[3])
    iqr_ylim=3*(sd(result$variates$X[,3],na.rm=TRUE)) #2.5*(summary(result$variates$X[,2])[5]-summary(result$variates$X[,2])[3])
   
  
  print(dataEllipse(x=result$variates$X[,1], y=result$variates$X[,3],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=NULL,col=col_per_group,pch=pch,xlab=paste("PLS1 (",r1[1],"% variation)",sep=""),ylab=paste("PLS3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=FALSE,cex.main=0.8,center.pch=NULL,add=FALSE,ylim=c(min(result$variates$X[,3])-iqr_ylim,max(result$variates$X[,3])+iqr_ylim),xlim=c(min(result$variates$X[,1])-iqr_xlim,max(result$variates$X[,1])+iqr_xlim)))
  
  
  print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
  

        #    print(plotIndiv(result, comp = c(2,3),ind.names = FALSE, group=samplelabels, cex = cex[1], pch = pch, ellipse=FALSE, ellipse.level = 0.95, X.label=paste("PLS2 (",r1[2],"% variation)",sep=""),Y.label=paste("PLS3 (",r1[3],"% variation)",sep=""),add.legend=TRUE))
        
         w <- 0.1
        par(omd=c(0, 1-w, 0, 1))
        #   print(dataEllipse(x=result$variates$X[,2], y=result$variates$X[,3],groups=as.factor(samplelabels),grid=FALSE,lwd=1,levels=c(ellipse.conf.level),pch=pch,col=col_per_group,center.pch=NULL,xlab=paste("PLS2 (",r1[2],"% variation)",sep=""),ylab=paste("PLS3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=TRUE,cex.main=0.7))
        #print(legend(legendlocation, l1, col = col_per_group,pch = pch_per_group, pt.cex = cex[1], title = "Class", cex=cex[1]))
       # iqr_xlim=2.5*(summary(result$variates$X[,2])[5]-summary(result$variates$X[,2])[3])
       # iqr_ylim=2.5*(summary(result$variates$X[,3])[5]-summary(result$variates$X[,3])[3])
       
         iqr_xlim=3*(sd(result$variates$X[,1],na.rm=TRUE)) #2.5*(summary(result$variates$X[,1])[5]-summary(result$variates$X[,1])[3])
        iqr_ylim=3*(sd(result$variates$X[,2],na.rm=TRUE)) #2.5*(summary(result$variates$X[,2])[5]-summary(result$variates$X[,2])[3])
   
        
        main_text1=paste(main_text,"\n with elliptical contours at 95% confidence interval",sep="")
        
        print(dataEllipse(x=result$variates$X[,2], y=result$variates$X[,3],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=c(ellipse.conf.level),col=col_per_group,pch=pch,xlab=paste("PLS2 (",r1[2],"% variation)",sep=""),ylab=paste("PLS3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text1,fill=TRUE,cex.main=0.8,center.pch=FALSE,ylim=c(min(result$variates$X[,3])-iqr_ylim,max(result$variates$X[,3])+iqr_ylim),xlim=c(min(result$variates$X[,2])-iqr_xlim,max(result$variates$X[,2])+iqr_xlim)))
        
        
        print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
        
        
        w <- 0.1
        par(omd=c(0, 1-w, 0, 1))
        
     #   iqr_xlim=2.5*(summary(result$variates$X[,2])[5]-summary(result$variates$X[,2])[3])
      #  iqr_ylim=2.5*(summary(result$variates$X[,3])[5]-summary(result$variates$X[,3])[3])
      
        iqr_xlim=3*(sd(result$variates$X[,2],na.rm=TRUE)) #2.5*(summary(result$variates$X[,1])[5]-summary(result$variates$X[,1])[3])
    iqr_ylim=3*(sd(result$variates$X[,3],na.rm=TRUE)) #2.5*(summary(result$variates$X[,2])[5]-summary(result$variates$X[,2])[3])
   
        
        print(dataEllipse(x=result$variates$X[,2], y=result$variates$X[,3],groups=as.factor(samplelabels),grid=FALSE,lwd=0.5,levels=NULL,col=col_per_group,pch=pch,xlab=paste("PLS2 (",r1[2],"% variation)",sep=""),ylab=paste("PLS3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=FALSE,cex.main=0.8,center.pch=NULL,add=FALSE,ylim=c(min(result$variates$X[,3])-iqr_ylim,max(result$variates$X[,3])+iqr_ylim),xlim=c(min(result$variates$X[,2])-iqr_xlim,max(result$variates$X[,2])+iqr_xlim)))
        
        
        
        # print(dataEllipse(x=result$variates$X[,2], y=result$variates$X[,3],groups=as.factor(samplelabels),grid=FALSE,lwd=1,levels=NULL,pch=pch,col=col_per_group,center.pch=NULL,xlab=paste("PLS2 (",r1[2],"% variation)",sep=""),ylab=paste("PLS3 (",r1[3],"% variation)",sep=""),group.labels="",main=main_text,fill=FALSE,cex.main=0.7,add=FALSE))
        #print(legend(legendlocation, l1, col = col_per_group,pch = pch_per_group, pt.cex = cex[1], title = "Class", cex=cex[1]))
        print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,l1, col = col_per_group,pch = pch_per_group, pt.cex = 0.6, title = "Class",cex=0.8))
        
    }

    

}


get_barplots<-function(feature_table_file,class_labels_file,X=NA,Y=NA,parentoutput_dir,newdevice=FALSE,ylabel="Intensity",bar.colors=NA,cex.val=0.6,barplot.col.opt="journal",error.bar=TRUE,barplot.xaxis="Factor2"){
    
    
    xaxis=barplot.xaxis
    if(is.na(X)==TRUE){
        data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE)
    }else{
        data_matrix<-X
        rm(X)
        
    }
    
    dir.create(parentoutput_dir)
    setwd(parentoutput_dir)
    
    
    mzlabels<-data_matrix[,1]
    
    timelabels<-data_matrix[,2]

    data_m<-data_matrix[,-c(1:2)]
    
    data_m<-as.matrix(data_m)
    
    col_samples<-TRUE
    
    if(is.na(Y)==TRUE){
        if(is.na(class_labels_file)==FALSE){
            classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
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
        
        tabbed_groups<-TRUE
        
    }else{
        
        Class<-classlabels[,2]
        
        tabbed_groups=FALSE
    }
    
  
    class_labels_levels<-levels(as.factor(Class))
    
    
    if(is.na(barplot.col.opt)==FALSE)
    {
        if(barplot.col.opt=="default"){
            
            col_vec<-c("#CC0000","#AAC000","blue","mediumpurple4","mediumpurple1","blueviolet","cornflowerblue","cyan4","skyblue",
            "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
            "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
            "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
            
        }else{
            if(barplot.col.opt=="topo"){
                #col_vec<-topo.colors(256) #length(class_labels_levels))
                
                #col_vec<-col_vec[seq(1,length(col_vec),)]
                
                col_vec <- topo.colors(length(class_labels_levels), alpha=alphacol)
            }else{
                if(barplot.col.opt=="heat"){
                    #col_vec<-heat.colors(256) #length(class_labels_levels))
                    
                    col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
                }else{
                    if(barplot.col.opt=="rainbow"){
                        #col_vec<-heat.colors(256) #length(class_labels_levels))
                        col_vec<-rainbow(length(class_labels_levels), start = 0, end = alphacol)
                        
                        #col_vec <- heat.colors(length(class_labels_levels), alpha=alphacol)
                    }else{
                        
                        if(barplot.col.opt=="terrain"){
                            
                            col_vec <- cm.colors(length(class_labels_levels), alpha=alphacol)
                        }else{
                            if(is.na(barplot.col.opt)==TRUE){
                                col_vec<-c("black")
                            }else{
                                
                                if(barplot.col.opt=="colorblind"){
                                    #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                                    # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                                    
                                    if(length(class_labels_levels)<9){
                                        
                                        col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7", "grey57")
                                        
                                    }else{
                                        
                                        #col_vec<-colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels))
                                        col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                        "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                                        
                                    }
                                    
                                    
                                }else{
                                    
                                    check_brewer<-grep(pattern="brewer",x=barplot.col.opt)
                                    
                                    if(length(check_brewer)>0){
                                        barplot.col.opt=gsub(x=barplot.col.opt,pattern="brewer.",replacement="")
                                        col_vec <- colorRampPalette(brewer.pal(10, barplot.col.opt))(length(class_labels_levels))
                                        
                                    }else{
                                        
                                        if(barplot.col.opt=="journal"){
                                            
                                            col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                            "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                            "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                            
                                            "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                            "#E64B3519","#4DBBD519","#631879E5","grey75")
                                            
                                        }else{
                                            col_vec <-barplot.col.opt
                                        }
                                        
                                    }
                                    
                                }
                            }
                        }
                        
                        
                    }
                    
                }
                
            }
        }
    }else{
        
        col_vec<-c("grey57")
        
    }

    
    mtcars<-t(data_m)
    Class<-as.character(Class)
    
  
    
    
    
    if(newdevice==TRUE){
    
    pdf("barplots.pdf")
    
    }
    
    par(mfrow=c(1,1),family="sans",cex=cex.val)
    # par(mfrow=c(par_rows,max_per_row))
    
    myData_list<-new("list")
    
   
    
   time_start<-Sys.time()
   save(list=ls(),file="debugbarplot.Rda")
   if(tabbed_groups==FALSE){
       
       mtcars<-cbind(Class,mtcars)
       mtcars<-as.data.frame(mtcars)
       mtcars <- do.call(data.frame, mtcars)
        
       mtcars_sum<-do.call(data.frame,aggregate(list(mtcars[-c(1)]),by=list(mtcars$Class),FUN = function(x) {x<-as.numeric(as.character(x));c(mean = mean(x), sd = sd(x),n = length(x),se=sd(x)/sqrt(length(x)))}))
       
       label_inc_list<-seq(2,dim(mtcars_sum)[2],4) #seq(1,dim(mtcars_sum)[2],3)
       ##save(list=ls(),file="debugbarplot.Rda")
       #  for(i in 2:dim(mtcars)[2]){
       myData_list<-lapply(seq(2,dim(mtcars_sum)[2],4),function(i){
           
           myData<-mtcars_sum[,c(1,i:(i+3))]
           colnames(myData) <- c("Class", "Intensity", "sd", "n", "se")
           
           get_label_ind<-which(label_inc_list==i)
           round_mzval<-sprintf("%.4f",mzlabels[get_label_ind])
           round_timeval<-sprintf("%.1f",timelabels[get_label_ind])
       
           mzlabel_cur<-paste("mz_time: ",round_mzval,"_",round_timeval,sep="")
         
           ymax = myData$Intensity + 1.96*myData$se
           
           ymin = myData$Intensity - 1.96*myData$se
           
           max_yval<-ceiling(max((myData$Intensity + (2.5*myData$se)),na.rm=TRUE)) #round(max(myData$Intensity+(4*myData$se),na.rm=TRUE))
           
           
           min_yval<-max(0,floor(min((myData$Intensity - (2.5*myData$se)),na.rm=TRUE)))
           
           below_zero_check<-which(ymin<0)
           if(length(below_zero_check)>0){
               
               ymin[below_zero_check]<-0
           }
           
           
           
           myData$Class<-as.factor(myData$Class)
           
           colnames(myData) <- c("Class", "Intensity", "sd", "n", "se")
           
           t1<-table(myData$Class)
           
           barplot.col.opt1=rep(col_vec[1:length(t1)],t1)
           
           
           
           w <- 0.1
           par(omd=c(0, 1-w, 0, 1))
           ylim=range(pretty(c(min_yval,max_yval)))
           rnum=length(unique(myData$Class))
           if(rnum<7)
           {
               
               space_vec=c(2,2,1.5,1,0.5,0.25)
               name_vec=unique(myData$Class)
               barCenters=barplot(myData$Intensity,ylab=ylabel,xlab="",main=mzlabel_cur,col=barplot.col.opt1,width=0.125,xlim=c(0,1),names.arg=name_vec[1:rnum],xpd=FALSE,ylim=ylim,space=space_vec[rnum])
            
           }else{
               barCenters=barplot(myData$Intensity,ylab=ylabel,xlab="",main=mzlabel_cur,col=barplot.col.opt1,space=0.1,names.arg=unique(myData$Class),xpd=FALSE,ylim=ylim)
           }
           
           if(error.bar==TRUE){
               arrows(barCenters, ymin,barCenters, ymax,angle=90,code=3,lty=1,lwd=1.25,length=0.05)
           }
           
           
           
           print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,unique(myData$Class), col = col_vec[1:length(t1)],pch = rep(19,length(col_vec[1:length(t1)])), pt.cex = 0.6, title = "Class",cex=0.8))
           
           return(myData)
           
       })
       
   }else{
       
       mtcars<-cbind(classlabels[,2:3],mtcars)
       mtcars<-as.data.frame(mtcars)
       cnames_mtcars<-colnames(mtcars)
       cnames_mtcars[1:2]<-c("Factor1","Factor2")
       colnames(mtcars)<-cnames_mtcars
       mtcars <- do.call(data.frame, mtcars)
       
       mtcars_sum<-do.call(data.frame,aggregate(list(mtcars[-c(1:2)]),by=list(Factor1=mtcars$Factor1,Factor2=mtcars$Factor2),FUN = function(x) {x<-as.numeric(as.character(x));c(mean = mean(x), sd = sd(x),n = length(x),se=sd(x)/sqrt(length(x)))}))
        
      
       
       
       label_inc_list<-seq(3,dim(mtcars_sum)[2],4) #seq(1,dim(mtcars_sum)[2],3)
       ##save(list=ls(),file="debugbarplot.Rda")
       #  for(i in 2:dim(mtcars)[2]){
       myData_list<-lapply(seq(3,dim(mtcars_sum)[2],4),function(i){
           
           myData<-mtcars_sum[,c(1:2,i:(i+3))]
           colnames(myData) <- c("Factor1", "Factor2", "Intensity", "sd", "n", "se")
           
           if(xaxis=="Factor1"){
               tabbedMeans <- tapply(myData$Intensity, list(myData$Factor2,myData$Factor1),
               function(x) c(x = x))
               
               tabbedSE <- tapply(myData$se, list(myData$Factor2,myData$Factor1),
               function(x) c(x = x))
           }else{
               
               if(xaxis=="Factor2"){
                   tabbedMeans <- tapply(myData$Intensity, list(myData$Factor1,myData$Factor2),
                   function(x) c(x = x))
                   
                   tabbedSE <- tapply(myData$se, list(myData$Factor1,myData$Factor2),
                   function(x) c(x = x))
               }
               
               
           }
           get_label_ind<-which(label_inc_list==i)
           round_mzval<-sprintf("%.4f",mzlabels[get_label_ind])
           round_timeval<-sprintf("%.1f",timelabels[get_label_ind])
           
           mzlabel_cur<-paste("mz_time: ",round_mzval,"_",round_timeval,sep="")
           
           ymax = tabbedMeans + 1.96*tabbedSE
           
           ymin = tabbedMeans - 1.96*tabbedSE
           
           max_yval<-ceiling(max((tabbedMeans + (2.5*tabbedSE)),na.rm=TRUE)) #round(max(myData$Intensity+(4*myData$se),na.rm=TRUE))
           
           min_yval<-max(0,floor(min((tabbedMeans - (2.5*tabbedSE)),na.rm=TRUE)))
           
           below_zero_check<-which(ymin<0)
           
           if(length(below_zero_check)>0){
               
               ymin[below_zero_check]<-0
           }
           
           if(xaxis=="Factor2"){
               
               t1<-table(factor(myData$Factor1))
               t2<-table(factor(myData$Factor2))
           }else{
               t1<-table(factor(myData$Factor2))
               t2<-table(factor(myData$Factor1))
           }
           
           barplot.col.opt1=rep(col_vec[1:length(t1)],length(t2))
           
           w <- 0.1
           par(omd=c(0, 1-w, 0, 1))
           ylim=range(pretty(c(min_yval,max_yval)))
           
          
          barCenters=barplot(tabbedMeans,beside=TRUE, ylab=ylabel,xlab="",main=mzlabel_cur,col=barplot.col.opt1,las=1,names.arg=unique(myData$Class),xpd=FALSE,ylim=ylim,border = "black")
          
          segments(barCenters, ymin, barCenters,ymax, lwd = 1.25)
          
           if(error.bar==TRUE){
               # arrows(barCenters, ymin,barCenters, ymax,angle=90,code=3,lty=1,lwd=1.25,length=0.05)
               arrows(barCenters, ymin, barCenters, ymax, lwd = 1.25, angle = 90, code = 3, length = 0.05)
           }
           
           if(xaxis=="Factor1"){
               print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,unique(myData$Factor2), col = col_vec[1:length(t1)],pch = rep(19,length(col_vec[1:length(t1)])), pt.cex = 0.6, title = "",cex=0.8))
           }else{
                print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,unique(myData$Factor1), col = col_vec[1:length(t1)],pch = rep(19,length(col_vec[1:length(t1)])), pt.cex = 0.6, title = "",cex=0.8))
               
           }
           return(myData)
           
       })
       
       
       
       
       
       
   }
   
    
    
    if(newdevice==TRUE){
        try(dev.off(),silent=TRUE)
    }
    
    suppressWarnings(dir.create("Tables"))
   
    #save(myData_list,file="Tables/barplots_data.Rda")

    
}



get_individualsampleplots<-function(feature_table_file,class_labels_file,X=NA,Y=NA,parentoutput_dir,newdevice=FALSE,ylabel="Intensity",bar.colors=NA,cex.val=0.75,sample.col.opt=c("grey57"),plottype="barplot"){
    
    if(is.na(X)==TRUE){
        data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE)
    }else{
        data_matrix<-X
        rm(X)
        
    }
    
     par(mfrow=c(1,1),family="sans",cex=cex.val)
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
            classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
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
                        
                        if(sample.col.opt=="colorblind"){
                            #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                            # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                            
                            if(length(class_labels_levels)<9){
                                
                                col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7", "grey57")
                                
                            }else{
                                
                                #col_vec<-colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels))
                                col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                                
                            }
                            
                            
                        }else{
                            
                            check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                            
                            if(length(check_brewer)>0){
                                
                                 sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                                col_vec <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(class_labels_levels))
                                
                            }else{
                                
                                if(sample.col.opt=="journal"){
                                    
                                    col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                    "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                    "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                    
                                    "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                    "#E64B3519","#4DBBD519","#631879E5","grey75")
                                    
                                }else{
                                    #col_vec <-rep(sample.col.opt,length(class_labels_levels))
                                    if(length(sample.col.opt)==1){
                                        col_vec <-rep(sample.col.opt,length(class_labels_levels))
                                    }else{
                                        
                                        colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                        
                                    }
                                }
                                
                            }
                            
                        }
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
    
    col_per_group<-{}
    pch_per_group<-{}
    
    l2<-levels(as.factor(Class))
    col1<-rep("red",length(Class))
    
    
    
    par_rows=2
    max_per_row=2
    
    # par(mfrow=c(par_rows,max_per_row))
    
    myData_list<-new("list")
    
    mtcars <- do.call(data.frame, mtcars)
    
    ##save(mtcars,file="mtcars.Rda")
    
    
    
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
        t1<-table(myData$Class)
        if(length(barplot.col.opt)<2){
            
            barplot.col.opt1=rep(barplot.col.opt,length(myData$Class))
        }else{
            
            
            if(length(barplot.col.opt)==length(levels(factor(myData$Class)))){
        
                    #barplot.col.opt1=rep(barplot.col.opt,t1)
                    
                    barplot.col.opt1=rep(barplot.col.opt[1:length(t1)],t1)
            
            
            }else{
                
                # print("Number of classes is greater than the length of the color vector. Using default colors.")
                #col_clust<-topo.colors(length(t1))
                #barplot.col.opt1=col1
                #t1<-table(mtcars$Class)
                #barplot.col.opt1=rep(barplot.col.opt,t1)
                barplot.col.opt1=rep(barplot.col.opt[1:length(t1)],t1)
            }
            
            
        }
        
        #print(barplot.col.opt1)
        #barplot.col.opt1=barplot.col.opt
        
        # print(barplot.col.opt1)
        
        myData$Class<-factor(myData$Class)
        
        Samples<-seq(1,nrow(myData))
        Var2<-rep(1,nrow(myData))
        myData<-cbind(myData,Samples,Var2)
        
        # #save(list=ls(),file="barplotdebug.Rda")
        
        if(plottype=="barplot"){
            #barplot(as.vector(myData[,2]),col=c(barplot.col.opt1),main=mzlabel_cur, ylab="Intensity",xlab="Sample") #,ylim=c(min(myData[,2])-1,max(myData[,2])+1),xpd=FALSE)
               
            
            w <- 0.1
            par(omd=c(0, 1-w, 0, 1))
             
	     plot(as.vector(myData[,2]),col=c(barplot.col.opt1),main=mzlabel_cur, ylab=ylabel,xlab="Sample",type="h",lwd=2,ylim=range(pretty(c(0,myData[,2]))))
            
          
            
            print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,levels(as.factor(mtcars$Class)), col = col_vec[1:length(t1)],pch = rep(19,length(col_vec[1:length(t1)])), pt.cex = 0.6, title = "Class",cex=0.8))
            
            
            if(FALSE){
             barplot1 <- ggplot(data=myData) + geom_bar(aes(x=Samples, y=Intensity, fill=Class), stat="identity") + facet_wrap(~Var2) + scale_fill_manual(values=unique(barplot.col.opt1)) + ggtitle(mzlabel_cur) + ylab(ylabel) + xlab("Samples") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.margin=unit(c(10,5,5,5),"mm"),
             strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
             strip.text = element_text(face="bold"))
            }
            
            #print(barplot1)
             
        }else{
            if(plottype=="point"){
                plot(as.vector(myData[,2]),col=c(barplot.col.opt1),main=mzlabel_cur, ylab=ylab,xlab="Samples",type="p")
            }else{
                
                if(plottype=="point_grouped"){
                    p <- ggplot(data = myData, aes(x = Class, y = Intensity,fill=Class)) + scale_fill_manual(values = c(barplot.col.opt1)) + ylab(ylabel)
                    print(p+geom_point(size=2)+ggtitle(mzlabel_cur)) #,aes(colour=barplot.col.opt))
                }
            }
        }
        
        
    })
    if(newdevice==TRUE){
        try(dev.off(),silent=TRUE)
    }
    
   
    
    
}




get_hca<-function(feature_table_file=NA,parentoutput_dir,class_labels_file=NA,X=NA,Y=NA,heatmap.col.opt="RdBu",cor.method="spearman",is.data.znorm=FALSE,analysismode="classification",
sample.col.opt="rainbow",plots.width=8,plots.height=8,plots.res=600, plots.type="cairo", alphacol=0.3, hca_type="two-way",newdevice=FALSE,input.type="intensity",mainlab="",cexRow=1, cexCol=1,plot.bycluster=FALSE,color.rows=FALSE,similarity.matrix="correlation",deepsplit=2,minclustsize=10,mergeCutHeight=0.1,num_nodes=2)
{
    
    h73<-get_hca_child(X=X,Y=Y,feature_table_file=feature_table_file,parentoutput_dir=parentoutput_dir,class_labels_file=class_labels_file,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=is.data.znorm,analysismode=analysismode,
    sample.col.opt=sample.col.opt,plots.width=plots.width,plots.height=plots.height,plots.res=plots.res, plots.type=plots.type, alphacol=alphacol, hca_type=hca_type,newdevice=newdevice,input.type=input.type,mainlab=mainlab,cexRow=cexRow, cexCol=cexCol,plot.bycluster=plot.bycluster,color.rows=color.rows,similarity.matrix,deepsplit,minclustsize,mergeCutHeight,num_nodes)
    return(h73)
}

get_hca_child<-function(feature_table_file,parentoutput_dir,class_labels_file,X=NA,Y=NA,heatmap.col.opt="RdBu",cor.method="spearman",is.data.znorm=FALSE,analysismode="classification",
sample.col.opt="rainbow",plots.width=8,plots.height=8,plots.res=600, plots.type="cairo", alphacol=0.3, hca_type,newdevice=FALSE,input.type="intensity",mainlab="",cexRow=1, cexCol=1,plot.bycluster=FALSE,color.rows=FALSE,similarity.matrix="correlation",deepsplit=2,minclustsize=10,mergeCutHeight=0.1,num_nodes=2)
{
    
    #print(dim(X))
    
    if(is.na(X)==TRUE){
    data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE)
    }else{
        data_matrix<-X
        rm(X)
        
    }
    
    dir.create(parentoutput_dir)
    setwd(parentoutput_dir)
    col_metabs=color.rows
    
    data_m<-data_matrix[,-c(1:2)]
    
    data_m<-as.matrix(data_m)
    
    col_samples<-TRUE
    
    suppressWarnings(dir.create("Tables"))
    
    if(is.na(Y)==TRUE){
    if(is.na(class_labels_file)==FALSE){
        classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
    }else{
        
        classlabels<-rep("classA",dim(data_m)[2])
        classlabels<-cbind(classlabels,classlabels)
        
        col_samples<-FALSE
        
    }
    }else{
        classlabels<-Y
        
    }
    
    
    #patientcolors<-rep("green",dim(data_m)[2])
    
    class_labels_levels<-levels(as.factor(classlabels[,2]))
    ordered_labels<-classlabels[,2]
    
    #class_label_alphabets<-c("A","B","C","D","E","F","G","H","I","J","K","L","M")
    class_label_alphabets<-class_labels_levels #paste("C",1:length(class_labels_levels),sep="")
    
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
                        if(sample.col.opt=="colorblind"){
                            #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                            # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                            
                            if(length(class_labels_levels)<9){
                                
                                col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7", "grey57")
                                
                            }else{
                                
                                # col_vec<-colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels))
                                col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                                
                            }
                            
                            
                        }else{
                            
                            print("CHECKING THIS")
                            print(sample.col.opt)
                            check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                            
                            if(length(check_brewer)>0){
                                
                                sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                                col_vec <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(class_labels_levels))
                                
                            }else{
                                
                                if(sample.col.opt=="journal"){
                                    
                                    col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                    "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                    "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                    
                                    "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                    "#E64B3519","#4DBBD519","#631879E5","grey75")
                                    
                                }else{
                                    #colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                    
                                    if(length(sample.col.opt)==1){
                                        col_vec <-rep(sample.col.opt,length(class_labels_levels))
                                    }else{
                                        
                                        colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                        
                                    }
                                }
                                
                            }
                            
                        }
                    }
                    
                    
                }
                
            }
            
        }
    }
    col_samples=FALSE
    
    print("CHECKING THIS 2")
    print(sample.col.opt)
    print(col_vec)
    if(analysismode=="classification")
    {
        
     
        sampleclass<-{}
        patientcolors<-rep("green",nrow(classlabels))
        #print(classlabels)
        classlabels<-as.data.frame(classlabels)
        f<-factor(classlabels[,1])
        
        col_samples=TRUE
        for(c in 1:length(class_labels_levels)){
            
            num_samps_group_cur=length(which(ordered_labels==class_labels_levels[c]))
            
            #classlabels<-c(classlabels,rep(paste("Class",class_label_alphabets,sep=""),num_samps_group_cur))
            #,rep("ClassB",num_samps_group[[2]]),rep("ClassC",num_samps_group[[3]]))
            sampleclass<-c(sampleclass,rep(paste("Class",class_label_alphabets[c],sep=""),num_samps_group_cur))
            
            #patientcolors <-c(patientcolors,rep(col_vec[c],num_samps_group_cur))
            patientcolors[which(ordered_labels==class_labels_levels[c])]<-col_vec[c]
            
        }
        
        
        
    }
    
    if(heatmap.col.opt=="RdBu"){
        
        heatmap.col.opt="redblue"
    }
    
    heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
    heatmap_cols<-rev(heatmap_cols)
    
    if(heatmap.col.opt=="topo"){
        heatmap_cols<-topo.colors(256)
         heatmap_cols<-rev(heatmap_cols)
    }else{
        if(heatmap.col.opt=="heat"){
            heatmap_cols<-heat.colors(256)
             heatmap_cols<-rev(heatmap_cols)
        }else{
            
            if(heatmap.col.opt=="yellowblue"){
                
                  heatmap_cols<-colorRampPalette(c("yellow","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                #heatmap_cols<-blue2yellow(256) #colorRampPalette(c("yellow","blue"))(256)
               heatmap_cols<-rev(heatmap_cols)
            }else{
                
                if(heatmap.col.opt=="redblue"){
                    
                    heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
                    heatmap_cols<-rev(heatmap_cols)
                }else{
                    
                    #my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                    if(heatmap.col.opt=="redyellowgreen"){
                        
                        heatmap_cols <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                        heatmap_cols<-rev(heatmap_cols)
                    }else{
                        if(heatmap.col.opt=="yellowwhiteblue"){
                            
                            heatmap_cols<-colorRampPalette(c("yellow2","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                            heatmap_cols<-rev(heatmap_cols)
                        }else{
                            
                            if(heatmap.col.opt=="redwhiteblue"){
                                
                                heatmap_cols<-colorRampPalette(c("red","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                                heatmap_cols<-rev(heatmap_cols)
                            }else{
                                
                                
                                
                                    heatmap_cols <- colorRampPalette(brewer.pal(10, heatmap.col.opt))(256)
                                    heatmap_cols<-rev(heatmap_cols)
                                
                            }

                        }
                        
                    }
                    
                }
                
            }
        }
        
    }
    
    # print(classlabels)
    #print(patientcolors)
    col_vec2 <- topo.colors(255, alpha=alphacol)
    
    if(input.type=="intensity"){
    
    
    if(similarity.matrix=="TOM"){
     
    
        data_m<-log2(data_m+1)
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        simmat<-WGCNA::cor(t(data_m),nThreads=num_nodes,method=cor.method,use = 'p')
        
        sft = try(pickSoftThreshold.fromSimilarity(similarity=simmat, powerVector = powers, verbose = 0),silent=TRUE)
        #power_val=sft$powerEstimate
        
        if(is(sft,"try-error")){
            power_val=6
        }else{
            power_val=sft$powerEstimate
        }
        
        if(is.na(power_val)==TRUE){
            power_val=6
        }
        
        
        ADJdataOne<-adjacency.fromSimilarity(similarity=simmat,power=power_val)
        
        
        dissTOMCormat=TOMdist(ADJdataOne) #(1-global_cor)
        hr = flashClust(as.dist(dissTOMCormat),method="complete");
        
        pdf("metabplot.pdf")
        plot(hr,labels=F,main="Dendrogram")
        dev.off()

         mycl_metabs <-cutreeDynamic(hr,distM= dissTOMCormat,deepSplit=deepsplit, minClusterSize=minclustsize, pamRespectsDendro = FALSE, pamStage=TRUE,verbose=0)
         
         m2=try(mergeCloseModules(t(data_m),colors=mycl_metabs,cutHeight=mergeCutheight),silent=TRUE)

        if(is(m2,"try-error")){
            
            mod_list<-mycl_metabs
        }else{
            
            mod_list<-as.numeric(m2$colors)
        }
        
        mycl_metabs<-mod_list
        
        
        ###samples
        simmat<-WGCNA::cor((data_m),nThreads=num_nodes,method=cor.method,use = 'p')
        
        sft = try(pickSoftThreshold.fromSimilarity(similarity=simmat, powerVector = powers, verbose = 0),silent=TRUE)
        #power_val=sft$powerEstimate
        
        if(is(sft,"try-error")){
            power_val=6
        }else{
            power_val=sft$powerEstimate
        }
        
        if(is.na(power_val)==TRUE){
            power_val=6
        }
        
        
        ADJdataOne<-adjacency.fromSimilarity(similarity=simmat,power=power_val)
        
        
        dissTOMCormat=TOMdist(ADJdataOne) #(1-global_cor)
        hc = flashClust(as.dist(dissTOMCormat),method="complete");
        
        pdf("sampleplot.pdf")
        plot(hr,labels=F,main="Dendrogram")
        dev.off()

        
        mycl_samples <-cutreeDynamic(hc,distM= dissTOMCormat,deepSplit=deepsplit, minClusterSize=minclustsize, pamRespectsDendro = FALSE, pamStage=TRUE,verbose=0)
        
        m3=try(mergeCloseModules((data_m),colors=mycl_samples,cutHeight=mergeCutheight),silent=TRUE)
        
        if(is(m3,"try-error")){
            
            mod_list2<-mycl_samples
        }else{
            
            mod_list2<-as.numeric(m3$colors)
        }
        
        mycl_samples<-mod_list2
        
         
     
    }else{
        hc <- try(hclust(as.dist(1-WGCNA::cor(data_m,method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #samples
        hr <- try(hclust(as.dist(1-WGCNA::cor(t(data_m),method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #metabolites
   
        mycl_samples <- cutree(hc, h=max(hc$height)/2)
        mycl_metabs <- cutree(hr, h=max(hr$height)/2)
    }
        
        
    }else{
        if(input.type=="correlation"){
            #hc <- try(hclust(as.dist(1-data_m)),silent=TRUE) #samples
            
             hc <- try(hclust(as.dist(1-data_m)),silent=TRUE) #samples
             hr <- try(hclust(as.dist(1-data_m)),silent=TRUE) #metabolites
             
             mycl_samples <- cutree(hc, h=max(hc$height)/2)
             mycl_metabs <- cutree(hr, h=max(hr$height)/2)
        }
        
    }
    
    
   
   mainlab1<-"" #paste("HCA using ",mainlab," selected features",sep="")
    
    if(is(hr,"try-error") || is(hc,"try-error")){
								
                                print("Hierarchical clustering can not be performed. ")
    }else{
        heatmap_file<-paste("heatmap.jpeg",sep="")
        
        if(newdevice==TRUE){
        png(heatmap_file,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        #png("test_png.png",res=600,width=8,height=8,units="in",type="cairo")
        }
        
        if(hca_type=="two-way"){
            
          
            
          
            if(col_metabs==TRUE){
                
                col_vec2 <- rainbow(length(unique(mycl_metabs)), alpha=alphacol)
                
                
                if(sample.col.opt=="topo"){
                    
                    col_vec2 <- topo.colors(length(unique(mycl_metabs)), alpha=alphacol)
                }else{
                    if(sample.col.opt=="heat"){
                  
                        col_vec2 <- heat.colors(length(unique(mycl_metabs)), alpha=alphacol)
                    }else{
                        if(sample.col.opt=="rainbow"){
                            col_vec2<-rainbow(length(unique(mycl_metabs)), start = 0, end = alphacol)
                            
                            
                        }else{
                            
                            if(sample.col.opt=="terrain"){
                                
                                col_vec2 <- cm.colors(length(unique(mycl_metabs)), alpha=alphacol)
                            }else{
                                
                                if(sample.col.opt=="colorblind"){
                                    #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                                    # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                                    
                                    if(length(unique(mycl_metabs))<9){
                                        
                                        col_vec2 <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7", "grey57")
                                        
                                    }else{
                                        
                                        #col_vec2<-colorRampPalette(brewer.pal(10, "RdBu"))(length(unique(mycl_metabs)))
                                        col_vec2<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                        "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                                        
                                    }
                                    
                                    
                                }else{
                                    
                                    check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                                    
                                    if(length(check_brewer)>0){
                                        
                                        sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                                        col_vec2 <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(unique(mycl_metabs)))
                                        
                                    }else{
                                        
                                        if(sample.col.opt=="journal"){
                                            
                                            col_vec2<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                            "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                            "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                            
                                            "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                            "#E64B3519","#4DBBD519","#631879E5","grey75")
                                            
                                        }else{
                                            col_vec2 <-sample.col.opt
                                        }
                                        
                                    }
                                    
                                }
                            }
                            
                            
                        }
                        
                    }
                    
                }
                
                rowcolors=col_vec2[as.numeric(mycl_metabs)]
            }else{
                rowcolors=NA #rep("",length(mycl_metabs))
            }
            
            
            if(plot.bycluster==TRUE){
                
                #col_vec2 <- rainbow(length(unique(mycl_samples)), alpha=alphacol)
                if(sample.col.opt=="topo"){
                    
                    col_vec2 <- topo.colors(length(unique(mycl_samples)), alpha=alphacol)
                }else{
                    if(sample.col.opt=="heat"){
                        
                        col_vec2 <- heat.colors(length(unique(mycl_samples)), alpha=alphacol)
                    }else{
                        if(sample.col.opt=="rainbow"){
                            col_vec2<-rainbow(length(unique(mycl_samples)), start = 0, end = alphacol)
                            
                            
                        }else{
                            
                            if(sample.col.opt=="terrain"){
                                
                                col_vec2 <- cm.colors(length(unique(mycl_samples)), alpha=alphacol)
                            }else{
                                
                                if(sample.col.opt=="colorblind"){
                                    #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                                    # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                                    
                                    if(length(unique(mycl_samples))<9){
                                        
                                        col_vec2 <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7", "grey57")
                                        
                                    }else{
                                        
                                        # col_vec2<-colorRampPalette(brewer.pal(10, "RdBu"))(length(unique(mycl_samples)))
                                        col_vec2<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                        "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                                        
                                    }
                                    
                                    
                                }else{
                                    
                                    check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                                    
                                    if(length(check_brewer)>0){
                                        
                                        sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                                        
                                        col_vec2 <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(unique(mycl_samples)))
                                        
                                    }else{
                                        
                                        if(sample.col.opt=="journal"){
                                            
                                            col_vec2<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                            "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                            "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                            
                                            "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                            "#E64B3519","#4DBBD519","#631879E5","grey75")
                                            
                                        }else{
                                            col_vec2 <-sample.col.opt
                                        }
                                        
                                    }
                                    
                                }
                            }
                            
                            
                        }
                        
                    }
                    
                }
                    print(table(mycl_samples))
                    patientcolors<-col_vec2[as.numeric(mycl_samples)]
                
            }
            
            # #save(list=ls(),file="debug.Rda")
            if(col_samples==FALSE){
                if(is.data.znorm==FALSE){
                    
                    if(is.na(rowcolors)==FALSE){
                        
                        w <- 0.1
                        par(omd=c(0, 1-w, 0, 1))
                        
                     
                        
                    h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,RowSideColors=rowcolors,labRow = FALSE, labCol = FALSE)
                    
                    print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                    
                    
                    
                    }else{
                         w <- 0.1
                        par(omd=c(0, 1-w, 0, 1))
                        
                        
                        
                        h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,labRow = FALSE, labCol = FALSE)
                  
                  print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                  
                    }
                    
                }else{
                    
                    if(is.na(rowcolors)==FALSE){
                         w <- 0.1
                        par(omd=c(0, 1-w, 0, 1))
                        
                        
                    h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,RowSideColors=rowcolors,labRow = FALSE, labCol = FALSE)
                   
                   print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                   
                   
                    }else{
                        w <- 0.1
                        par(omd=c(0, 1-w, 0, 1))
                        
                        
                        h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,labRow = FALSE, labCol = FALSE)
                        
                        print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                        
                  
                    }
                }
                
            }else{

            
            if(is.data.znorm==FALSE){
                
                if(is.na(rowcolors)==FALSE){
                   w <- 0.1
                    par(omd=c(0, 1-w, 0, 1))
                    
                h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol, xlab="Samples",ylab="mzfeatures", main=mainlab1, ColSideColors=patientcolors,RowSideColors=rowcolors,labRow = FALSE, labCol = FALSE)
                
                print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                
                }else{
                    
                      w <- 0.1
                    par(omd=c(0, 1-w, 0, 1))
                    
                     h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol, xlab="Samples",ylab="mzfeatures", main=mainlab1, ColSideColors=patientcolors,labRow = FALSE, labCol = FALSE)
                     
                     print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                     
                }
                
            }else{
                
                if(is.na(rowcolors)==FALSE){
                    
                     w <- 0.1
                    par(omd=c(0, 1-w, 0, 1))
                    
                h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1, ColSideColors=patientcolors,RowSideColors=rowcolors,labRow = FALSE, labCol = FALSE)
                
                print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                
                }else{
                    
                     w <- 0.1
                    par(omd=c(0, 1-w, 0, 1))
                    
                    h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1, ColSideColors=patientcolors,labRow = FALSE, labCol = FALSE)
                    
                    print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                    
                }
            }
            }
   

            if(newdevice==TRUE){
            try(dev.off(),silent=TRUE)
            }
            
            
           
            ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],data_matrix[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])
        }
        else{
            if(hca_type=="one-way"){
            hc<-seq(1,dim(data_m)[2])
            
            
            # mycl_metabs <- cutree(hr, h=max(hr$height)/2)
            
            if(col_metabs==TRUE){
                
                #col_vec2 <- rainbow(length(unique(mycl_metabs)), alpha=alphacol)
                
                if(sample.col.opt=="topo"){
                    
                    col_vec2 <- topo.colors(length(unique(mycl_metabs)), alpha=alphacol)
                }else{
                    if(sample.col.opt=="heat"){
                        
                        col_vec2 <- heat.colors(length(unique(mycl_metabs)), alpha=alphacol)
                    }else{
                        if(sample.col.opt=="rainbow"){
                            col_vec2<-rainbow(length(unique(mycl_metabs)), start = 0, end = alphacol)
                            
                            
                        }else{
                            
                            if(sample.col.opt=="terrain"){
                                
                                col_vec2 <- cm.colors(length(unique(mycl_metabs)), alpha=alphacol)
                            }else{
                                
                                if(sample.col.opt=="colorblind"){
                                    #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                                    # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                                    
                                    if(length(unique(mycl_metabs))<9){
                                        
                                        col_vec2 <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7", "grey57")
                                        
                                    }else{
                                        
                                        #col_vec2<-colorRampPalette(brewer.pal(10, "RdBu"))(length(unique(mycl_metabs)))
                                        col_vec2<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                        "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                                        
                                    }
                                    
                                    
                                }else{
                                    
                                    check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                                    
                                    if(length(check_brewer)>0){
                                        
                                        sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                                        col_vec2 <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(unique(mycl_metabs)))
                                        
                                    }else{
                                        
                                        if(sample.col.opt=="journal"){
                                            
                                            col_vec2<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                            "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                            "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                            
                                            "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                            "#E64B3519","#4DBBD519","#631879E5","grey75")
                                            
                                        }else{
                                            col_vec2 <-sample.col.opt
                                        }
                                        
                                    }
                                    
                                }
                            }
                            
                            
                        }
                        
                    }
                    
                }
                
                
                rowcolors=col_vec2[as.numeric(mycl_metabs)]
            }else{
                rowcolors=NA
            }
            
           
            
            if(col_samples==FALSE){
                if(is.data.znorm==FALSE){
                    
                    if(is.na(rowcolors)==FALSE){
                        
                        w <- 0.1
                        par(omd=c(0, 1-w, 0, 1))
                        
                    h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,RowSideColors=rowcolors,labRow = FALSE, labCol = FALSE)
                    
                    print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                    
                    }else{
                        
                          w <- 0.1
                        par(omd=c(0, 1-w, 0, 1))
                        
                        h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,labRow = FALSE, labCol = FALSE)
                        
                        print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                        
                        
                    }
                }else{
                    
                    if(is.na(rowcolors)==FALSE){
                        
                          w <- 0.1
                        par(omd=c(0, 1-w, 0, 1))
                        
                    h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,RowSideColors=rowcolors,labRow = FALSE, labCol = FALSE)
                    
                    print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                    
                    }else{
                        
                        w <- 0.1
                        par(omd=c(0, 1-w, 0, 1))
                        
                        h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,labRow = FALSE, labCol = FALSE)
                        
                        print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                        
                        
                    }
                }
                
            }else{
            if(is.data.znorm==FALSE){
                
                if(is.na(rowcolors)==FALSE){
                    
                     w <- 0.1
                    par(omd=c(0, 1-w, 0, 1))
                    
                h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,dendrogram = c("row"),ColSideColors=patientcolors,RowSideColors=rowcolors,labRow = FALSE, labCol = FALSE)
                
                print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                
                }else{
                     w <- 0.1
                    par(omd=c(0, 1-w, 0, 1))
                    
                    h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,dendrogram = c("row"),ColSideColors=patientcolors,labRow = FALSE, labCol = FALSE)
                    
                    print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                    
                }
            }else{
                
                if(is.na(rowcolors)==FALSE){
                     w <- 0.1
                    par(omd=c(0, 1-w, 0, 1))
                    
                h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,dendrogram = c("row"),ColSideColors=patientcolors,RowSideColors=rowcolors,labRow = FALSE, labCol = FALSE)
                
                print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                
                }else{
                     w <- 0.1
                    par(omd=c(0, 1-w, 0, 1))
                    
                    h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,dendrogram = c("row"),ColSideColors=patientcolors,labRow = FALSE, labCol = FALSE)
                    
                    print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                    
                    
                }
            }
            }
           
            if(newdevice==TRUE){
                try(dev.off(),silent=TRUE)
            }
            
            mycl_samples<-seq(1,dim(data_m)[2])
            #mycl_samples <- cutree(hc, h=max(hc$height)/2)
            #mycl_metabs <- cutree(hr, h=max(hr$height)/2)
            
            ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],data_matrix[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])
            
            }else{
                if(hca_type=="samples"){
                hr<-seq(1,dim(data_m)[1])
              
                
                
            
                if(col_samples==FALSE){
                    if(is.data.znorm==FALSE){
                        
                        w <- 0.1
                        par(omd=c(0, 1-w, 0, 1))
                        h73<-heatmap.2(data_m, Rowv=NULL, Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,labRow = FALSE, labCol = FALSE)
                        
                        print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                        
                    }else{
                        w <- 0.1
                        par(omd=c(0, 1-w, 0, 1))
                        h73<-heatmap.2(data_m, Rowv=NULL, Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol, xlab="Samples",ylab="mzfeatures", main=mainlab1,labRow = FALSE, labCol = FALSE)
                        
                        print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                        
                    }
                    
                }else{
                    if(is.data.znorm==FALSE){
                        
                          w <- 0.1
                        par(omd=c(0, 1-w, 0, 1))
                        
                        
                        h73<-heatmap.2(data_m, Rowv=NULL, Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,dendrogram = c("col"),ColSideColors=patientcolors,labRow = FALSE, labCol = FALSE)
                        
                        print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                        
                    }else{
                        w <- 0.1
                        par(omd=c(0, 1-w, 0, 1))
                        
                        
                        h73<-heatmap.2(data_m, Rowv=NULL, Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,dendrogram = c("col"),ColSideColors=patientcolors,labRow = FALSE, labCol = FALSE)
                        
                        print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                        
                    }
                }
                
             
                
                if(newdevice==TRUE){
                   try(dev.off(),silent=TRUE)
                }
                
                #mycl_samples<-seq(1,dim(data_m)[2])
                # mycl_samples <- cutree(hc, h=max(hc$height)/2)
                mycl_metabs <- seq(1,dim(data_m)[1]) #cutree(hr, h=max(hr$height)/2)
                
                ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],data_matrix[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])
                
                }else{
                    
                    
                    #only plot heatmap; no clustering
                    hr<-seq(1,dim(data_m)[1])
                    hc<-seq(1,dim(data_m)[2])
                    
                    
                    if(col_samples==FALSE){
                        if(is.data.znorm==FALSE){
                            w <- 0.1
                            par(omd=c(0, 1-w, 0, 1))
                            
                            h73<-heatmap.2(data_m, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,dendrogram = c("none"),labRow = FALSE, labCol = FALSE)
                            
                            print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                            
                        }else{
                            
                            w <- 0.1
                            par(omd=c(0, 1-w, 0, 1))
                            
                            h73<-heatmap.2(data_m, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,dendrogram = c("none"),labRow = FALSE, labCol = FALSE)
                            
                            print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                            
                        }
                        
                    }else{
                        if(is.data.znorm==FALSE){
                            
                            w <- 0.1
                            par(omd=c(0, 1-w, 0, 1))
                            
                            h73<-heatmap.2(data_m, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,ColSideColors=patientcolors,dendrogram = c("none"),labRow = FALSE, labCol = FALSE)
                        
                            print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                        
                        
                        }else{
                            
                            w <- 0.1
                            par(omd=c(0, 1-w, 0, 1))
                            
                            h73<-heatmap.2(data_m, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=cexRow, cexCol=cexCol,xlab="Samples",ylab="mzfeatures", main=mainlab1,ColSideColors=patientcolors,dendrogram = c("none"),labRow = FALSE, labCol = FALSE)
                            
                            print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                            
                        }
                    }
                    
                    #  par(xpd=TRUE)
                    #legend("bottomleft",legend=levels(ordered_labels),text.col=unique(patientcolors),pch=13,cex=0.4)
                    #par(xpd=FALSE)
                    
                    if(newdevice==TRUE){
                        try(dev.off(),silent=TRUE)
                    }
                    
                    mycl_samples<-seq(1,dim(data_m)[2])
                    #mycl_samples <- cutree(hc, h=max(hc$height)/2)
                    mycl_metabs <- seq(1,dim(data_m)[1]) #cutree(hr, h=max(hr$height)/2)
                    
                    ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],data_matrix[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])
                    
                    
                }
            }
        }
        
        cnames1<-colnames(ord_data)
        cnames1[1]<-"mz_cluster_label"
        colnames(ord_data)<-cnames1
        fname1<-paste("Tables/Clustering_based_sorted_intensity_using_",mainlab,"features.txt",sep="")
        write.table(ord_data,file=fname1,sep="\t",row.names=FALSE)
        
        fname2<-paste("Tables/Sample_clusterlabels_using_",mainlab,"features.txt",sep="")
        
        sample_clust_num<-mycl_samples[h73$colInd]
        temp1<-classlabels[h73$colInd,1]
        temp2<-classlabels
        
        #print(head(temp2))
        #temp1<-as.data.frame(temp1)
        
        #print(dim(temp1))
        match_ind<-match(temp1,temp2[,1])
        
        temp3<-temp2[match_ind,]
        
        #print(head(temp3))
        temp4<-cbind(temp1,temp3,sample_clust_num)
       
       #write.table(temp4,file="s1.txt",sep="\t",row.names=FALSE)
    #   print(head(temp1))
        print(head(temp4))
        
        rnames1<-rownames(temp4)
        #temp4<-cbind(rnames1,temp4)
        temp4<-as.data.frame(temp4)
        temp4<-temp4[,-c(1)]
        # print(temp4[,1:4])
        
        
        
        
        if(analysismode=="regression"){
            

            
            temp3<-temp4 #[,-c(1)]
            temp3<-as.data.frame(temp3)
            temp3<-apply(temp3,2,as.numeric)
            
            
            temp_vec<-as.vector(temp3[,2])
            
            
            
            names(temp_vec)<-as.character(temp4[,1])
            
            
            # if(output.device.type!="pdf"){
                
             if(newdevice==TRUE){
                
                temp_filename_1<-"Figures/Barplot_dependent_variable_ordered_by_HCA.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
            }
            
            
            #print(temp_vec)
            #tiff("Barplot_sample_cluster_ymat.tiff", width=plots.width,height=plots.height,res=plots.res, compression="lzw")
            barplot(temp_vec,col="brown",ylab="Y",cex.axis=0.5,cex.names=0.5,main="Dependent variable levels in samples; \n ordered based on hierarchical clustering")
            #dev.off()
            
            
           
            if(newdevice==TRUE){
                try(dev.off(),silent=TRUE)
            }
            
            
            
            
            
        }
        
        write.table(temp4,file=fname2,sep="\t",row.names=FALSE)
        
        
        
        fname3<-paste("Tables/Metabolite_clusterlabels_for_",mainlab,"features.txt",sep="")
        
        mycl_metabs_ord<-mycl_metabs[rev(h73$rowInd)]
        write.table(mycl_metabs_ord,file=fname3,sep="\t",row.names=TRUE)
        
    }
    return(h73)
}


data_preprocess<-function(Xmat=NA,Ymat=NA,feature_table_file,parentoutput_dir,class_labels_file,num_replicates=3,feat.filt.thresh=NA,summarize.replicates=TRUE,summary.method="mean",
all.missing.thresh=0.5,group.missing.thresh=0.7,
log2transform=TRUE,medcenter=TRUE,znormtransform=FALSE,quantile_norm=TRUE,lowess_norm=FALSE,madscaling=FALSE,missing.val=0,samplermindex=NA, rep.max.missing.thresh=0.5,summary.na.replacement="zeros",featselmethod=NA,TIC_norm=FALSE,pairedanalysis=FALSE){
    
    options(warn=-1)
    
    #read file; First row is column headers
    if(is.na(Xmat==TRUE)){
        data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE)
    }else{
        data_matrix<-Xmat
        #rm(Xmat)
    }

    #print("signal filter threshold ")
    #print(group.missing.thresh)
    
    print("missing val is")
    print(missing.val)
    
    if(is.na(all.missing.thresh)==TRUE){
        
        all.missing.thresh=(-1)
    }
    
    if(is.na(samplermindex)==FALSE){
        data_matrix<-data_matrix[,-c(samplermindex)]
    }
    
    #use only unique records
    data_matrix<-unique(data_matrix)
    
    if(is.na(missing.val)==FALSE){
        
        print("Replacing missing values with NAs.")
        data_matrix<-replace(as.matrix(data_matrix),which(data_matrix==missing.val),NA)
    }
    
    # print(data_matrix[1:10,1:5])

    #print("dim of original data matrix")
    #print(dim(data_matrix))
    data_matrix_orig<-data_matrix
    
    
    snames<-colnames(data_matrix)
    
    
    
    
    
    dir.create(parentoutput_dir,showWarnings=FALSE)
    parentoutput_dir<-paste(parentoutput_dir,"/Stage1/",sep="")
    
    dir.create(parentoutput_dir,showWarnings=FALSE)
    fheader="transformed_log2fc_threshold_"
    setwd(parentoutput_dir)
    
    data_m<-as.matrix(data_matrix[,-c(1:2)])
    
    if(is.na(Xmat)==FALSE){
        
        #   write.table(Xmat,file="organized_featuretableA.txt",sep="\t",row.names=TRUE)
       
        
    }
    
    if(is.na(Ymat)==FALSE){
        # write.table(Ymat,file="organized_classlabelsA.txt",sep="\t",row.names=FALSE)
        
    }
    
    #Step 2) Average replicates
    if(summarize.replicates==TRUE)
    {
        if(num_replicates>1)
        {
            
            data_m<-getSumreplicates(data_matrix,alignment.tool="apLCMS",numreplicates=num_replicates,numcluster=10,rep.max.missing.thresh=rep.max.missing.thresh,summary.method=summary.method,summary.na.replacement, missing.val=missing.val)
            
            #data_m<-round(data_m,3)
            
            data_m<-replace(data_m,which(is.na(data_m)==TRUE),missing.val)
            
            if(summary.method=="mean"){
                print("Replicate averaging done")
                filename<-paste("Rawdata_averaged.txt",sep="")
            }else{
                if(summary.method=="median"){
                    print("Replicate median summarization done")
                    filename<-paste("Rawdata_median_summarized.txt",sep="")
                }
                
            }
            
            data_m_prenorm<-cbind(data_matrix[,c(1:2)],data_m)
            
            write.table(data_m_prenorm, file=filename,sep="\t",row.names=FALSE)
            
            data_matrix<-cbind(data_matrix[,c(1:2)],data_m)
            #num_samps_group[[1]]=(1/num_replicates)*num_samps_group[[1]]
            #num_samps_group[[2]]=(1/num_replicates)*num_samps_group[[2]]
        }
    }
    
    
    data_matrix<-cbind(data_matrix[,c(1:2)],data_m)
    
    data_matrix_orig<-data_matrix
    data_subjects<-data_m
    
    ordered_labels={}
    
    num_samps_group<-new("list")
    
    if(is.na(class_labels_file)==FALSE)
    {
        
        print("Class labels file:")
        print(class_labels_file)
        
        data_matrix={}
        
        if(is.na(Ymat)==TRUE){
            classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
        }else{
            classlabels<-Ymat
        }
        
        #class_labels_sampnames<-classlabels[,1]
        #data_matrix_sampnames<-colnames(data_m)
        
        #classlabels<-classlabels[match(class_labels_sampnames,data_matrix_sampnames),]
        
        classlabels<-as.data.frame(classlabels)
        if(pairedanalysis==TRUE){
             classlabels<-classlabels[,-c(2)]
            
        }
        print(head(classlabels))
       
        
        cnames1<-colnames(classlabels)
        cnames1[1]<-c("SampleID")
        cnames1[2]<-c("Class")
        
        colnames(classlabels)<-cnames1 #c("SampleID","Class")
        f1<-table(classlabels$SampleID)
        
        
        classlabels<-as.data.frame(classlabels)
        classlabels<-classlabels[seq(1,dim(classlabels)[1],num_replicates),]
        #print(classlabels)
        class_labels_levels<-levels(as.factor(classlabels[,2]))
        bad_rows<-which(class_labels_levels=="")
        if(length(bad_rows)>0){
            class_labels_levels<-class_labels_levels[-bad_rows]
        }
        
        for(c in 1:length(class_labels_levels))
        {
            
            if(c>1){
                data_matrix<-cbind(data_matrix,data_subjects[,which(classlabels[,2]==class_labels_levels[c])])
            }else{
                data_matrix<-data_subjects[,which(classlabels[,2]==class_labels_levels[c])]
            }
            classlabels_index<-which(classlabels[,2]==class_labels_levels[c])
            ordered_labels<-c(ordered_labels,as.character(classlabels[classlabels_index,2]))
            num_samps_group[[c]]<-length(classlabels_index)
            
        }
        
        #colnames(data_matrix)<-as.character(ordered_labels)
        data_matrix<-cbind(data_matrix_orig[,c(1:2)],data_matrix)
        data_m<-as.matrix(data_matrix[,-c(1:2)])
        
        
    }else
    {
        if(is.na(Ymat)==TRUE)
        {
            classlabels<-rep("A",dim(data_m)[2])
            classlabels<-as.data.frame(classlabels)
            ordered_labels<-classlabels
            num_samps_group[[1]]<-dim(data_m)[2]
            class_labels_levels<-c("A")
            data_m<-as.matrix(data_matrix[,-c(1:2)])
            
        }else{
            classlabels<-Ymat
            classlabels<-as.data.frame(classlabels)
            
            if(pairedanalysis==TRUE){
                classlabels<-classlabels[,-c(2)]
                
            }
            print(head(classlabels))
            cnames1<-colnames(classlabels)
            cnames1[1]<-c("SampleID")
            cnames1[2]<-c("Class")
            
            colnames(classlabels)<-cnames1
            #colnames(classlabels)<-c("SampleID","Class")
            f1<-table(classlabels$SampleID)
            
            
            classlabels<-as.data.frame(classlabels)
            classlabels<-classlabels[seq(1,dim(classlabels)[1],num_replicates),]
            #print(classlabels)
            class_labels_levels<-levels(as.factor(classlabels[,2]))
            bad_rows<-which(class_labels_levels=="")
            if(length(bad_rows)>0){
                class_labels_levels<-class_labels_levels[-bad_rows]
            }
            
            for(c in 1:length(class_labels_levels))
            {
                #if(c>1){
                #data_matrix<-cbind(data_matrix,data_subjects[,which(classlabels[,2]==class_labels_levels[c])])
                #}else{
                #	data_matrix<-data_subjects[,which(classlabels[,2]==class_labels_levels[c])]
                #}
                classlabels_index<-which(classlabels[,2]==class_labels_levels[c])
                ordered_labels<-c(ordered_labels,as.character(classlabels[classlabels_index,2]))
                num_samps_group[[c]]<-length(classlabels_index)
                
            }
            
            #colnames(data_matrix)<-as.character(ordered_labels)
            #data_matrix<-cbind(data_matrix_orig[,c(1:2)],data_matrix)
            #data_m<-as.matrix(data_matrix[,-c(1:2)])
        }
        
        
        
    }
    
    #  #save(classlabels,file="classlabels1.Rda")
    ##save(num_samps_group,file="num_samps_group.Rda")
    ##save(ordered_labels,file="ordered_labels1.Rda")
    
    rnames_xmat<-colnames(data_matrix[,-c(1:2)])
    
    if(is.na(classlabels)==FALSE){
    
    
       #print("here")
     rnames_ymat<-as.character(classlabels[,1])
    }else{
        rnames_ymat<-rnames_xmat
        
    }
    
    for(ind1 in 1:length(rnames_xmat)){
    check_ylabel<-regexpr(rnames_ymat[ind1],pattern="^[0-9]*",perl=TRUE)
    check_xlabel<-regexpr(rnames_xmat[ind1],pattern="^X[0-9]*",perl=TRUE)
    
    if(length(check_ylabel)>0 && length(check_xlabel)>0){
        if(attr(check_ylabel,"match.length")>0 && attr(check_xlabel,"match.length")>0){
            
            rnames_ymat<-paste("X",rnames_ymat,sep="") #gsub(rnames_ymat,pattern="\\.[0-9]*",replacement="")
            
            
        }
    }
    
    }
    
    match_names<-match(rnames_xmat,rnames_ymat)
    
    bad_colnames<-length(which(is.na(match_names)==TRUE))
    
    classlabels<-classlabels[match_names,]

    #Step 3a) Remove features if signal is not detected in at least x% of all samples
    ##################################################################################
    metab_zeros={}
    data_clean<-{}
    clean_metabs<-{}
   

    if(is.na(all.missing.thresh)==FALSE)
    {
        
        total_sigs<-apply(data_m,1,function(x){
            if(is.na(missing.val)==FALSE){return(length(which(x>missing.val)))
            }else{
                return(length(which(is.na(x)==FALSE)))
            }})
        
        
        
        total_sig_thresh<-dim(data_m)[2]*all.missing.thresh
        
        total_good_metabs<-which(total_sigs>total_sig_thresh)
        
    }
    
    #remove bad features based on all missing values criteria
    if(length(total_good_metabs)>0){
        data_m<-data_m[total_good_metabs,]
        data_matrix<-data_matrix[total_good_metabs,]
        #print(paste("Dimension of data matrix after overall ",all.missing.thresh,"% signal threshold filtering",sep=""))
        print(paste("Dimension of data matrix after using overall ",100*all.missing.thresh, "% signal criteria for filtering:"),sep="")
        print(dim(data_matrix))
    }else{
        stop(paste("None of the metabolites have signal in ",all.missing.thresh*100, "% of samples",sep=""))
    }
    
    
    #Step 3b) Find features for which the signal is not detected in at least x% of samples in either of the groups
    
    
    data_m<-data_matrix[,-c(1:2)]
    
    if(is.na(group.missing.thresh)==FALSE)
    {
        
        if(length(class_labels_levels)==0)
        {
            
            sig_thresh_groupA<-group.missing.thresh*num_samps_group[[1]]
            sig_thresh_groupB<-group.missing.thresh*num_samps_group[[2]]
            
            for(metab_num in 1:dim(data_matrix)[1])
            {
                #print(missing.val)
                if(is.na(missing.val)==FALSE){
                    
                    num_sigsA<-length(which(data_m[metab_num,1:num_samps_group[[1]]]>missing.val))
                    
                    
                    num_sigsB<-length(which(data_m[metab_num,(num_samps_group[[1]]+1):(num_samps_group[[1]]+num_samps_group[[2]])]>missing.val))
                    
                }else{
                    
                    num_sigsA<-length(which(is.na(data_m[metab_num,1:num_samps_group[[1]]])==FALSE))
                    num_sigsB<-length(which(is.na(data_m[metab_num,(num_samps_group[[1]]+1):(num_samps_group[[1]]+num_samps_group[[2]])])==FALSE))
                    
                }
                
                if((num_sigsA>=sig_thresh_groupA) || (num_sigsB>=sig_thresh_groupB))
                {
                    clean_metabs<-c(clean_metabs,metab_num)
                }
                
            }
            
            #print(length(clean_metabs))
        }else{
            if(length(class_labels_levels)==0){
                
                
                sig_thresh_groupA<-group.missing.thresh*num_samps_group[[1]]
                sig_thresh_groupB<-group.missing.thresh*num_samps_group[[2]]
                sig_thresh_groupC<-group.missing.thresh*num_samps_group[[3]]
                
                for(metab_num in 1:dim(data_matrix)[1])
                {
                    if(is.na(missing.val)==FALSE){
                        num_sigsA<-length(which(data_m[metab_num,1:num_samps_group[[1]]]>missing.val))
                        num_sigsB<-length(which(data_m[metab_num,(num_samps_group[[1]]+1):(num_samps_group[[1]]+num_samps_group[[2]])]>missing.val))
                        num_sigsC<-length(which(data_m[metab_num,(num_samps_group[[1]]+num_samps_group[[2]]+1):(num_samps_group[[1]]+num_samps_group[[2]]+num_samps_group[[3]])]>missing.val))
                    }else{
                        
                        num_sigsA<-length(which(is.na(data_m[metab_num,1:num_samps_group[[1]]])==FALSE))
                        num_sigsB<-length(which(is.na(data_m[metab_num,(num_samps_group[[1]]+1):(num_samps_group[[1]]+num_samps_group[[2]])])==FALSE))
                        num_sigsC<-length(which(is.na(data_m[metab_num,(num_samps_group[[1]]+num_samps_group[[2]]+1):(num_samps_group[[1]]+num_samps_group[[2]]+num_samps_group[[3]])])==FALSE))
                        
                    }
                    
                    if((num_sigsA>=sig_thresh_groupA) || (num_sigsB>=sig_thresh_groupB) || (num_sigsC>=sig_thresh_groupC))
                    {
                        clean_metabs<-c(clean_metabs,metab_num)
                    }
                    
                }
            }else{
                if(length(class_labels_levels)>1){
                    

                   

                    for(metab_num in 1:dim(data_m)[1])
                    {
                        for(c in 1:length(class_labels_levels)){

				classlabels_index<-which(classlabels[,2]==class_labels_levels[c])
				templabels<-classlabels[,2]
                            if(is.na(missing.val)==FALSE){
                                
                                num_cursig<-length(which(data_m[metab_num,which(templabels==class_labels_levels[c])]>missing.val))
                                
                            }else{
                                
                                num_cursig<-length(which(is.na(data_m[metab_num,which(templabels==class_labels_levels[c])])==FALSE))
                                
                                
                            }
                            
                            sig_thresh_cur<-length(which(templabels==class_labels_levels[c]))*group.missing.thresh
                            if(num_cursig>=sig_thresh_cur)
                            {
                                clean_metabs<-c(clean_metabs,metab_num)
                                break   #for(i in 1:4){if(i==3){break}else{print(i)}}
                                
                            }
                            
                        }
                    }
                }
                else{
                    
                    
                    
                    if(length(class_labels_levels)==1){
                        num_samps_group[[1]]<-num_samps_group[[1]]
                        
                        
                        sig_thresh_groupA<-group.missing.thresh*num_samps_group[[1]]
                        
                        
                        for(metab_num in 1:dim(data_matrix)[1])
                        {
                            if(is.na(missing.val)==FALSE){
                                num_sigsA<-length(which(data_m[metab_num,1:num_samps_group[[1]]]>missing.val))
                                
                            }else{
                                
                                num_sigsA<-length(which(is.na(data_m[metab_num,1:num_samps_group[[1]]])==FALSE))
                            }
                            
                            if((num_sigsA>=sig_thresh_groupA) )
                            {
                                clean_metabs<-c(clean_metabs,metab_num)
                            }
                            
                        }
                    }
                }
            }
        }
        
        
        
    }
    ####################################################################################
    
    #Step 4) Replace missing values
    if(summarize.replicates==TRUE)
    {
        
        {
            
            if(is.na(missing.val)==FALSE){
                
                print("Replacing missing values with NAs.")
                data_m<-replace(as.matrix(data_m),which(data_m==missing.val),NA)
            }

            
            if(summary.na.replacement=="zeros"){
                data_m<-replace(data_m,which(is.na(data_m)==TRUE),0)
            }else{
                if(summary.na.replacement=="halfsamplemin"){
                    data_m<-apply(data_m,2,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-(min(x,na.rm=TRUE)/2)}; return(x)})
                }else{
                    
                    if(summary.na.replacement=="halfdatamin"){
                        
                        
                        min_val<-min(data_m,na.rm=TRUE)*0.5
                        data_m<-replace(data_m,which(is.na(data_m)==TRUE),min_val)
                        
                        #data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-min(x,na.rm=TRUE)/2}; return(x)})
                    }else{
                        if(summary.na.replacement=="halffeaturemin"){
                            data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-(min(x,na.rm=TRUE)/2)}; return(x)})
                            data_m<-t(data_m)
                        }else{
                            
                            
                            if(summary.na.replacement=="bpca"){
                                library(pcaMethods)
                                
                                print(Sys.time())
                                pc1 <- pcaMethods::pca(t(data_m), method="bpca", nPcs=3,scale="uv")
                               
                                data_m<-pcaMethods::completeObs(pc1)
                                print(Sys.time())
                                try(detach("package:pcaMethods",unload=TRUE),silent=TRUE)
                                
                                data_m<-t(data_m)
                                
                            }
                            
                            
                            
                        }
                    }
                }
                
                
            }
        }
    }else
    {
        data_m<-data_matrix[,-c(1:2)]
        
        if(is.na(missing.val)==FALSE){
            
            print("Replacing missing values with NAs.")
            data_m<-replace(as.matrix(data_m),which(data_m==missing.val),NA)
        }
        
        if(summary.na.replacement=="zeros"){
            data_m<-replace(data_m,which(is.na(data_m)==TRUE),0)
        }else{
            if(summary.na.replacement=="halfsamplemin"){
                data_m<-apply(data_m,2,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-(min(x,na.rm=TRUE)/2)}; return(x)})
            }else{
                
                if(summary.na.replacement=="halfdatamin"){
                    
                    
                    min_val<-min(data_m,na.rm=TRUE)*0.5
                    data_m<-replace(data_m,which(is.na(data_m)==TRUE),min_val)
                    
                    #data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-min(x,na.rm=TRUE)/2}; return(x)})
                }else{
                    if(summary.na.replacement=="halffeaturemin"){
                        data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-(min(x,na.rm=TRUE)/2)}; return(x)})
                        data_m<-t(data_m)
                    }else{
                        
                        if(summary.na.replacement=="bpca"){
                            
                            pc1 <- pca(t(data_m), method="bpca", nPcs=2)
                            
                            data_m<-completeObs(pc1)
                            data_m<-t(data_m)
                        }
                        
                    }
                }
            }
            
            
        }
        
    }
    
    
    
    #group-wise missing values
    if(length(clean_metabs)>0)
    {
        data_m<-data_m[clean_metabs,]
        data_matrix<-data_matrix[clean_metabs,]
        
        print(paste("Dimension of data matrix after using group-wise (Factor 1) ",100*group.missing.thresh, "% signal criteria for filtering:"),sep="")
        print(dim(data_matrix))
        
    }
    
    data_m<-as.data.frame(data_m)
    data_matrix<-as.data.frame(data_matrix)
    
    ##save(data_matrix,file="data_matrix.Rda")
    # #save(data_m,file="data_m.Rda")
    
    
    #    data_matrix<-cbind(data_matrix[,c(1:2)],data_m)
    #write.table(data_matrix,file="pretransformation.txt",sep="\t",row.names=FALSE)
    ####################################################################
    #Step 4) Data transformation and normalization
    
    data_m_prescaling<-data_m
    
    
    

    
    if(TIC_norm==TRUE){
        
        ##save(data_m,file="data_m_raw.Rda")
        sum_int<-apply(data_m,2,function(x){sum(x,na.rm=TRUE)})
        
        data_m<-sweep(data_m,2,sum_int,'/')
        
       
  
        #if(log2transform==TRUE){
            
        #data_m<-10^3*(data_m)
        #}else
        
        {
         data_m<-10^10*(data_m)
        }
        
         
    }
    
      if(log2transform==TRUE)
    {
        data_m<-log2(data_m+1)
        
        # print("log scale")
        #print(head(data_m))
    }
    
    # if(RUVrand_norm==TRUE){
        
        #    temp1=t(data_m)
        #data_m<-NormalizeRUVRand(Y=temp1,ctl,k=NULL,lambda=NULL,plotk=TRUE)
        #}
	
   if(quantile_norm==TRUE)
    {
        data_m<-normalizeQuantiles(data_m)
        #print("quant norm")
        #print(head(data_m))
        
        
    }
    
    if(lowess_norm==TRUE)
    {
        data_m<-normalizeCyclicLoess(data_m)
        #print("lowess")
    }
    
  
    
    if(medcenter==TRUE)
    {
        colmedians=apply(data_m,1,function(x){median(x,na.rm=TRUE)})
        data_m=sweep(data_m,1,colmedians)
        
        
    }
    if(znormtransform==TRUE)
    {
        data_m<-scale(t(data_m))
        data_m<-t(data_m)
    }
    
    
    if(madscaling==TRUE)
    {
        colmedians=apply(data_m,2,function(x){median(x,na.rm=TRUE)})
        
        Y=sweep(data_m,2,colmedians)
        mad<-apply(abs(Y),2,function(x){median(x,na.rm=TRUE)})
        const<-prod(mad)^(1/length(mad))
        scale.normalized<-sweep(data_m,2,const/mad,"*")
        data_m<-scale.normalized
    }
    
    
  
    data_matrix<-cbind(data_matrix[,c(1:2)],data_m)
    
    #print(dim(data_matrix))
    #print(dim(data_m))
    
    data_m<-as.data.frame(data_m)
    
    num_rows<-dim(data_m)[1]
    num_columns<-dim(data_m)[2]
    
    #print("num rows is ")
    #print(num_rows)
    #for apLCMS:
    rnames<-paste("mzid_",seq(1,num_rows),sep="")
    rownames(data_m)=rnames
    
    mzid_mzrt<-data_matrix[,c(1:2)]
    colnames(mzid_mzrt)<-c("mz","time")
    rownames(mzid_mzrt)=rnames
    write.table(mzid_mzrt, file="mzid_mzrt.txt",sep="\t",row.names=FALSE)
    
    
    
    filename<-paste("ordered_classlabels_file.txt",sep="")
    write.table(classlabels, file=filename,sep="\t",row.names=FALSE)
    
    filename<-paste("Normalized_sigthreshfilt_averaged_data.txt",sep="")
    data_matrix<-cbind(data_matrix[,c(1:2)],data_m)
    write.table(data_matrix, file=filename,sep="\t",row.names=FALSE)
    data_matrix_prescaling<-cbind(data_matrix[,c(1:2)],data_m_prescaling)
    return(list(data_matrix_afternorm_scaling=data_matrix,data_matrix_prescaling=data_matrix_prescaling,classlabels=classlabels))
    #return(data_matrix)
}




replace_outliers<-function(cdata,replace.by.NA=FALSE){
data_sum<-summary(cdata,na.rm=TRUE)
iqr_value<-data_sum[5]-data_sum[2]
upper_limit<-data_sum[5]+1.5*iqr_value
lower_limit<-data_sum[2]-(1.5*iqr_value)


if(is.na(iqr_value)==FALSE){
if(iqr_value>0){
if(length(which(cdata<lower_limit)==TRUE)>0){
    if(replace.by.NA==TRUE){
        cdata[which(cdata<lower_limit)]<-NA
    }else{
        cdata[which(cdata<lower_limit)]<-min(cdata[which(cdata>lower_limit & cdata<upper_limit)],na.rm=TRUE)
        
    }
}

if(length(which(cdata>upper_limit)==TRUE)>0){
    
        if(replace.by.NA==TRUE){
            cdata[which(cdata>upper_limit)]<-NA
        }else{
    
            cdata[which(cdata>upper_limit)]<-max(cdata[which(cdata>lower_limit & cdata<upper_limit)],na.rm=TRUE)
        }
    }
}
}

return(cdata)
}
 
diffexp<-function(Xmat=NA,Ymat=NA,feature_table_file,parentoutput_dir,class_labels_file,num_replicates=3,summarize.replicates=TRUE,summary.method="mean",summary.na.replacement="zeros",missing.val=0,rep.max.missing.thresh=0.3,
 all.missing.thresh=0.1,group.missing.thresh=0.7,input.intensity.scale="raw",
log2transform=TRUE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=TRUE,lowess_norm=FALSE,madscaling=FALSE,TIC_norm=FALSE,rsd.filt.list=1,
pairedanalysis=FALSE,featselmethod=c("limma","pls"),fdrthresh=0.05,fdrmethod="BH",cor.method="spearman",networktype="complete",abs.cor.thresh=0.4,cor.fdrthresh=0.05,kfold=10,
pred.eval.method="BER",globalcor=TRUE,
target.metab.file=NA,target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=100, 
numtrees=20000,analysismode="classification",net_node_colors=c("green","red"), net_legend=TRUE,
svm_kernel="radial",heatmap.col.opt="redblue",manhattanplot.col.opt=c("darkblue","red3"),boxplot.col.opt=c("white"),barplot.col.opt=c("journal"),sample.col.opt="journal",lineplot.col.opt="black",hca_type="two-way",pls_vip_thresh=2,num_nodes=2,max_varsel=100,
pls_ncomp=5,pca.stage2.eval=FALSE,scoreplot_legend=TRUE,pca.global.eval=TRUE,rocfeatlist=seq(2,6,1),rocfeatincrement=TRUE,rocclassifier="svm",foldchangethresh=1,wgcnarsdthresh=20,WGCNAmodules=FALSE,
optselect=TRUE,max_comp_sel=1,saveRda=FALSE,legendlocation="topleft",pca.cex.val=4,
pca.ellipse=FALSE,ellipse.conf.level=0.95,pls.permut.count=1000,svm.acc.tolerance=5,limmadecideTests=TRUE,pls.vip.selection="max",globalclustering=FALSE,plots.res=600,plots.width=8,plots.height=8,plots.type="cairo",
output.device.type="pdf",pvalue.thresh=0.05,individualsampleplot.col.opt="journal",pamr.threshold.select.max=FALSE,aggregation.method="RankAggreg",aggregation.max.iter=1000,mars.gcv.thresh=1,error.bar=TRUE,cex.plots=0.7,modeltype="RI",barplot.xaxis="Factor2",lineplot.lty.option=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),...)
{
    
    time_start<-Sys.time()
    


    print("**")
    print("**")
    
    #print(paste("g is ",group.missing.thresh,sep=""))
    
    degree_rank_method="diffK"
    feat.filt.thresh=NA
    feat_weight=1
    samplermindex=NA
    pcacenter=TRUE
    pcascale=TRUE
    alphacol=0.3
  
    
    # print(is(group.missing.thresh<0.8))
    
    if(is.na(group.missing.thresh)==FALSE){
        if(group.missing.thresh<0.8){
            
            
            print("********************************************************************************")
            print("***** WARNING: group.missing.thresh is set to below 0.8. This can lead to false significance for class or group-wise comparisons.****")
            print("********************************************************************************")
            print("**")
            print("**")
            
        }
    }
    options(warn=-1)
    
    if(input.intensity.scale=="raw" || input.intensity.scale=="log2"){
        
        print("##################################################################################")
        print("Note 1: The order of samples should be same in the feature table and classlabels file")
        print(paste("Note 2: Treating input intensities as ",input.intensity.scale," values",sep=""))
        
    }else{
        
        stop("Input intensities should either be at raw or log2 scale")
    }
    
    
    suppressWarnings(suppressWarnings(sink(file=NULL)))
    x<-date()
    x<-strsplit(x,split=" ")
    
    #x<-gsub(x,pattern=":",replacement="_")
    targeted_feat_raw<-{}
    logistic_reg=FALSE
    x1<-unlist(x)
    x1<-gsub(x1,pattern=":",replacement="_")
    #fname<-paste(x1[2:5],collapse="_")
    
    #fname<-paste(x1[2:3],x1[5],x1[4],collapse="_")
    
    fname<-paste(x1[2:3],collapse="")
    
    #fname<-paste(fname,x1[6],sep="")
    x1[4]<-gsub(x1[4],pattern=":",replacement="_")
    fname<-paste(fname,x1[5],sep="")
    fname<-paste(fname,x1[4],sep="_")
    
    
    dir.create(parentoutput_dir)
    setwd(parentoutput_dir)
    
    
    #fname<-paste(parentoutput_dir,"/Log",fname,".txt",sep="")
    
    fname<-paste(parentoutput_dir,"/Log.txt",sep="")
    
    
    if(is.na(foldchangethresh)==FALSE){
     if(log2transform==TRUE && znormtransform==TRUE){
        
           stop("Both log2transform and znormtransform can not be true if foldchangethresh is not equal to NA.")
        }
    }
    
    if(featselmethod=="limma2way")
    {
        
        print("Note 3: lm2wayanova option is recommended for greater than 2x2 designs and this includes post-hoc comparisons")
        
        
    }
    
   
    if(featselmethod=="lm1wayanovarepeat" | featselmethod=="spls1wayrepeat" | featselmethod=="limma1wayrepeat")
    {
        print("Note 3: Class labels format should be: Sample ID, Subject, Time. lm1wayanovarepeat is based on the nlme::lme() function with post-hoc Tukey HSD test.")
    }else{
        
        if(featselmethod=="lm2wayanovarepeat" | featselmethod=="spls2wayrepeat" | featselmethod=="limma2wayrepeat")
        {
            print("Note 3: Class labels format should be: Sample ID, Subject, Factor, Time. lm2wayanovarepeat is based on the nlme::lme() funciton with post-hoc Tukey HSD test. ")
        }
        
    }
    
    
    print("##############################Starting processing now################################")
    print(paste("**Program is running. Please check the logfile for runtime status: ",fname,"**",sep=""))
    
    fname_params<-paste(parentoutput_dir,"/InputParameters.csv",sep="")
    #sink(fname_params)
    # #save(list=ls(),file="cur.Rda")
    c1<-"feature_table_file:"
    c2<-feature_table_file
    #c2<-rbind(c2,feature_table_file)
    c1<-rbind(c1,"parentoutput_dir:")
    c2<-rbind(c2,parentoutput_dir)
    c1<-rbind(c1,"class_labels_file:")
    c2<-rbind(c2,class_labels_file)
    
    c1<-rbind(c1,"num_replicates:")
    c2<-rbind(c2,num_replicates)
    c1<-rbind(c1,"summarize.replicates:")
    c2<-rbind(c2,summarize.replicates)
    c1<-rbind(c1,"summary.method:")
    c2<-rbind(c2,summary.method)
    c1<-rbind(c1,"summary.na.replacement:")
    c2<-rbind(c2,summary.na.replacement)
    c1<-rbind(c1,"rep.max.missing.thresh:")
    c2<-rbind(c2,rep.max.missing.thresh)
    c1<-rbind(c1,"all.missing.thresh:")
    c2<-rbind(c2,all.missing.thresh)
    c1<-rbind(c1,"group.missing.thresh:")
    c2<-rbind(c2,group.missing.thresh)
    c1<-rbind(c1,"input.intensity.scale:")
    c2<-rbind(c2,input.intensity.scale)
    c1<-rbind(c1,"log2transform:")
    c2<-rbind(c2,log2transform)
    c1<-rbind(c1,"medcenter:")
    c2<-rbind(c2,medcenter)
    c1<-rbind(c1,"znormtransform:")
    c2<-rbind(c2,znormtransform)
    c1<-rbind(c1,"quantile_norm:")
    c2<-rbind(c2,quantile_norm)
  
    c1<-rbind(c1,"TIC_norm:")
    c2<-rbind(c2,TIC_norm)
    c1<-rbind(c1,"lowess_norm:")
    c2<-rbind(c2,lowess_norm)
    c1<-rbind(c1,"madscaling:")
    c2<-rbind(c2,madscaling)
    c1<-rbind(c1,"rsd.filt.list:")
    c2<-rbind(c2,rsd.filt.list)
    c1<-rbind(c1,"pairedanalysis:")
    c2<-rbind(c2,pairedanalysis)
    c1<-rbind(c1,"featselmethod:")
    
    c2<-rbind(c2,paste(featselmethod,collapse=";"))
    
    c1<-rbind(c1,"pvalue.thresh:")
    c2<-rbind(c2,pvalue.thresh)
    c1<-rbind(c1,"fdrthresh:")
    c2<-rbind(c2,fdrthresh)
    c1<-rbind(c1,"fdrmethod:")
    c2<-rbind(c2,fdrmethod)
    c1<-rbind(c1,"cor.method:")
    c2<-rbind(c2,cor.method)
    c1<-rbind(c1,"abs.cor.thresh:")
    c2<-rbind(c2,abs.cor.thresh)
    c1<-rbind(c1,"cor.fdrthresh:")
    c2<-rbind(c2,cor.fdrthresh)
    c1<-rbind(c1,"kfold:")
    c2<-rbind(c2,kfold)
    c1<-rbind(c1,"globalcor:")
    c2<-rbind(c2,globalcor)
    c1<-rbind(c1,"target.metab.file:")
    c2<-rbind(c2,target.metab.file)
    c1<-rbind(c1,"target.mzmatch.diff:")
    c2<-rbind(c2,target.mzmatch.diff)
    c1<-rbind(c1,"target.rtmatch.diff:")
    c2<-rbind(c2,target.rtmatch.diff)
    c1<-rbind(c1,"max.cor.num:")
    c2<-rbind(c2,max.cor.num)
    c1<-rbind(c1,"missing.val:")
    c2<-rbind(c2,missing.val)
    c1<-rbind(c1,"networktype:")
    c2<-rbind(c2,networktype)
    c1<-rbind(c1,"samplermindex:")
    c2<-rbind(c2,samplermindex)
    c1<-rbind(c1,"numtrees:")
    c2<-rbind(c2,numtrees)
    c1<-rbind(c1,"analysismode:")
    c2<-rbind(c2,analysismode)
    c1<-rbind(c1,"net_node_colors:")
    c2<-rbind(c2,net_node_colors)
    c1<-rbind(c1,"net_legend:")
    c2<-rbind(c2,net_legend)
    c1<-rbind(c1,"heatmap.col.opt:")
    c2<-rbind(c2,heatmap.col.opt)
    c1<-rbind(c1,"manhattanplot.col.opt:")
    c2<-rbind(c2,paste(manhattanplot.col.opt,collapse=";"))
    c1<-rbind(c1,"boxplot.col.opt:")
    c2<-rbind(c2,boxplot.col.opt)
    
    c1<-rbind(c1,"barplot.col.opt:")
    c2<-rbind(c2,barplot.col.opt)
    
    c1<-rbind(c1,"sample.col.opt:")
    c2<-rbind(c2,sample.col.opt)
    c1<-rbind(c1,"alphacol:")
    c2<-rbind(c2,alphacol)
    c1<-rbind(c1,"pls_vip_thresh:")
    c2<-rbind(c2,pls_vip_thresh)
    c1<-rbind(c1,"alphacol:")
    c2<-rbind(c2,alphacol)
    c1<-rbind(c1,"max_varsel:")
    c2<-rbind(c2,max_varsel)
    c1<-rbind(c1,"pls_ncomp:")
    c2<-rbind(c2,pls_ncomp)
    c1<-rbind(c1,"pcacenter:")
    c2<-rbind(c2,pcacenter)
    c1<-rbind(c1,"pcascale:")
    c2<-rbind(c2,pcascale)
    c1<-rbind(c1,"pred.eval.method:")
    c2<-rbind(c2,pred.eval.method)
    c1<-rbind(c1,"rocfeatlist:")
    c2<-rbind(c2,paste(rocfeatlist,collapse=";"))
    c1<-rbind(c1,"rocfeatincrement:")
    c2<-rbind(c2,rocfeatincrement)
    c1<-rbind(c1,"rocclassifier:")
    c2<-rbind(c2,rocclassifier)
    c1<-rbind(c1,"foldchangethresh:")
    c2<-rbind(c2,foldchangethresh)
    c1<-rbind(c1,"wgcnarsdthresh:")
    c2<-rbind(c2,wgcnarsdthresh)
    c1<-rbind(c1,"WGCNAmodules:")
    c2<-rbind(c2,WGCNAmodules)
    c1<-rbind(c1,"optselect:")
    c2<-rbind(c2,optselect)
    c1<-rbind(c1,"max_comp_sel:")
    c2<-rbind(c2,max_comp_sel)
    c1<-rbind(c1,"saveRda:")
    c2<-rbind(c2,saveRda)
    c1<-rbind(c1,"pca.cex.val:")
    c2<-rbind(c2,pca.cex.val)
    c1<-rbind(c1,"pls.permut.count:")
    c2<-rbind(c2,pls.permut.count)
    c1<-rbind(c1,"pca.ellipse:")
    c2<-rbind(c2,pca.ellipse)
    c1<-rbind(c1,"ellipse.conf.level:")
    c2<-rbind(c2,ellipse.conf.level)
    c1<-rbind(c1,"legendlocation:")
    c2<-rbind(c2,legendlocation)
    c1<-rbind(c1,"svm.acc.tolerance:")
    c2<-rbind(c2,svm.acc.tolerance)
    c1<-rbind(c1,"limmadecideTests:")
    c2<-rbind(c2,limmadecideTests)
    c1<-rbind(c1,"pls.vip.selection:")
    c2<-rbind(c2,pls.vip.selection)
    c1<-rbind(c1,"globalclustering:")
    c2<-rbind(c2,globalclustering)
    c1<-rbind(c1,"plots.res:")
    c2<-rbind(c2,plots.res)
    c1<-rbind(c1,"plots.width:")
    c2<-rbind(c2,plots.width)
    c1<-rbind(c1,"plots.height:")
    c2<-rbind(c2,plots.height)
    c1<-rbind(c1,"plots.type:")
    c2<-rbind(c2,plots.type)
    c1<-rbind(c1,"output.device.type:")
    c2<-rbind(c2,output.device.type)
    c1<-rbind(c1,"pamr.threshold.select.max:")
    c2<-rbind(c2,pamr.threshold.select.max)
    c1<-rbind(c1,"aggregation.method:")
    c2<-rbind(c2,aggregation.method)
    c1<-rbind(c1,"mars.gcv.thresh")
    c2<-rbind(c2,mars.gcv.thresh)
    c1<-rbind(c1,"error.bar")
    c2<-rbind(c2,error.bar)
    

    c1<-cbind(c1,c2)
    c1<-as.data.frame(c1)

    colnames(c1)<-c("InputParameter:","Value")
    write.csv(c1,file=fname_params,row.names=FALSE)
    rm(c1)
    

    
    
    sink(fname)
    print(sessionInfo())
    
if(length(featselmethod)>1){
        
        if(length(rsd.filt.list)>1){
            
            print("Warning: only one RSD threshold allowed for multiple feature selection methods. Only the first RSD threshold will be used.")
            rsd.filt.list=rsd.filt.list[1]
            
        }
        
        consensus_res<-{}
        diffexp.res<-new("list")
        consensus_analysis=TRUE
        ranked_list<-{}
        common_feats<-{}
for(i in 1:length(featselmethod))
{
    
    if(featselmethod[i]=="rfesvm" && analysismode=="regression"){
        try(dev.off(),silent=TRUE)
        next;
    }
    
    
                    outloc<-paste(parentoutput_dir,featselmethod[i],sep="/")
           suppressWarnings(diffexp.res[[i]]<-diffexp.child(Xmat,Ymat,feature_table_file,parentoutput_dir,class_labels_file,num_replicates,feat.filt.thresh,summarize.replicates,summary.method,summary.na.replacement,missing.val,rep.max.missing.thresh,
             all.missing.thresh,group.missing.thresh,input.intensity.scale,
            log2transform,medcenter,znormtransform,quantile_norm,lowess_norm,madscaling,TIC_norm,rsd.filt.list,
            pairedanalysis,featselmethod[i],fdrthresh,fdrmethod,cor.method,networktype,abs.cor.thresh,cor.fdrthresh,kfold,pred.eval.method,feat_weight,globalcor,
            target.metab.file,target.mzmatch.diff,target.rtmatch.diff,max.cor.num,samplermindex,pcacenter,pcascale,
            numtrees,analysismode,net_node_colors,net_legend,svm_kernel,heatmap.col.opt,manhattanplot.col.opt,boxplot.col.opt,barplot.col.opt,sample.col.opt,lineplot.col.opt,hca_type,alphacol,pls_vip_thresh,num_nodes,max_varsel, pls_ncomp=pls_ncomp,pca.stage2.eval=pca.stage2.eval,scoreplot_legend=scoreplot_legend,pca.global.eval=pca.global.eval,rocfeatlist=rocfeatlist,rocfeatincrement=rocfeatincrement,rocclassifier=rocclassifier,foldchangethresh=foldchangethresh,wgcnarsdthresh=wgcnarsdthresh,WGCNAmodules=WGCNAmodules,
            optselect=optselect,max_comp_sel=max_comp_sel,saveRda=saveRda,legendlocation=legendlocation,degree_rank_method=degree_rank_method,pca.cex.val=pca.cex.val,pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,pls.permut.count=pls.permut.count,svm.acc.tolerance=svm.acc.tolerance,limmadecideTests=limmadecideTests,pls.vip.selection=pls.vip.selection,globalclustering=globalclustering,plots.res=plots.res,plots.width=plots.width,plots.height=plots.height,plots.type=plots.type,output.device.type=output.device.type,pvalue.thresh,individualsampleplot.col.opt,pamr.threshold.select.max,mars.gcv.thresh,error.bar,cex.plots,modeltype,barplot.xaxis,lineplot.lty.option)
            

            
            )
            
if(is.na(aggregation.method)==TRUE){
 
    consensus_analysis=FALSE
    
}

##save(diffexp.res,file="diffexp.res.Rda")
            
if(length(diffexp.res[[i]]$all_metabs)>0 && consensus_analysis==TRUE){
    

        diffexp.res[[i]]$all_metabs<-diffexp.res[[i]]$all_metabs[order(diffexp.res[[i]]$all_metabs$mz),]
 
        mz_rt_all<-paste(diffexp.res[[i]]$all_metabs$mz,"_",diffexp.res[[i]]$all_metabs$time,sep="")
        tname<-paste("mz_rt_",i,".Rda",sep="")
        
        check_dup_index<-which(duplicated(mz_rt_all)==TRUE)
        
        tname<-paste("rank_",i,".Rda",sep="")
        tvec=diffexp.res[[i]]$all_metabs$diffexp_rank
        ##save(tvec,file=tname)
        
        if(length(tvec)<1){
            
            try(dev.off(),silent=TRUE)
            next;
        }
        if(length(check_dup_index)>0){
            mz_rt_all<-mz_rt_all[-check_dup_index]
            tvec<-tvec[-check_dup_index]
            diffexp.res[[i]]$all_metabs<-diffexp.res[[i]]$all_metabs[-check_dup_index,]
        }
        
        cur_data_res<-diffexp.res[[i]]$all_metabs
        
if(i==1){
    
    #sort(rankingCriteria, index.return = TRUE)$ix
        ranked_list<-mz_rt_all[sort(tvec,index.return=TRUE)$ix]
}

        if(i>1){
            mz_rt_1<-paste(diffexp.res[[(i-1)]]$diffexp_metabs$mz,"_",diffexp.res[[(i-1)]]$diffexp_metabs$time,sep="")
            mz_rt_2<-paste(diffexp.res[[i]]$diffexp_metabs$mz,"_",diffexp.res[[i]]$diffexp_metabs$time,sep="")
            cnamesd1<-colnames(diffexp.res[[(i-1)]]$diffexp_metabs)
            time_ind<-which(cnamesd1=="time")
            mz_ind<-which(cnamesd1=="mz")
            common_feat_ind<-which(mz_rt_1%in%mz_rt_2)
            common_feat_ind2<-which(mz_rt_2%in%mz_rt_1)
            
            
            ranked_list<-rbind(ranked_list,mz_rt_all[sort(tvec,index.return=TRUE)$ix])
            
            if(length(common_feat_ind)<1){
             ermsg<-paste("No common selected features found between ",featselmethod[(i-1)]," and ",featselmethod[i],sep="")
             #stop(ermsg)
             if(aggregation.method=="Consensus"){
                 consensus_analysis=FALSE
             }
            }else{
            mat_A<-diffexp.res[[(i-1)]]$diffexp_metabs[common_feat_ind,c(1:(mz_ind+1))]
            mat_B<-diffexp.res[[i]]$diffexp_metabs[common_feat_ind2,]
            mz_rt_1<-paste(mat_A$mz,"_",mat_A$time,sep="")
            mz_rt_2<-paste(mat_B$mz,"_",mat_B$time,sep="")
            m1<-match(mz_rt_1,mz_rt_2)
            cnamesd1<-colnames(mat_A)
            time_ind<-which(cnamesd1=="time")
            mz_ind<-which(cnamesd1=="mz")
            #mat_A<-mat_A[,c(1:(mz_ind-1))]
            
            mat_A<-mat_A[,c(mz_ind:time_ind)]
            cnamesd2<-colnames(mat_B)
            time_ind<-which(cnamesd2=="time")
            mz_ind<-which(cnamesd2=="mz")
            
            mat_B<-mat_B[,-c(1:(time_ind))]
            common_feats<-cbind(mat_A,mat_B[m1,]) #merge(mat_B,mat_A,by.x="mz_rt_2",by.y="mz_rt_1")
            print(paste("Number of common sig feats between ",featselmethod[(i-1)]," and ",featselmethod[i],sep=""))
            print(dim(common_feats))
            
            }
        }
}
}

if(length(ranked_list)<1){
    consensus_analysis=FALSE
    
}

if(consensus_analysis==TRUE){
    
    print("################")
    print(paste("Aggregating results from different methods using ",aggregation.method," aggregation method.",sep=""))


        print("################")
    
    if(aggregation.method=="RankAggreg" | aggregation.method=="RankAggregGA"){
       
    if(max_varsel>dim(ranked_list)[2]){
        max_varsel=round(dim(ranked_list)[2]*0.3)
        
        
    }
   
   if(aggregation.method=="RankAggreg"){
    r1<-RankAggreg(x=ranked_list,k=max_varsel,verbose=TRUE,distance="Spearman",method="CE",maxIter=aggregation.max.iter)
   }else{
       
       r1<-RankAggreg(x=ranked_list,k=max_varsel,verbose=TRUE,distance="Spearman",method="GA",maxIter=aggregation.max.iter)
   }
   
    common_row_index<-which(mz_rt_all%in%r1$top.list)
    
    
    common_feats<-cur_data_res[common_row_index,]
  
    cnamesd1<-colnames(common_feats)
    time_ind<-which(cnamesd1=="time")
    mz_ind<-which(cnamesd1=="mz")

    
    }
    
    cnames_1<-try(colnames(common_feats),silent=TRUE)
    
    
    if(consensus_analysis==FALSE | is(cnames_1,"try-error") | length(common_feats)<1){
        
        print("Skipping aggregation.")
    }else{
       
    
    bad_colind<-grep(tolower(cnames_1),pattern="rank")
    
    bad_colind_2<-grep(tolower(cnames_1),pattern="max.fold.change.log2")
    
    if(length(bad_colind_2)>1){
        bad_colind_2<-bad_colind_2[-c(1)]
    }
    
    bad_colind<-c(bad_colind,bad_colind_2)
    
    
    if(nrow(common_feats)>0){
        
        if(length(bad_colind)>0){
            common_feats<-common_feats[,-c(bad_colind)]
        }
    }


        common_feats<-unique(common_feats)
        print("Dimension of aggregated feature table of selected features:")
        print(dim(common_feats))
        
        num_common_feats<-dim(common_feats)[1]
        
        if(num_common_feats<1){
            
            stop("No common features found.")
        }
        ##save(common_feats,file="common_feats.Rda")
        
        
        cnamesd1<-colnames(common_feats)
        time_ind<-which(cnamesd1=="time")
        mz_ind<-which(cnamesd1=="mz")
        
 
        #Xmat<-common_feats[,-c(1:2)]
        mz<-common_feats[,mz_ind]
        time<-common_feats[,time_ind]
        rnames1<-paste(mz,time,sep="_")
        
        Xmat<-cbind(mz,time,common_feats[,-c(1:time_ind)])
        
        rownames(Xmat)<-rnames1
        
        if(max_varsel>nrow(Xmat)){
            
            max_varsel<-nrow(Xmat)
        }

        #Ymat<-cbind(colnames(demetabs_res$norm_data[,-c(1:2)]),diffexp.res[i]$classlabels)
        Ymat<-diffexp.res[[i]]$classlabels
        
        if(nrow(Xmat)>1){
            
        outloc<-paste(parentoutput_dir,"AggregatedResults/",sep="/")
        
        if(log2transform==TRUE){
            
            input.intensity.scale="log2"
        }else{
            input.intensity.scale="raw"
        }
        
        dir.create(outloc)
        setwd(outloc)
        
        dir.create("Figures")
        dir.create("Tables")
        
        write.table(Xmat,file="Tables/Aggregated_selected_features.txt",sep="\t",row.names=FALSE)
        
        subdata<-t(Xmat[,-c(1:2)])
        classlabels<-Ymat[,2]
        classlabels<-as.data.frame(classlabels)
        if(analysismode=="classification"){
        svm_model<-try(svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
        #
        
        
        
                if(is(svm_model,"try-error")){
                    kfold_acc<-NA
                    
                }else{
                    kfold_acc<-svm_model$avg_acc
                }
        }else{
            
            svm_model_reg<-try(svm(x=subdata,y=(classlabels[,1]),type="eps",cross=kfold),silent=TRUE)
            
            if(is(svm_model_reg,"try-error")){
                kfold_acc<-NA
                
            }else{
            
                kfold_acc<-svm_model_reg$tot.MSE
            }
        }
        numcores<-num_nodes #round(detectCores()*0.6)
        
        kfold_acc_rand<-{}
        #for(p1 in 1:100){
         cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
         clusterEvalQ(cl,library(e1071))
         clusterEvalQ(cl,library(pROC))
          clusterEvalQ(cl,library(ROCR))
         clusterExport(cl,"svm_cv",envir = .GlobalEnv)
         clusterExport(cl,"svm",envir = .GlobalEnv)
         if(analysismode=="classification"){
         kfold_acc_rand<-parLapply(cl,1:100,function(p1){
            
            sample_ind<-sample(1:dim(classlabels)[1],size=dim(classlabels)[1])
            classlabels_rand<-classlabels[sample_ind,]
            classlabels_rand<-as.data.frame(classlabels_rand)
            svm_model<-try(svm_cv(v=kfold,x=subdata,y=classlabels_rand,kname=svm_kernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
            #
            
            if(is(svm_model,"try-error")){
                # kfold_acc_rand<-c(kfold_acc_rand,NA)
                return(NA)
            }else{
                #kfold_acc_rand<-c(kfold_acc_rand,svm_model$avg_acc)
                return(svm_model$avg_acc)
            }
        })
         
        }else{
            
            kfold_acc_rand<-parLapply(cl,1:100,function(p1){
                
                sample_ind<-sample(1:dim(classlabels)[1],size=dim(classlabels)[1])
                classlabels_rand<-classlabels[sample_ind,]
                classlabels_rand<-as.data.frame(classlabels_rand)
               
                svm_model_reg<-try(svm(x=subdata,y=(classlabels[,1]),type="eps",cross=kfold),silent=TRUE)
                
                if(is(svm_model_reg,"try-error")){
                    # kfold_acc_rand<-c(kfold_acc_rand,NA)
                    return(NA)
                }else{
                    #kfold_acc_rand<-c(kfold_acc_rand,svm_model$avg_acc)
                    return(svm_model_reg$tot.MSE)
                }
            })
        }
         stopCluster(cl)
        
        kfold_acc_rand<-unlist(kfold_acc_rand)
        
        kfold_acc_rand<-mean(kfold_acc_rand,na.rm=TRUE)
        
        num_common_feats<-dim(common_feats)[1]
        summary_res<-cbind(num_common_feats,kfold_acc,kfold_acc_rand)
        colnames(summary_res)<-c("Number of selected features after aggregation",paste(pred.eval.method,"-accuracy",sep=""),paste(pred.eval.method," permuted accuracy",sep=""))

        
        file_name<-paste("../Results_summary_aggregated.txt",sep="")
        write.table(summary_res,file=file_name,sep="\t",row.names=FALSE)
        
        if(output.device.type=="pdf"){
        pdf("Figures/Aggregatedresults.pdf")
        }
        
        if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/HCA_aggregated_selectedfeats.png"
            
            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        }
        
        #   print("here")
        #print(head(Ymat))
        ##save(Xmat,file="Xmat.Rda")
        ##save(Ymat,file="Ymat.Rda")
        
        g1<-get_hca(feature_table_file=NA,parentoutput_dir=outloc,class_labels_file=NA,X=Xmat,Y=Ymat,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode=analysismode,
        sample.col.opt=sample.col.opt,plots.width=plots.width,plots.height=plots.height,plots.res=plots.res, plots.type=plots.type, alphacol=0.3, hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="",cexRow=0.4,cexCol=0.4)
        
        if(output.device.type!="pdf"){
            
            try(dev.off(),silent=TRUE)
        }


if(analysismode=="classification")
{
        best_subset<-{}
        best_acc<-0
        
        xvec<-{}
        yvec<-{}
        for(i in 2:max_varsel){
            
            subdata<-t(Xmat[1:i,-c(1:2)])
            svm_model<-try(svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
            #svm_model<-svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95
            
            if(is(svm_model,"try-error")){
                
                svm_model<-NA
                
            }else{
                
                xvec<-c(xvec,i)
                yvec<-c(yvec,svm_model$avg_acc)
                if(svm_model$avg_acc>best_acc){
                    
                    best_acc<-svm_model$avg_acc
                    best_subset<-seq(1,i)
                    
                    
                }
                
                if(svm_model$avg_acc<best_acc){
                    
                    diff_acc<-best_acc-svm_model$avg_acc
                    if(diff_acc>50){
                        
                        break;
                        
                    }
                    
                }
            }
            
            
        }
        
        if(pred.eval.method=="CV"){
            ylab_text=paste(pred.eval.method," accuracy (%)",sep="")
            
        }else{
            if(pred.eval.method=="BER"){
                ylab_text=paste("Balanced accuracy"," (%)",sep="")
            }else{
                
                ylab_text=paste("AUC"," (%)",sep="")
            }
        }
        
        
        if(length(yvec)>0){
            
            msg1<-paste("k-fold CV classification accuracy based on forward selection of\n aggregated features ordered by ",featselmethod[1],sep="")
            
            
            if(output.device.type!="pdf"){
                
                temp_filename_1<-"Figures/kfold_forward_selection_aggregated_selectedfeats.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
            }
            
            
            plot(x=xvec,y=yvec,main=msg1,xlab="Feature index",ylab=ylab_text,type="b",col="brown")
            
            
            if(output.device.type!="pdf"){
                
                try(dev.off(),silent=TRUE)
            }


            cv_mat<-cbind(xvec,yvec)
            colnames(cv_mat)<-c("Feature Index",ylab_text)
            
            write.table(cv_mat,file="Tables/aggregated_kfold_cv_mat.txt",sep="\t")
        }
        
       
        
        if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/Boxplots_aggregated_selectedfeats.png"
            
            #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
            pdf(temp_filename_1)
            
        }
        
        if(log2transform==TRUE || input.intensity.scale=="log2"){
            
            if(znormtransform==TRUE){
                ylab_text_2="scale normalized"
            }else{
                if(quantile_norm==TRUE){
                    
                    ylab_text_2="quantile normalized"
                }else{
                    
                    ylab_text_2=""
                }
            }
            ylab_text=paste("log2 intensity ",ylab_text_2,sep="")
        }else{
            if(znormtransform==TRUE){
                ylab_text_2="scale normalized"
            }else{
                if(quantile_norm==TRUE){
                    
                    ylab_text_2="quantile normalized"
                }else{
                    ylab_text_2=""
                }
            }
            ylab_text=paste("Raw intensity ",ylab_text_2,sep="")
        }
        
        # par(mfrow=c(2,2))
        par(mfrow=c(1,1),family="sans",cex=cex.plots)
        get_boxplots(X=Xmat,Y=Ymat,parentoutput_dir=outloc,sample.col.opt=sample.col.opt,boxplot.col.opt=boxplot.col.opt, alphacol=0.3,newdevice=FALSE,cex=0.6,ylabel=ylab_text)
        
       
        
        if(output.device.type!="pdf"){
            
            try(dev.off(),silent=TRUE)
        }
        
      
         
         if(output.device.type!="pdf"){
             
             temp_filename_1<-"Figures/Barplots_aggregated_selectedfeats.pdf"
             
             #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
             pdf(temp_filename_1)
         }
         # par(mfrow=c(2,2))
         
         
         
           par(mfrow=c(1,1),family="sans",cex=cex.plots)
           
           get_barplots(feature_table_file=NA,class_labels_file=NA,X=Xmat,Y=Ymat,parentoutput_dir=outloc,newdevice=FALSE,ylabel=ylab_text,cex.val=cex.plots,barplot.col.opt=barplot.col.opt,error.bar=error.bar,barplot.xaxis=barplot.xaxis)

if(output.device.type!="pdf"){
    
    try(dev.off(),silent=TRUE)
}


if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/Individual_sample_plots_aggregated_selectedfeats.png"
    
    #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
    pdf(temp_filename_1)
}

par(mfrow=c(1,1),family="sans",cex=cex.plots)
get_individualsampleplots(feature_table_file=NA,class_labels_file=NA,X=Xmat,Y=Ymat,parentoutput_dir=outloc,newdevice=FALSE,
ylabel=ylab_text,cex.val=cex.plots,sample.col.opt=individualsampleplot.col.opt)
            
            
              if(output.device.type!="pdf"){
        
                try(dev.off(),silent=TRUE)
              }
 
if(pairedanalysis==TRUE)
 {
     
    
      if(output.device.type!="pdf"){
          
          temp_filename_1<-"Figures/Lineplots_aggregated_selectedfeats.png"
          pdf(temp_filename_1)
          #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
      }
      
      par(mfrow=c(1,1),family="sans",cex=cex.plots)
     classlabels_orig<-read.table("../Stage2/classlabels_orig.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE,check.names = FALSE, quote="")
     classlabels_orig<-classlabels_orig[,-c(1)]
     
     ##save(list=ls(),file="debuga_lineplots.Rda")
      try(get_lineplots(X=Xmat,Y=classlabels_orig,feature_table_file=NA,parentoutput_dir=getwd(),class_labels_file=NA,lineplot.col.opt=lineplot.col.opt, alphacol=0.3,col_vec=col_vec,pairedanalysis=pairedanalysis,pca.cex.val=pca.cex.val,legendlocation=legendlocation,pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,filename="selected",ylabel=ylab_text,error.bar=error.bar,cex.val=cex.plots,lineplot.lty.option=lineplot.lty.option),silent=TRUE)
 }
 
    if(output.device.type!="pdf"){
        
        try(dev.off(),silent=TRUE)
    }
 }
 
 
 if(output.device.type=="pdf"){
        try(dev.off(),silent=TRUE)
 }
        }
    }
print("#######################")
print("#")
print("#")
print("Program ended successfully.")
print("#")
print("#")
print("#######################")

suppressWarnings(sink(file=NULL))

#print("###############################")
#print("###############################")
#print("###############################")
time_end<-Sys.time()

time_taken_panda<-round(time_end-time_start,2)


#print("*********")
print(paste("**Program ended successfully in ",time_taken_panda," ",units(time_taken_panda),". Please see the ReadMe.txt file for description of output files and folders.**", sep=""))

#print("*********")
#print(paste("All result files are in the specified output location: ",parentoutput_dir,sep=""))
#print("*********")
#print("There will be a sub-folder for each step: Stage 1: pre-processing, Stage 2: statistical analysis (e.g. limma, PLS), and Stages 3 and 4: network analysis (Global and/or Targeted).")
#print("*********")
#print("")
#print("*********")
# print("Enjoy!")
print("##############################END################################")




s1<-"Stage 1 results: Preprocessing (Normalization, transformation)"
s2<-"Stage 2 results: Feature selection & evaluation results for each feature selection method. Description of files and folder: A) *selected.features.final.txt file includes final table of selected features; B) Figures subfolder: includes Manhattan plots, boxplots, barplots, and other graphics, and C) Tables sub-folder: includes data files with feature selection results for all features, PCA (and PLS) scores and loadings, HCA clusters, k-fold CV results, and files for generating boxplots and barplots."
s3<-"Stage 3 results: Correlation based network analysis"
s4<-"Stage 4 results: Correlation based targeted network analysis"
s5<-"Consensus results: HCA, k-fold CV, boxplots, and barplots for aggregated selected features"
sm<-rbind(s1,s2,s3,s4,s5)
setwd(parentoutput_dir)
write.table(sm,file="ReadMe.txt",sep="\t",row.names=FALSE)
return(list("individual.featsel.res"=diffexp.res,"aggregated.res"=common_feats))
}else{
    
    print("#######################")
    print("#")
    print("#")
    print("Program ended successfully.")
    print("#")
    print("#")
    print("#######################")
    
    suppressWarnings(sink(file=NULL))
    
    # print("###############################")
    #print("###############################")
    #print("###############################")
    time_end<-Sys.time()
    
    time_taken_panda<-round(time_end-time_start,2)
    
print(paste("**Program ended successfully in ",time_taken_panda," ",units(time_taken_panda),". Please see the ReadMe.txt file for description of output files and folders.**", sep=""))
   
    
    # print(paste("*******Program ended successfully in ",round(time_taken_panda,2)," minutes*******", sep=""))
    
    #  print("*     *")
    #  print("Consensus analysis could not be performed as not all features were selected by all feature selection methods.")
    #print("*     *")
    #print(paste("All result files are in the specified output location: ",parentoutput_dir,sep=""))
    # print("*     *")
    #print("There will be a sub-folder for each step: Stage 1: pre-processing, Stage 2: statistical analysis (e.g. limma, PLS), and Stages 3 and 4: network analysis (Global and/or Targeted).")
    #  print("*     *")
    
    #print("Please see the ReadMe.txt file for more information.")
    # print("*     *")
    # print("Enjoy!")
    print("##############################END################################")
    
    
    
    s1<-"Stage 1 results: Preprocessing (Normalization, transformation)"
        s2<-"Stage 2 results: Feature selection & evaluation results for each feature selection method. Description of files and folder: A) *selected.features.final.txt file includes final table of selected features; B) Figures subfolder: includes Manhattan plots, boxplots, barplots, and other graphics, and C) Tables sub-folder: includes data files with feature selection results for all features, PCA (and PLS) scores and loadings, HCA clusters, k-fold CV results, and files for generating boxplots and barplots."
    s3<-"Stage 3 results: Correlation based network analysis"
    s4<-"Stage 4 results: Correlation based targeted network analysis"
    
    sm<-rbind(s1,s2,s3,s4)
    setwd(parentoutput_dir)
    write.table(sm,file="ReadMe.txt",sep="\t",row.names=FALSE)
    return(list("individual.featsel.res"=diffexp.res))
    
}


    
    }else{
        
        suppressWarnings(
        diffexp.res<-diffexp.child(Xmat,Ymat,feature_table_file,parentoutput_dir,class_labels_file,num_replicates,feat.filt.thresh,summarize.replicates,summary.method,summary.na.replacement,missing.val,rep.max.missing.thresh,
        all.missing.thresh,group.missing.thresh,input.intensity.scale,
        log2transform,medcenter,znormtransform,quantile_norm,lowess_norm,madscaling,TIC_norm,rsd.filt.list,
        pairedanalysis,featselmethod,fdrthresh,fdrmethod,cor.method,networktype,abs.cor.thresh,cor.fdrthresh,kfold,pred.eval.method,feat_weight,globalcor,
        target.metab.file,target.mzmatch.diff,target.rtmatch.diff,max.cor.num,samplermindex,pcacenter,pcascale,
        numtrees,analysismode,net_node_colors,net_legend,svm_kernel,heatmap.col.opt,manhattanplot.col.opt,boxplot.col.opt,barplot.col.opt,sample.col.opt,lineplot.col.opt,hca_type,alphacol,pls_vip_thresh,num_nodes,max_varsel, pls_ncomp,pca.stage2.eval,scoreplot_legend,pca.global.eval,rocfeatlist,rocfeatincrement,rocclassifier,foldchangethresh,wgcnarsdthresh,WGCNAmodules,
        optselect,max_comp_sel,saveRda,legendlocation,degree_rank_method,pca.cex.val,pca.ellipse,ellipse.conf.level,pls.permut.count,svm.acc.tolerance,limmadecideTests,pls.vip.selection,globalclustering,plots.res,plots.width,plots.height,plots.type,output.device.type,pvalue.thresh,individualsampleplot.col.opt,pamr.threshold.select.max,mars.gcv.thresh,error.bar,cex.plots,modeltype,barplot.xaxis,lineplot.lty.option)
        )
        
        time_end<-Sys.time()
        
        time_taken_panda<-round(time_end-time_start,2)
        
        
        print("#######################")
        print("#")
        print("#")
        print("Program ended successfully.")
        print("#")
        print("#")
        print("#######################")
        suppressWarnings(sink(file=NULL))
        
print(paste("**Program ended successfully in ",time_taken_panda," ",units(time_taken_panda),". Please see the ReadMe.txt file for description of output files and folders.**", sep=""))
        
        #     print("###############################")
        #print("###############################")
        #print("###############################")
        #print("*     *")
        
        # print(paste("***Program ended successfully in ",round(time_taken_panda,2)," minutes***", sep=""))
        #  print(paste("*******Program ended successfully in ",round(time_taken_panda,2)," minutes*******", sep=""))
        
        # print("*     *")
        #    print(paste("All result files are in the specified output location: ",parentoutput_dir,sep=""))
        # print("*     *")
        #   print("There will be a sub-folder for each step: Stage 1: pre-processing, Stage 2: statistical analysis (e.g. limma, PLS), and Stages 3 and 4: network analysis (Global and/or Targeted).")
        # print("*     *")
        
        #   print("Please see the ReadMe.txt file for more information.")
        #   print("*     *")
        # print("Enjoy!")
        print("##############################END################################")
        
        
        s1<-"Stage 1 results: Preprocessing (Normalization, transformation)"
        s2<-"Stage 2 results: Feature selection & evaluation results for each feature selection method. Description of files and folder: A) *selected.features.final.txt file includes final table of selected features; B) Figures subfolder: includes Manhattan plots, boxplots, barplots, and other graphics, and C) Tables sub-folder: includes data files with feature selection results for all features, PCA (and PLS) scores and loadings, HCA clusters, k-fold CV results, and files for generating boxplots and barplots."
        s3<-"Stage 3 results: Correlation based network analysis"
        s4<-"Stage 4 results: Correlation based targeted network analysis"
        sm<-rbind(s1,s2,s3,s4)
        setwd(parentoutput_dir)
        write.table(sm,file="ReadMe.txt",sep="\t",row.names=FALSE)
        return(diffexp.res)
        
        
    }
}
	
diffexp.child<-function(Xmat,Ymat,feature_table_file,parentoutput_dir,class_labels_file,num_replicates,feat.filt.thresh,summarize.replicates,summary.method,
summary.na.replacement,missing.val,rep.max.missing.thresh,
 all.missing.thresh,group.missing.thresh,input.intensity.scale,
log2transform,medcenter,znormtransform,quantile_norm,lowess_norm,madscaling,TIC_norm,rsd.filt.list,
pairedanalysis,featselmethod,fdrthresh,fdrmethod,cor.method,networktype,abs.cor.thresh,cor.fdrthresh,kfold,pred.eval.method,feat_weight,globalcor,
target.metab.file,target.mzmatch.diff,target.rtmatch.diff,max.cor.num, samplermindex,pcacenter,pcascale,
numtrees,analysismode,net_node_colors,net_legend,svm_kernel,heatmap.col.opt,manhattanplot.col.opt,boxplot.col.opt,barplot.col.opt,sample.col.opt,lineplot.col.opt,hca_type,alphacol,pls_vip_thresh,num_nodes,max_varsel,pls_ncomp,pca.stage2.eval,scoreplot_legend,pca.global.eval,rocfeatlist,rocfeatincrement,rocclassifier,foldchangethresh,wgcnarsdthresh,WGCNAmodules,optselect,max_comp_sel,saveRda,legendlocation,degree_rank_method,pca.cex.val,pca.ellipse,ellipse.conf.level,pls.permut.count,svm.acc.tolerance,limmadecideTests,pls.vip.selection,globalclustering,plots.res,plots.width,plots.height,plots.type,output.device.type,pvalue.thresh,individualsampleplot.col.opt,pamr.threshold.select.max,mars.gcv.thresh,error.bar,cex.plots,modeltype,barplot.xaxis,lineplot.lty.option)
{
    
   
   
		#############
        
        remove_firstrun=FALSE #TRUE or FALSE
        run_number=1
        minmaxtransform=FALSE
        pca.CV=TRUE
        max_rf_var=5000
        
        logistic_reg=FALSE
        poisson_reg=FALSE
        goodfeats_allfields={}
        mwan_fdr={}
        targetedan_fdr={}
        data_m_fc_withfeats={}
        classlabels_orig={}
        
        if(input.intensity.scale=="log2"){
            
            log2transform=FALSE
        }
        featselmethod<-unique(featselmethod)
        parentfeatselmethod<-featselmethod
        rfconditional=FALSE
        
        print("############################")
        print("Feature selection method:")
      
        print(featselmethod)
         print("############################")
        if(featselmethod=="rf" | featselmethod=="RF"){
            
            featselmethod="RF"
            
            rfconditional=FALSE
        }else{
            
            if(featselmethod=="rfconditional" | featselmethod=="RFconditional" | featselmethod=="RFcond" | featselmethod=="rfcond"){
                
                featselmethod="RF"
                
                rfconditional=TRUE
            }
        }
        
        
        if(featselmethod=="rf"){
            
            featselmethod="RF"
        }else{
            if(featselmethod=="mars"){
            
                featselmethod="MARS"
            }
        }
		log2.fold.change.thresh_list<-rsd.filt.list
		if(featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma2wayrepeat" | featselmethod=="limma1wayrepeat"){
			
			if(analysismode=="regression"){
				
				stop("Invalid analysis mode. Please set analysismode=\"classification\".")
			}
		
		print("##############Level 1: Using LIMMA function to find differentially expressed metabolites###########")
		}else{
			if(featselmethod=="RF"){
			
			print("##############Level 1: Using random forest function to find discriminatory metabolites###########")
			#log2.fold.change.thresh_list<-c(0)
			
			}else{
					if(featselmethod=="RFcond"){
			
						print("##############Level 1: Using conditional random forest function to find discriminatory metabolites###########")
						#stop("Please use \"limma\", \"RF\", or \"MARS\".")
			
					}else{
					if(featselmethod=="MARS"){
			
						print("##############Level 1: Using MARS to find discriminatory metabolites###########")
							#log2.fold.change.thresh_list<-c(0)
					}else{
						
							if(featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="poissonreg" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="rfesvm" | featselmethod=="wilcox" | featselmethod=="ttest" | featselmethod=="pamr" | featselmethod=="ttestrepeat" | featselmethod=="wilcoxrepeat"){
								print("##########Level 1: Finding discriminatory metabolites ###########")

							if(featselmethod=="logitreg"){

								featselmethod="lmreg"
								logistic_reg=TRUE
                                poisson_reg=FALSE
                            }else{
                                
                                if(featselmethod=="poissonreg"){
                                        poisson_reg=TRUE
                                        featselmethod="lmreg"
                                        logistic_reg=FALSE
                                }else{
                                logistic_reg=FALSE
                                poisson_reg=FALSE
                                }
                                
                            }
								}else{
									
									if(featselmethod=="pls" | featselmethod=="o1pls" | featselmethod=="o2pls" | featselmethod=="spls" | featselmethod=="spls1wayrepeat" | featselmethod=="spls2wayrepeat" | featselmethod=="pls2way" | featselmethod=="spls2way" | featselmethod=="o1spls" | featselmethod=="o2spls"){
										
										print("##########Level 1: Finding discriminatory metabolites ###########")
										
									}else{
                                        
                                        
                                        stop("Invalid featselmethod specified.")
									}
									
									}
						
						#stop("Invalid featselmethod specified. Please use \"limma\", \"RF\", or \"MARS\".")
						}
			
					}
			
			}
		
		}
		####################################################################################
		
		if(is.na(Xmat)==TRUE){
		X<-read.table(feature_table_file,sep="\t",header=TRUE)
		
			
            #print(X[1:3,])
		
		#X<-X[order(X$mz),]
		
        X[,1]<-round(X[,1],5)
        X[,2]<-round(X[,2],2)
        
        
		Xmat<-t(X[,-c(1:2)])
		
		Xmat<-as.data.frame(Xmat)
		
        #print(Xmat[1:3,1:5])
		
				}else{
					X<-Xmat
                    
                    X[,1]<-round(X[,1],5)
                    X[,2]<-round(X[,2],2)
                    
					Xmat<-t(X[,-c(1:2)])
		
					Xmat<-as.data.frame(Xmat)
		}
		
        
		
        ##save(Xmat,file="Xmat.Rda")
						
	if(analysismode=="regression")
	{
	
					#log2.fold.change.thresh_list<-c(0)
					
				print("Performing regression analysis")
				if(is.na(Ymat)==TRUE){
								classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
						
					}else{
					classlabels<-Ymat
					
						}
                    
                    classlabels_orig<-classlabels
                    classlabels_sub<-classlabels
						class_labels_levels<-c("A")
									if(featselmethod=="lmregrepeat"){
				     				colnames(classlabels)<-c("SampleID","SubjectNum",paste("Response",seq(1,dim(factor_inf)[2]),sep=""))
											
											
											#Xmat<-chocolate[,1]
											Xmat_temp<-Xmat #t(Xmat)
											Xmat_temp<-cbind(classlabels,Xmat_temp)
											
											#Xmat_temp<-Xmat_temp[order(Xmat_temp[,3],Xmat_temp[,2]),]

											cnames<-colnames(Xmat_temp)
											
											factor_lastcol<-grep("^Response", cnames)
											
											classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
																						
											subject_inf<-classlabels[,2]
											classlabels<-classlabels[,-c(2)]
											
											Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
											
									}
							

					classlabels<-as.data.frame(classlabels)
					classlabels_response_mat<-classlabels[,-c(1)]
					classlabels_response_mat<-as.data.frame(classlabels_response_mat)
					Ymat<-classlabels
					Ymat<-as.data.frame(Ymat)
                    
                    #print(Xmat[1:3,1:5])
                    
					
									
                                    
                                    rnames_xmat<-rownames(Xmat)
                                    rnames_ymat<-as.character(Ymat[,1])
                                    
                                    
                                    if(length(which(duplicated(rnames_ymat)==TRUE))>0){
                                        
                                        stop("Duplicate sample IDs are not allowed. Please represent replicates by _1,_2,_3.")
                                    }
                                    
                                    check_ylabel<-regexpr(rnames_ymat[1],pattern="^[0-9]*",perl=TRUE)
                                    check_xlabel<-regexpr(rnames_xmat[1],pattern="^X[0-9]*",perl=TRUE)
                                   
                                   #print(check_xlabel)
                                   #print(check_ylabel)
                                    if(length(check_ylabel)>0 && length(check_xlabel)>0){
                                        if(attr(check_ylabel,"match.length")>0 && attr(check_xlabel,"match.length")>0){
                                            
                                            rnames_ymat<-paste("X",rnames_ymat,sep="") #gsub(rnames_ymat,pattern="\\.[0-9]*",replacement="")
                                            
                                            
                                        }
                                    }
                                    
                                    match_names<-match(rnames_xmat,rnames_ymat)
                                    
                                    bad_colnames<-length(which(is.na(match_names)==TRUE))
                                   
                                    #if(is.na()==TRUE){
                                    
                                    if(bad_colnames>0){
                                        print("Sample names do not match between feature table and class labels files.\n Please try replacing any \"-\" with \".\" in sample names.")
                                        print("Sample names in feature table")
                                        print(head(rnames_xmat))
                                        print("Sample names in classlabels file")
                                        
                                        print(head(rnames_ymat))
                                        stop("Sample names do not match between feature table and class labels files.\n Please try replacing any \"-\" with \".\" in sample names. Please try again.")
                                    }
                                   
                                   Xmat<-t(Xmat)
                                   Xmat<-cbind(X[,c(1:2)],Xmat)
                                   Xmat<-as.data.frame(Xmat)
                                   # print(Xmat[1:3,1:5])
                                   #TIC_norm=FALSE
                        if(is.na(all(diff(match(rnames_xmat,rnames_ymat))))==FALSE){
                                    if(all(diff(match(rnames_xmat,rnames_ymat)) > 0)==TRUE){
                                        
                                        data_matrix<-data_preprocess(Xmat,Ymat,feature_table_file,parentoutput_dir,class_labels_file=NA,num_replicates=num_replicates,feat.filt.thresh,summarize.replicates,
                                        summary.method, all.missing.thresh,group.missing.thresh=NA,
                                        log2transform,medcenter,znormtransform,quantile_norm,lowess_norm,madscaling,missing.val, samplermindex,rep.max.missing.thresh,summary.na.replacement,featselmethod,TIC_norm)
                                    }
                                    }else{
                                        
                                        #print(diff(match(rnames_xmat,rnames_ymat)))
                                        stop("Orders of feature table and classlabels do not match")
                                    }
                                    

                                    #			print("Before preprocessing")
                                    #print("Ymat")
                                    #	print(Ymat[1:3,])
                                    #	print("Xmat:")
                                    #	print(Xmat[1:3,1:5])
			
					
						
						write.table(classlabels,file="ordered_response_matrix.txt",sep="\t",row.names=FALSE)
		
						
				}else{
							if(analysismode=="classification")
							{
		
				classlabels_sub<-NA
                
                if(pairedanalysis==TRUE){
                    
                    if(featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma1way" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="MARS" | featselmethod=="RF" | featselmethod=="pls" | featselmethod=="o1pls" | featselmethod=="o2pls" | featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="spls" | featselmethod=="pls2way" | featselmethod=="spls2way" | featselmethod=="o1spls" | featselmethod=="o2spls" |
                        featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="rfesvm" | featselmethod=="pamr" | featselmethod=="poissonreg")
                    {
                        # print("here")
                            
                        stop(paste(featselmethod," for paired analysis, featselmethod should be: limma1wayrepeat, limma2wayrepeat, lm1wayanovarepeat, lm2wayanovarepeat, spls1wayrepeat, or spls2wayrepeat",sep=""))
                        
                    }
                    
                    
                }
		
      
								if(is.na(Ymat)==TRUE){
											classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
                                            Ymat<-classlabels
									
								}else{
								classlabels<-Ymat
								
									}
                                
                                print(paste("Number of samples in class labels file:",dim(Ymat)[1],sep=""))
                                print(paste("Number of samples in feature table:",dim(Xmat)[1],sep=""))
                                
                                if(dim(Ymat)[1]!=(dim(Xmat)[1]))
                                {
                                    
                                    stop("Number of samples are different in feature table and class labels file.")
                                }
                                
                                
                                if(fdrmethod=="none"){
                                    
                                    fdrthresh=pvalue.thresh
                                }
									
		if(featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma2wayrepeat" | featselmethod=="limma1way" | featselmethod=="limma1wayrepeat" | featselmethod=="MARS" | featselmethod=="RF" | featselmethod=="pls" | featselmethod=="o1pls" | featselmethod=="o2pls" | featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="spls" | featselmethod=="pls1wayrepeat" | featselmethod=="spls1wayrepeat" | featselmethod=="pls2wayrepeat" | featselmethod=="spls2wayrepeat" | featselmethod=="pls2way" | featselmethod=="spls2way" | featselmethod=="o1spls" | featselmethod=="o2spls" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="rfesvm" | featselmethod=="wilcox" | featselmethod=="ttest" | featselmethod=="pamr" | featselmethod=="ttestrepeat" | featselmethod=="poissonreg" | featselmethod=="wilcoxrepeat")
		{
			#analysismode="classification"

                #if(is.na(Ymat)==TRUE)
				{
						#classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)

						if(analysismode=="classification"){
						if(featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="poissonreg")
						{
                        				
								levels_classA<-levels(factor(classlabels[,2]))		
								for(l1 in levels_classA){
									g1<-grep(x=l1,pattern="[0-9]")
									
									if(length(g1)>0){
										#stop("Class labels or factor levels should not have any numbers.")
									}
								
								}	

						}else{
								for(c1 in 2:dim(classlabels)[2]){
								#classlabels[,c1]<-paste("x",classlabels[,c1],sep="")	

									levels_classA<-levels(factor(classlabels[,c1]))		
									for(l1 in levels_classA){
										g1<-grep(x=l1,pattern="[0-9]")
										
										if(length(g1)>0){
										#stop("Class labels or factor levels should not have any numbers.")
										}
								
								}	
	
								}
						}
						}

						
						classlabels_orig<-classlabels
                        
                        
						if(featselmethod=="limma1way"){
							
							featselmethod="limma"
						}
                        
                       
                        
						
						# | featselmethod=="limma1wayrepeat"
							if(featselmethod=="limma" | featselmethod=="limma1way" | featselmethod=="MARS" | featselmethod=="RF" | featselmethod=="pls" | featselmethod=="o1pls" | featselmethod=="o2pls" | featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="spls" | featselmethod=="o1spls" | featselmethod=="o2spls" | featselmethod=="rfesvm" | featselmethod=="pamr" | featselmethod=="poissonreg")
							{
								
								
								if(featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="poissonreg")
								{
									factor_inf<-classlabels[,-c(1)]
									factor_inf<-as.data.frame(factor_inf)
                                    #print(factor_inf)
												
								classlabels_orig<-colnames(classlabels[,-c(1)])
								colnames(classlabels)<-c("SampleID",paste("Factor",seq(1,dim(factor_inf)[2]),sep=""))
												
												#print(classlabels)
											
											
												
								Xmat_temp<-Xmat #t(Xmat)
								
								#print(Xmat_temp[1:2,1:3])
								Xmat_temp<-cbind(classlabels,Xmat_temp)
								#print("here")				
								
								Xmat_temp<-Xmat_temp[order(Xmat_temp[,2]),]
								
								cnames<-colnames(Xmat_temp)
								
								factor_lastcol<-grep("^Factor", cnames)
								
								classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
								classlabels_xyplots<-classlabels
								
								Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
								
                                #print(Xmat_temp[1:2,1:3])
							
								classlabels_response_mat<-classlabels[,-c(1)]
								
                                classlabels<-as.data.frame(classlabels)
								classlabels_response_mat<-classlabels[,-c(1)]
								
                                classlabels_response_mat<-as.data.frame(classlabels_response_mat)
									
									colnames(classlabels_response_mat)<-as.character(classlabels_orig)
									Ymat<-classlabels
                                    
                                    classlabels_orig<-classlabels
									
								}else
								{
								
								
							
								if(dim(classlabels)[2]>2){
									if(pairedanalysis==FALSE){	
									print("Invalid classlabels file format. Correct format: \nColumnA: SampleID\nColumnB: Class")
									print("Using the second column as Class")
									classlabels<-classlabels[,c(1:2)]
									}
								}

						if(analysismode=="classification")
						{
									factor_inf<-classlabels[,-c(1)]
									factor_inf<-as.data.frame(factor_inf)
												
								colnames(classlabels)<-c("SampleID",paste("Factor",seq(1,dim(factor_inf)[2]),sep=""))
												
								Xmat_temp<-Xmat #t(Xmat)
								
								Xmat_temp<-cbind(classlabels,Xmat_temp)
								#print("here")				
								Xmat_temp<-Xmat_temp[order(Xmat_temp[,2]),]
								
								cnames<-colnames(Xmat_temp)
								
								factor_lastcol<-grep("^Factor", cnames)
								
								classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
								

								Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
								classlabels_xyplots<-classlabels
                                
                                classlabels_sub<-classlabels[,-c(1)]
                                #classlabels_orig<-classlabels
                                
                                # print("DOING DOING")
                        
                        }
								
								classlabels_response_mat<-classlabels[,-c(1)]
								
								
									classlabels<-as.data.frame(classlabels)
					classlabels_response_mat<-classlabels[,-c(1)]
					classlabels_response_mat<-as.data.frame(classlabels_response_mat)
					
					#classlabels[,1]<-as.factor(classlabels[,1])
					Ymat<-classlabels
                    
                    classlabels_orig<-classlabels
								
                                }
					
					#print("here 2")
					
							}
							
								if(featselmethod=="limma1wayrepeat"){
											factor_inf<-classlabels[,-c(1:2)]
											factor_inf<-as.data.frame(factor_inf)
											
                                            # print("here")
											colnames(classlabels)<-c("SampleID","SubjectNum",paste("Factor",seq(1,length(factor_inf)),sep=""))
											
											
											#Xmat<-chocolate[,1]
											Xmat_temp<-Xmat #t(Xmat)
											Xmat_temp<-cbind(classlabels,Xmat_temp)
											
											Xmat_temp<-Xmat_temp[order(Xmat_temp[,3],Xmat_temp[,2]),]
											
											cnames<-colnames(Xmat_temp)
											
											factor_lastcol<-grep("^Factor", cnames)
											
											classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
											
											subject_inf<-classlabels[,2]
											classlabels_sub<-classlabels[,-c(1)]
											subject_inf<-subject_inf[seq(1,dim(classlabels)[1],num_replicates)]
											classlabels<-classlabels[,-c(2)]
                                            
                                            classlabels_xyplots<-classlabels
											
											Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]

											classlabels_response_mat<-classlabels[,-c(1)]
											
												classlabels<-as.data.frame(classlabels)
					classlabels_response_mat<-classlabels[,-c(1)]
					classlabels_response_mat<-as.data.frame(classlabels_response_mat)
					Ymat<-classlabels
                    
                    # classlabels_orig<-classlabels
                    #	print("classlabels is")
                    #						print(classlabels)
			
										if(featselmethod=="limma1wayrepeat"){	
										featselmethod="limma"

										}else{

											if(featselmethod=="spls1wayrepeat"){
												featselmethod="spls"
											}else{
												if(featselmethod=="pls1wayrepeat"){
                                                                                                featselmethod="pls"
                                                                                        }
											}
										}
										pairedanalysis = TRUE
	
										}
						
							
							
							if(featselmethod=="limma2way"){
								
									factor_inf<-classlabels[,-c(1)]
									factor_inf<-as.data.frame(factor_inf)
												
								colnames(classlabels)<-c("SampleID",paste("Factor",seq(1,dim(factor_inf)[2]),sep=""))
												
								Xmat_temp<-Xmat #t(Xmat)
								
                                ##save(Xmat,file="Xmat.Rda")
                                
                                ##save(classlabels,file="Xmat_classlabels.Rda")
                                
								if(dim(classlabels)[2]>2){
												levels_classA<-levels(factor(classlabels[,2]))
												print(head(classlabels))
												print(levels_classA)
										
												#if(FALSE)
												{
													if(length(levels_classA)>2){
												
                                                stop("Factor 1 can only have two levels/categories. Factor 2 can have upto 6 levels. \nPlease rearrange the factors in your classlabels file. Or use lm2wayanova option.")
													#	classtemp<-classlabels[,3]
													#	classlabels[,3]<-classlabels[,2]
													#	classlabels[,2]<-classtemp
													
                                                    #    stop("Please select lm2wayanova option for greater than 2x2 designs.")
													}
													
													levels_classA<-levels(factor(classlabels[,2]))
												}
												
												levels_classB<-levels(factor(classlabels[,3]))
												print(levels_classB)
												if(length(levels_classB)>7){
                                                    #stop("Only one of the factors can have more than 2 levels/categories. \nPlease rearrange the factors in your classlabels file or use lm2wayanova.")
											
                                                    stop("Please select lm2wayanova option for greater than 2x7 designs.")
												}	
																		
												
												Xmat_temp<-cbind(classlabels,Xmat_temp)
								
                                
                                # #save(Xmat_temp,file="Xmat_temp.Rda")
                                
                                #		print(Xmat_temp[1:10,1:10])
												
								Xmat_temp<-Xmat_temp[order(Xmat_temp[,2],Xmat_temp[,3]),]
								
                                #		print(Xmat_temp[1:10,1:10])
								cnames<-colnames(Xmat_temp)
								
								factor_lastcol<-grep("^Factor", cnames)
								
								classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
								
								Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
								classlabels_sub<-classlabels[,-c(1)]
                                
                                classlabels_response_mat<-classlabels[,-c(1)]
                                classlabels<-as.data.frame(classlabels)
                                
                                classlabels_response_mat<-as.data.frame(classlabels_response_mat)
												
												levels_classA<-levels(factor(classlabels[,2]))
												
												levels_classB<-levels(factor(classlabels[,3]))
												print(paste("Factor 1 levels: ",paste(levels_classA,collapse=","),sep=""))
												
												print(paste("Factor 2 levels: ",paste(levels_classB,collapse=","),sep=""))
										classlabels_class<-as.factor(classlabels[,2]):as.factor(classlabels[,3])
										classtable1<-table(classlabels[,2],classlabels[,3])
										
                                        #print("classtable1 is")
                                        #print(classtable1)
										classlabels_xyplots<-classlabels
                                        #classlabels_orig<-classlabels
										classlabels<-cbind(as.character(classlabels[,1]),as.character(classlabels_class))
                                        #classlabels_response_mat<-classlabels[,-c(1)]
											classlabels<-as.data.frame(classlabels)
                                            #classlabels_response_mat<-classlabels[,-c(1)]
                                            #classlabels_response_mat<-as.data.frame(classlabels_response_mat)
					Ymat<-classlabels
				
                
                #classlabels_orig<-classlabels


										
								}
								else{
									stop("Only one factor specificied in the class labels file.")			
								}
					
							}
							
							
							
									if(featselmethod=="limma2wayrepeat"){
											factor_inf<-classlabels[,-c(1:2)]
											factor_inf<-as.data.frame(factor_inf)
                                            
                                            
											
											colnames(classlabels)<-c("SampleID","SubjectNum",paste("Factor",seq(1,dim(factor_inf)[2]),sep=""))
											
											Xmat_temp<-Xmat
											if(dim(classlabels)[2]>2)
											{
												
												levels_classA<-levels(factor(classlabels[,3]))
											
												if(length(levels_classA)>2){
											
								#stop("Factor 1 can only have two levels/categories. Factor 2 can have upto 6 levels. \nPlease rearrange the factors in your classlabels file.")
												#	classtemp<-classlabels[,3]
												#	classlabels[,3]<-classlabels[,4]
												#	classlabels[,4]<-classtemp
												}
												
												levels_classA<-levels(factor(classlabels[,3]))
												
												if(length(levels_classA)>2){
                                                    #stop("Only one of the factors can have more than 2 levels/categories. \nPlease rearrange the factors in your classlabels file or use lm2wayanovarepeat.")
                                                    #stop("Please select lm2wayanova or lm2wayanovarepeat option for greater than 2x2 designs.")
                                                    stop("Factor 1 can only have two levels/categories. Factor 2 can have upto 6 levels. \nPlease rearrange the factors in your classlabels file. Or use lm2wayanova option.")
                                                }
												
												levels_classB<-levels(factor(classlabels[,4]))
												if(length(levels_classB)>7){
											#stop("Only one of the factors can have more than 2 levels/categories. \nPlease rearrange the factors in your classlabels file or use lm2wayanova.")
											
                                            stop("Please select lm2wayanovarepeat option for greater than 2x7 designs.")
												}							
											
											#Xmat<-chocolate[,1]
											 #t(Xmat)
											
											Xmat_temp<-cbind(classlabels,Xmat_temp)
											
                                            # #save(Xmat_temp,file="Xmat_temp.Rda")
                                             
                                             
											#print("here")
											
										#	print(Xmat_temp[1:2,1:4])
										#	print(dim(Xmat_temp))
											
											Xmat_temp<-Xmat_temp[order(Xmat_temp[,3],Xmat_temp[,4],Xmat_temp[,2]),]
											
										#	print(dim(Xmat_temp))
										#	print(Xmat_temp[1:2,1:4])
											#Xmat_temp<-t(Xmat_temp)
											
											cnames<-colnames(Xmat_temp)
											
											factor_lastcol<-grep("^Factor", cnames)
											
											classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
											
                                            classlabels_sub<-classlabels[,-c(1)]
											subject_inf<-classlabels[,2]
											classlabels<-classlabels[,-c(2)]
											
                                            classlabels_response_mat<-classlabels[,-c(1)]
                                            classlabels<-as.data.frame(classlabels)
                                            
                                            classlabels_response_mat<-as.data.frame(classlabels_response_mat)
                                            
                                            classlabels_xyplots<-classlabels
											subject_inf<-subject_inf[seq(1,dim(classlabels)[1],num_replicates)]
											
                                            #write.table(classlabels,file="organized_classlabelsA1.txt",sep="\t",row.names=FALSE)
											Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
											
											
                                            #write.table(Xmat_temp,file="organized_featuretableA1.txt",sep="\t",row.names=TRUE)




												
												levels_classA<-levels(factor(classlabels[,2]))
												
												levels_classB<-levels(factor(classlabels[,3]))
												print(paste("Factor 1 levels: ",paste(levels_classA,collapse=","),sep=""))
												
												print(paste("Factor 2 levels: ",paste(levels_classB,collapse=","),sep=""))
										classlabels_class<-as.factor(classlabels[,2]):as.factor(classlabels[,3])
										classtable1<-table(classlabels[,2],classlabels[,3])
										
                                        #classlabels_orig<-classlabels
										classlabels<-cbind(as.character(classlabels[,1]),as.character(classlabels_class))
										
					Ymat<-classlabels
			
            
            print("Class labels file limma2wayrep:")
              print(head(classlabels))
					#rownames(Xmat)<-as.character(classlabels[,1])
					
					#write.table(classlabels,file="organized_classlabels.txt",sep="\t",row.names=FALSE)
					
					Xmat1<-cbind(classlabels,Xmat)
                    #write.table(Xmat1,file="organized_featuretable.txt",sep="\t",row.names=TRUE)

						featselmethod="limma2way"
						pairedanalysis = TRUE
										
								}
								else{
									print(head(classlabels))
									stop("Only one factor specificied in the class labels file.")			
								}
										}
								
						
						
					}
  
				classlabels<-as.data.frame(classlabels)
				
				
	
					
				
			
			
				
				if(featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="pls2way" | featselmethod=="spls2way" | featselmethod=="wilcox" | featselmethod=="ttest"){
					
					analysismode="classification"

					classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
					
         

							#cnames[2]<-"Factor1"
														
							cnames<-colnames(classlabels)
							
							factor_inf<-classlabels[,-c(1)]
							factor_inf<-as.data.frame(factor_inf)
							
							colnames(classlabels)<-c("SampleID",paste("Factor",seq(1,dim(factor_inf)[2]),sep=""))
							
							analysismode="classification"

							Xmat_temp<-Xmat #t(Xmat)
							Xmat_temp<-cbind(classlabels,Xmat_temp)
							
                            
                            # #save(Xmat_temp,file="Xmat_temp.Rda")
                             
							if(featselmethod=="lm1wayanova" | featselmethod=="wilcox" | featselmethod=="ttest"){
						
	
								Xmat_temp<-Xmat #t(Xmat)
        							classlabels<-classlabels[,c(1:2)]
		                                                Xmat_temp<-cbind(classlabels,Xmat_temp)
	
								Xmat_temp<-Xmat_temp[order(Xmat_temp[,2]),]

								levels_classA<-levels(factor(classlabels[,2]))
								print(paste("Factor 1 levels: ",paste(levels_classA,collapse=","),sep=""))
							}else{
									if(featselmethod=="lm2wayanova" | featselmethod=="pls2way" | featselmethod=="spls2way"){
										
										Xmat_temp<-Xmat_temp[order(Xmat_temp[,2],Xmat_temp[,3]),]

                                        levels_classA<-levels(factor(classlabels[,2]))
                                        levels_classB<-levels(factor(classlabels[,3]))
										
                                        print(paste("Factor 1 levels: ",paste(levels_classA,collapse=","),sep=""))
										print(paste("Factor 2 levels: ",paste(levels_classB,collapse=","),sep=""))
								
									}
								}
							cnames<-colnames(Xmat_temp)
							
							factor_lastcol<-grep("^Factor", cnames)
							
							classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
							classlabels_sub<-classlabels[,-c(1)]

							classlabels_response_mat<-classlabels[,-c(1)]
                            
                           
							Ymat<-classlabels
                            
                             classlabels_orig<-classlabels
							
							Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
							
							if(featselmethod=="lm2wayanova" | featselmethod=="pls2way" | featselmethod=="spls2way"){
								
																
								classlabels_class<-as.factor(classlabels[,2]):as.factor(classlabels[,3])
										classtable1<-table(classlabels[,2],classlabels[,3])
										
                                        classlabels_xyplots<-classlabels
                                        #classlabels_orig<-classlabels
                                        
                                        # classlabels_orig<-classlabels_orig[seq(1,dim(classlabels)[1],num_replicates),]
										classlabels<-cbind(as.character(classlabels[,1]),as.character(classlabels_class))
										Ymat<-classlabels
										if(featselmethod=="pls2way"){
											featselmethod="pls"
										}else{

											if(featselmethod=="spls2way"){
                                                                                        featselmethod="spls"
                                                                                	}
										}
								
							}

						
                        # write.table(classlabels,file="organized_classlabelsB.txt",sep="\t",row.names=FALSE)
											Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
											
											
                                            #write.table(Xmat_temp,file="organized_featuretableA.txt",sep="\t",row.names=TRUE)
										#write.table(classlabels,file="organized_classlabelsA.txt",sep="\t",row.names=FALSE)

					
				}
				
						if(featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="pls1wayrepeat" | featselmethod=="spls1wayrepeat" | featselmethod=="pls2wayrepeat" | featselmethod=="spls2wayrepeat" | featselmethod=="ttestrepeat" | featselmethod=="wilcoxrepeat"){
							
                            #analysismode="classification"
                            pairedanalysis=TRUE
							classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)
					
							#cnames[2]<-"Factor1"
							#classlabels<-classlabels[,-c(2)]
							
							cnames<-colnames(classlabels)
							
							factor_inf<-classlabels[,-c(1:2)]
							factor_inf<-as.data.frame(factor_inf)
							
							colnames(classlabels)<-c("SampleID","SubjectNum",paste("Factor",seq(1,dim(factor_inf)[2]),sep=""))
							
							classlabels_orig<-classlabels
							#Xmat<-chocolate[,1]
							Xmat_temp<-Xmat #t(Xmat)
							Xmat_temp<-cbind(classlabels,Xmat_temp)
		
							pairedanalysis=TRUE
							if(featselmethod=="lm1wayanovarepeat" | featselmethod=="pls1wayrepeat" | featselmethod=="spls1wayrepeat" | featselmethod=="ttestrepeat" | featselmethod=="wilcoxrepeat"){
								
										Xmat_temp<-Xmat_temp[order(Xmat_temp[,3],Xmat_temp[,2]),]
							
										cnames<-colnames(Xmat_temp)
							
										factor_lastcol<-grep("^Factor", cnames)
							
										classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
							
										subject_inf<-classlabels[,2]
							
										subject_inf<-subject_inf[seq(1,dim(classlabels)[1],num_replicates)]
										
										classlabels_response_mat<-classlabels[,-c(1:2)]
										
										
                                        
                                        # classlabels_orig<-classlabels
										classlabels_sub<-classlabels[,-c(1)]
										levels_classA<-levels(factor(classlabels[,3]))
										print(paste("Factor 1 levels: ",paste(levels_classA,collapse=","),sep=""))
										classlabels<-classlabels[,-c(2)]
                                        
                                        
                                        classlabels_xyplots<-classlabels
                                        
                                        
                                        # classlabels<-classlabels[seq(1,dim(classlabels)[1],num_replicates),]
										Ymat<-classlabels
                                        print("Class labels file:")
                                        print(dim(classlabels))
                                        #					print(head(classlabels))
                                        #	print(head(classlabels_sub))
                                        
                                        Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
                                        
                                        
                                        
                                        # write.table(Xmat_temp,file="organized_featuretableA.txt",sep="\t",row.names=FALSE)

##save(Ymat,file="Ymat.Rda")
#                                       #save(Xmat,file="Xmat.Rda")

										if(featselmethod=="spls1wayrepeat"){
                                                                                                featselmethod="spls"
												
                                                                                        }else{
                                                                                                if(featselmethod=="pls1wayrepeat"){
                                                                                                featselmethod="pls"
                                                                                        }														
											}
                                                                                        
                                                                                        
                                                                                        if(featselmethod=="wilcoxrepeat"){
                                                                                            
                                                                                            featselmethod=="wilcox"
                                                                                            pairedanalysis=TRUE
                                                                                        }
                                                                                        
                                                                                        if(featselmethod=="ttestrepeat"){
                                                                                            
                                                                                            featselmethod=="ttest"
                                                                                            pairedanalysis=TRUE
                                                                                        }
							}
									if(featselmethod=="lm2wayanovarepeat" | featselmethod=="pls2wayrepeat" | featselmethod=="spls2wayrepeat"){
								
                                
										Xmat_temp<-Xmat_temp[order(Xmat_temp[,3],Xmat_temp[,4],Xmat_temp[,2]),]
                                      	cnames<-colnames(Xmat_temp)
							
										factor_lastcol<-grep("^Factor", cnames)
							
										classlabels<-Xmat_temp[,c(1:factor_lastcol[length(factor_lastcol)])]
                                        classlabels_sub<-classlabels[,-c(1)]
                                        
										subject_inf<-classlabels[,2]
																			
										subject_inf<-subject_inf[seq(1,dim(classlabels)[1],num_replicates)]
										classlabels_response_mat<-classlabels[,-c(1:2)]
										
										Ymat<-classlabels
                                        
                                        
                                        
                                        classlabels_xyplots<-classlabels[,-c(2)]


									    levels_classA<-levels(factor(classlabels[,3]))
										print(paste("Factor 1 levels: ",paste(levels_classA,collapse=","),sep=""))
									    levels_classB<-levels(factor(classlabels[,4]))
													
										print(paste("Factor 2 levels: ",paste(levels_classB,collapse=","),sep=""))
										Ymat<-classlabels
										
                                        
														#if(featselmethod!="spls2wayrepeat")
														{	
														classlabels<-classlabels[,-c(2)]
                                                        #classlabels_class<-as.factor(classlabels[,2]):as.factor(classlabels[,3])
                                                        
                                                        classlabels_class<-paste(classlabels[,2],":",classlabels[,3],sep="")
                                                        
														classtable1<-table(classlabels[,2],classlabels[,3])
														
                                                        #classlabels_orig<-classlabels
														
														classlabels<-cbind(as.character(classlabels[,1]),as.character(classlabels_class))
														}
														
														Ymat<-classlabels
			
				
                
                                            # write.table(classlabels,file="organized_classlabelsA1.txt",sep="\t",row.names=FALSE)
											Xmat<-Xmat_temp[,-c(1:factor_lastcol[length(factor_lastcol)])]
											
                                            
                                            #write.table(Xmat_temp,file="organized_featuretableA.txt",sep="\t",row.names=FALSE)
                                            #write.table(Xmat,file="organized_featuretableB1.txt",sep="\t",row.names=FALSE)
											pairedanalysis=TRUE
											if(featselmethod=="spls2wayrepeat"){
                                                                                                featselmethod="spls"

												}
									}
                                    
                                    if(featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="poissonreg"){
                                        
                                        classlabels<-as.data.frame(classlabels)
                                        classlabels_response_mat<-classlabels[,-c(1)]
                                        classlabels_response_mat<-as.data.frame(classlabels_response_mat)
                                        
                                        classlabels_response_mat[,1]<-as.factor(classlabels_response_mat)
                                        Ymat<-classlabels
                                        
                                        #classlabels_orig<-classlabels
                                        Ymat<-as.data.frame(Ymat)
                                        Xmat<-t(Xmat)
                                        Xmat<-cbind(X[,c(1:2)],Xmat)
                                        Xmat<-as.data.frame(Xmat)
                                        classlabels_xyplots<-classlabels
                                        
                                    }



						}
				
        }
				
                rnames_xmat<-rownames(Xmat)
                rnames_ymat<-as.character(Ymat[,1])
                
                
                if(length(which(duplicated(rnames_ymat)==TRUE))>0){
                    
                    stop("Duplicate sample IDs are not allowed. Please represent replicates by _1,_2,_3.")
                }
                
                    check_ylabel<-regexpr(rnames_ymat[1],pattern="^[0-9]*",perl=TRUE)
                    check_xlabel<-regexpr(rnames_xmat[1],pattern="^X[0-9]*",perl=TRUE)
                   
                   if(length(check_ylabel)>0 && length(check_xlabel)>0){
                        if(attr(check_ylabel,"match.length")>0 && attr(check_xlabel,"match.length")>0){
                
                            rnames_ymat<-paste("X",rnames_ymat,sep="") #gsub(rnames_ymat,pattern="\\.[0-9]*",replacement="")
                
                
                        }
                   }
                
                
                
			Xmat<-t(Xmat)
            
            
          
            colnames(Xmat)<-as.character(Ymat[,1])
            
			Xmat<-cbind(X[,c(1:2)],Xmat)
						
						Xmat<-as.data.frame(Xmat)
						Ymat<-as.data.frame(Ymat)
						
					
            
            match_names<-match(rnames_xmat,rnames_ymat)
            
            bad_colnames<-length(which(is.na(match_names)==TRUE))
           
	    #print(match_names) 
            #if(is.na()==TRUE){
            
            if(bad_colnames>0){
            print("Sample names do not match between feature table and class labels files.\n Please try replacing any \"-\" with \".\" in sample names.")
            print("Sample names in feature table")
            print(head(rnames_xmat))
            print("Sample names in classlabels file")
            
            print(head(rnames_ymat))
            stop("Sample names do not match between feature table and class labels files.\n Please try replacing any \"-\" with \".\" in sample names. Please try again.")
                            }
            #TIC_norm=FALSE
                            if(is.na(all(diff(match(rnames_xmat,rnames_ymat))))==FALSE){
                                if(all(diff(match(rnames_xmat,rnames_ymat)) > 0)==TRUE){
                                    
                                    data_matrix<-data_preprocess(Xmat=Xmat,Ymat=Ymat,feature_table_file=NA,parentoutput_dir,class_labels_file=NA,num_replicates=num_replicates,feat.filt.thresh,summarize.replicates,summary.method,all.missing.thresh,group.missing.thresh,
                                    log2transform,medcenter,znormtransform,quantile_norm,lowess_norm,madscaling,missing.val, samplermindex,rep.max.missing.thresh,summary.na.replacement,featselmethod,TIC_norm)
                                }else{

                                #print(diff(match(rnames_xmat,rnames_ymat)))
                                stop("Orders of feature table and classlabels do not match")
                            }


				
                            }else{
                                
                                #print(diff(match(rnames_xmat,rnames_ymat)))
                                stop("Orders of feature table and classlabels do not match")
                            }
            
           
		
		if(FALSE){
		data_matrix<-data_preprocess(Xmat,Ymat,
feature_table_file,
parentoutput_dir="C:/Users/kuppal2/Documents/Projects/EGCG_pos//xmsPANDA_preprocess3/",
class_labels_file=NA,num_replicates=1,feat.filt.thresh=NA,summarize.replicates=TRUE,
	summary.method="mean", all.missing.thresh=0.5,group.missing.thresh=0.5,
	 log2transform =FALSE, medcenter=FALSE, znormtransform = FALSE, 
    quantile_norm = FALSE, lowess_norm = FALSE, madscaling = FALSE, 
missing.val=0, samplermindex=NA,rep.max.missing.thresh=0.5,summary.na.replacement="zeros")

}
		}else{
			stop("Invalid value for analysismode parameter. Please use regression or classification.")
			}
	
	
				}


	data_matrix_beforescaling<-data_matrix$data_matrix_prescaling
  
    data_matrix_beforescaling<-as.data.frame( data_matrix_beforescaling)
   data_matrix<-data_matrix$data_matrix_afternorm_scaling
   

   
	data_m<-data_matrix[,-c(1:2)]
	

	
	Xmat<-data_matrix
	
		
			#classlabels<-as.data.frame(classlabels)
			if(dim(classlabels)[2]<2){
				
				stop("The class labels/response matrix should have two columns: SampleID, Class/Response. Please see the example.")
			}
			
			
			classlabelsA<-classlabels
			
			classlabels<-classlabels[seq(1,dim(classlabels)[1],num_replicates),]
		
        #if(dim(classlabels_orig)==TRUE){

            if(length(ncol(classlabels_sub))<1){
                classlabels_sub<-classlabels_sub[seq(1,length(classlabels_sub),num_replicates)]
			}else{
			classlabels_sub<-classlabels_sub[seq(1,dim(classlabels_sub)[1],num_replicates),]
            			}


	    
            classlabels_orig<-classlabels_orig[seq(1,dim(classlabels_orig)[1],num_replicates),]

			classlabels_response_mat<-as.data.frame(classlabels_response_mat)
			
			classlabels_response_mat<-classlabels_response_mat[seq(1,dim(classlabels_response_mat)[1],num_replicates),]
			
			
            
			class_labels_levels_main<-c("S")
			Ymat<-classlabels
			
          
			rnames1<-as.character(Ymat[,1])
			rnames2<-as.character(classlabels_orig[,1])
            
			sorted_index<-{}
			for(i in 1:length(rnames1)){


				sorted_index<-c(sorted_index,grep(x=rnames2,pattern=paste("^",rnames1[i],"$",sep="")))

			}
			classlabels_orig<-classlabels_orig[sorted_index,]

            #write.table(classlabels_response_mat,file="original_classlabelsB.txt",sep="\t",row.names=TRUE)
			classlabelsA<-classlabels
			
	
			if(length(which(duplicated(classlabels)==TRUE))>0){
				rownames(classlabels)<-paste("S",seq(1,dim(classlabels)[1]),sep="")
			}else{
				rownames(classlabels)<-as.character(classlabels[,1])
			}#as.character(classlabels[,1])
			#print(classlabels)
			#print(classlabels[1:10,])
            
 #           #save(classlabels,file="classlabels.Rda")
  #          #save(classlabels_orig,file="classlabels_orig.Rda")
   #         #save(classlabels_response_mat,file="classlabels_response_mat.Rda")
            
            if(pairedanalysis==TRUE){
                
                #save(subject_inf,file="subjectinf.Rda")
            }
            
		if(analysismode=="classification")
		{	
			
						
									class_labels_levels<-levels(as.factor(classlabels[,2]))
									
									print("Using the following class labels")			
									print(class_labels_levels)
                                    
                                    class_labels_levels_main<-class_labels_levels
									
									class_labels_levels<-unique(class_labels_levels)
						
									
									bad_rows<-which(class_labels_levels=="")
									if(length(bad_rows)>0){
									class_labels_levels<-class_labels_levels[-bad_rows]
									}
							ordered_labels={}
							num_samps_group<-new("list")
							num_samps_group[[1]]<-0
							groupwiseindex<-new("list")
							groupwiseindex[[1]]<-0
						
							for(c in 1:length(class_labels_levels))
							{
								
								classlabels_index<-which(classlabels[,2]==class_labels_levels[c])
								ordered_labels<-c(ordered_labels,as.character(classlabels[classlabels_index,2]))
								num_samps_group[[c]]<-length(classlabels_index)
								groupwiseindex[[c]]<-classlabels_index
							}
						
							Ymatorig<-classlabels

							##save(classlabels,file="classlabels_1.Rda")
							##save(class_labels_levels,file="class_labels_levels.Rda")

                            class_label_alphabets<-c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
                            if(length(class_labels_levels)>25){	
	  		     	class_label_alphabets<-paste("C",1:length(class_labels_levels),sep="")
				}		
							classlabels<-{}
							
							if(length(class_labels_levels)==2){
								#num_samps_group[[1]]=length(which(ordered_labels==class_labels_levels[1]))
								#num_samps_group[[2]]=length(which(ordered_labels==class_labels_levels[2]))
								class_label_A<-class_labels_levels[[1]]
								class_label_B<-class_labels_levels[[2]]
								classlabels<-c(rep("ClassA",num_samps_group[[1]]),rep("ClassB",num_samps_group[[2]]))
							}else{
								if(length(class_labels_levels)==3){
									
									class_label_A<-class_labels_levels[[1]]
								class_label_B<-class_labels_levels[[2]]
								class_label_C<-class_labels_levels[[3]]
									classlabels<-c(rep("ClassA",num_samps_group[[1]]),rep("ClassB",num_samps_group[[2]]),rep("ClassC",num_samps_group[[3]]))
								}else{
									
									for(c in 1:length(class_labels_levels)){
										
										num_samps_group_cur=length(which(Ymatorig[,2]==class_labels_levels[c]))
										
										classlabels<-c(classlabels,rep(paste("Class",class_label_alphabets[c],sep=""),num_samps_group_cur))
										#,rep("ClassB",num_samps_group[[2]]),rep("ClassC",num_samps_group[[3]]))
										
									}
									
									}
							}
	
	####################################################################################
	#print(head(data_m))
	
	snames<-colnames(data_m)
	
	Ymat<-as.data.frame(classlabels)
	m1<-match(snames,Ymat[,1])
	#Ymat<-Ymat[m1,]
	
	data_temp<-data_matrix_beforescaling[,-c(1:2)]
    
    
    cl<-makeSOCKcluster(num_nodes)
    
    mean_overall<-apply(data_temp,1,do_mean)
    
    #clusterExport(cl,"do_mean")
    #mean_overall<-parApply(cl,data_temp,1,do_mean)
    
    #stopCluster(cl)
	
    #mean_overall<-unlist(mean_overall)
    #print("mean overall")
    #print(summary(mean_overall))
	bad_feat<-which(mean_overall==0)

	if(length(bad_feat)>0){

		data_matrix_beforescaling<-data_matrix_beforescaling[-bad_feat,]
		data_m<-data_m[-bad_feat,]
		data_matrix<-data_matrix[-bad_feat,]
	
	} 
	
	
	#Step 5) RSD/CV calculation
	
	
						
	
	}else{
		
		classlabels<-(classlabels[,-c(1)])
		}
		
        #	print("######classlabels#########")
        #print(classlabels)
	
    
    class_labels_levels_new<-levels(classlabels)
    
    if(analysismode=="classification"){
    test_classlabels<-cbind(class_labels_levels_main,class_labels_levels_new)
    }
    
 

	
	#print("here 2")
	######################################################################################

#Step 6) Log2 mean fold change criteria from 0 to 1 with step of 0.1
feat_eval<-{}
feat_sigfdrthresh<-{}
feat_sigfdrthresh_cv<-{}
feat_sigfdrthresh_permut<-{}

permut_acc<-{}
feat_sigfdrthresh<-rep(0,length(log2.fold.change.thresh_list))
feat_sigfdrthresh_cv<-rep(NA,length(log2.fold.change.thresh_list))

feat_sigfdrthresh_permut<-rep(NA,length(log2.fold.change.thresh_list))
res_score_vec<-rep(0,length(log2.fold.change.thresh_list))
#feat_eval<-seq(0,1,0.1)

if(analysismode=="classification"){
	
	best_cv_res<-(-1)*10^30
}else{
	best_cv_res<-(1)*10^30
	
}


best_feats<-{}

goodfeats<-{}
mwan_fdr<-{}
targetedan_fdr<-{}
best_limma_res<-{}
best_acc<-{}
termA<-{}

fheader="transformed_log2fc_threshold_"


X<-t(data_m)

X<-replace(as.matrix(X),which(is.na(X)==TRUE),0)




   
  # rm(pcaMethods)
  #try(detach("package:pcaMethods",unload=TRUE),silent=TRUE)
   library(mixOmics)
 

#if(pca.global.eval==TRUE)


#"#CC0000","#AAC000",
# col_vec<-c("blue","mediumpurple4","mediumpurple1","blueviolet","cornflowerblue","cyan4","skyblue",
  #  "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
   # "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
    #"aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
    
  
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
                                            
                                            if(sample.col.opt=="colorblind"){
                                                #col_vec <-c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
                                                # col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","black")
                                                
                                                if(length(class_labels_levels)<9){
                                                    
                                                    col_vec <- c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7", "grey57")
                                                    
                                                }else{
                                                    
                                                    #col_vec<-colorRampPalette(brewer.pal(10, "RdBu"))(length(class_labels_levels))
                                                    col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35B2", "#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2",
                                                    "#374E55B2","#DF8F44B2","#00A1D5B2","#B24745B2","#79AF97B2","#6A6599B2","#80796BB2","#0073C2B2","#EFC000B2", "#868686B2","#CD534CB2","#7AA6DCB2","#003C67B2","grey57")
                                                    
                                                }
                                                
                                                
                                            }else{
                                                
                                                check_brewer<-grep(pattern="brewer",x=sample.col.opt)
                                                
                                                if(length(check_brewer)>0){
                                                    
                                                    sample.col.opt_temp=gsub(x=sample.col.opt,pattern="brewer.",replacement="")
                                                    col_vec <- colorRampPalette(brewer.pal(10, sample.col.opt_temp))(length(class_labels_levels))
                                                    
                                                }else{
                                                    
                                                    if(sample.col.opt=="journal"){
                                                        
                                                        col_vec<-c("#0072B2", "#E69F00", "#009E73", "gold1", "#56B4E9", "#D55E00", "#CC79A7","#E64B35FF","#3C5488FF","#F39B7FFF",
                                                        "#8491B4FF","#91D1C2FF","#DC0000FF","#B09C85FF","#5F559BFF",
                                                        "#808180FF","#20854EFF","#FFDC91FF","#B24745FF",
                                                        
                                                        "#374E55FF","#8F7700FF","#5050FFFF","#6BD76BFF",
                                                        "#E64B3519","#4DBBD519","#631879E5","grey75")
                                                        
                                                    }else{
                                                        #colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                                        if(length(sample.col.opt)==1){
                                                            col_vec <-rep(sample.col.opt,length(class_labels_levels))
                                                        }else{
                                                            
                                                            colfunc <-colorRampPalette(sample.col.opt);col_vec<-colfunc(length(class_labels_levels))
                                                            
                                                        }
                                                    }
                                                    
                                                }
                                                
                                            }
                                        }
									
											
									}
							
						}
			
		}	
}
#pca_col_vec<-col_vec

pca_col_vec<-c("mediumpurple4","mediumpurple1","blueviolet","darkblue","blue","cornflowerblue","cyan4","skyblue",
"darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
"red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
"aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")


if(is.na(individualsampleplot.col.opt)==TRUE){
    
    individualsampleplot.col.opt=col_vec
}




#cl<-makeSOCKcluster(num_nodes)
#feat_sds<-parApply(cl,data_m,1,sd)
feat_sds<-apply(data_m,1,sd)
			
#stopCluster(cl)
		
bad_sd_ind<-c(which(feat_sds==0),which(is.na(feat_sds)==TRUE))

bad_sd_ind<-unique(bad_sd_ind)

if(length(bad_sd_ind)>0){

	data_matrix<-data_matrix[-c(bad_sd_ind),]
	
	data_m<-data_m[-c(bad_sd_ind),]
                
        data_matrix_beforescaling<-data_matrix_beforescaling[-c(bad_sd_ind),]
}

data_temp<-data_matrix_beforescaling[,-c(1:2)]





#cl<-makeSOCKcluster(num_nodes)

#clusterExport(cl,"do_rsd")
			#feat_rsds<-parApply(cl,data_temp,1,do_rsd)
			feat_rsds<-apply(data_temp,1,do_rsd)
			#stopCluster(cl)
			
		sum_rsd<-summary(feat_rsds,na.rm=TRUE)
		max_rsd<-max(feat_rsds,na.rm=TRUE)
		max_rsd<-round(max_rsd,2)
		
			print("Summary of RSD across all features:")
			print(sum_rsd)

if(log2.fold.change.thresh_list[length(log2.fold.change.thresh_list)]>max_rsd){
stop(paste("The maximum relative standard deviation threshold in rsd.filt.list should be below ",max_rsd,sep=""))
}

parent_data_m<-round(data_m,2)
	

res_score<-0
#best_cv_res<-0
		best_feats<-{}

		best_acc<-0
		best_limma_res<-{}
		best_logfc_ind<-1

		

output_dir1<-paste(parentoutput_dir,"/Stage2/",sep="")
		dir.create(output_dir1,showWarnings=FALSE)

			setwd(output_dir1)

classlabels_sub_parent<-classlabels_sub
classlabels_orig_parent<-classlabels_orig

#write.table(classlabels_orig,file="classlabels.txt",sep="\t",row.names=FALSE)
classlabels_response_mat_parent<-classlabels_response_mat

rocfeatlist<-rocfeatlist+1

if(pairedanalysis==TRUE){
    #print(subject_inf)
write.table(subject_inf,file="subject_inf.txt",sep="\t")
paireddesign=subject_inf
}else{
    paireddesign=NA
    
    
}
write.table(classlabels_orig,file="classlabels_orig.txt",sep="\t")
#write.table(classlabels,file="classlabels.txt",sep="\t")
#write.table(classlabels_response_mat,file="classlabels_response_mat.txt",sep="\t")

if(is.na(max_varsel)==TRUE){
    
    max_varsel=dim(data_m)[1]
}

for(lf in 1:length(log2.fold.change.thresh_list))
{
		
        allmetabs_res<-{}
            classlabels_response_mat<-classlabels_response_mat_parent
			classlabels_sub<-classlabels_sub_parent
			classlabels_orig<-classlabels_orig_parent
            
			setwd(parentoutput_dir)
			log2.fold.change.thresh=log2.fold.change.thresh_list[lf]
			
			
			#print("filter is")
			#print(log2.fold.change.thresh)
            if(logistic_reg==TRUE){
			output_dir<-paste(output_dir1,"logitreg","signalthresh",group.missing.thresh,"RSD",log2.fold.change.thresh,"/",sep="")
            }else{
                
               
                if(poisson_reg==TRUE){
                    output_dir<-paste(output_dir1,"poissonreg","signalthresh",group.missing.thresh,"RSD",log2.fold.change.thresh,"/",sep="")
                }else{
                output_dir<-paste(output_dir1,featselmethod,"signalthresh",group.missing.thresh,"RSD",log2.fold.change.thresh,"/",sep="")
                }
                
            }
            
            dir.create(output_dir,showWarnings=FALSE)

			setwd(output_dir)
            
            dir.create("Figures")
             
            dir.create("Tables")
	    
	    data_m<-parent_data_m
			
            #print("dim of data_m")
            #print(dim(data_m))
			
			pdf_fname<-paste("Figures/Results_RSD",log2.fold.change.thresh,".pdf",sep="")
			
			#zip_fname<-paste("Results_RSD",log2.fold.change.thresh,".zip",sep="")
			
            if(output.device.type=="pdf"){
                pdf(pdf_fname,width=10,height=10)
            }
			
			if(analysismode=="classification" | analysismode=="regression"){
				
				print(paste("Performing RSD filtering using ",log2.fold.change.thresh, " as threshold",sep=""))
				if(log2.fold.change.thresh>=0){

					if(log2.fold.change.thresh==0){
						log2.fold.change.thresh=0.001
					}
                    
                    
					#good_metabs<-which(abs(mean_groups)>log2.fold.change.thresh)
					abs_feat_rsds<-abs(feat_rsds)
					
					good_metabs<-which(abs_feat_rsds>log2.fold.change.thresh)
					
					#print("length of good_metabs")
					#print(good_metabs)
					
				}else{
					good_metabs<-seq(1,dim(data_m)[1])
					}
				if(length(good_metabs)>0){
						
						data_m_fc<-data_m[good_metabs,]

						data_m_fc_withfeats<-data_matrix[good_metabs,c(1:2)]
				
						data_matrix_beforescaling_rsd<-data_matrix_beforescaling[good_metabs,]
				
				}else{
					#data_m_fc<-data_m
					#data_m_fc_withfeats<-data_matrix[,c(1:2)]

					stop(paste("Please decrease the maximum relative standard deviation (rsd.filt.thresh) threshold to ",max_rsd,sep=""))
					
					}
			}else{
				
				data_m_fc<-data_m
				data_m_fc_withfeats<-data_matrix[,c(1:2)]
			}
			
			#print("here 3")
			
            if(log2transform==TRUE || input.intensity.scale=="log2"){
                
                if(znormtransform==TRUE){
                    ylab_text_2="scale normalized"
                }else{
                    if(quantile_norm==TRUE){
                        
                        ylab_text_2="quantile normalized"
                    }else{
                        ylab_text_2=""
                    }
                }
                ylab_text=paste("log2 intensity ",ylab_text_2,sep="")
            }else{
                if(znormtransform==TRUE){
                    ylab_text_2="scale normalized"
                }else{
                    if(quantile_norm==TRUE){
                        
                        ylab_text_2="quantile normalized"
                    }else{
                        ylab_text_2=""
                    }
                }
                ylab_text=paste("Raw intensity ",ylab_text_2,sep="")
            }
            
            print(cex.plots)
            
            if(dim(data_m_fc)[2]>50){
                
                if(output.device.type!="pdf"){
                    
                    temp_filename_1<-"Figures/SampleIntensityDistribution.png"
                    
                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                
                }
                
                size_num<-min(100,dim(data_m_fc)[2])
                
                par(mfrow=c(1,1),family="sans",cex=cex.plots)
                samp_index<-sample(x=1:dim(data_m_fc)[2],size=size_num)
                try(boxplot(data_m_fc[,samp_index],main="Intensity distribution across samples after preprocessing",xlab="Samples",ylab=ylab_text,col=boxplot.col.opt),silent=TRUE)
                
                if(output.device.type!="pdf"){
                    
                    try(dev.off(),silent=TRUE)
                }
                
            }else{
                
                if(output.device.type!="pdf"){
                    
                    temp_filename_1<-"Figures/SampleIntensityDistribution.png"
                    
                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                }
                par(mfrow=c(1,1),family="sans",cex=cex.plots)
                try(boxplot(data_m_fc,main="Intensity distribution across samples after preprocessing",xlab="Samples",ylab=ylab_text,col=boxplot.col.opt),silent=TRUE)
                if(output.device.type!="pdf"){
                    
                    try(dev.off(),silent=TRUE)
                }
                
            }
           
            
            
			data_m_fc_withfeats<-cbind(data_m_fc_withfeats,data_m_fc)
			
			
			feat_eval[lf]<-0
			res_score_vec[lf]<-0
			#feat_sigfdrthresh_cv[lf]<-0
			
			filename<-paste(fheader,log2.fold.change.thresh,".txt",sep="")
			#write.table(data_m_fc_withfeats, file=filename,sep="\t",row.names=FALSE)

			
			if(length(data_m_fc)>=dim(parent_data_m)[2])
			{
					
					
					if(dim(data_m_fc)[1]>0){
					
					

						feat_eval[lf]<-dim(data_m_fc)[1]
						
						# col_vec<-c("#CC0000","#AAC000","blue","mediumpurple4","mediumpurple1","blueviolet","darkblue","blue","cornflowerblue","cyan4","skyblue",
					    #"darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
					    #"red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
					    #"aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
					    
						if(analysismode=="classification")
						{
                                                    
                                                    sampleclass<-{}
                                                    patientcolors<-{}
                                                    #print(classlabels)
                                                    classlabels<-as.data.frame(classlabels)
                                                    f<-factor(classlabels[,1])
                                                    
                                                    for(c in 1:length(class_labels_levels)){
                                                        
                                                                num_samps_group_cur=length(which(ordered_labels==class_labels_levels[c]))
                                                                
                                                                #classlabels<-c(classlabels,rep(paste("Class",class_label_alphabets,sep=""),num_samps_group_cur))
                                                                #,rep("ClassB",num_samps_group[[2]]),rep("ClassC",num_samps_group[[3]]))
                                                                sampleclass<-c(sampleclass,rep(paste("Class",class_label_alphabets[c],sep=""),num_samps_group_cur))

                                                                patientcolors <-c(patientcolors,rep(col_vec[c],num_samps_group_cur))
                                                            }




                                                                                   library(pcaMethods)

                                #p1<-pcaMethods::pca(data_m_fc,method="rnipals",center=TRUE,scale="uv",cv="q2",nPcs=3)

                                tempX<-t(data_m_fc)



#  p1<-pcaMethods::pca(tempX,method="rnipals",center=TRUE,scale="uv",cv="q2",nPcs=10)



                                if(output.device.type!="pdf"){

                                temp_filename_2<-"Figures/PCAdiagnostics_allfeats.png"

				# png(temp_filename_2,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                }



                                if(output.device.type!="pdf"){

# dev.off()
                                }

                                 try(detach("package:pcaMethods",unload=TRUE),silent=TRUE)



                                if(dim(classlabels)[2]>2){
                                    classgroup<-paste(classlabels[,1],":",classlabels[,2],sep="") #classlabels[,1]:classlabels[,2]
                                }else{
                                    
                                    classgroup<-classlabels
                                }

                                classlabels_orig<-classlabels_orig_parent
                              
                                if(pairedanalysis==TRUE){

                                        classlabels_orig<-classlabels_orig[,-c(2)]


                                }else{

                                        if(featselmethod=="lmreg" || featselmethod=="logitreg" ||  featselmethod=="poissonreg"){
                                            classlabels_orig<-classlabels_orig[,c(1:2)]
                                            classlabels_orig<-as.data.frame(classlabels_orig)
                                        }
                                }



                                if(output.device.type!="pdf"){

                                temp_filename_1<-"Figures/PCAplots_allfeats.pdf"

                                #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")

                                pdf(temp_filename_1)
                                }


                            #save(list=ls(),file="pcaplotsall.Rda")
                               
                               
                               if(FALSE)
                               { try(get_pcascoredistplots(X=data_m_fc_withfeats,Y=classlabels_orig,feature_table_file=NA,parentoutput_dir=getwd(),class_labels_file=NA,sample.col.opt=sample.col.opt,
                                  plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,col_vec=col_vec,pairedanalysis=pairedanalysis,pca.cex.val=pca.cex.val,legendlocation=legendlocation, pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,filename="all",paireddesign=paireddesign),silent=TRUE)
                               }
                              
                                      get_pcascoredistplots(X=data_m_fc_withfeats,Y=classlabels_orig,feature_table_file=NA,parentoutput_dir=getwd(),class_labels_file=NA,sample.col.opt=sample.col.opt,
                   plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,col_vec=col_vec,pairedanalysis=pairedanalysis,pca.cex.val=pca.cex.val,legendlocation=legendlocation,
                   pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,filename="all",paireddesign=paireddesign,lineplot.col.opt=lineplot.col.opt,lineplot.lty.option=lineplot.lty.option)
                   
                   
                                  
                                if(output.device.type!="pdf"){

                                    try(dev.off(),silent=TRUE)
                                }


                                classlabels_orig<-classlabels_orig_parent
                   }else{
                       #regression
                            tempgroup<-rep("A",dim(data_m_fc)[2])
                            col_vec1<-rep("black",dim(data_m_fc)[2])
                            class_labels_levels_main1<-c("A")
                            
                            
                            
                                                     try(get_pca(X=data_m_fc,samplelabels=tempgroup,legendlocation=legendlocation,filename="all",ncomp=3,center=TRUE,scale=TRUE,legendcex=0.5,outloc=getwd(),col_vec=col_vec1,sample.col.opt=sample.col.opt,alphacol=0.3,class_levels=NA,pca.cex.val=pca.cex.val,pca.ellipse=FALSE,paireddesign=paireddesign),silent=TRUE)
                            
                            
                            
                        }


if(featselmethod=="pamr"){
    
    #print("HERE")
    ##save(data_m_fc,classlabels,file="pamdebug.Rda")
    
    pamr.res<-do_pamr(X=data_m_fc,Y=classlabels,fdrthresh=fdrthresh,nperms=pls.permut.count,pamr.threshold.select.max=pamr.threshold.select.max,kfold=kfold)
    
    goodip<-pamr.res$feature.list
    
    pamr.threshold_value<-pamr.res$threshold_value
   
    feature_rowindex<-seq(1,nrow(data_m_fc))
    
    discore<-rep(0,nrow(data_m_fc))
    
    discore_all<-pamr.res$max.discore.allfeats
    
    if(is.na(goodip)==FALSE){
        discore[goodip]<-pamr.res$max.discore.sigfeats
    
        sel.diffdrthresh<-feature_rowindex%in%goodip
        
        max_absolute_standardized_centroids_thresh0<-pamr.res$max.discore.allfeats[goodip]
        
        selected_id_withmztime<-cbind(data_m_fc_withfeats[goodip,c(1:2)],pamr.res$pam_toplist,max_absolute_standardized_centroids_thresh0)
        #save(pamr.res,file="pamr.res.Rda")
        write.csv(selected_id_withmztime,file="dscores.selectedfeats.csv",row.names=FALSE)
        
        rank_vec<-rank(-discore_all)
        
        max_absolute_standardized_centroids_thresh0<-pamr.res$max.discore.allfeats
        
        data_limma_fdrall_withfeats<-cbind(max_absolute_standardized_centroids_thresh0,data_m_fc_withfeats)
        write.table(data_limma_fdrall_withfeats,file="Tables/pamr_ranked_feature_table.txt",sep="\t",row.names=FALSE)
        
	
    }else{
        goodip<-{}
        sel.diffdrthresh<-rep(FALSE,length(feature_rowindex))
    }
    
    rank_vec<-rank(-discore_all)
    
    pamr_ythresh<-pamr.res$max.discore.all.thresh-0.00000001


    
}


if(featselmethod=="rfesvm"){
 
    svm_classlabels<-classlabels[,1]
    
    if(analysismode=="classification"){
    svm_classlabels<-as.data.frame(svm_classlabels)
    }
    
  
    if(length(class_labels_levels)<3){
        featureRankedList = diffexpsvmrfe(x=t(data_m_fc),y=svm_classlabels,svmkernel=svm_kernel)
        
        #best_subset<-featureRankedList$best_subset
        
    }else{
        
        featureRankedList = diffexpsvmrfemulticlass(x=t(data_m_fc),y=svm_classlabels,svmkernel=svm_kernel)
        
        
    }
    

    rank_vec<-seq(1,dim(data_m_fc_withfeats)[1])
    goodip<-featureRankedList[1:max_varsel]
    #dtemp1<-data_m_fc_withfeats[goodip,]

    sel.diffdrthresh<-rank_vec%in%goodip
    
   
    rank_vec<-sort(featureRankedList,index.return=TRUE)$ix
    
    data_limma_fdrall_withfeats<-cbind(rank_vec,data_m_fc_withfeats)
    
    
}

						
								
										if(featselmethod=="limma" | featselmethod=="limma1way")
										{
														
												design <- model.matrix(~ -1+f)
												 colnames(design) <- levels(f)

											
											options(digit=3)
											parameterNames<-colnames(design)
											if(dim(design)[2]==2){
											cont.matrix <- makeContrasts(Grp1vs2="ClassA-ClassB",levels=design)
											
											design_mat_names<-c(paste(class_labels_levels[1],"vs",class_labels_levels[2],sep=""))
												
											}else{
											if(dim(design)[2]==3){
											
											cont.matrix <- makeContrasts(Grp1vs2="ClassA-ClassB",Grp1vs3="ClassA-ClassC",Grp2vs3="ClassB-ClassC",levels=design)
											
											design_mat_names<-c(paste(class_labels_levels[1],"vs",class_labels_levels[2],sep=""),
												paste(class_labels_levels[1],"vs",class_labels_levels[3],sep=""),
												paste(class_labels_levels[2],"vs",class_labels_levels[3],sep=""))
												
											}else{
												if(dim(design)[2]==4){
											
												cont.matrix <- makeContrasts(Grp1vs2="ClassA-ClassB",Grp1vs3="ClassA-ClassC",Grp1vs4="ClassA-ClassD",Grp2vs3="ClassB-ClassC",Grp2vs4="ClassB-ClassD",
												Grp3vs4="ClassC-ClassD",levels=design)
												
												design_mat_names<-c(paste(class_labels_levels[1],"vs",class_labels_levels[2],sep=""),
												paste(class_labels_levels[1],"vs",class_labels_levels[3],sep=""),
												paste(class_labels_levels[1],"vs",class_labels_levels[4],sep=""),
												paste(class_labels_levels[2],"vs",class_labels_levels[3],sep=""),
												paste(class_labels_levels[2],"vs",class_labels_levels[4],sep=""),
												paste(class_labels_levels[3],"vs",class_labels_levels[4],sep=""))
												
												}else{
													if(dim(design)[2]==5){
											
												cont.matrix <- makeContrasts(Grp1vs2="ClassA-ClassB",Grp1vs3="ClassA-ClassC",Grp1vs4="ClassA-ClassD",Grp1vs5="ClassA-ClassE",
												Grp2vs3="ClassB-ClassC",Grp2vs4="ClassB-ClassD",Grp2vs5="ClassB-ClassE",
												Grp3vs4="ClassC-ClassD",Grp3vs5="ClassC-ClassE",
												Grp4vs5="ClassD-ClassE",levels=design)
												
												design_mat_names<-c(paste(class_labels_levels[1],"vs",class_labels_levels[2],sep=""),
												paste(class_labels_levels[1],"vs",class_labels_levels[3],sep=""),
												paste(class_labels_levels[1],"vs",class_labels_levels[4],sep=""),
												paste(class_labels_levels[1],"vs",class_labels_levels[5],sep=""),
												paste(class_labels_levels[2],"vs",class_labels_levels[3],sep=""),
												paste(class_labels_levels[2],"vs",class_labels_levels[4],sep=""),
												paste(class_labels_levels[2],"vs",class_labels_levels[5],sep=""),
												paste(class_labels_levels[3],"vs",class_labels_levels[4],sep=""),
												paste(class_labels_levels[3],"vs",class_labels_levels[5],sep=""),
												paste(class_labels_levels[4],"vs",class_labels_levels[5],sep="")
												)
												
												}else{
													if(dim(design)[2]==6){
											
												cont.matrix <- makeContrasts(Grp1vs2="ClassA-ClassB",Grp1vs3="ClassA-ClassC",Grp1vs4="ClassA-ClassD",Grp1vs5="ClassA-ClassE",Grp1vs6="ClassA-ClassF",
												Grp2vs3="ClassB-ClassC",Grp2vs4="ClassB-ClassD",Grp2vs5="ClassB-ClassE",Grp2vs6="ClassB-ClassF",
												Grp3vs4="ClassC-ClassD",Grp3vs5="ClassC-ClassE",Grp3vs6="ClassC-ClassF",
												Grp4vs5="ClassD-ClassE",Grp4vs6="ClassD-ClassF",
												Grp5vs6="ClassE-ClassF",levels=design)
												
												design_mat_names<-c(paste(class_labels_levels[1],"vs",class_labels_levels[2],sep=""),
												paste(class_labels_levels[1],"vs",class_labels_levels[3],sep=""),
												paste(class_labels_levels[1],"vs",class_labels_levels[4],sep=""),
												paste(class_labels_levels[1],"vs",class_labels_levels[5],sep=""),
												paste(class_labels_levels[1],"vs",class_labels_levels[6],sep=""),
												paste(class_labels_levels[2],"vs",class_labels_levels[3],sep=""),
												paste(class_labels_levels[2],"vs",class_labels_levels[4],sep=""),
												paste(class_labels_levels[2],"vs",class_labels_levels[5],sep=""),
												paste(class_labels_levels[2],"vs",class_labels_levels[6],sep=""),
												paste(class_labels_levels[3],"vs",class_labels_levels[4],sep=""),
												paste(class_labels_levels[3],"vs",class_labels_levels[5],sep=""),
												paste(class_labels_levels[3],"vs",class_labels_levels[6],sep=""),
												paste(class_labels_levels[4],"vs",class_labels_levels[5],sep=""),
												paste(class_labels_levels[4],"vs",class_labels_levels[6],sep=""),
												paste(class_labels_levels[5],"vs",class_labels_levels[6],sep="")
												)
												
												}else{
                                                    
                                                        if(dim(design)[2]==7){
                                                            
                                                            cont.matrix <- makeContrasts(Grp1vs2="ClassA-ClassB",Grp1vs3="ClassA-ClassC",Grp1vs4="ClassA-ClassD",Grp1vs5="ClassA-ClassE",Grp1vs6="ClassA-ClassF",Grp1vs7="ClassA-ClassG",
                                                            Grp2vs3="ClassB-ClassC",Grp2vs4="ClassB-ClassD",Grp2vs5="ClassB-ClassE",Grp2vs6="ClassB-ClassF",Grp2vs7="ClassB-ClassG",
                                                            Grp3vs4="ClassC-ClassD",Grp3vs5="ClassC-ClassE",Grp3vs6="ClassC-ClassF",Grp3vs7="ClassC-ClassG",
                                                            Grp4vs5="ClassD-ClassE",Grp4vs6="ClassD-ClassF",Grp4vs7="ClassD-ClassG",
                                                            Grp5vs6="ClassE-ClassF",Grp5vs7="ClassE-ClassG",Grp6vs7="ClassF-ClassG",levels=design)
                                                            
                                                            design_mat_names<-c(paste(class_labels_levels[1],"vs",class_labels_levels[2],sep=""),
                                                            paste(class_labels_levels[1],"vs",class_labels_levels[3],sep=""),
                                                            paste(class_labels_levels[1],"vs",class_labels_levels[4],sep=""),
                                                            paste(class_labels_levels[1],"vs",class_labels_levels[5],sep=""),
                                                            paste(class_labels_levels[1],"vs",class_labels_levels[6],sep=""),
                                                            paste(class_labels_levels[1],"vs",class_labels_levels[7],sep=""),
                                                            paste(class_labels_levels[2],"vs",class_labels_levels[3],sep=""),
                                                            paste(class_labels_levels[2],"vs",class_labels_levels[4],sep=""),
                                                            paste(class_labels_levels[2],"vs",class_labels_levels[5],sep=""),
                                                            paste(class_labels_levels[2],"vs",class_labels_levels[6],sep=""),
                                                            paste(class_labels_levels[2],"vs",class_labels_levels[7],sep=""),
                                                            paste(class_labels_levels[3],"vs",class_labels_levels[4],sep=""),
                                                            paste(class_labels_levels[3],"vs",class_labels_levels[5],sep=""),
                                                            paste(class_labels_levels[3],"vs",class_labels_levels[6],sep=""),
                                                            paste(class_labels_levels[3],"vs",class_labels_levels[7],sep=""),
                                                            paste(class_labels_levels[4],"vs",class_labels_levels[5],sep=""),
                                                            paste(class_labels_levels[4],"vs",class_labels_levels[6],sep=""),
                                                            paste(class_labels_levels[4],"vs",class_labels_levels[7],sep=""),
                                                            paste(class_labels_levels[5],"vs",class_labels_levels[6],sep=""),
                                                            paste(class_labels_levels[5],"vs",class_labels_levels[7],sep=""),
                                                            paste(class_labels_levels[6],"vs",class_labels_levels[7],sep="")
                                                            )
                                                            
                                                        
                                                    }else{
													
                                                    if(dim(design)[2]==8){
                                                        
                                                        cont.matrix <- makeContrasts(Grp1vs2="ClassA-ClassB",Grp1vs3="ClassA-ClassC",Grp1vs4="ClassA-ClassD",Grp1vs5="ClassA-ClassE",Grp1vs6="ClassA-ClassF",Grp1vs7="ClassA-ClassG",Grp1vs8="ClassA-ClassH",
                                                        Grp2vs3="ClassB-ClassC",Grp2vs4="ClassB-ClassD",Grp2vs5="ClassB-ClassE",Grp2vs6="ClassB-ClassF",Grp2vs7="ClassB-ClassG",Grp2vs8="ClassB-ClassH",
                                                        Grp3vs4="ClassC-ClassD",Grp3vs5="ClassC-ClassE",Grp3vs6="ClassC-ClassF",Grp3vs7="ClassC-ClassG",Grp3vs8="ClassC-ClassH",
                                                        Grp4vs5="ClassD-ClassE",Grp4vs6="ClassD-ClassF",Grp4vs7="ClassD-ClassG",Grp4vs8="ClassD-ClassH",
                                                        Grp5vs6="ClassE-ClassF",Grp5vs7="ClassE-ClassG",Grp5vs8="ClassE-ClassH",
                                                        Grp6vs7="ClassF-ClassG",Grp6vs8="ClassF-ClassH",Grp7vs8="ClassG-ClassH",levels=design)
                                                        
                                                        design_mat_names<-c(paste(class_labels_levels[1],"vs",class_labels_levels[2],sep=""),
                                                        paste(class_labels_levels[1],"vs",class_labels_levels[3],sep=""),
                                                        paste(class_labels_levels[1],"vs",class_labels_levels[4],sep=""),
                                                        paste(class_labels_levels[1],"vs",class_labels_levels[5],sep=""),
                                                        paste(class_labels_levels[1],"vs",class_labels_levels[6],sep=""),
                                                        paste(class_labels_levels[1],"vs",class_labels_levels[7],sep=""),
                                                        paste(class_labels_levels[1],"vs",class_labels_levels[8],sep=""),
                                                        paste(class_labels_levels[2],"vs",class_labels_levels[3],sep=""),
                                                        paste(class_labels_levels[2],"vs",class_labels_levels[4],sep=""),
                                                        paste(class_labels_levels[2],"vs",class_labels_levels[5],sep=""),
                                                        paste(class_labels_levels[2],"vs",class_labels_levels[6],sep=""),
                                                        paste(class_labels_levels[2],"vs",class_labels_levels[7],sep=""),
                                                        paste(class_labels_levels[2],"vs",class_labels_levels[8],sep=""),
                                                        paste(class_labels_levels[3],"vs",class_labels_levels[4],sep=""),
                                                        paste(class_labels_levels[3],"vs",class_labels_levels[5],sep=""),
                                                        paste(class_labels_levels[3],"vs",class_labels_levels[6],sep=""),
                                                        paste(class_labels_levels[3],"vs",class_labels_levels[7],sep=""),
                                                        paste(class_labels_levels[3],"vs",class_labels_levels[8],sep=""),
                                                        paste(class_labels_levels[4],"vs",class_labels_levels[5],sep=""),
                                                        paste(class_labels_levels[4],"vs",class_labels_levels[6],sep=""),
                                                        paste(class_labels_levels[4],"vs",class_labels_levels[7],sep=""),
                                                        paste(class_labels_levels[4],"vs",class_labels_levels[8],sep=""),
                                                        paste(class_labels_levels[5],"vs",class_labels_levels[6],sep=""),
                                                        paste(class_labels_levels[5],"vs",class_labels_levels[7],sep=""),
                                                        paste(class_labels_levels[5],"vs",class_labels_levels[8],sep=""),
                                                        paste(class_labels_levels[6],"vs",class_labels_levels[7],sep=""),
                                                        paste(class_labels_levels[6],"vs",class_labels_levels[8],sep=""),
                                                        paste(class_labels_levels[7],"vs",class_labels_levels[8],sep="")
                                                        )
                                                        
                                                        
                                                    }
                                                    else{
                                                        if(dim(design)[2]==9){
                                                            
                                                            cont.matrix <- makeContrasts(Grp1vs2="ClassA-ClassB",Grp1vs3="ClassA-ClassC",Grp1vs4="ClassA-ClassD",Grp1vs5="ClassA-ClassE",Grp1vs6="ClassA-ClassF",Grp1vs7="ClassA-ClassG",Grp1vs8="ClassA-ClassH",Grp1vs9="ClassA-ClassI",
                                                            Grp2vs3="ClassB-ClassC",Grp2vs4="ClassB-ClassD",Grp2vs5="ClassB-ClassE",Grp2vs6="ClassB-ClassF",Grp2vs7="ClassB-ClassG",Grp2vs8="ClassB-ClassH",Grp2vs9="ClassB-ClassI",
                                                            Grp3vs4="ClassC-ClassD",Grp3vs5="ClassC-ClassE",Grp3vs6="ClassC-ClassF",Grp3vs7="ClassC-ClassG",Grp3vs8="ClassC-ClassH",Grp3vs9="ClassC-ClassI",
                                                            Grp4vs5="ClassD-ClassE",Grp4vs6="ClassD-ClassF",Grp4vs7="ClassD-ClassG",Grp4vs8="ClassD-ClassH",Grp4vs9="ClassD-ClassI",
                                                            Grp5vs6="ClassE-ClassF",Grp5vs7="ClassE-ClassG",Grp5vs8="ClassE-ClassH",Grp5vs9="ClassE-ClassI",
                                                            Grp6vs7="ClassF-ClassG",Grp6vs8="ClassF-ClassH",Grp6vs9="ClassF-ClassI",Grp7vs8="ClassG-ClassH",Grp7vs9="ClassG-ClassI",Grp8vs9="ClassH-ClassI",levels=design)
                                                            
                                                            design_mat_names<-c(paste(class_labels_levels[1],"vs",class_labels_levels[2],sep=""),
                                                            paste(class_labels_levels[1],"vs",class_labels_levels[3],sep=""),
                                                            paste(class_labels_levels[1],"vs",class_labels_levels[4],sep=""),
                                                            paste(class_labels_levels[1],"vs",class_labels_levels[5],sep=""),
                                                            paste(class_labels_levels[1],"vs",class_labels_levels[6],sep=""),
                                                            paste(class_labels_levels[1],"vs",class_labels_levels[7],sep=""),
                                                            paste(class_labels_levels[1],"vs",class_labels_levels[8],sep=""),
                                                            paste(class_labels_levels[1],"vs",class_labels_levels[9],sep=""),
                                                            paste(class_labels_levels[2],"vs",class_labels_levels[3],sep=""),
                                                            paste(class_labels_levels[2],"vs",class_labels_levels[4],sep=""),
                                                            paste(class_labels_levels[2],"vs",class_labels_levels[5],sep=""),
                                                            paste(class_labels_levels[2],"vs",class_labels_levels[6],sep=""),
                                                            paste(class_labels_levels[2],"vs",class_labels_levels[7],sep=""),
                                                            paste(class_labels_levels[2],"vs",class_labels_levels[8],sep=""),
                                                            paste(class_labels_levels[2],"vs",class_labels_levels[9],sep=""),
                                                            paste(class_labels_levels[3],"vs",class_labels_levels[4],sep=""),
                                                            paste(class_labels_levels[3],"vs",class_labels_levels[5],sep=""),
                                                            paste(class_labels_levels[3],"vs",class_labels_levels[6],sep=""),
                                                            paste(class_labels_levels[3],"vs",class_labels_levels[7],sep=""),
                                                            paste(class_labels_levels[3],"vs",class_labels_levels[8],sep=""),
                                                            paste(class_labels_levels[3],"vs",class_labels_levels[9],sep=""),
                                                            paste(class_labels_levels[4],"vs",class_labels_levels[5],sep=""),
                                                            paste(class_labels_levels[4],"vs",class_labels_levels[6],sep=""),
                                                            paste(class_labels_levels[4],"vs",class_labels_levels[7],sep=""),
                                                            paste(class_labels_levels[4],"vs",class_labels_levels[8],sep=""),
                                                            paste(class_labels_levels[4],"vs",class_labels_levels[9],sep=""),
                                                            paste(class_labels_levels[5],"vs",class_labels_levels[6],sep=""),
                                                            paste(class_labels_levels[5],"vs",class_labels_levels[7],sep=""),
                                                            paste(class_labels_levels[5],"vs",class_labels_levels[8],sep=""),
                                                            paste(class_labels_levels[5],"vs",class_labels_levels[9],sep=""),
                                                            paste(class_labels_levels[6],"vs",class_labels_levels[7],sep=""),
                                                            paste(class_labels_levels[6],"vs",class_labels_levels[8],sep=""),
                                                            paste(class_labels_levels[6],"vs",class_labels_levels[9],sep=""),
                                                            paste(class_labels_levels[7],"vs",class_labels_levels[8],sep=""),
                                                            paste(class_labels_levels[7],"vs",class_labels_levels[9],sep=""),
                                                            paste(class_labels_levels[8],"vs",class_labels_levels[9],sep="")
                                                            )
                                                            
                                                            
                                                        }
                                                        else{
                                                            stop("Please use other featselmethod options for more than 7 classes.")
                                                        }
                                                        }
                                                    }
                                                    
                                                    }
												}
												}
												}
											}
											
											# Ordinary fit
											
											
											if(pairedanalysis==TRUE){
											
											
												f1<-{}
												for(c in 1:length(class_labels_levels)){
												
													f1<-c(f1,seq(1,num_samps_group[[c]]))
													
															
															
												}
												
												print("Paired samples order")
												
												
												f1<-subject_inf
                                                						
												print(subject_inf)
												print("Design matrix")
												print(design)
                                                print(cont.matrix)
												#print(data_m[1:10,1:10])
											
                                                ##save(list=ls(),file="limma.Rda")
												##save(f,file="f.Rda")
												##save(design,file="design.Rda")
												##save(data_m_fc,file="data_m_fc.Rda")
												##save(subject_inf,file="subject_inf.Rda")
												
												corfit<-duplicateCorrelation(data_m_fc,design=design,block=subject_inf,ndups=1)
												
                                               

												fit<-lmFit(data_m_fc,design,block=f1,cor=corfit$consensus)

											}else{
											
                                            ##save(list=ls(),file="d.Rda")
												fit <- lmFit(data_m_fc,design)
											
												#fit<-treat(fit,lfc=lf)
												
											}
											
											
											#print(data_m_fc[1:3,])
											fit2  <- contrasts.fit(fit, cont.matrix)

											fit2 <- eBayes(fit2)
											#as.data.frame(fit2[1:10,])

											# Various ways of summarising or plotting the results
											#topTable(fit,coef=2)
                                            
                                        
                                            #write.table(t1,file="topTable_limma.txt",sep="\t")
                                            

											if(dim(design)[2]>2){
												pvalues<-fit2$F.p.value
												p.value<-fit2$F.p.value
												
												}else{
											pvalues<-fit2$p.value
											p.value<-fit2$p.value
											}
											
											if(fdrmethod=="BH"){
											fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
											}else{
													if(fdrmethod=="ST"){
                                                        #fdr_adjust_pvalue<-qvalue(pvalues)
                                                        #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                        
                                                        fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
                                                        
                                                        if(is(fdr_adjust_pvalue,"try-error")){
                                                            
                                                            fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
                                                        }
                                                        
                                                        fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                        
                                                        
													}else{
															if(fdrmethod=="Strimmer"){
																pdf("fdrtool.pdf")
																#par_rows=1
																#par(mfrow=c(par_rows,1))
															fdr_adjust_pvalue<-fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE)
															fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
																try(dev.off(),silent=TRUE)
															}else{
																if(fdrmethod=="none"){
																		fdr_adjust_pvalue<-pvalues
                                                                        #fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
																}else{
																	if(fdrmethod=="BY"){
                                                                        fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                                                                    }else{
                                                                        if(fdrmethod=="bonferroni"){
                                                                            fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                                                                        }
                                                                    }
                                                                    
                                                                    
																	}
															}
														}
													
												}
												
												if(dim(design)[2]<3){
								
										if(fdrmethod=="none"){
										filename<-"limma_pvalall_withfeats.txt"
										}else{
										filename<-"limma_fdrall_withfeats.txt"
										}
										cnames_tab<-colnames(data_m_fc_withfeats)
										cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
										
                                        
										data_limma_fdrall_withfeats<-cbind(p.value,fdr_adjust_pvalue,data_m_fc_withfeats)

										colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
                                        
                                        pvalues<-p.value
										
                                        #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
                                        #write.table(data_limma_fdrall_withfeats,file=filename,sep="\t",row.names=FALSE)
                                        
                                      
                                        final.pvalues<-pvalues
                                        sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
                                        
                                       
                                        
                                        goodip<-which(sel.diffdrthresh==TRUE)
                                        d4<-as.data.frame(data_limma_fdrall_withfeats)
                                        logp<-(-1)*log((d4[,1]+(10^-20)),10)
                                        
                                        #tiff("pval_dist.tiff",compression="lzw")
                                        #hist(d4[,1],xlab="p",main="Distribution of p-values")
                                        #dev.off()

								
								}else{
									
									
									adjusted.P.value<-fdr_adjust_pvalue
                                    if(limmadecideTests==TRUE){
                                    				results2<-decideTests(fit2,method="nestedF",adjust.method="BH",p.value=fdrthresh)
									#tiff("comparison_contrast_overlap.tiff",width=plots.width,height=plots.height,res=plots.res, compression="lzw")
                                    if(length(class_labels_levels)<4){
                                        
                                        if(output.device.type!="pdf"){
                                            
                                            temp_filename_5<-"Figures/LIMMA_venn_diagram.png"
                                            
                                            png(temp_filename_5,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                        }
                                       

                                            vennDiagram(results2,cex=0.8)
                                            
                                            if(output.device.type!="pdf"){
                                                
                                                try(dev.off(),silent=TRUE)
                                            }

                                    }
                                    }else{
                                    #dev.off()
                                     
                                     results2<-fit2$p.value
                                    }
									 cnames_tab<-colnames(data_m_fc_withfeats)
									 cnames_tab2<-design_mat_names #colnames(results2)
									 
									cnames_tab<-c("P.value","adjusted.P.value",cnames_tab2,cnames_tab)

									 data_limma_fdrall_withfeats<-cbind(p.value,adjusted.P.value,results2,data_m_fc_withfeats)
										data_limma_fdrall_withfeats<-as.data.frame(data_limma_fdrall_withfeats)
										filename<-"Tables/limma_posthoc1wayanova_results.txt"
										colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
								
								#data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
                                write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
								
                                data_limma_fdrall_withfeats<-cbind(p.value,adjusted.P.value,data_m_fc_withfeats)


									if(fdrmethod=="none"){
									filename<-paste("limma_posthoc1wayanova_pval",fdrthresh,"_results.txt",sep="")
									}else{
									filename<-paste("limma_posthoc1wayanova_fdr",fdrthresh,"_results.txt",sep="")
									}
								if(length(which(data_limma_fdrall_withfeats$adjusted.P.value<fdrthresh & data_limma_fdrall_withfeats$p.value<pvalue.thresh))>0){
								data_limma_sig_withfeats<-data_limma_fdrall_withfeats[data_limma_fdrall_withfeats$adjusted.P.value<fdrthresh & data_limma_fdrall_withfeats$p.value<pvalue.thresh,]
                                #write.table(data_limma_sig_withfeats, file=filename,sep="\t",row.names=FALSE)
								}
								
								# data_limma_fdrall_withfeats<-cbind(p.value,adjusted.p.value,data_m_fc_withfeats)
								 
								 data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
								 final.pvalues<-pvalues
								 
									 cnames_tab<-colnames(data_m_fc_withfeats)
													cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
													colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
								 }
								 
                                 #pvalues<-data_limma_fdrall_withfeats$p.value
                                 
                                 #final.pvalues<-pvalues
                                 
								# print("checking here")
								sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
								
								goodip<-which(sel.diffdrthresh==TRUE)
							   d4<-as.data.frame(data_limma_fdrall_withfeats)
							   logp<-(-1)*log((d4[,1]+(10^-20)),10)
								
								#tiff("pval_dist.tiff",compression="lzw")
                                #hist(d4[,1],xlab="p",main="Distribution of p-values")
								#dev.off()
								
									
	}
											
									
											if(featselmethod=="limma2way")
											{
											#design <- cbind(Grp1vs2=c(rep(1,num_samps_group[[1]]),rep(0,num_samps_group[[2]])),Grp2vs1=c(rep(0,num_samps_group[[1]]),rep(1,num_samps_group[[2]])))
									 #print("here")
									 design <- model.matrix(~ -1+f)
									
									# print(data_m_fc[1:4,])
									 #colnames(design) <- levels(f)
									colnames(design)<-levels(factor(sampleclass))
								
								
								options(digit=3)
								parameterNames<-colnames(design)
								
								print("Design matrix")
								print(design)
								
								if(pairedanalysis==TRUE)
								{

									
                                    print("Paired design is")
									#print(f1)
                                    print(subject_inf)
									
									f1<-subject_inf
									
									#print(data_m_fc[1:10,1:10])
									
									
								}
								
								if(dim(design)[2]==4){
								
								#cont.matrix <- makeContrasts(Grp1vs2="ClassA-ClassB",Grp1vs3="ClassC-ClassD",Grp2vs3=("ClassA-ClassB")-("ClassC-ClassD"),levels=design)
								#cont.matrix <- makeContrasts(Grp1vs2=ClassA-ClassB,Grp1vs3=ClassC-ClassD,Grp2vs3=(ClassA-ClassB)-(ClassC-ClassD),Grp3vs4=ClassA-ClassC,Group2vs4=ClassB-ClassD,levels=design)
								
								cont.matrix <- makeContrasts(Factor1=(ClassA+ClassB)-(ClassC+ClassD),Factor2=(ClassA+ClassC)-(ClassB+ClassD),Factor1x2=(ClassA-ClassB)-(ClassC-ClassD),levels=design)
								if(pairedanalysis==TRUE){
								
									class_table_facts<-table(classlabels)
									
									#f1<-c(seq(1,num_samps_group[[1]]),seq(1,num_samps_group[[2]]),seq(1,num_samps_group[[1]]),seq(1,num_samps_group[[2]]))
									
								
									corfit<-duplicateCorrelation(data_m_fc,design=design,block=subject_inf,ndups=1)
									
                                    #print(f1)
									

									fit<-lmFit(data_m_fc,design,block=f1,cor=corfit$consensus)

								}else{
								
								
									fit <- lmFit(data_m_fc,design)
								}

								
								}else{
									
									if(dim(design)[2]==6){
								
								#cont.matrix <- makeContrasts(Grp1vs2="ClassA-ClassB",Grp1vs3="ClassC-ClassD",Grp2vs3=("ClassA-ClassB")-("ClassC-ClassD"),levels=design)
								#cont.matrix <- makeContrasts(Grp1vs2=ClassA-ClassB,Grp1vs3=ClassC-ClassD,Grp2vs3=(ClassA-ClassB)-(ClassC-ClassD),Grp3vs4=ClassA-ClassC,Group2vs4=ClassB-ClassD,levels=design)
								
								cont.matrix <- makeContrasts(Factor1=(ClassA+ClassB+ClassC)-(ClassD+ClassE+ClassF),Factor2=(ClassA+ClassD)-(ClassB+ClassE)-(ClassC+ClassF),Factor1x2=(ClassA-ClassB-ClassC)-(ClassD-ClassE-ClassF),levels=design)
								
								if(pairedanalysis==TRUE){
								
									class_table_facts<-table(classlabels)
									
									#f1<-c(seq(1,num_samps_group[[1]]),seq(1,num_samps_group[[2]]),seq(1,num_samps_group[[1]]),seq(1,num_samps_group[[2]]))
								
									corfit<-duplicateCorrelation(data_m_fc,design=design,block=subject_inf,ndups=1)
									
									

									fit<-lmFit(data_m_fc,design,block=f1,cor=corfit$consensus)

								}else{
								
								
									fit <- lmFit(data_m_fc,design)
								}
								
								
								}else{
									if(dim(design)[2]==8){
								
								#cont.matrix <- makeContrasts(Grp1vs2="ClassA-ClassB",Grp1vs3="ClassC-ClassD",Grp2vs3=("ClassA-ClassB")-("ClassC-ClassD"),levels=design)
								#cont.matrix <- makeContrasts(Grp1vs2=ClassA-ClassB,Grp1vs3=ClassC-ClassD,Grp2vs3=(ClassA-ClassB)-(ClassC-ClassD),Grp3vs4=ClassA-ClassC,Group2vs4=ClassB-ClassD,levels=design)
								
								cont.matrix <- makeContrasts(Factor1=(ClassA+ClassB+ClassC+ClassD)-(ClassE+ClassF+ClassG+ClassH),Factor2=(ClassA+ClassE)-(ClassB+ClassF)-(ClassC+ClassG)-(ClassD+ClassH),Factor1x2=(ClassA-ClassB-ClassC-ClassD)-(ClassE-ClassF-ClassG-ClassH),levels=design)

								if(pairedanalysis==TRUE){
								
									class_table_facts<-table(classlabels)
									
									#f1<-c(seq(1,num_samps_group[[1]]),seq(1,num_samps_group[[2]]),seq(1,num_samps_group[[1]]),seq(1,num_samps_group[[2]]))
									
								
									corfit<-duplicateCorrelation(data_m_fc,design=design,block=subject_inf,ndups=1)
									
									

									fit<-lmFit(data_m_fc,design,block=f1,cor=corfit$consensus)

								}else{
								
								
									fit <- lmFit(data_m_fc,design)
								}



								}else{
									
										if(dim(design)[2]==10){
								
								#cont.matrix <- makeContrasts(Grp1vs2="ClassA-ClassB",Grp1vs3="ClassC-ClassD",Grp2vs3=("ClassA-ClassB")-("ClassC-ClassD"),levels=design)
								#cont.matrix <- makeContrasts(Grp1vs2=ClassA-ClassB,Grp1vs3=ClassC-ClassD,Grp2vs3=(ClassA-ClassB)-(ClassC-ClassD),Grp3vs4=ClassA-ClassC,Group2vs4=ClassB-ClassD,levels=design)
								
								cont.matrix <- makeContrasts(Factor1=(ClassA+ClassB+ClassC+ClassD+ClassE)-(ClassF+ClassG+ClassH+ClassI+ClassJ),Factor2=(ClassA+ClassF)-(ClassB+ClassG)-(ClassC+ClassH)-(ClassD+ClassI)-(ClassE+ClassJ),Factor1x2=(ClassA-ClassB-ClassC-ClassD-ClassE)-(ClassF-ClassG-ClassH-ClassI-ClassJ),levels=design)
								
								if(pairedanalysis==TRUE){
								
									class_table_facts<-table(classlabels)
									
									#f1<-c(seq(1,num_samps_group[[1]]),seq(1,num_samps_group[[2]]),seq(1,num_samps_group[[1]]),seq(1,num_samps_group[[2]]))
									
								
								
									corfit<-duplicateCorrelation(data_m_fc,design=design,block=subject_inf,ndups=1)
									
									

									fit<-lmFit(data_m_fc,design,block=f1,cor=corfit$consensus)

								}else{
								
								
									fit <- lmFit(data_m_fc,design)
								}
								
								}else{
									
									if(dim(design)[2]==12){
								
								#cont.matrix <- makeContrasts(Grp1vs2="ClassA-ClassB",Grp1vs3="ClassC-ClassD",Grp2vs3=("ClassA-ClassB")-("ClassC-ClassD"),levels=design)
								#cont.matrix <- makeContrasts(Grp1vs2=ClassA-ClassB,Grp1vs3=ClassC-ClassD,Grp2vs3=(ClassA-ClassB)-(ClassC-ClassD),Grp3vs4=ClassA-ClassC,Group2vs4=ClassB-ClassD,levels=design)
								
								cont.matrix <- makeContrasts(Factor1=(ClassA+ClassB+ClassC+ClassD+ClassE+ClassF)-(ClassG+ClassH+ClassI+ClassJ+ClassK+ClassL),Factor2=(ClassA+ClassG)-(ClassB+ClassH)-(ClassC+ClassI)-(ClassD+ClassJ)-(ClassE+ClassK)-(ClassF-ClassL),Factor1x2=(ClassA-ClassB-ClassC-ClassD-ClassE-ClassF)-(ClassG-ClassH-ClassI-ClassJ-ClassK-ClassL),levels=design)

								if(pairedanalysis==TRUE){
								
									class_table_facts<-table(classlabels)
									
									#f1<-c(seq(1,num_samps_group[[1]]),seq(1,num_samps_group[[2]]),seq(1,num_samps_group[[1]]),seq(1,num_samps_group[[2]]))
									
									
								
									corfit<-duplicateCorrelation(data_m_fc,design=design,block=subject_inf,ndups=1)
									
									

									fit<-lmFit(data_m_fc,design,block=f1,cor=corfit$consensus)

								}else{
								
								
									fit <- lmFit(data_m_fc,design)
								}

                                    }else{
                                     
                                     cont.matrix <- makeContrasts(Factor1=(ClassA+ClassB+ClassC+ClassD+ClassE+ClassF+ClassG)-(ClassH+ClassI+ClassJ+ClassK+ClassL+ClassM+ClassN),Factor2=(ClassA+ClassH)-(ClassB+ClassI)-(ClassC+ClassJ)-(ClassD+ClassK)-(ClassE+ClassL)-(ClassF+ClassM)-(ClassG+ClassN),Factor1x2=(ClassA-ClassB-ClassC-ClassD-ClassE-ClassF-ClassG)-(ClassH-ClassI-ClassJ-ClassK-ClassL-ClassM-ClassN),levels=design)
                                     
                                     if(pairedanalysis==TRUE){
                                         
                                         class_table_facts<-table(classlabels)
                                         
                                         #f1<-c(seq(1,num_samps_group[[1]]),seq(1,num_samps_group[[2]]),seq(1,num_samps_group[[1]]),seq(1,num_samps_group[[2]]))
                                         
                                         
                                         
                                         corfit<-duplicateCorrelation(data_m_fc,design=design,block=subject_inf,ndups=1)
                                         
                                         
                                         
                                         fit<-lmFit(data_m_fc,design,block=f1,cor=corfit$consensus)
                                         
                                     }else{
                                         
                                         
                                         fit <- lmFit(data_m_fc,design)
                                     }
                                     
                                    }
									
									}

									
									}
									
									
									}
									
									#stop("Please use \"RF\" or \"MARS\" as featselmethod for more than 3 classes.")
									}
								
                                
                                #   print("Limma classlabels")
                                
                                #print(head(classlabels))
                                # print(head(classlabels_orig))
                                
                                fit2<-contrasts.fit(fit,cont.matrix)
                                
                                fit2<-eBayes(fit2)
                                
                                #fit3$F.p.value[1:10]
                                
                               # t1<-topTable(fit2,coef=1:3,n=dim(data_m_fc)[1])
                                #	write.table(t1,file="topTable_limma.txt",sep="\t")
								
								# Ordinary fit
								
								
                                
									
                                #fit2  <- contrasts.fit(fit, cont.matrix)

                                #fit2 <- eBayes(fit2)
								#as.data.frame(fit2[1:10,])

								# Various ways of summarising or plotting the results
								#topTable(fit,coef=2)

												if(dim(design)[2]>2){
													pvalues<-fit2$F.p.value
													p.value<-fit2$F.p.value
													
													 
													 
													}else{
												pvalues<-fit2$p.value
												p.value<-fit2$p.value
												}
																
												if(fdrmethod=="BH"){
												fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
												}else{
														if(fdrmethod=="ST"){
                                                            #fdr_adjust_pvalue<-qvalue(pvalues)
                                                            #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                            
                                                            fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
                                                            
                                                            if(is(fdr_adjust_pvalue,"try-error")){
                                                                
                                                                fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
                                                            }
                                                            
                                                            fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                            
                                                            
														}else{
																if(fdrmethod=="Strimmer"){
																	pdf("fdrtool.pdf")
																	#par_rows=1
																	#par(mfrow=c(par_rows,1))
																fdr_adjust_pvalue<-fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE)
																fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
																	try(dev.off(),silent=TRUE)
																}else{
																	if(fdrmethod=="none"){
																		#	fdr_adjust_pvalue<-pvalues
																		fdr_adjust_pvalue<-p.adjust(pvalues,method="none")	
																	}else{
																		if(fdrmethod=="BY"){
											fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                                                                        }else{
                                                                            if(fdrmethod=="bonferroni"){
                                                                                fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                                                                            }
                                                                        }
																		}
																}
															}
														
												}

								#print("Doing this:")
								
								adjusted.p.value<-fdr_adjust_pvalue
								
								data_limma_fdrall_withfeats<-cbind(p.value,adjusted.p.value,data_m_fc_withfeats)
								
                                if(limmadecideTests==TRUE){
                                    results2<-decideTests(fit2,method="nestedF",adjust.method="BH",p.value=fdrthresh)
                                    #tiff("comparison_contrast_overlap.tiff",width=plots.width,height=plots.height,res=plots.res, compression="lzw")
                                    
                                
                                    if(length(class_labels_levels)<4){
                                        
                                        
                                        if(output.device.type!="pdf"){
                                            
                                            temp_filename_1<-"Figures/LIMMA_venn_diagram.png"
                                            
                                            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                        }
                                      
                                           vennDiagram(results2,cex=0.8)
                                           
                                            if(output.device.type!="pdf"){
                                               
                                              try(dev.off(),silent=TRUE)
                                           }

                                        
                                    }
                                }else{
                                    #dev.off()
                                    
                                    results2<-fit2$p.value
                                }
                                
                                
									#tiff("comparison_contrast_overlap.tiff",width=plots.width,height=plots.height,res=plots.res, compression="lzw")
									
									 #dev.off()
                                     
                                     #results2<-fit2$p.value
									 
									 cnames_tab<-colnames(data_m_fc_withfeats)
									 cnames_tab2<-colnames(results2)
									cnames_tab<-c("P.value","adjusted.P.value",cnames_tab2,cnames_tab)

									 data_limma_fdrall_withfeats<-cbind(p.value,adjusted.p.value,results2,data_m_fc_withfeats)
                                
                               
                                
                                classlabels_orig<-as.data.frame(classlabels_orig)
                                
                                if(limmadecideTests==TRUE){
                                
                                
                                X1=data_m_fc_withfeats[which(data_limma_fdrall_withfeats[,3]!=0),]
                                X2=data_m_fc_withfeats[which(data_limma_fdrall_withfeats[,4]!=0),]
                                 X3=data_m_fc_withfeats[which(data_limma_fdrall_withfeats[,5]!=0),]
                                 
                                }else{
                                    X1=data_m_fc_withfeats[which(data_limma_fdrall_withfeats[,3]<fdrthresh),]
                                    X2=data_m_fc_withfeats[which(data_limma_fdrall_withfeats[,4]<fdrthresh),]
                                    X3=data_m_fc_withfeats[which(data_limma_fdrall_withfeats[,5]<fdrthresh),]
                                    
                                }
                                if(pairedanalysis==FALSE){
                                    Y1=classlabels_orig[,c(1:2)] #,classlabels_orig[,2])
                                }else{
                                    Y1=classlabels_orig[,c(1,3)] #cbind(classlabels_orig[,1],classlabels_orig[,3])
                                    
                                }
                                Y1<-as.data.frame(Y1)
                                
                                
                                
                                
                                if(output.device.type!="pdf"){
                                    
                                    temp_filename_1<-"Figures/HCA_Factor1selectedfeats.png"
                                    
                                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                }
                                
                                #save(list=ls(),file="hca_factor1.Rda")
                                                                get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X1,Y=Y1,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                                sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="Factor1")
                                
                               if(output.device.type!="pdf"){
                                    
                                    try(dev.off(),silent=TRUE)
                                }

                                
                                #Y2=cbind(classlabels_orig[,1],classlabels_orig[,2])
                                
                                if(pairedanalysis==FALSE){
                                    Y2=classlabels_orig[,c(1,3)] #cbind(classlabels_orig[,1],classlabels_orig[,3])
                                }else{
                                    Y2=classlabels_orig[,c(1,4)] #cbind(classlabels_orig[,1],classlabels_orig[,4])
                                    
                                }
                                
                                Y2<-as.data.frame(Y2)
                                
                                if(output.device.type!="pdf"){
                                    
                                    temp_filename_1<-"Figures/HCA_Factor2selectedfeats.png"
                                    
                                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                }
                               
                                                                get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X2,Y=Y2,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                                sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="Factor2")
                                
                                 if(output.device.type!="pdf"){
                                    
                                    try(dev.off(),silent=TRUE)
                                }

                                
                                if(pairedanalysis==FALSE){
                                    
                                    
                                    class_interact<-paste(classlabels_orig[,2],":",classlabels_orig[,3],sep="") #classlabels_orig[,2]:classlabels_orig[,3]

                                }else{
                                    class_interact<-paste(classlabels_orig[,3],":",classlabels_orig[,4],sep="") #classlabels_orig[,3]:classlabels_orig[,4]
                                    
                                }
                                
                                
                               
                                Y3=cbind(classlabels_orig[,1],class_interact)
                                Y3<-as.data.frame(Y3)
                                
                                # print(Y3)
                                #print(classlabels_orig)
                                #print(dim(X3))
                              
                              ##save(classlabels_orig,file="classlabels_orig.Rda")
                              
                              if(output.device.type!="pdf"){
                                  
                                  temp_filename_1<-"Figures/HCA_Factor1xFactor2selectedfeats.png"
                                  
                                  png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                              }
                                 get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X3,Y=Y3,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                                sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="Factor1 x Factor2")
                                
                                 if(output.device.type!="pdf"){
                                    
                                    try(dev.off(),silent=TRUE)
                                }


                                
										filename<-"Tables/limma_2wayanova_posthocresults.txt"
										
										colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
								
								#data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
								write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)

								 data_limma_fdrall_withfeats<-cbind(p.value,adjusted.p.value,data_m_fc_withfeats)
								 
								 		# data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
								 		
								 					 cnames_tab<-colnames(data_m_fc_withfeats)
													cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
													colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
                                    #write.table(data_limma_fdrall_withfeats,file="Limma_posthoc2wayanova_results.txt",sep="\t",row.names=FALSE)
								#print("checking here")
                                pvalues<-p.value
                                
                                final.pvalues<-pvalues
                                
								sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
								
								goodip<-which(sel.diffdrthresh==TRUE)
							   d4<-as.data.frame(data_limma_fdrall_withfeats)
							   logp<-(-1)*log((d4[,1]+(10^-20)),10)
								
								
								
								results2<-decideTests(fit2,method="nestedF",adjust.method=fdrmethod,p.value=fdrthresh)
									
									
									
								}
									
									
					
                        if(featselmethod=="RF")
                        {
							
							maxint<-apply(data_m_fc,1,max)
							
							
                            data_m_fc_withfeats<-as.data.frame(data_m_fc_withfeats)
                            
                            data_m_fc<-as.data.frame(data_m_fc)
                            #write.table(classlabels,file="classlabels_rf.txt",sep="\t",row.names=FALSE)
                  
                            ##save(data_m_fc,classlabels,numtrees,analysismode,file="rfdebug.Rda")
                            
                            rf_classlabels<-classlabels[,1]
                            
                            if(rfconditional==TRUE){
                                
                                print("Performing random forest analysis using the cforest and Boruta functions")
                                
                            rfcondres1<-do_rf_conditional(X=data_m_fc,rf_classlabels,ntrees=numtrees,analysismode) #,silent=TRUE)
                            filename<-"RFconditional_VIM_withBoruta_allfeats.txt"
                            }else{
                                
                                    print("Performing random forest analysis using the randomForest and Boruta functions")
                                
                                    #rfcondres1<-do_rf(X=data_m_fc,rf_classlabels,ntrees=numtrees,analysismode)
                                    filename<-"RF_VIM_withBoruta_allfeats.txt"
                                
                            }
                            
                                temp1<-cbind(rf_classlabels,t(data_m_fc))
                                temp1<-as.data.frame(temp1)
                              #  #save(temp1,file="temp1.Rda")
                                set.seed(290875)
                                Boruta(rf_classlabels~.,data=temp1,maxRuns=1000)->Bor.test
                                
                                
                                mean_imphistory<-apply(Bor.test$ImpHistory,2,mean)
                  
                                varimp_res2<-mean_imphistory[1:(length(mean_imphistory)-3)] #rep(0,length(Bor.test$finalDecision))
                                
                                temp_decision_vec<-as.character(Bor.test$finalDecision)
                                
                                if(length(which(temp_decision_vec=="Rejected"))>0){
                                    temp_decision_vec<-replace(temp_decision_vec,which(temp_decision_vec=="Rejected"),0)
                                }
                                
                                if(length(which(temp_decision_vec=="Tentative"))>0){
                                    temp_decision_vec<-replace(temp_decision_vec,which(temp_decision_vec=="Tentative"),1)
                                }
                                if(length(which(temp_decision_vec=="Confirmed"))>0){
                                    temp_decision_vec<-replace(temp_decision_vec,which(temp_decision_vec=="Confirmed"),2)
                                }
                                temp_decision_vec<-as.numeric(as.character(temp_decision_vec))
                                
                                
                                #vim_matrix<-Bor.test$ImpHistory
                                
                                
                                
                                #varimp_res2<-temp_decision_vec
                                #write.table(varimp_res2, file=filename,sep="\t",row.names=FALSE)
                                ##save(rfcondres1,file="rfcondres1.Rda")
                                
                                #varimp_res2<-rfcondres1$rf_varimp
                                                    
                                                    if(length(which(temp_decision_vec<2))>0){
                                                        varimp_res2[which(temp_decision_vec<2)]<-0
                                                    }
                                                   
                                                   varimp_res3<-cbind(data_m_fc_withfeats[,c(1:2)],varimp_res2)
                                                   
						   filename<-paste("Tables/",filename,sep="")
                                                   write.table(varimp_res3, file=filename,sep="\t",row.names=TRUE)
                                                
                                            
                                            
                                                                        goodip<-which(varimp_res2>0)
                                                                        
                                                                        if(length(goodip)<1){
                                                                            print("No features were selected using the selection criteria.")
                                                                        }
                                                                        var_names<-paste(sprintf("%.3f",data_m_fc_withfeats[,1]),sprintf("%.1f",data_m_fc_withfeats[,2]),sep="_")
                                                                        
                                                                        names(varimp_res2)<-as.character(var_names)
                                                                        sel.diffdrthresh<-varimp_res2>0
                                                                        
                                                                        if(length(which(sel.diffdrthresh==TRUE))<1){
                                                                            print("No features were selected using the selection criteria")
                                                                        }
                                                                        
                                                   

						num_var_rf<-length(which(sel.diffdrthresh==TRUE))
						
						if(num_var_rf>10){
							
							num_var_rf=10
						}
						sorted_varimp_res<-varimp_res2[order(varimp_res2,decreasing=TRUE)[1:(num_var_rf)]]

						sorted_varimp_res<-sort(sorted_varimp_res)
		
						barplot_text=paste("Variable Importance measures (top ",length(sorted_varimp_res)," shown)\n",sep="")
						
				        if(output.device.type!="pdf"){
                            
                            temp_filename_1<-"Figures/RF_selectfeats_VIMbarplot.png"
                            
                            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                        }
                        

						barplot(sorted_varimp_res, xlab="Selected features", main=barplot_text,cex.axis=0.5,cex.names=0.4, ylab="VIM",range(pretty(c(0,sorted_varimp_res))),space=0.1)
					
                            if(output.device.type!="pdf"){
                                
                                try(dev.off(),silent=TRUE)
                            }




                    rank_num<-rank(-varimp_res2)
                
					data_limma_fdrall_withfeats<-cbind(varimp_res2,rank_num,data_m_fc_withfeats)
					
					cnames_tab<-colnames(data_m_fc_withfeats)
					cnames_tab<-c("VIM","Rank",cnames_tab)
					
                    goodip<-which(sel.diffdrthresh==TRUE)
					    
					feat_sigfdrthresh[lf]<-length(which(sel.diffdrthresh==TRUE))
	
					colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
					
					#data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
                    #write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
					

						}
							
						if(featselmethod=="MARS"){
							
						mars_classlabels<-classlabels[,1]
						marsres1<-do_mars(X=data_m_fc,mars_classlabels, analysismode,kfold)
						
						varimp_marsres1<-marsres1$mars_varimp
						
						mars_mznames<-rownames(varimp_marsres1)
						
						
						all_names<-paste("mz",seq(1,dim(data_m_fc)[1]),sep="")
					
						com1<-match(all_names,mars_mznames)
					

						filename<-"MARS_variable_importance.txt"
						
                        
                        if(is.na(max_varsel)==FALSE){
                            
                            if(max_varsel>dim(data_m_fc)[1]){
                                max_varsel=dim(data_m_fc)[1]
                            }
                            varimp_res2<-varimp_marsres1[,4]
                            
                            #sort by VIM; and keep the top max_varsel scores
                            sorted_varimp_res<-varimp_res2[order(varimp_res2,decreasing=TRUE)[1:(max_varsel)]]
                            
                            #get the minimum VIM from the top max_varsel scores
                            min_thresh<-min(sorted_varimp_res[which(sorted_varimp_res>=mars.gcv.thresh)],na.rm=TRUE)
                            
                
                            row_num_vec<-seq(1,length(varimp_res2))
                            
                            #only use the top max_varsel scores
                            #goodip<-order(varimp_res2,decreasing=TRUE)[1:(max_varsel)]
                            #sel.diffdrthresh<-row_num_vec%in%goodip
                            
                            #use a threshold of mars.gcv.thresh
                            sel.diffdrthresh<-varimp_marsres1[,4]>=min_thresh
                            
                            goodip<-which(sel.diffdrthresh==TRUE)
                            
                        }else{
                        
                        #use a threshold of mars.gcv.thresh
						sel.diffdrthresh<-varimp_marsres1[,4]>=mars.gcv.thresh
						
                        goodip<-which(sel.diffdrthresh==TRUE)
                        }
                        
                        
                        num_var_rf<-length(which(sel.diffdrthresh==TRUE))
                        
                        if(num_var_rf>10){
                            
                            num_var_rf=10
                        }
                        sorted_varimp_res<-varimp_res2[order(varimp_res2,decreasing=TRUE)[1:(num_var_rf)]]
                        
                        sorted_varimp_res<-sort(sorted_varimp_res)
                        
                        barplot_text=paste("Generalized cross validation (top ",length(sorted_varimp_res)," shown)\n",sep="")
                        
                        if(output.device.type!="pdf"){
                            
                            temp_filename_1<-"Figures/MARS_selectfeats_GCVbarplot.png"
                            
                            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                        }
                        
                        
                        barplot(sorted_varimp_res, xlab="Selected features", main=barplot_text,cex.axis=0.5,cex.names=0.4, ylab="GCV",range(pretty(c(0,sorted_varimp_res))),space=0.1)
                        
                        if(output.device.type!="pdf"){
                            
                            try(dev.off(),silent=TRUE)
                        }

						
						data_limma_fdrall_withfeats<-cbind(varimp_marsres1[,c(4,6)],data_m_fc_withfeats)
					
						cnames_tab<-colnames(data_m_fc_withfeats)
						cnames_tab<-c("GCV importance","RSS importance",cnames_tab)
						feat_sigfdrthresh[lf]<-length(which(sel.diffdrthresh==TRUE))
						
						colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
					
					
						goodip<-which(sel.diffdrthresh==TRUE)
						
					
						
						}
						
								if(featselmethod=="pls" | featselmethod=="o1pls" | featselmethod=="o2pls" | featselmethod=="spls" | featselmethod=="o1spls" | featselmethod=="o2spls")
								{
									
                                    

									classlabels<-as.data.frame(classlabels)
									
							
                                    if(is.na(max_comp_sel)==TRUE){
                                            max_comp_sel=pls_ncomp
                                    }
                
                                    rand_pls_sel<-{} #new("list")
									if(featselmethod=="spls" | featselmethod=="o1spls" | featselmethod=="o2spls"){
								
                                
                                            if(featselmethod=="o1spls"){
                                    
                                                featselmethod="o1pls"
                                                
                                            }else{
                                    
                                                if(featselmethod=="o2spls"){
                                                        featselmethod="o2pls"
                                                    }
                                    
                                            }
                                
									if(pairedanalysis==TRUE){
                                        
									classlabels_temp<-cbind(classlabels_sub[,2],classlabels)


                                        set.seed(999)

                                        plsres1<-do_plsda(X=data_m_fc,Y=classlabels_sub,oscmode=featselmethod,numcomp=pls_ncomp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=TRUE,analysismode,sample.col.opt=sample.col.opt,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=optselect,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,output.device.type=output.device.type,plots.res=plots.res,plots.width=plots.width,plots.height=plots.height,plots.type=plots.type)

                                        if (is(plsres1, "try-error")){
                                            print(paste("sPLS could not be performed at RSD threshold: ",log2.fold.change.thresh,sep=""))
                                            #break;
                                        }

                                        opt_comp<-plsres1$opt_comp
                                                                            #for(randindex in 1:100)

                                       
									if(is.na(pls.permut.count)==FALSE){
                                        set.seed(999)
                                        seedvec<-runif(pls.permut.count,10,10*pls.permut.count)
                                        
                                        
                                        
                                        if(pls.permut.count>0){
                                            
                                            cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                                            clusterEvalQ(cl,library(plsgenomics))
                                            clusterEvalQ(cl,library(dplyr))
                                            
                                            clusterEvalQ(cl,library(plyr))
                                            clusterExport(cl,"pls.lda.cv",envir = .GlobalEnv)
                                            clusterExport(cl,"plsda_cv",envir = .GlobalEnv)
                                            #clusterExport(cl,"%>%",envir = .GlobalEnv) #%>%
                                            clusterExport(cl,"do_plsda_rand",envir = .GlobalEnv)
                                            clusterEvalQ(cl,library(mixOmics))
                                            clusterEvalQ(cl,library(pls))
                                            
                                            
                                            rand_pls_sel<-parLapply(cl,1:pls.permut.count,function(x)
                                            {
                                                
                                                set.seed(seedvec[x])


                                        plsresrand<-do_plsda_rand(X=data_m_fc,Y=classlabels_sub[sample(x=seq(1,dim(classlabels_sub)[1]),size=dim(classlabels_sub)[1]),],oscmode=featselmethod,numcomp=opt_comp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=TRUE,analysismode,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=FALSE,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,plotindiv=FALSE) #,silent=TRUE)

										#rand_pls_sel<-cbind(rand_pls_sel,plsresrand$vip_res[,1])
										if (is(plsresrand, "try-error")){
                                          
											return(rep(0,dim(data_m_fc)[1]))
											}else{
											return(plsresrand$vip_res[,1])
										}
									})
                                        }
									}
                                    
                                    
									
									}else{	
                                        #plsres1<-try(do_plsda(X=data_m_fc,Y=classlabels,oscmode=featselmethod,numcomp=pls_ncomp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=TRUE,analysismode,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=optselect,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,pls.vip.selection=pls.vip.selection),silent=TRUE)
                                    
                                    set.seed(999)
                                    plsres1<-do_plsda(X=data_m_fc,Y=classlabels,oscmode=featselmethod,numcomp=pls_ncomp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=TRUE,analysismode,sample.col.opt=sample.col.opt,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=optselect,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,pls.vip.selection=pls.vip.selection,output.device.type=output.device.type,plots.res=plots.res,plots.width=plots.width,plots.height=plots.height,plots.type=plots.type)
                                    
                                    opt_comp<-plsres1$opt_comp
                                    
                                    if (is(plsres1, "try-error")){
                                        print(paste("sPLS could not be performed at RSD threshold: ",log2.fold.change.thresh,sep=""))
                                        break;
                                    }
											#for(randindex in 1:100)
									if(is.na(pls.permut.count)==FALSE){
                                        
                                        set.seed(999)
                                        seedvec<-runif(pls.permut.count,10,10*pls.permut.count)
                                        
                                        
                                        
                                        if(pls.permut.count>0){
                                            
                                            cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                                            clusterEvalQ(cl,library(plsgenomics))
                                            clusterEvalQ(cl,library(dplyr))
                                            
                                            clusterEvalQ(cl,library(plyr))
                                            clusterExport(cl,"pls.lda.cv",envir = .GlobalEnv)
                                            clusterExport(cl,"plsda_cv",envir = .GlobalEnv)
                                            #clusterExport(cl,"%>%",envir = .GlobalEnv) #%>%
                                            clusterExport(cl,"do_plsda_rand",envir = .GlobalEnv)
                                            clusterEvalQ(cl,library(mixOmics))
                                            clusterEvalQ(cl,library(pls))
                                            
                                            
                                            rand_pls_sel<-parLapply(cl,1:pls.permut.count,function(x)
                                            {
                                                
                                                set.seed(seedvec[x])


plsresrand<-do_plsda_rand(X=data_m_fc,Y=classlabels[sample(x=seq(1,dim(classlabels)[1]),size=dim(classlabels)[1]),],oscmode=featselmethod,numcomp=opt_comp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=TRUE,analysismode,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=FALSE,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,plotindiv=FALSE)
                                                                                
                                                                                #rand_pls_sel<-cbind(rand_pls_sel,plsresrand$vip_res[,1])
                                                                        	#return(plsresrand$vip_res[,1])		
										if (is(plsresrand, "try-error")){
                                            
									                return(rep(0,dim(data_m_fc)[1]))
                                                                                        }else{
                                                                                        return(plsresrand$vip_res[,1])
                                                                                }
									})
                                            
                                            stopCluster(cl)
                                            
									}
                                    }
                                    
									}
									pls_vip_thresh<-0
									
									if (is(plsres1, "try-error")){
										print(paste("sPLS could not be performed at RSD threshold: ",log2.fold.change.thresh,sep=""))
                                        break;
									}else{	
									opt_comp<-plsres1$opt_comp
									}

									}else{
                                        #PLS
											if(pairedanalysis==TRUE){
                                                                        classlabels_temp<-cbind(classlabels_sub[,2],classlabels)
                                                                        plsres1<-do_plsda(X=data_m_fc,Y=classlabels_temp,oscmode=featselmethod,numcomp=pls_ncomp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=FALSE,analysismode=analysismode,vip.thresh=pls_vip_thresh,sample.col.opt=sample.col.opt,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=optselect,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,pls.vip.selection=pls.vip.selection,output.device.type=output.device.type,plots.res=plots.res,plots.width=plots.width,plots.height=plots.height,plots.type=plots.type)
                                                                        
                                                                        if (is(plsres1, "try-error")){
                                                                            print(paste("PLS could not be performed at RSD threshold: ",log2.fold.change.thresh,sep=""))
                                                                            break;
                                                                        }else{
                                                                            opt_comp<-plsres1$opt_comp
                                                                        }

                                            }else{
                                                
											plsres1<-do_plsda(X=data_m_fc,Y=classlabels,oscmode=featselmethod,numcomp=pls_ncomp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=FALSE,analysismode=analysismode,vip.thresh=pls_vip_thresh,sample.col.opt=sample.col.opt,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=optselect,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,pls.vip.selection=pls.vip.selection,output.device.type=output.device.type,plots.res=plots.res,plots.width=plots.width,plots.height=plots.height,plots.type=plots.type)
                                            
                                            if (is(plsres1, "try-error")){
                                                print(paste("PLS could not be performed at RSD threshold: ",log2.fold.change.thresh,sep=""))
                                                break;
                                            }else{
                                                opt_comp<-plsres1$opt_comp
                                            }
                                        				#for(randindex in 1:100){
									if(is.na(pls.permut.count)==FALSE){
                                        
                                        
                                        set.seed(999)
                                        seedvec<-runif(pls.permut.count,10,10*pls.permut.count)
                                        
                                       
                                        
                                        if(pls.permut.count>0){
                                            
                                            cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                                            clusterEvalQ(cl,library(plsgenomics))
                                            clusterEvalQ(cl,library(dplyr))
                                            
                                            clusterEvalQ(cl,library(plyr))
                                            clusterExport(cl,"pls.lda.cv",envir = .GlobalEnv)
                                            clusterExport(cl,"plsda_cv",envir = .GlobalEnv)
                                            #clusterExport(cl,"%>%",envir = .GlobalEnv) #%>%
                                            clusterExport(cl,"do_plsda_rand",envir = .GlobalEnv)
                                            clusterEvalQ(cl,library(mixOmics))
                                            clusterEvalQ(cl,library(pls))
                                            
                                            #here
                                                        rand_pls_sel<-parLapply(cl,1:pls.permut.count,function(x)
                                                        {

                                                                set.seed(seedvec[x])
                                                                #t1fname<-paste("ranpls",x,".Rda",sep="")
                                                                ##save(list=ls(),file=t1fname)
                                                                print(paste("PLSDA permutation number: ",x,sep=""))

                                                                plsresrand<-do_plsda_rand(X=data_m_fc,Y=classlabels[sample(x=seq(1,dim(classlabels)[1]),size=dim(classlabels)[1]),],oscmode=featselmethod,numcomp=opt_comp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=FALSE,analysismode,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=FALSE,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,plotindiv=FALSE) #,silent=TRUE)

                                                                                if (is(plsresrand, "try-error")){
                                                                                    
                                                                                    
                                                                                    return(1)
                                                                                }else{
                                                                                    return(plsresrand$vip_res[,1])
                                                                                    
                                                                                }
                                                                                
										
                                                       })
                                                        
                                                        stopCluster(cl)
                                        }
                                                        ##save(rand_pls_sel,file="rand_pls_sel1.Rda")
                                                        
									}
										
									}
											opt_comp<-plsres1$opt_comp
                                    }
										
									if(length(plsres1$bad_variables)>0){
										
										data_m_fc_withfeats<-data_m_fc_withfeats[-c(plsres1$bad_variables),]
											data_m_fc<-data_m_fc[-c(plsres1$bad_variables),]
									}
																			

print("pls.permut.count")
print(pls.permut.count)
					
									if(is.na(pls.permut.count)==FALSE){
                                        
                                         if(pls.permut.count>0){
                                            
                                         #save(rand_pls_sel,file="rand_pls_sel.Rda")
                                        
                                        #rand_pls_sel<-ldply(rand_pls_sel,rbind)  #unlist(rand_pls_sel)
                                        rand_pls_sel<-as.data.frame(rand_pls_sel)
                                        rand_pls_sel<-t(rand_pls_sel)
                                        rand_pls_sel<-as.data.frame(rand_pls_sel)
                                        
                                                    if(featselmethod=="spls"){
                                                     
                                                            rand_pls_sel[rand_pls_sel!=0]<-1
                                                    }else{
                                                        
                                                        rand_pls_sel[rand_pls_sel<pls_vip_thresh]<-0
                                                        rand_pls_sel[rand_pls_sel>=pls_vip_thresh]<-1
                                                        
                                                    }
                                    
                                                    ##save(rand_pls_sel,file="rand_pls_sel2.Rda")
                                                    rand_pls_sel_prob<-apply(rand_pls_sel,2,sum)/pls.permut.count
                                                    #rand_pls_sel_fdr<-p.adjust(rand_pls_sel_prob,method=fdrmethod)
                                                    pvalues<-rand_pls_sel_prob
                                                    if(fdrmethod=="BH"){
                                                        fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
                                                    }else{
                                                        if(fdrmethod=="ST"){
                                                            #fdr_adjust_pvalue<-qvalue(pvalues)
                                                            #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                            
                                                            fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
                                                            
                                                            if(is(fdr_adjust_pvalue,"try-error")){
                                                                
                                                                fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
                                                            }
                                                            
                                                            fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                            
                                                            
                                                        }else{
                                                            if(fdrmethod=="Strimmer"){
                                                                pdf("fdrtool.pdf")
                                                                #par_rows=1
                                                                #par(mfrow=c(par_rows,1))
                                                                fdr_adjust_pvalue<-fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE)
                                                                fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
                                                                try(dev.off(),silent=TRUE)
                                                            }else{
                                                                if(fdrmethod=="none"){
                                                                    fdr_adjust_pvalue<-pvalues
                                                                    #fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
                                                                }else{
                                                                    if(fdrmethod=="BY"){
                                                                        fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                                                                    }else{
                                                                        if(fdrmethod=="bonferroni"){
                                                                            fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                                                                        }
                                                                    }
                                                                    
                                                                }
                                                            }
                                                        }
                                                        
                                                    }
                                    
                                    rand_pls_sel_fdr<-fdr_adjust_pvalue
                                  
									
									vip_res<-cbind(data_m_fc_withfeats[,c(1:2)],plsres1$vip_res,rand_pls_sel_prob,rand_pls_sel_fdr)
                                    
                                         }else{
                                             vip_res<-cbind(data_m_fc_withfeats[,c(1:2)],plsres1$vip_res)
                                             rand_pls_sel_fdr<-rep(0,dim(data_m_fc_withfeats[,c(1:2)])[1])
                                             rand_pls_sel_prob<-rep(0,dim(data_m_fc_withfeats[,c(1:2)])[1])
                                         }
									}else{
										vip_res<-cbind(data_m_fc_withfeats[,c(1:2)],plsres1$vip_res)
										rand_pls_sel_fdr<-rep(0,dim(data_m_fc_withfeats[,c(1:2)])[1])
                                        rand_pls_sel_prob<-rep(0,dim(data_m_fc_withfeats[,c(1:2)])[1])
									}

									write.table(vip_res,file="Tables/vip_res.txt",sep="\t",row.names=FALSE)
									
									
					#				write.table(r2_q2_valid_res,file="pls_r2_q2_res.txt",sep="\t",row.names=TRUE)
									
									
									varimp_plsres1<-plsres1$selected_variables
						
					opt_comp<-plsres1$opt_comp				
                    if(max_comp_sel>opt_comp){
                        
                        max_comp_sel<-opt_comp
                    }
                    
                    #	print("opt comp is")
                    #print(opt_comp)
									if(featselmethod=="spls"){
										
										cnames_tab<-colnames(data_m_fc_withfeats)
										cnames_tab<-c("Loading (absolute)","Rank",cnames_tab)
										
									
										if(opt_comp>1){
                                            
                                            #abs
											vip_res1<-abs(plsres1$vip_res)
											
                                            if(max_comp_sel>1){
											vip_res1<-apply(vip_res1[,c(1:max_comp_sel)],1,max)
											
                                            }else{

												vip_res1<-vip_res1[,c(1)]
											}		
										}else{

											vip_res1<-abs(plsres1$vip_res)
										}	
									
												pls_vip<-vip_res1 #(plsres1$vip_res)
									
											
											#based on loadings for sPLS
											sel.diffdrthresh<-pls_vip!=0 & rand_pls_sel_fdr<fdrthresh & rand_pls_sel_prob<pvalue.thresh


                                            goodip<-which(sel.diffdrthresh==TRUE)
                                            
                                            
                                            
                                            
                                            
										
										}else{
											
												cnames_tab<-colnames(data_m_fc_withfeats)
												cnames_tab<-c("VIP","Rank",cnames_tab)
												
						
                        if(max_comp_sel>opt_comp){
                            
                            max_comp_sel<-opt_comp
                        }
														
									
												#pls_vip<-plsres1$vip_res[,c(1:max_comp_sel)]
									
                                   
                                                if(opt_comp>1){
                                                                                        vip_res1<-(plsres1$vip_res)
                                                                                        	if(max_comp_sel>1){
                                                                                
                                                                                if(pls.vip.selection=="mean"){
                                                                                vip_res1<-apply(vip_res1[,c(1:max_comp_sel)],1,mean)
                                                                                
                                                                                }else{
                                                                                    
                                                                                     vip_res1<-apply(vip_res1[,c(1:max_comp_sel)],1,max)
                                                                                }
                                                                                       		 }else{

                                                                                                	vip_res1<-vip_res1[,c(1)]
                                                                                        		}
                                                    }else{

                                                        vip_res1<-plsres1$vip_res
                                                    }


                                    
                                    #vip_res1<-plsres1$vip_res
                                    pls_vip<-vip_res1
											
                                            
                                           
                                            
                                            #pls
                                                sel.diffdrthresh<-pls_vip>=pls_vip_thresh & rand_pls_sel_fdr<fdrthresh & rand_pls_sel_prob<pvalue.thresh
                                                
                                                goodip<-which(sel.diffdrthresh==TRUE)
									
									}
									
									rank_vec<-order(pls_vip,decreasing=TRUE)
										rank_vec2<-seq(1,length(rank_vec))

				ranked_vec<-pls_vip[rank_vec]
				rank_num<-match(pls_vip,ranked_vec)


									
						data_limma_fdrall_withfeats<-cbind(pls_vip,rank_num,data_m_fc_withfeats)
					
					
						feat_sigfdrthresh[lf]<-length(which(sel.diffdrthresh==TRUE)) #length(plsres1$selected_variables) #length(which(sel.diffdrthresh==TRUE))
						
						filename<-paste("Tables/",featselmethod,"_variable_importance.txt",sep="")

						colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
					
						#data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
                        write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
					
							
								
                                
								}
								
								#stop("Please choose limma, RF, RFcond, or MARS for featselmethod.")
								if(featselmethod=="lmreg" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat"| featselmethod=="logitreg" | featselmethod=="wilcox" | featselmethod=="ttest" | featselmethod=="ttestrepeat" |  featselmethod=="poissonreg" | featselmethod=="wilcoxrepeat")
								{
								pvalues<-{}
								
								classlabels_response_mat<-as.data.frame(classlabels_response_mat)
					
                    if(featselmethod=="ttestrepeat"){
                        featselmethod="ttest"
                        pairedanalysis=TRUE
                    }
								          
                                          if(featselmethod=="wilcoxrepeat"){
                                              
                                              featselmethod="wilcox"
                                              pairedanalysis=TRUE
                                          }
								
										if(featselmethod=="lm1wayanova")
										{
										
											print("Performing one-way ANOVA analysis")
											
											#print(dim(data_m_fc))
											#print(dim(classlabels_response_mat))
											#print(dim(classlabels))
											
											#data_mat_anova<-cbind(t(data_m_fc),classlabels_response_mat)
                                            
                                            
                                            numcores<-round(detectCores()*0.6)
                                            
                                            cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                                            
                                            clusterExport(cl,"diffexponewayanova",envir = .GlobalEnv)
                                            
                                             clusterExport(cl,"anova",envir = .GlobalEnv)
                                             
                                             
                                             clusterExport(cl,"TukeyHSD",envir = .GlobalEnv)
                                             
                                              clusterExport(cl,"aov",envir = .GlobalEnv)
                                            


                                                res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat){
														xvec<-x	
														
																												
														data_mat_anova<-cbind(xvec,classlabels_response_mat)
														
														data_mat_anova<-as.data.frame(data_mat_anova)
														cnames<-colnames(data_mat_anova)
														
														cnames[1]<-"Response"
														
														colnames(data_mat_anova)<-c("Response","Factor1")
														
														
														
														data_mat_anova$Factor1<-as.factor(data_mat_anova$Factor1)
														
														anova_res<-diffexponewayanova(dataA=data_mat_anova)

															
															
														return(anova_res)
													},classlabels_response_mat)
														
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
														
														#print(summary(pvalues))
														
													if(fdrmethod=="BH"){
															fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
													}else{
												if(fdrmethod=="ST"){
                                                    #fdr_adjust_pvalue<-qvalue(pvalues)
                                                    #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                    
                                                    fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
                                                    
                                                    if(is(fdr_adjust_pvalue,"try-error")){
                                                        
                                                        fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
                                                    }
                                                    
                                                    fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                    
                                                    
												}else{
														if(fdrmethod=="Strimmer"){
															pdf("fdrtool.pdf")
															#par_rows=1
															#par(mfrow=c(par_rows,1))
														fdr_adjust_pvalue<-fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE)
														fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
															try(dev.off(),silent=TRUE)
														}else{
															if(fdrmethod=="none"){
																	#fdr_adjust_pvalue<-pvalues
																	fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
															}else{
																if(fdrmethod=="BY"){
																fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                                                                }else{
                                                                    if(fdrmethod=="bonferroni"){
                                                                        fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                                                                    }
                                                                }
																}
														}
													}
												
											}
											
										if(fdrmethod=="none"){
										filename<-"lm1wayanova_pvalall_posthoc.txt"
										}else{
										filename<-"lm1wayanova_fdrall_posthoc.txt"
										}
										cnames_tab<-colnames(data_m_fc_withfeats)
									
					
										posthoc_names<-colnames(posthoc_pval_mat)
										if(length(posthoc_names)<1){

											posthoc_names<-c("Factor1vs2")
										}																

										cnames_tab<-c("P.value","adjusted.P.value",posthoc_names,cnames_tab)
										
										#cnames_tab<-c("P.value","adjusted.P.value","posthoc.pvalue",cnames_tab)
										
										pvalues<-as.data.frame(pvalues)
										
                                        #pvalues<-t(pvalues)
										
                                        pvalues<-as.data.frame(pvalues)
                                        
                                        final.pvalues<-pvalues
                                        #final.pvalues<-pvalues
										
										data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,posthoc_pval_mat,data_m_fc_withfeats)

										colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
										
										#data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
										filename<-paste("Tables/",filename,sep="")
                                        write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
											
											
									       data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)

										}
										
                                   
                                   
                                        if(featselmethod=="ttest" && pairedanalysis==TRUE)
                                        {
                                            
                                            print("Performing paired t-test analysis")
                                            
                                            #print(dim(data_m_fc))
                                            #print(dim(classlabels_response_mat))
                                            #print(dim(classlabels))
                                            
                                            #data_mat_anova<-cbind(t(data_m_fc),classlabels_response_mat)
                                            
                                            numcores<-round(detectCores()*0.5)
                                            
                                            cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                                            
                                            clusterExport(cl,"t.test",envir = .GlobalEnv)
                                            
                                            
                                            res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat){
                                                
                                                xvec<-x
                                                
                                                
                                                data_mat_anova<-cbind(xvec,classlabels_response_mat)
                                                
                                                data_mat_anova<-as.data.frame(data_mat_anova)
                                                cnames<-colnames(data_mat_anova)
                                                
                                                cnames[1]<-"Response"
                                                
                                                colnames(data_mat_anova)<-c("Response","Factor1")
                                                
                                                #print(data_mat_anova)
                                                
                                                data_mat_anova$Factor1<-as.factor(data_mat_anova$Factor1)
                                                
                                                #anova_res<-diffexponewayanova(dataA=data_mat_anova)
                                                
                                                x1<-data_mat_anova$Response[which(data_mat_anova$Factor1==class_labels_levels[1])]
                                                
                                                y1<-data_mat_anova$Response[which(data_mat_anova$Factor1==class_labels_levels[2])]
                                                
                                                w1<-t.test(x=x1,y=y1,alternative="two.sided",paired=TRUE)
                                                
                                                return(w1$p.value)
                                            },classlabels_response_mat)
                                            
                                            stopCluster(cl)
                                            
                                            main_pval_mat<-{}
                                            
                                            posthoc_pval_mat<-{}
                                            pvalues<-{}
                                            
                                            
                                            
                                            pvalues<-unlist(res1)
                                            
                                            
                                            #print(summary(pvalues))
                                            
                                            if(fdrmethod=="BH"){
                                                fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
                                            }else{
                                                if(fdrmethod=="ST"){
                                                    #fdr_adjust_pvalue<-qvalue(pvalues)
                                                    #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                    
                                                    fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
                                                    
                                                    if(is(fdr_adjust_pvalue,"try-error")){
                                                        
                                                        fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
                                                    }
                                                    
                                                    fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                    
                                                    
                                                }else{
                                                    if(fdrmethod=="Strimmer"){
                                                        pdf("fdrtool.pdf")
                                                        #par_rows=1
                                                        #par(mfrow=c(par_rows,1))
                                                        fdr_adjust_pvalue<-fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE)
                                                        fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
                                                        try(dev.off(),silent=TRUE)
                                                    }else{
                                                        if(fdrmethod=="none"){
                                                            #fdr_adjust_pvalue<-pvalues
                                                            fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
                                                        }else{
                                                            if(fdrmethod=="BY"){
                                                                fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                                                            }else{
                                                                if(fdrmethod=="bonferroni"){
                                                                    fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                                
                                            }
                                            
                                            if(fdrmethod=="none"){
                                                filename<-"pairedttest_pvalall_withfeats.txt"
                                            }else{
                                                filename<-"pairedttest_fdrall_withfeats.txt"
                                            }
                                            cnames_tab<-colnames(data_m_fc_withfeats)
                                            
                                            
                                            posthoc_names<-colnames(posthoc_pval_mat)
                                            if(length(posthoc_names)<1){
                                                
                                                posthoc_names<-c("Factor1vs2")
                                            }
                                            
                                            cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
                                            
                                            #cnames_tab<-c("P.value","adjusted.P.value","posthoc.pvalue",cnames_tab)
                                            
                                            pvalues<-as.data.frame(pvalues)
                                            
                                            #pvalues<-t(pvalues)
                                            # print(dim(pvalues))
                                            #print(dim(data_m_fc_withfeats))
                                            
                                            final.pvalues<-pvalues
                                            sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
                                            
                                            
                                            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
                                            
                                            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
                                            
                                            #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
                                            #  write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
                                            
                                            
                                            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
                                            
                                        }
                                        
                                        if(featselmethod=="ttest" && pairedanalysis==FALSE)
                                            {
                                                
                                                print("Performing t-test analysis")
                                                
                                                #print(dim(data_m_fc))
                                                #print(dim(classlabels_response_mat))
                                                #print(dim(classlabels))
                                                
                                                #data_mat_anova<-cbind(t(data_m_fc),classlabels_response_mat)
                                                
                                                numcores<-round(detectCores()*0.5)
                                                
                                                cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                                                
                                                clusterExport(cl,"t.test",envir = .GlobalEnv)
                                                
                                                
                                                res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat){
                                                    
                                                    xvec<-x
                                                    
                                                    
                                                    data_mat_anova<-cbind(xvec,classlabels_response_mat)
                                                    
                                                    data_mat_anova<-as.data.frame(data_mat_anova)
                                                    cnames<-colnames(data_mat_anova)
                                                    
                                                    cnames[1]<-"Response"
                                                    
                                                    colnames(data_mat_anova)<-c("Response","Factor1")
                                                    
                                                    #print(data_mat_anova)
                                                    
                                                    data_mat_anova$Factor1<-as.factor(data_mat_anova$Factor1)
                                                    
                                                    #anova_res<-diffexponewayanova(dataA=data_mat_anova)
                                                    
                                                    x1<-data_mat_anova$Response[which(data_mat_anova$Factor1==class_labels_levels[1])]
                                                    
                                                    y1<-data_mat_anova$Response[which(data_mat_anova$Factor1==class_labels_levels[2])]
                                                    
                                                    w1<-t.test(x=x1,y=y1,alternative="two.sided")
                                                    
                                                    return(w1$p.value)
                                                },classlabels_response_mat)
                                                
                                                stopCluster(cl)
                                                
                                                main_pval_mat<-{}
                                                
                                                posthoc_pval_mat<-{}
                                                pvalues<-{}
                                                
                                                
                                                
                                                pvalues<-unlist(res1)
                                                
                                                
                                                #print(summary(pvalues))
                                                
                                                if(fdrmethod=="BH"){
                                                    fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
                                                }else{
                                                    if(fdrmethod=="ST"){
                                                        #fdr_adjust_pvalue<-qvalue(pvalues)
                                                        #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                        
                                                        fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
                                                        
                                                        if(is(fdr_adjust_pvalue,"try-error")){
                                                            
                                                            fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
                                                        }
                                                        
                                                        fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                        
                                                        
                                                    }else{
                                                        if(fdrmethod=="Strimmer"){
                                                            pdf("fdrtool.pdf")
                                                            #par_rows=1
                                                            #par(mfrow=c(par_rows,1))
                                                            fdr_adjust_pvalue<-fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE)
                                                            fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
                                                            try(dev.off(),silent=TRUE)
                                                        }else{
                                                            if(fdrmethod=="none"){
                                                                #fdr_adjust_pvalue<-pvalues
                                                                fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
                                                            }else{
                                                                if(fdrmethod=="BY"){
                                                                    fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                                                                }else{
                                                                    if(fdrmethod=="bonferroni"){
                                                                        fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                    
                                                }
                                                
                                                if(fdrmethod=="none"){
                                                    filename<-"ttest_pvalall_withfeats.txt"
                                                }else{
                                                    filename<-"ttest_fdrall_withfeats.txt"
                                                }
                                                cnames_tab<-colnames(data_m_fc_withfeats)
                                                
                                                
                                                posthoc_names<-colnames(posthoc_pval_mat)
                                                if(length(posthoc_names)<1){
                                                    
                                                    posthoc_names<-c("Factor1vs2")
                                                }
                                                
                                                cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
                                                
                                                #cnames_tab<-c("P.value","adjusted.P.value","posthoc.pvalue",cnames_tab)
                                                
                                                pvalues<-as.data.frame(pvalues)
                                                
                                                #pvalues<-t(pvalues)
                                                # print(dim(pvalues))
                                                #print(dim(data_m_fc_withfeats))
                                                
                                                final.pvalues<-pvalues
                                                sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
                                                
                                                
                                                data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
                                                
                                                colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
                                                
                                                #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
                                                #  write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
                                                
                                                
                                                data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
                                                
                                            }
                                            
                                        
                                        if(featselmethod=="wilcox")
                                        {
                                            
                                            print("Performing Wilcox rank-sum analysis")
                                            
                                            #print(dim(data_m_fc))
                                            #print(dim(classlabels_response_mat))
                                            #print(dim(classlabels))
                                            
                                            #data_mat_anova<-cbind(t(data_m_fc),classlabels_response_mat)
                                            
                                            numcores<-round(detectCores()*0.5)
                                            
                                            cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                                            
                                            clusterExport(cl,"wilcox.test",envir = .GlobalEnv)
                                            
                                           
                                           res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat){
                                                
                                                xvec<-x
                                                
                                                
                                                data_mat_anova<-cbind(xvec,classlabels_response_mat)
                                                
                                                data_mat_anova<-as.data.frame(data_mat_anova)
                                                cnames<-colnames(data_mat_anova)
                                                
                                                cnames[1]<-"Response"
                                                
                                                colnames(data_mat_anova)<-c("Response","Factor1")
                                                
                                                #print(data_mat_anova)
                                                
                                                data_mat_anova$Factor1<-as.factor(data_mat_anova$Factor1)
                                                
                                                #anova_res<-diffexponewayanova(dataA=data_mat_anova)
                                                
                                                x1<-data_mat_anova$Response[which(data_mat_anova$Factor1==class_labels_levels[1])]
                                                
                                                y1<-data_mat_anova$Response[which(data_mat_anova$Factor1==class_labels_levels[2])]
                                                
                                                w1<-wilcox.test(x=x1,y=y1,alternative="two.sided")
                                                
                                                return(w1$p.value)
                                            },classlabels_response_mat)
                                           
                                            stopCluster(cl)
                                            
                                            main_pval_mat<-{}
                                            
                                            posthoc_pval_mat<-{}
                                            pvalues<-{}
                                            
                                           
                                            
                                            pvalues<-unlist(res1)
                                            
                                            
                                            #print(summary(pvalues))
                                            
                                            if(fdrmethod=="BH"){
                                                fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
                                            }else{
                                                if(fdrmethod=="ST"){
                                                    #fdr_adjust_pvalue<-qvalue(pvalues)
                                                    #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                    
                                                    fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
                                                    
                                                    if(is(fdr_adjust_pvalue,"try-error")){
                                                        
                                                        fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
                                                    }
                                                    
                                                    fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                    
                                                    
                                                }else{
                                                    if(fdrmethod=="Strimmer"){
                                                        pdf("fdrtool.pdf")
                                                        #par_rows=1
                                                        #par(mfrow=c(par_rows,1))
                                                        fdr_adjust_pvalue<-fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE)
                                                        fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
                                                        try(dev.off(),silent=TRUE)
                                                    }else{
                                                        if(fdrmethod=="none"){
                                                            #fdr_adjust_pvalue<-pvalues
                                                            fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
                                                        }else{
                                                            if(fdrmethod=="BY"){
                                                                fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                                                            }else{
                                                                if(fdrmethod=="bonferroni"){
                                                                    fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                                
                                            }
                                            
                                            if(fdrmethod=="none"){
                                                filename<-"wilcox_pvalall_withfeats.txt"
                                            }else{
                                                filename<-"wilcox_fdrall_withfeats.txt"
                                            }
                                            cnames_tab<-colnames(data_m_fc_withfeats)
                                            
                                            
                                            posthoc_names<-colnames(posthoc_pval_mat)
                                            if(length(posthoc_names)<1){
                                                
                                                posthoc_names<-c("Factor1vs2")
                                            }																
                                            
                                            cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
                                            
                                            #cnames_tab<-c("P.value","adjusted.P.value","posthoc.pvalue",cnames_tab)
                                            
                                            pvalues<-as.data.frame(pvalues)
                                            
                                           
                                            final.pvalues<-pvalues
                                            sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
                                            
                                            
                                            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
                                            
                                            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
                                            
                                            #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
                                            #  write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
                                            
                                            
                                            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
                                            
                                        }
                                        
                                        
										if(featselmethod=="lmreg")
										{
										
                                            if(logistic_reg==TRUE){
                                               
                                               
                                               if(length(levels(classlabels_response_mat[,1]))>2){
                                                   
                                                    print("More than 2 classes found. Skipping logistic regression analysis.")
                                                   next;
                                               }
                                               
                                                print("Performing logistic regression analysis:")
                                               
                                                classlabels_response_mat[,1]<-as.numeric((classlabels_response_mat[,1]))-1
                                                
                                                fileheader="logitreg"
                                            }else{
                                                
                                                if(poisson_reg==TRUE){
                                                    
                                                    
                                                    print("Performing poisson regression analysis")
                                                    fileheader="poissonreg"
                                                    classlabels_response_mat[,1]<-as.numeric((classlabels_response_mat[,1]))
                                                    
                                                }else{
                                                    print("Performing linear regression analysis:")
                                                    fileheader="lmreg"
                                                }
                                            }


   
   numcores<-num_nodes #round(detectCores()*0.5)
   
   cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
   
   clusterExport(cl,"diffexplmreg",envir = .GlobalEnv)
  clusterExport(cl,"lm",envir = .GlobalEnv)
  clusterExport(cl,"glm",envir = .GlobalEnv)
   clusterExport(cl,"summary",envir = .GlobalEnv)
   clusterExport(cl,"anova",envir = .GlobalEnv)

											#data_mat_anova<-cbind(t(data_m_fc),classlabels_response_mat)
													res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat,logistic_reg,poisson_reg){
												
														xvec<-x
														
                                                        
                                                        if(dim(classlabels_response_mat)[2]>1){
                                                            
                                                            #for(cnum in 2:dim(classlabels_response_mat)[2]){
                                                                
                                                                # classlabels_response_mat[,cnum]<-as.numeric(classlabels_response_mat[,cnum])
                                                            
                                                            #}
                                                        }
                                                        
														data_mat_anova<-cbind(xvec,classlabels_response_mat)
														
														cnames<-colnames(data_mat_anova)
														cnames[1]<-"Response"
														
														colnames(data_mat_anova)<-cnames
														
														
														anova_res<-diffexplmreg(dataA=data_mat_anova,logistic_reg,poisson_reg)
														
														return(anova_res)
													},classlabels_response_mat,logistic_reg,poisson_reg)
														
													stopCluster(cl)
													main_pval_mat<-{}
														
													posthoc_pval_mat<-{}
													pvalues<-{}
													
													all_inf_mat<-{}

													for(i in 1:length(res1)){
														
														main_pval_mat<-rbind(main_pval_mat,res1[[i]]$mainpvalues)
														pvalues<-c(pvalues,res1[[i]]$mainpvalues[1])
														#posthoc_pval_mat<-rbind(posthoc_pval_mat,res1[[i]]$posthocfactor1)
														
														cur_pvals<-t(res1[[i]]$mainpvalues)
														cur_est<-t(res1[[i]]$estimates)
														cur_stderr<-t(res1[[i]]$stderr)
														cur_tstat<-t(res1[[i]]$statistic)
														
														#cur_pvals<-as.data.frame(cur_pvals)
														
														cur_res<-cbind(cur_pvals,cur_est,cur_stderr,cur_tstat)
														
														all_inf_mat<-rbind(all_inf_mat,cur_res)
														
													
													}
													
                                                    
                                                    
                                                    cnames_1<-c(paste("P.value_",colnames(cur_pvals),sep=""),paste("Estimate_",colnames(cur_pvals),sep=""),paste("StdError_var_",colnames(cur_pvals),sep=""),paste("t-statistic_",colnames(cur_pvals),sep=""))
                                                    
                                                    
													
													#	print("here after lm reg")
														
														#print(summary(pvalues))
														
													if(fdrmethod=="BH"){
															fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
													}else{
												if(fdrmethod=="ST"){
                                                    #fdr_adjust_pvalue<-qvalue(pvalues)
                                                    #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                    
                                                    fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
                                                    
                                                    if(is(fdr_adjust_pvalue,"try-error")){
                                                        
                                                        fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
                                                    }
                                                    
                                                    fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                    
                                                    
												}else{
														if(fdrmethod=="Strimmer"){
															pdf("fdrtool.pdf")
															#par_rows=1
															#par(mfrow=c(par_rows,1))
														fdr_adjust_pvalue<-fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE)
														fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
															try(dev.off(),silent=TRUE)
														}else{
															if(fdrmethod=="none"){
																	#fdr_adjust_pvalue<-pvalues
																	fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
															}else{
																if(fdrmethod=="BY"){
																fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                                                                }else{
                                                                    if(fdrmethod=="bonferroni"){
                                                                        fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                                                                    }
                                                                }
																}
														}
													}
												
												
												
												}
												
												
												if(fdrmethod=="none"){
										filename<-paste(fileheader,"_pvalall_withfeats.txt",sep="")
                                        
										}else{
										filename<-paste(fileheader,"_fdrall_withfeats.txt",sep="")
										}
										cnames_tab<-colnames(data_m_fc_withfeats)
										cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
										
										pvalues<-as.data.frame(pvalues)
										
										#pvalues<-t(pvalues)
										#print(dim(pvalues))
										#print(dim(data_m_fc_withfeats))
										
                                       
                                        
										
                                        final.pvalues<-pvalues
										sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
										
                                       
										data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)

										colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
										
										#data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
                                        #write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
											
										
                                        if(analysismode=="regression"){
                                            
                                        filename<-paste(fileheader,"_results_allfeatures.txt",sep="")
					filename<-paste("Tables/",filename,sep="")

                                        write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
                                        }
                                        
                                        filename<-paste(fileheader,"_pval_coef_stderr.txt",sep="")
                                        


										data_allinf_withfeats<-cbind(all_inf_mat,data_m_fc_withfeats)
                                        filename<-paste("Tables/",filename,sep="")
                                     #   write.table(data_allinf_withfeats, file=filename,sep="\t",row.names=FALSE)


										cnames_tab<-colnames(data_m_fc_withfeats)
										
                                      
                                        cnames_tab<-c(cnames_1,cnames_tab)
                                        
                                           class_column_names<-colnames(classlabels_response_mat)
                                           

										colnames(data_allinf_withfeats)<-as.character(cnames_tab)
										
										
										write.table(data_allinf_withfeats, file=filename,sep="\t",row.names=FALSE)
											
											
											
										}
											
											
											
													if(featselmethod=="lm2wayanova")
													{
										
														print("Performing two-way ANOVA analysis with Tukey post hoc comparisons")
											
                                            #print(dim(data_m_fc))
                                            #			print(dim(classlabels_response_mat))
		
        
        numcores<-num_nodes #round(detectCores()*0.5)
        
        cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
        
        
                                                            clusterExport(cl,"diffexplmtwowayanova",envir = .GlobalEnv)
                                                            
                                                            clusterExport(cl,"TukeyHSD",envir = .GlobalEnv)
                                                            clusterExport(cl,"aov",envir = .GlobalEnv)
                                                            clusterExport(cl,"anova",envir = .GlobalEnv)


                                                        #res1<-apply(data_m_fc,1,function(x){
                                                            
                                                            res1<-parRapply(cl,data_m_fc,function(x,classlabels_response_mat){
												
														xvec<-x
														
														colnames(classlabels_response_mat)<-paste("Factor",seq(1,dim(classlabels_response_mat)[2]),sep="")
														
														data_mat_anova<-cbind(xvec,classlabels_response_mat)
                                                        #print("2way anova")
                                                        #	print(data_mat_anova[1:2,])
														cnames<-colnames(data_mat_anova)
														cnames[1]<-"Response"
														
														colnames(data_mat_anova)<-cnames
													
														##save(data_mat_anova,file="data_mat_anova.Rda")
	
														#diffexplmtwowayanova
														anova_res<-diffexplmtwowayanova(dataA=data_mat_anova)
														
														
														return(anova_res)
														},classlabels_response_mat)
														
														
                                                        stopCluster(cl)
                                                        #	print("done")

##save(res1,file="res1.Rda")
            
            main_pval_mat<-{}
																		posthoc_pval_mat<-{}
																		pvalues1<-{}
				pvalues2<-{}
				pvalues3<-{}
                
                
																		for(i in 1:length(res1)){
																			
																			#print(i)
																			#print(res1[[i]]$mainpvalues)
																			#print(res1[[i]]$posthoc)
																			main_pval_mat<-rbind(main_pval_mat,res1[[i]]$mainpvalues)
																			pvalues1<-c(pvalues1,as.vector(res1[[i]]$mainpvalues[1]))
																			pvalues2<-c(pvalues2,as.vector(res1[[i]]$mainpvalues[2]))
																			pvalues3<-c(pvalues3,as.vector(res1[[i]]$mainpvalues[3]))
																			posthoc_pval_mat<-rbind(posthoc_pval_mat,res1[[i]]$posthoc)
																		}
                                                                        twoanova_res<-cbind(data_m_fc_withfeats[,c(1:2)],main_pval_mat,posthoc_pval_mat)
                                                                        
                                                                         write.table(twoanova_res,file="Tables/twoanova_with_posthoc_pvalues.txt",sep="\t",row.names=FALSE)
                                                                         pvalues1<-main_pval_mat[,1]
                                                                         pvalues2<-main_pval_mat[,2]
                                                                         pvalues3<-main_pval_mat[,3]
                                                                         
                                                                        if(fdrmethod=="none"){
                                                                            fdr_adjust_pvalue1<-p.adjust(pvalues1,method="none")
                                                                            fdr_adjust_pvalue2<-p.adjust(pvalues2,method="none")
                                                                            fdr_adjust_pvalue3<-p.adjust(pvalues3,method="none")
                                                                        }
                                                                
																		if(fdrmethod=="BH"){
																				fdr_adjust_pvalue1<-p.adjust(pvalues1,method="BH")
																				fdr_adjust_pvalue2<-p.adjust(pvalues2,method="BH")
																				fdr_adjust_pvalue3<-p.adjust(pvalues3,method="BH")
													}else{
												if(fdrmethod=="ST"){
                                                   
                                                   
                                                   fdr_adjust_pvalue1<-try(qvalue(pvalues1),silent=TRUE)
                                                   
                                                   if(is(fdr_adjust_pvalue1,"try-error")){
                                                       
                                                       fdr_adjust_pvalue1<-qvalue(pvalues1,lambda=max(pvalues1,na.rm=TRUE))
                                                   }
                                                   
                                                   fdr_adjust_pvalue1<-fdr_adjust_pvalue1$qvalues
                                                   
                                                   fdr_adjust_pvalue2<-try(qvalue(pvalues2),silent=TRUE)
                                                   
                                                   if(is(fdr_adjust_pvalue2,"try-error")){
                                                       
                                                       fdr_adjust_pvalue2<-qvalue(pvalues2,lambda=max(pvalues2,na.rm=TRUE))
                                                   }
                                                   
                                                   fdr_adjust_pvalue2<-fdr_adjust_pvalue2$qvalues
                                                   
                                                   fdr_adjust_pvalue3<-try(qvalue(pvalues3),silent=TRUE)
                                                   
                                                   if(is(fdr_adjust_pvalue3,"try-error")){
                                                       
                                                       fdr_adjust_pvalue3<-qvalue(pvalues3,lambda=max(pvalues3,na.rm=TRUE))
                                                   }
                                                   
                                                   fdr_adjust_pvalue3<-fdr_adjust_pvalue3$qvalues
                                                   
                                                   
												}else{
														if(fdrmethod=="Strimmer"){
															pdf("fdrtool.pdf")
															#par_rows=1
															#par(mfrow=c(par_rows,1))
														fdr_adjust_pvalue1<-fdrtool(as.vector(pvalues1),statistic="pvalue",verbose=FALSE)
														fdr_adjust_pvalue1<-fdr_adjust_pvalue1$qval
														
														fdr_adjust_pvalue2<-fdrtool(as.vector(pvalues2),statistic="pvalue",verbose=FALSE)
														fdr_adjust_pvalue2<-fdr_adjust_pvalue2$qval
														
														fdr_adjust_pvalue3<-fdrtool(as.vector(pvalues3),statistic="pvalue",verbose=FALSE)
														fdr_adjust_pvalue3<-fdr_adjust_pvalue3$qval
															try(dev.off(),silent=TRUE)
														}else{
															if(fdrmethod=="none"){
                                                                fdr_adjust_pvalue1<-p.adjust(pvalues1,method="none")
                                                                fdr_adjust_pvalue2<-p.adjust(pvalues2,method="none")
                                                                fdr_adjust_pvalue3<-p.adjust(pvalues3,method="none")

																	
															}else{
																if(fdrmethod=="BY"){
											fdr_adjust_pvalue1<-p.adjust(pvalues1,method="BY")
											fdr_adjust_pvalue2<-p.adjust(pvalues2,method="BY")
											fdr_adjust_pvalue3<-p.adjust(pvalues3,method="BY")

											
                                                                }else{
                                                                    if(fdrmethod=="bonferroni"){
                                                                        # fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                                                                        fdr_adjust_pvalue1<-p.adjust(pvalues1,method="bonferroni")
                                                                        fdr_adjust_pvalue2<-p.adjust(pvalues2,method="bonferroni")
                                                                        fdr_adjust_pvalue3<-p.adjust(pvalues3,method="bonferroni")
                                                                    }
                                                                }
                                                                
																}
														}
													}
										
											
											
													}
													
													if(fdrmethod=="none"){
										filename<-"lm2wayanova_pvalall_withfeats.txt"
										}else{
										filename<-"lm2wayanova_fdrall_withfeats.txt"
										}
										cnames_tab<-colnames(data_m_fc_withfeats)
										
										posthoc_names<-colnames(posthoc_pval_mat)
										
                                        cnames_tab<-c("Factor1.P.value","Factor1.adjusted.P.value","Factor2.P.value","Factor2.adjusted.P.value","Interact.P.value","Interact.adjusted.P.value",posthoc_names,cnames_tab)
										
                                        if(FALSE)
                                        {
                                        pvalues1<-as.data.frame(pvalues1)
                                        pvalues1<-t(pvalues1)
                                        fdr_adjust_pvalue1<-as.data.frame(fdr_adjust_pvalue1)
                                         pvalues2<-as.data.frame(pvalues2)
                                        pvalues2<-t(pvalues2)
                                        fdr_adjust_pvalue2<-as.data.frame(fdr_adjust_pvalue2)
                                         pvalues3<-as.data.frame(pvalues3)
                                        pvalues3<-t(pvalues3)
                                        fdr_adjust_pvalue3<-as.data.frame(fdr_adjust_pvalue3)
                                        posthoc_pval_mat<-as.data.frame(posthoc_pval_mat)
                                        }
                                        
                                        # #save(data_m_fc_withfeats,file="data_m_fc_withfeats.Rda")
                                        data_limma_fdrall_withfeats<-cbind(pvalues1,fdr_adjust_pvalue1,pvalues2,fdr_adjust_pvalue2,pvalues3,fdr_adjust_pvalue3,posthoc_pval_mat,data_m_fc_withfeats)

fdr_adjust_pvalue<-cbind(fdr_adjust_pvalue1,fdr_adjust_pvalue2,fdr_adjust_pvalue3)
fdr_adjust_pvalue<-apply(fdr_adjust_pvalue,1,function(x){min(x,na.rm=TRUE)})


                                                                                colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)

                                                                                #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
										
										filename<-paste("Tables/",filename,sep="")
                                                                                write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)

 
                                        fdr_matrix<-cbind(fdr_adjust_pvalue1,fdr_adjust_pvalue2,fdr_adjust_pvalue3)
                                        
                                        fdr_matrix<-as.data.frame(fdr_matrix)
                                        
                                        fdr_adjust_pvalue_all<-apply(fdr_matrix,1,function(x){return(min(x,na.rm=TRUE)[1])})
                                        
                                        pvalues<-cbind(pvalues1,pvalues2,pvalues3)
                                        pvalues<-apply(pvalues,1,function(x){min(x,na.rm=TRUE)[1]})
                                        
										#pvalues1<-t(pvalues1)
										
										#print("here")
										pvalues1<-as.data.frame(pvalues1)
										pvalues1<-t(pvalues1)
										#print(dim(pvalues1))
										
										#pvalues2<-t(pvalues2)
										pvalues2<-as.data.frame(pvalues2)
										pvalues2<-t(pvalues2)
										
										#pvalues3<-t(pvalues3)
										pvalues3<-as.data.frame(pvalues3)
										pvalues3<-t(pvalues3)
										
										#pvalues<-t(pvalues)
										#print(dim(pvalues1))
										#print(dim(pvalues2))
										#print(dim(pvalues3))
										#print(dim(data_m_fc_withfeats))
                                        
                                        final.pvalues<-pvalues
                                        
                                        
										
										sel.diffdrthresh<-fdr_adjust_pvalue_all<fdrthresh & final.pvalues<pvalue.thresh
                                        
                                        if(length(which(fdr_adjust_pvalue1<fdrthresh))>0){
                                            
                                            
                                           
                                            
                                        X1=data_m_fc_withfeats[which(fdr_adjust_pvalue1<fdrthresh),]
                                        Y1=cbind(classlabels_orig[,1],as.character(classlabels_response_mat[,1]))
                                        Y1<-as.data.frame(Y1)

                                    #save(classlabels_orig,file="classlabels_orig.Rda")
                                    #save(classlabels_response_mat,file="classlabels_response_mat.Rda")
                                    #save(Y1,file="Y1.Rda")


if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/HCA_Factor1selectedfeats.png"
    
    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
}

get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X1,Y=Y1,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                                        sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="Factor1")
                                        
                                        if(output.device.type!="pdf"){
                                            
                                           try(dev.off(),silent=TRUE)
                                        }
                                        }else{
                                            
                                            print("No significant features for Factor 1.")
                                        }

if(length(which(fdr_adjust_pvalue2<fdrthresh))>0){
                                        X2=data_m_fc_withfeats[which(fdr_adjust_pvalue2<fdrthresh),]
                                        
                                        
                                        Y2=cbind(classlabels_orig[,1],as.character(classlabels_response_mat[,2]))
                                        Y2<-as.data.frame(Y2)
                                        
                                        #save(classlabels_orig,file="classlabels_orig.Rda")
                                        #save(classlabels_response_mat,file="classlabels_response_mat.Rda")
                                        #save(Y2,file="Y2.Rda")
                                        
                                        if(output.device.type!="pdf"){
                                            
                                            temp_filename_1<-"Figures/HCA_Factor2selectedfeats.png"
                                            
                                            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                        }
                                        
                                        
                                        get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X2,Y=Y2,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                                        sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="Factor2")
                                        

                                        if(output.device.type!="pdf"){
                                            
                                            try(dev.off(),silent=TRUE)
                                        }
                                        
}else{
    
    print("No significant features for Factor 2.")
}
                                        class_interact<-paste(classlabels_response_mat[,1],":",classlabels_response_mat[,2],sep="")
                                        #classlabels_response_mat[,1]:classlabels_response_mat[,2]

if(length(which(fdr_adjust_pvalue3<fdrthresh))>0){
                                        X3=data_m_fc_withfeats[which(fdr_adjust_pvalue3<fdrthresh),]
                                        Y3=cbind(classlabels_orig[,1],class_interact)
                                        Y3<-as.data.frame(Y3)
                                        
                                        
                                        if(output.device.type!="pdf"){
                                            
                                            temp_filename_1<-"Figures/HCA_Factor1xFactor2selectedfeats.png"
                                            
                                            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                        }
                                        
                                        get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X3,Y=Y3,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                                        sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="Factor1 x Factor2")
                                        

                                        if(output.device.type!="pdf"){
                                            
                                            try(dev.off(),silent=TRUE)
                                        }


}else{
    
    print("No significant features for the interaction.")
}



                                    data_limma_fdrall_withfeats<-cbind(final.pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
                                    
                                    cnames_tab<-colnames(data_m_fc_withfeats)
                                    cnames_tab<-c("P.value.Min(Factor1,Factor2,Interaction)","adjusted.P.value.Min(Factor1,Factor2,Interaction)",cnames_tab)
                                    colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)


                                    
                                    #filename2<-"test2.txt"
                                    #data_limma_fdrsig_withfeats<-data_limma_fdrall_withfeats[sel.diffdrthresh==TRUE,]
                                    #write.table(data_limma_fdrsig_withfeats, file=filename2,sep="\t",row.names=FALSE)
                                    
                                    fdr_adjust_pvalue<-fdr_adjust_pvalue_all
                                    
									}
												
														if(featselmethod=="lm1wayanovarepeat"){
														
														print("Performing one-way ANOVA with repeated measurements analysis using nlme::lme()")
														
                                                        
                                                        
                                                        
                                                        
                                                        numcores<-num_nodes #round(detectCores()*0.5)
                                                        
                                                        cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                                                        
                                                        clusterExport(cl,"diffexplmonewayanovarepeat",envir = .GlobalEnv)
                                                        clusterEvalQ(cl,library(nlme))
                                                        clusterEvalQ(cl,library(multcomp))
                                                        clusterEvalQ(cl,library(lsmeans))
                                                        clusterExport(cl,"lme",envir = .GlobalEnv)
                                                        clusterExport(cl,"interaction",envir = .GlobalEnv)
                                                        clusterExport(cl,"anova",envir = .GlobalEnv)
                                                        #clusterExport(cl,"classlabels_response_mat",envir = .GlobalEnv)
                                                        #clusterExport(cl,"subject_inf",envir = .GlobalEnv)
                                                        
                                                        #res1<-apply(data_m_fc,1,function(x){
                                                        
                                                        res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat,subject_inf,modeltype){
                                                            
                                                            #res1<-apply(data_m_fc,1,function(x){
												
														xvec<-x
														
														colnames(classlabels_response_mat)<-paste("Factor",seq(1,dim(classlabels_response_mat)[2]),sep="")
														
														data_mat_anova<-cbind(xvec,classlabels_response_mat)
														
														cnames<-colnames(data_mat_anova)
														cnames[1]<-"Response"
														
														colnames(data_mat_anova)<-cnames
														
														anova_res<-diffexplmonewayanovarepeat(dataA=data_mat_anova,subject_inf=subject_inf,modeltype=modeltype)
														
														return(anova_res)
                                                        
                                                        
														},classlabels_response_mat,subject_inf,modeltype)
                                                        
                                                        stopCluster(cl)
                                                        
														main_pval_mat<-{}
														pvalues<-{}

                                                         #save(res1,file="lmres1.Rda")
														
																		posthoc_pval_mat<-{}
                                                                        
                                                                        bad_lm1feats<-{}

																		for(i in 1:length(res1)){
																			
                                                                            #print(i)
                                                                            #print(res1)
																			
                                                                            if(is.na(res1[[i]]$mainpvalues)==FALSE){
                                                                            
																			main_pval_mat<-rbind(main_pval_mat,res1[[i]]$mainpvalues)
																			pvalues<-c(pvalues,res1[[i]]$mainpvalues[1])
																			posthoc_pval_mat<-rbind(posthoc_pval_mat,res1[[i]]$posthoc)
                                                                           
                                                                            
                                                                            }else{
                                                                                
                                                                                 bad_lm1feats<-c(bad_lm1feats,i)
                                                                                
                                                                            }
																		}
																 
                                                                         pvalues<-unlist(pvalues)
																		if(fdrmethod=="BH"){
																				fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
													}else{
												if(fdrmethod=="ST"){
                                                    #fdr_adjust_pvalue<-qvalue(pvalues)
                                                    #fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                    
                                                   
                                                   
                                                  
                                                    fdr_adjust_pvalue<-try(qvalue(pvalues),silent=TRUE)
                                                    
                                                    if(is(fdr_adjust_pvalue,"try-error")){
                                                        
                                                        fdr_adjust_pvalue<-qvalue(pvalues,lambda=max(pvalues,na.rm=TRUE))
                                                    }
                                                    
                                                    fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
                                                    
                                                    
												}else{
														if(fdrmethod=="Strimmer"){
															pdf("fdrtool.pdf")
															#par_rows=1
															#par(mfrow=c(par_rows,1))
														fdr_adjust_pvalue<-fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE)
														fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
															try(dev.off(),silent=TRUE)
														}else{
															if(fdrmethod=="none"){
																	fdr_adjust_pvalue<-pvalues
																	
															}else{
																if(fdrmethod=="BY"){
											fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                                                                }else{
                                                                    if(fdrmethod=="bonferroni"){
                                                                         fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                                                                        
                                                                    }
                                                                }
																}
														}
													}
										
											
											
													}
													if(fdrmethod=="none"){
										filename<-"lm1wayanovarepeat_pvalall_withfeats.txt"
										}else{
										filename<-"lm1wayanovarepeat_fdrall_withfeats.txt"
										}
										cnames_tab<-colnames(data_m_fc_withfeats)
										
										posthoc_names<-colnames(posthoc_pval_mat)
										
										cnames_tab<-c("P.value","adjusted.P.value",posthoc_names,cnames_tab)
										
										#cnames_tab<-c("P.value","adjusted.P.value","posthoc.pvalue",cnames_tab)
										
										#pvalues<-t(pvalues)
										pvalues<-as.data.frame(pvalues)
                                        
                                        final.pvalues<-pvalues
										sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
                                        #pvalues<-t(pvalues)
										#print(dim(pvalues))
										#print(dim(data_m_fc_withfeats))
                                        
                                        if(length(bad_lm1feats)>0){
                                        data_m_fc_withfeats<-data_m_fc_withfeats[-c(bad_lm1feats),]
                                        
                                        data_m_fc<-data_m_fc[-c(bad_lm1feats),]
                                        }
                                        
										data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,posthoc_pval_mat,data_m_fc_withfeats)
                                        

										
										
										colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
										
										#data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
										
										filename<-paste("Tables/",filename,sep="")
                                        write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
											
											data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)

		
									}
													
																if(featselmethod=="lm2wayanovarepeat"){
														
														print("Performing two-way ANOVA with repeated measurements analysis using nlme::lme()")
														
                                                        
                                                        numcores<-num_nodes #round(detectCores()*0.5)
                                                        
                                                        cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                                                       
                                                        clusterExport(cl,"diffexplmtwowayanovarepeat",envir = .GlobalEnv)
                                                        clusterEvalQ(cl,library(nlme))
                                                        clusterEvalQ(cl,library(multcomp))
                                                        clusterEvalQ(cl,library(lsmeans))
                                                        clusterExport(cl,"lme",envir = .GlobalEnv)
                                                        clusterExport(cl,"interaction",envir = .GlobalEnv)
                                                        clusterExport(cl,"anova",envir = .GlobalEnv)
                                                        
                                                        #clusterExport(cl,"classlabels_response_mat",envir = .GlobalEnv)
                                                        #clusterExport(cl,"subject_inf",envir = .GlobalEnv)
                                                        
                                                        #res1<-apply(data_m_fc,1,function(x){

#	print(dim(data_m_fc))
#	print(dim(classlabels_response_mat))
                                                        
                                                        res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat,subject_inf,modeltype){
                                                            
                                                            #  res1<-apply(data_m_fc,1,function(x){
												
                                                # #save(classlabels_response_mat,file="classlabels_response_mat.Rda")
                                                        
                                                        #       #save(subject_inf,file="subject_inf.Rda")

                                                        xvec<-x
														
                                                        ##save(xvec,file="xvec.Rda")
														colnames(classlabels_response_mat)<-paste("Factor",seq(1,dim(classlabels_response_mat)[2]),sep="")
														
														data_mat_anova<-cbind(xvec,classlabels_response_mat)
														
														cnames<-colnames(data_mat_anova)
														cnames[1]<-"Response"
														
														colnames(data_mat_anova)<-cnames
														
														#print(subject_inf)
														#print(dim(data_mat_anova))
														
														
														subject_inf<-as.data.frame(subject_inf)
														#print(dim(subject_inf))
														
														anova_res<-diffexplmtwowayanovarepeat(dataA=data_mat_anova,subject_inf=subject_inf[,1],modeltype=modeltype)
														
														return(anova_res)
														},classlabels_response_mat,subject_inf,modeltype)
                                                        
                                                        
														main_pval_mat<-{}
														
                                                        stopCluster(cl)
																		posthoc_pval_mat<-{}

#print(head(res1))
#	print("here")

															pvalues<-{}


                                                                bad_lm1feats<-{}

#save(res1,file="res1.Rda")
																		for(i in 1:length(res1)){
																			
                                                                            if(is.na(res1[[i]]$mainpvalues)==FALSE){
																			main_pval_mat<-rbind(main_pval_mat,res1[[i]]$mainpvalues)
																			pvalues<-c(pvalues,res1[[i]]$mainpvalues[1])
																			posthoc_pval_mat<-rbind(posthoc_pval_mat,res1[[i]]$posthoc)
                                                                            }else{
                                                                                
                                                                                bad_lm1feats<-c(bad_lm1feats,i)
                                                                                
                                                                               
                                                                            }
																		}
                                                                        
                                                                        if(length(bad_lm1feats)>0){
                                                                            
                                                                            data_m_fc_withfeats<-data_m_fc_withfeats[-c(bad_lm1feats),]
                                                                            
                                                                            data_m_fc<-data_m_fc[-c(bad_lm1feats),]
                                                                        }
                                                                        twoanovarepeat_res<-cbind(data_m_fc_withfeats[,c(1:2)],main_pval_mat,posthoc_pval_mat)
																
                                                                write.table(twoanovarepeat_res,file="Tables/2wayanovarepeat_with_posthoc_pvalues.txt",sep="\t",row.names=FALSE)
                                                                

pvalues1<-main_pval_mat[,1]
pvalues2<-main_pval_mat[,2]
pvalues3<-main_pval_mat[,3]

twoanova_res<-cbind(data_m_fc_withfeats[,c(1:2)],main_pval_mat,posthoc_pval_mat)

# write.table(twoanova_res,file="twoanova_with_posthoc_pvalues.txt",sep="\t",row.names=FALSE)

if(fdrmethod=="none"){
    fdr_adjust_pvalue1<-p.adjust(pvalues1,method="none")
    fdr_adjust_pvalue2<-p.adjust(pvalues2,method="none")
    fdr_adjust_pvalue3<-p.adjust(pvalues3,method="none")
}

if(fdrmethod=="BH"){
    fdr_adjust_pvalue1<-p.adjust(pvalues1,method="BH")
    fdr_adjust_pvalue2<-p.adjust(pvalues2,method="BH")
    fdr_adjust_pvalue3<-p.adjust(pvalues3,method="BH")
}else{
    if(fdrmethod=="ST"){
        #print(head(pvalues1))
        #print(head(pvalues2))
        #print(head(pvalues3))
        #print(summary(pvalues1))
        #print(summary(pvalues2))
        #print(summary(pvalues3))
        fdr_adjust_pvalue1<-try(qvalue(pvalues1),silent=TRUE)
        fdr_adjust_pvalue2<-try(qvalue(pvalues2),silent=TRUE)
        fdr_adjust_pvalue3<-try(qvalue(pvalues3),silent=TRUE)
        
        if(is(fdr_adjust_pvalue1,"try-error")){
            
            fdr_adjust_pvalue1<-qvalue(pvalues1,lambda=max(pvalues1,na.rm=TRUE))
        }
        
        if(is(fdr_adjust_pvalue2,"try-error")){
            fdr_adjust_pvalue2<-qvalue(pvalues2,lambda=max(pvalues2,na.rm=TRUE))
        }
        
        if(is(fdr_adjust_pvalue3,"try-error")){
            fdr_adjust_pvalue3<-qvalue(pvalues3,lambda=max(pvalues3,na.rm=TRUE))
        }
        
        
        fdr_adjust_pvalue1<-fdr_adjust_pvalue1$qvalues
        
        fdr_adjust_pvalue2<-fdr_adjust_pvalue2$qvalues
        
        fdr_adjust_pvalue3<-fdr_adjust_pvalue3$qvalues
    }else{
        if(fdrmethod=="Strimmer"){
            pdf("fdrtool.pdf")
            #par_rows=1
            #par(mfrow=c(par_rows,1))
            fdr_adjust_pvalue1<-fdrtool(as.vector(pvalues1),statistic="pvalue",verbose=FALSE)
            fdr_adjust_pvalue1<-fdr_adjust_pvalue1$qval
            
            fdr_adjust_pvalue2<-fdrtool(as.vector(pvalues2),statistic="pvalue",verbose=FALSE)
            fdr_adjust_pvalue2<-fdr_adjust_pvalue2$qval
            
            fdr_adjust_pvalue3<-fdrtool(as.vector(pvalues3),statistic="pvalue",verbose=FALSE)
            fdr_adjust_pvalue3<-fdr_adjust_pvalue3$qval
            try(dev.off(),silent=TRUE)
        }else{
            if(fdrmethod=="none"){
                fdr_adjust_pvalue1<-p.adjust(pvalues1,method="none")
                fdr_adjust_pvalue2<-p.adjust(pvalues2,method="none")
                fdr_adjust_pvalue3<-p.adjust(pvalues3,method="none")
                
                
            }else{
																if(fdrmethod=="BY"){
                                                                    fdr_adjust_pvalue1<-p.adjust(pvalues1,method="BY")
                                                                    fdr_adjust_pvalue2<-p.adjust(pvalues2,method="BY")
                                                                    fdr_adjust_pvalue3<-p.adjust(pvalues3,method="BY")
                                                                    
                                                                    
                                                                }else{
                                                                    if(fdrmethod=="bonferroni"){
                                                                        # fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                                                                        fdr_adjust_pvalue1<-p.adjust(pvalues1,method="bonferroni")
                                                                        fdr_adjust_pvalue2<-p.adjust(pvalues2,method="bonferroni")
                                                                        fdr_adjust_pvalue3<-p.adjust(pvalues3,method="bonferroni")
                                                                    }
                                                                }
            }
        }
    }
    
    
    
}

if(fdrmethod=="none"){
    filename<-"lm2wayanova_pvalall_withfeats.txt"
}else{
    filename<-"lm2wayanova_fdrall_withfeats.txt"
}
cnames_tab<-colnames(data_m_fc_withfeats)

posthoc_names<-colnames(posthoc_pval_mat)

#
cnames_tab<-c("Factor1.P.value","Factor1.adjusted.P.value","Factor2.P.value","Factor2.adjusted.P.value","Interact.P.value","Interact.adjusted.P.value",posthoc_names,cnames_tab)


fdr_matrix<-cbind(fdr_adjust_pvalue1,fdr_adjust_pvalue2,fdr_adjust_pvalue3)

fdr_matrix<-as.data.frame(fdr_matrix)

fdr_adjust_pvalue_all<-apply(fdr_matrix,1,function(x){return(min(x,na.rm=TRUE))})

pvalues_all<-cbind(pvalues1,pvalues2,pvalues3)
pvalue_matrix<-as.data.frame(pvalues_all)

pvalue_all<-apply(pvalue_matrix,1,function(x){return(min(x,na.rm=TRUE)[1])})


#pvalues1<-t(pvalues1)

#print("here")
#pvalues1<-as.data.frame(pvalues1)
#pvalues1<-t(pvalues1)
#print(dim(pvalues1))

#pvalues2<-t(pvalues2)
#pvalues2<-as.data.frame(pvalues2)
#pvalues2<-t(pvalues2)

#pvalues3<-t(pvalues3)
#pvalues3<-as.data.frame(pvalues3)
#pvalues3<-t(pvalues3)

#pvalues<-t(pvalues)
#print(dim(pvalues1))
#print(dim(pvalues2))
#print(dim(pvalues3))
#print(dim(data_m_fc_withfeats))

pvalues<-pvalue_all

final.pvalues<-pvalues
sel.diffdrthresh<-fdr_adjust_pvalue_all<fdrthresh & final.pvalues<pvalue.thresh


if(length(which(fdr_adjust_pvalue1<fdrthresh))>0){
X1=data_m_fc_withfeats[which(fdr_adjust_pvalue1<fdrthresh),]
Y1=cbind(classlabels_orig[,1],as.character(classlabels_response_mat[,1]))
Y1<-as.data.frame(Y1)

#save(classlabels_orig,file="classlabels_orig.Rda")
#save(classlabels_response_mat,file="classlabels_response_mat.Rda")


if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/HCA_Factor1selectedfeats.png"
    
    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
}


get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X1,Y=Y1,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="Factor 1")

if(output.device.type!="pdf"){
    
    try(dev.off(),silent=TRUE)
}
}else{
    print("No significant features for Factor 1.")
    
}

if(length(which(fdr_adjust_pvalue2<fdrthresh))>0){

X2=data_m_fc_withfeats[which(fdr_adjust_pvalue2<fdrthresh),]
Y2=cbind(classlabels_orig[,1],as.character(classlabels_response_mat[,2]))
Y2<-as.data.frame(Y2)


if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/HCA_Factor2selectedfeats.png"
    
    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
}


get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X2,Y=Y2,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="Factor 2")

if(output.device.type!="pdf"){
    
    try(dev.off(),silent=TRUE)
}

}else{
    
    print("No significant features for Factor 2.")
}

class_interact<-paste(classlabels_response_mat[,1],":",classlabels_response_mat[,2],sep="") #classlabels_response_mat[,1]:classlabels_response_mat[,2]

if(length(which(fdr_adjust_pvalue3<fdrthresh))>0){
    
X3=data_m_fc_withfeats[which(fdr_adjust_pvalue3<fdrthresh),]
Y3=cbind(classlabels_orig[,1],class_interact)
Y3<-as.data.frame(Y3)


if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/HCA_Factor1xFactor2selectedfeats.png"
    
    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
}

get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X3,Y=Y3,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type=hca_type,newdevice=FALSE,input.type="intensity",mainlab="Factor 1 x Factor 2")


if(output.device.type!="pdf"){
    
    try(dev.off(),silent=TRUE)
}
}else{
    print("No significant features for Factor 1x2 interaction.")
}


#data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,posthoc_pval_mat,data_m_fc_withfeats)


#
data_limma_fdrall_withfeats<-cbind(pvalues1,fdr_adjust_pvalue1,pvalues2,fdr_adjust_pvalue2,pvalues3,fdr_adjust_pvalue3,posthoc_pval_mat,data_m_fc_withfeats)

fdr_adjust_pvalue<-cbind(fdr_adjust_pvalue1,fdr_adjust_pvalue2,fdr_adjust_pvalue3)
fdr_adjust_pvalue<-apply(fdr_adjust_pvalue,1,function(x){min(x,na.rm=TRUE)})

colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)

#data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)

data_limma_fdrall_withfeats<-cbind(final.pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)

cnames_tab<-colnames(data_m_fc_withfeats)
cnames_tab<-c("P.value.Min(Factor1,Factor2,Interaction)","adjusted.P.value.Min(Factor1,Factor2,Interaction)",cnames_tab)
colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)

#filename2<-"test2.txt"
#data_limma_fdrsig_withfeats<-data_limma_fdrall_withfeats[sel.diffdrthresh==TRUE,]
#write.table(data_limma_fdrsig_withfeats, file=filename2,sep="\t",row.names=FALSE)

fdr_adjust_pvalue<-fdr_adjust_pvalue_all

													}
													
													
													
												
								
								
									
								}
							
							
							
						
							
						
						
							
						
						
					    
					    
					    
					  if(featselmethod=="lmreg" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat"
					  			| featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="logitreg" | featselmethod=="limma2wayrepeat" | featselmethod=="wilcox" | featselmethod=="ttest" |  featselmethod=="poissonreg")
					  {



                      sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh



					   goodip<-which(sel.diffdrthresh==TRUE)
                       

							classlabels<-as.data.frame(classlabels)
						
						
						
								
								if(featselmethod=="limma2way"){
								vennDiagram(results2,cex=0.8)
								}
									
                                    #print(summary(fdr_adjust_pvalue))
                                    #print(summary(final.pvalues))
						
					 }
					
					pred_acc<-0 #("NA")
					
                    
                    
					#print("here")
					feat_sigfdrthresh[lf]<-length(goodip) #which(sel.diffdrthresh==TRUE))
					if(kfold>dim(data_m_fc)[2]){
					kfold=dim(data_m_fc)[2]
					}
					if(analysismode=="classification"){
					
					
					
					#print("classification")
						
					if(length(goodip)>0 & dim(data_m_fc)[2]>=kfold){

                        # print("Time 1")
                        #print(Sys.time())
                
					Targetvar<-classlabels[,1]
                                        dataA<-cbind(Targetvar,t(data_m_fc))
                                        dataA<-as.data.frame(dataA)
                                        #df.summary <- dataA %>% group_by(Targetvar) %>%  summarize_all(funs(mean))

                                        # df.summary <- dataA %>% group_by(Targetvar) %>%  summarize_all(funs(mean))
                                       
                                       df.summary <-aggregate(x=dataA,by=list(dataA$Targetvar),mean)
                                       
                                       df2<-as.data.frame(df.summary[,-c(1:2)])

                                    #  #save(df2,file="df2.Rda")
                                  #    #save(dataA,file="dataA.Rda")
                                   #   #save(Targetvar,file="Targetvar.Rda")
                                        
                                        if(log2transform==TRUE || input.intensity.scale=="log2"){

                                            # foldchangeres<-apply(df2,2,dist)
                                        
                                        #  if(length(nrow(foldchangeres))>0){
                                            
                                            #   foldchangeres<-apply(foldchangeres,2,function(x)
                                            #   {
                                                
                                                #   max_ind<-which(x==max(abs(x)))[1];
                                                #return(x[max_ind])
                                                
                                                #}
                                                #)
                                            
                                            #}
                                         cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                                        foldchangeres<-parApply(cl,df2,2,function(x){
                                            
                                           
                                            res<-lapply(1:length(x),function(i){
                                                return((x[i]-x[-i]))
                                                
                                            })
                                            res<-unlist(res)
                                            
                                            tempres<-abs(res)
                                            res_ind<-which(tempres==max(tempres,na.rm=TRUE))
                                            return(res[res_ind[1]])
                                            
                                        })
                                        
                                        stopCluster(cl)
                                        
                                        
                                        print("Using log2 fold change threshold of")
                                        print(foldchangethresh)


                                        }else{
                                            
                                             #raw intensities
                                            if(znormtransform==FALSE)
                                            {
                                                #   foldchangeres<-apply(log2(df2+1),2,function(x){res<-{};for(i in 1:length(x)){res<-c(res,(x[i]-x[-i]));};tempres<-abs(res);res_ind<-which(tempres==max(tempres,na.rm=TRUE));return(res[res_ind[1]]);})
                                                
                                                if(FALSE){
                                                foldchangeres<-apply(log2(df2+1),2,dist)
                                                
                                                if(length(nrow(foldchangeres))>0){
                                                    
                                                    foldchangeres<-apply(foldchangeres,2,function(x)
                                                    {
                                                        
                                                        max_ind<-which(x==max(abs(x)))[1];
                                                        return(x[max_ind])
                                                        
                                                    }
                                                    )
                                                    
                                                }
                                                
                                                }
                                                
                                                cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                                                foldchangeres<-parApply(cl,log2(df2+1),2,function(x){
                                                    
                                                    
                                                    res<-lapply(1:length(x),function(i){
                                                        return((x[i]-x[-i]))
                                                        
                                                    })
                                                    res<-unlist(res)
                                                    
                                                    tempres<-abs(res)
                                                    res_ind<-which(tempres==max(tempres,na.rm=TRUE))
                                                    return(res[res_ind[1]])
                                                    
                                                })
                                                
                                                stopCluster(cl)
                                                
                                                
                                                foldchangethresh=foldchangethresh
                                                print("Using raw fold change threshold of")
                                                print(foldchangethresh)

                                            }else{
                                                
                                                #  foldchangeres<-apply(df2,2,function(x){res<-{};for(i in 1:length(x)){res<-c(res,(x[i]-(x[-i])));};tempres<-abs(res);res_ind<-which(tempres==max(tempres,na.rm=TRUE));return(res[res_ind[1]]);})
                                                
                                                if(FALSE){
                                                foldchangeres<-apply(df2,2,dist)
                                                
                                                if(length(nrow(foldchangeres))>0){
                                                    
                                                    foldchangeres<-apply(foldchangeres,2,function(x)
                                                    {
                                                        
                                                        max_ind<-which(x==max(abs(x)))[1];
                                                        return(x[max_ind])
                                                        
                                                    }
                                                    )
                                                    
                                                }
                                                
                                               }
                                                
                                                cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                                                foldchangeres<-parApply(cl,df2,2,function(x){
                                                    
                                                    
                                                    res<-lapply(1:length(x),function(i){
                                                        return((x[i]-x[-i]))
                                                        
                                                    })
                                                    res<-unlist(res)
                                                    
                                                    tempres<-abs(res)
                                                    res_ind<-which(tempres==max(tempres,na.rm=TRUE))
                                                    return(res[res_ind[1]])
                                                    
                                                })
                                                
                                                stopCluster(cl)
                                                
                                                #print(summary(foldchangeres))
                                               
                                               
                                                #foldchangethresh=2^foldchangethresh
                                                print("Using Z-score change threshold of")
                                                print(foldchangethresh)
                                                
                                            }
                        
						
                                        }
 

	                                       

					if(length(class_labels_levels)==2){

							
							zvec=foldchangeres
					}else{
	
							zvec=NA
	
							if(featselmethod=="lmreg" && analysismode=="regression"){
                                
                            cnames_matrix<-colnames(data_limma_fdrall_withfeats)
                            cnames_colindex<-grep("Estimate_",cnames_matrix)
                            
                           
                            
							zvec<-data_limma_fdrall_withfeats[,c(cnames_colindex[1])]
							}


					}


					maxfoldchange<-foldchangeres
                    
                    goodipfoldchange<-which(abs(maxfoldchange)>foldchangethresh)
                    
                    #if(FALSE)
                    {
                        if(input.intensity.scale=="raw" && log2transform==FALSE && znormtransform==FALSE){
                            
                            foldchangeres<-2^((foldchangeres))
                            
                            
                            
                            
                            
                        }
                    }
                    
                    maxfoldchange<-foldchangeres

					roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10))
                    {
                            if(length(x) != 1) stop("'x' must be of length 1")
                            10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
                    }
                    
					d4<-as.data.frame(data_limma_fdrall_withfeats)
                    
                    max_mz_val<-roundUpNice(max(d4$mz)[1])
                    max_time_val<-roundUpNice(max(d4$time)[1])
                    
                    x1increment=round_any(max_mz_val/10,10,f=floor)
                    
                    
                    x2increment=round_any(max_time_val/10,10,f=floor)

					 if(featselmethod=="lmreg" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat"
					  			| featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="logitreg" | featselmethod=="limma2wayrepeat" | featselmethod=="wilcox" | featselmethod=="ttest" |  featselmethod=="poissonreg")
					  {



                        # print("Plotting manhattan plots")
					  
                      sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh


					  			 goodip<-which(sel.diffdrthresh==TRUE)
                       

							classlabels<-as.data.frame(classlabels)
						
						
							   
							   logp<-(-1)*log((d4[,1]+(10^-20)),10)
                               
                              
                                if(fdrmethod=="none"){
							   ythresh<-(-1)*log10(pvalue.thresh)
                               }else{
                                   
                                      ythresh<-min(logp[goodip],na.rm=TRUE)
                                }
			maintext1="Type 1 manhattan plot (-logp vs mz) \n m/z features above the dashed horizontal line meet the selection criteria"
			maintext2="Type 2 manhattan plot (-logp vs time) \n m/z features above the dashed horizontal line meet the selection criteria"



			if(is.na(zvec[1])==FALSE){
				 maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
							maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
			}

    yvec_val=logp
    ylabel="(-)log10p"
    yincrement=1
    y2thresh=(-1)*log10(0.05)
    

    ##save(list=c("d4","yvec_val","ythresh","zvec","x1increment","yincrement","maintext1","x2increment","maintext2","ylabel","y2thresh"),file="manhattanplot_objects.Rda")
    
    
    if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
        
        png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
    }
    
    # get_manhattanplots(xvec=d4$mz,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab=ylabel,xincrement=x1increment,yincrement=yincrement,maintext=maintext1,col_seq=c("black"),y2thresh=y2thresh,colorvec=manhattanplot.col.opt)
    
    ##save(list=ls(),file="m1.Rda")
       try(get_manhattanplots(xvec=d4$mz,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab=ylabel,xincrement=x1increment,yincrement=yincrement,maintext=maintext1,col_seq=c("black"),y2thresh=y2thresh,colorvec=manhattanplot.col.opt),silent=TRUE)


    if(output.device.type!="pdf"){
    
        try(dev.off(),silent=TRUE)
    }

if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
    
    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
}

try(get_manhattanplots(xvec=d4$time,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="Retention time",ylab="-log10p",xincrement=x2increment,yincrement=1,maintext=maintext2,col_seq=c("black"),y2thresh=y2thresh,colorvec=manhattanplot.col.opt),silent=TRUE)

if(output.device.type!="pdf"){
    
    try(dev.off(),silent=TRUE)
}

if(length(class_labels_levels)==2){
if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/VolcanoPlot.png"
    
    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
}

maintext1="Volcano plot (-logp vs log2(fold change)) \n colored m/z features meet the selection criteria"
if(is.na(zvec[1])==FALSE){
    maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
    maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
}


        try(get_volcanoplots(xvec=maxfoldchange,yvec=logp,up_or_down=zvec,ythresh=ythresh,y2thresh=y2thresh,xthresh=foldchangethresh,maintext=maintext1,ylab="-log10(p-value)",xlab="log2(fold change)"),silent=TRUE)

if(output.device.type!="pdf"){
    
    try(dev.off(),silent=TRUE)
}
}

				}else{
							  
							if(featselmethod=="pls" | featselmethod=="o1pls"){

                                
                                # print("Time 2")
                                #print(Sys.time())
                                
				maintext1="Type 1 manhattan plot (VIP vs mz) \n m/z features above the dashed horizontal line meet the selection criteria"
                maintext2="Type 2 manhattan plot (VIP vs time) \n m/z features above the dashed horizontal line meet the selection criteria"

                if(is.na(zvec[1])==FALSE){
                     maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                                maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                }


            yvec_val<-data_limma_fdrall_withfeats[,1]
            ythresh=pls_vip_thresh
            vip_res<-as.data.frame(vip_res)
            bad.feature.index={}
            if(is.na(pls.permut.count)==FALSE){
                #yvec_val[which(vip_res$rand_pls_sel_prob>=pvalue.thresh | vip_res$rand_pls_sel_fdr>=fdrthresh)]<-0 #(ythresh)*0.5
                bad.feature.index=which(vip_res$rand_pls_sel_prob>=pvalue.thresh | vip_res$rand_pls_sel_fdr>=fdrthresh)
            }
            ylabel="VIP"
            yincrement=0.5
            y2thresh=NA
            
            # #save(list=ls(),file="manhattandebug.Rda")
            
                if(output.device.type!="pdf"){
                    
                    temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                
                }
                
                try(get_manhattanplots(xvec=d4$mz,yvec=yvec_val,ythresh=pls_vip_thresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="VIP",xincrement=x1increment,yincrement=0.5,maintext=maintext1,col_seq=c("black"),colorvec=manhattanplot.col.opt,bad.feature.index=bad.feature.index),silent=TRUE)

                if(output.device.type!="pdf"){
                    
                    try(dev.off(),silent=TRUE)
                }
                
                            if(output.device.type!="pdf"){
                                
                                temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                                 png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                            }
                    
                    try(get_manhattanplots(xvec=d4$time,yvec=yvec_val,ythresh=pls_vip_thresh,up_or_down=zvec,xlab="Retention time",ylab="VIP",xincrement=x2increment,yincrement=0.5,maintext=maintext2,col_seq=c("black"),colorvec=manhattanplot.col.opt,bad.feature.index=bad.feature.index),silent=TRUE)
							
                            if(output.device.type!="pdf"){
                                
                                try(dev.off(),silent=TRUE)
                            }

                            if(length(class_labels_levels)==2){
                                
                                if(output.device.type!="pdf"){
                                    temp_filename_1<-"Figures/VolcanoPlot_VIP_vs_foldchange.png"
                                    
                                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                }
                                
                                maintext1="Volcano plot (VIP vs log2(fold change)) \n colored m/z features meet the selection criteria"
                                
                                maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                                
                                maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                           
                           #  #save(list=ls(),file="volcanodebug.Rda")
                           try(get_volcanoplots(xvec=maxfoldchange,yvec=yvec_val,up_or_down=maxfoldchange,ythresh=ythresh,xthresh=foldchangethresh,maintext=maintext1,ylab="VIP",xlab="log2(fold change)",bad.feature.index=bad.feature.index),silent=TRUE)
                            
                                if(output.device.type!="pdf"){
                                    
                                    try(dev.off(),silent=TRUE)
                                }
                            }
                            
                            
                    }else{
										if(featselmethod=="spls" | featselmethod=="o1spls"){
											

maintext1="Type 1 manhattan plot (|loading| vs mz) \n m/z features with non-zero loadings meet the selection criteria"
			maintext2="Type 2 manhattan plot (|loading| vs time) \n m/z features with non-zero loadings meet the selection criteria"
	if(is.na(zvec[1])==FALSE){
				 maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
							maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
			}					
    yvec_val<-data_limma_fdrall_withfeats[,1]
    vip_res<-as.data.frame(vip_res)
    
    bad.feature.index={}
    if(is.na(pls.permut.count)==FALSE){
        # yvec_val[which(vip_res$rand_pls_sel_prob>=pvalue.thresh | vip_res$rand_pls_sel_fdr>=fdrthresh)]<-0
        bad.feature.index=which(vip_res$rand_pls_sel_prob>=pvalue.thresh | vip_res$rand_pls_sel_fdr>=fdrthresh)
    }
    ythresh=0
ylabel="Loading (absolute)"
yincrement=0.1
y2thresh=NA
##save(list=c("d4","yvec_val","ythresh","zvec","x1increment","yincrement","maintext1","x2increment","maintext2","ylabel","y2thresh"),file="manhattanplot_objects.Rda")

if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
    
    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
}

						try(get_manhattanplots(xvec=d4$mz,yvec=yvec_val,ythresh=0,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="Loading (absolute)",xincrement=x1increment,yincrement=0.1,maintext=maintext1,col_seq=c("black"),colorvec=manhattanplot.col.opt,bad.feature.index=bad.feature.index),silent=TRUE)
							
                            if(output.device.type!="pdf"){
                                
                                try(dev.off(),silent=TRUE)
                            }
                            
                            if(output.device.type!="pdf"){
                                
                                temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                                
                                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                            }
                            
                            
                            try(get_manhattanplots(xvec=d4$time,yvec=yvec_val,ythresh=0,up_or_down=zvec,xlab="Retention time",ylab="Loading (absolute)",xincrement=x2increment,yincrement=0.1,maintext=maintext2,col_seq=c("black"),colorvec=manhattanplot.col.opt,bad.feature.index=bad.feature.index),silent=TRUE)
										
                                        
                                        if(output.device.type!="pdf"){
                                            
                                            try(dev.off(),silent=TRUE)
                                        }
                                        
                                        #volcanoplot
                                        if(length(class_labels_levels)==2){
                                            if(output.device.type!="pdf"){
                                                
                                                temp_filename_1<-"Figures/VolcanoPlot_Loading_vs_foldchange.png"
                                                
                                                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                            }
                                            
                                            maintext1="Volcano plot (absolute) Loading vs log2(fold change)) \n colored m/z features meet the selection criteria"
                                            
                                            maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                                            maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                                            
                                            
                                            
                                            try(get_volcanoplots(xvec=maxfoldchange,yvec=yvec_val,up_or_down=maxfoldchange,ythresh=ythresh,xthresh=foldchangethresh,maintext=maintext1,ylab="(absolute) Loading",xlab="log2(fold change)",yincrement=0.1,bad.feature.index=bad.feature.index),silent=TRUE)
                                            
                                            if(output.device.type!="pdf"){
                                                
                                                try(dev.off(),silent=TRUE)
                                            }
                                        }
                                        
                                        
                                        }else{
                                            
                                            if(featselmethod=="pamr"){
                                                
                                                
                                                
                                                
                                                maintext1="Type 1 manhattan plot (max |standardized centroids (d-statistic)| vs mz) \n m/z features with above the horizontal line meet the selection criteria"
                                                maintext2="Type 2 manhattan plot (max |standardized centroids (d-statistic)| vs time) \n m/z features with above the horizontal line meet the selection criteria"
                                                if(is.na(zvec[1])==FALSE){
                                                    maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                                                    maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                                                }
                                                
                                                
                                                yvec_val<-data_limma_fdrall_withfeats[,1]
                                                
                                                ##error point
                                                #vip_res<-as.data.frame(vip_res)
                                                
                                                discore<-as.data.frame(discore)
                                                
                                                bad.feature.index={}
                                                if(is.na(pls.permut.count)==FALSE){
                                                    # yvec_val[which(vip_res$rand_pls_sel_prob>=pvalue.thresh | vip_res$rand_pls_sel_fdr>=fdrthresh)]<-0
                                                    #   bad.feature.index=which(vip_res$rand_pls_sel_prob>=pvalue.thresh | vip_res$rand_pls_sel_fdr>=fdrthresh)
                                                }
                                                ythresh=pamr_ythresh
                                                ylabel="d-statistic (absolute)"
                                                yincrement=0.1
                                                y2thresh=NA
                                                ##save(list=c("d4","yvec_val","ythresh","zvec","x1increment","yincrement","maintext1","x2increment","maintext2","ylabel","y2thresh"),file="manhattanplot_objects.Rda")
                                                
                                                if(output.device.type!="pdf"){
                                                    
                                                    temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
                                                    
                                                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                                }
                                                
                                                try(get_manhattanplots(xvec=d4$mz,yvec=yvec_val,ythresh=pamr_ythresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="d-statistic (absolute) at threshold=0",xincrement=x1increment,yincrement=0.1,maintext=maintext1,col_seq=c("black"),colorvec=manhattanplot.col.opt,bad.feature.index=NA),silent=TRUE)
                                                
                                                if(output.device.type!="pdf"){
                                                    
                                                    try(dev.off(),silent=TRUE)
                                                }
                                                
                                                if(output.device.type!="pdf"){
                                                    
                                                    temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                                                    
                                                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                                }
                                                
                                                
                                                try(get_manhattanplots(xvec=d4$time,yvec=yvec_val,ythresh=pamr_ythresh,up_or_down=zvec,xlab="Retention time",ylab="d-statistic (absolute) at threshold=0",xincrement=x2increment,yincrement=0.1,maintext=maintext2,col_seq=c("black"),colorvec=manhattanplot.col.opt,bad.feature.index=NA),silent=TRUE)
                                                
                                                
                                                if(output.device.type!="pdf"){
                                                    
                                                    try(dev.off(),silent=TRUE)
                                                }
                                                
                                                #volcanoplot
                                                if(length(class_labels_levels)==2){
                                                    if(output.device.type!="pdf"){
                                                        
                                                        temp_filename_1<-"Figures/VolcanoPlot_Dstatistic_vs_foldchange.png"
                                                        
                                                        png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                                    }
                                                    
                                                    maintext1="Volcano plot (absolute) max standardized centroid (d-statistic) vs log2(fold change)) \n colored m/z features meet the selection criteria"
                                                    
                                                    maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                                                    maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
                                                    
                                                    
                                                    
                                                    try(get_volcanoplots(xvec=maxfoldchange,yvec=yvec_val,up_or_down=maxfoldchange,ythresh=pamr_ythresh,xthresh=foldchangethresh,maintext=maintext1,ylab="(absolute) d-statistic at threshold=0",xlab="log2(fold change)",yincrement=0.1,bad.feature.index=NA),silent=TRUE)
                                                    
                                                    if(output.device.type!="pdf"){
                                                        
                                                        try(dev.off(),silent=TRUE)
                                                    }
                                                }
                                                
                                                
                                            }
                                            
                                        }
							}
								
					}
				
                
                

                                        goodip<-intersect(goodip,goodipfoldchange)


                                        dataA<-cbind(maxfoldchange,data_m_fc_withfeats)
                                        #write.table(dataA,file="foldchange.txt",sep="\t",row.names=FALSE)
						
                        goodfeats_allfields<-{}
                        
                        
                        if(length(goodip)>0){
					feat_sigfdrthresh[lf]<-length(goodip)
					subdata<-t(data_m_fc[goodip,])
					
					data_minval<-min(parent_data_m[,-c(1:2)],na.rm=TRUE)*0.5
                    
                     #svm_model<-svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95)
                    svm_model<-try(svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
					
					classlabels<-as.data.frame(classlabels)
                    



					if(is(svm_model,"try-error")){
						print("SVM could not be performed. Please try lowering the kfold or set kfold=total number of samples for Leave-one-out CV. Skipping to the next step.")
						termA<-(-1)
						pred_acc<-termA
                        permut_acc<-(-1)
					}else{
					
						
						pred_acc<-svm_model$avg_acc
					
                            print("Accuracy is:")
                            print(pred_acc)
                            
                            print("Calculating permuted CV accuracy")
                            
                            permut_acc<-{}
                            #permut_acc<-lapply(1:100,function(j){
                            numcores<-num_nodes #round(detectCores()*0.5)
                            cl <- parallel::makeCluster(getOption("cl.cores", num_nodes))
                            clusterEvalQ(cl,library(e1071))
                            clusterEvalQ(cl,library(pROC))
                            clusterEvalQ(cl,library(ROCR))
                            clusterExport(cl,"svm_cv",envir = .GlobalEnv)
                            permut_acc<-parLapply(cl,1:100,function(p1){
                                
                                rand_order<-sample(1:dim(classlabels)[1],size=dim(classlabels)[1])
                                classlabels_permut<-classlabels[rand_order,]
                                classlabels_permut<-as.data.frame(classlabels_permut)
                                svm_permut_res<-try(svm_cv(v=kfold,x=subdata,y=classlabels_permut,kname=svm_kernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
                                
                            
                                if(is(svm_permut_res,"try-error")){
                                    
                                    cur_perm_acc<-NA
                                }else{
                                    cur_perm_acc<-svm_permut_res$avg_acc
                                }
                                return(cur_perm_acc)
                            })
                            
                            stopCluster(cl)
                            
                            permut_acc<-unlist(permut_acc)
                            permut_acc<-mean(permut_acc,na.rm=TRUE)
                            permut_acc<-round(permut_acc,2)
                            
                            print("mean Permuted accuracy is:")
                            print(permut_acc)

					}
					



					termA<-100*pred_acc
					exp_fp<-1
					
                    best_feats<-goodip
                    
					
					
					if(featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma2wayrepeat" | featselmethod=="lmreg" | featselmethod=="logitreg"
| featselmethod=="lm2wayanova" | featselmethod=="lm1wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="wilcox" | featselmethod=="ttest" |  featselmethod=="poissonreg")
					{
                        if(fdrmethod=="none"){
                            exp_fp<-(dim(data_m_fc)[1]*fdrthresh)+1
                        }else{
						exp_fp<-(feat_sigfdrthresh[lf]*fdrthresh)+1
                        }
                    }
					
                    
                    
					termB<-(dim(parent_data_m)[1]*dim(parent_data_m)[1])/(dim(data_m_fc)[1]*dim(data_m_fc)[1]*100)
					
					
					res_score<-(100*(termA-permut_acc))-(feat_weight*termB*exp_fp)
					res_score<-round(res_score,2)
                
                    print("Debug variables")
                    print(termA)
                    print(permut_acc)
                    print(feat_weight)
                    print(exp_fp)
                    
                    print(termB)
                    print(dim(parent_data_m))
                    print(dim(data_m_fc))
                    
                    print("Res score")
                    print(res_score)
                    
                    print("best_cv_res")
                    print(best_cv_res)
                    print("End debug")
                    
						if(lf==0)
						{
						best_logfc_ind<-lf
						
						best_feats<-goodip
						best_cv_res<-res_score
						best_acc<-pred_acc
						best_limma_res<-data_limma_fdrall_withfeats[goodip,] #[sel.diffdrthresh==TRUE,]
						
                       
                        
						}else{
						
					if(res_score>best_cv_res){
				
						best_logfc_ind<-lf
						
						best_feats<-goodip
						best_cv_res<-res_score
						best_acc<-pred_acc
						best_limma_res<-data_limma_fdrall_withfeats[goodip,] #[sel.diffdrthresh==TRUE,]
						}
						
					}
					
					res_score_vec[lf]<-res_score
					
                        }else{
                            
                            print("No features meet the fold change criteria.")
                            
                            if(input.intensity.scale=="raw" && log2transform==FALSE){
                                
                                
                                
                                max.fold.change.log2<-maxfoldchange
                                data_limma_fdrall_withfeats_2<-cbind(max.fold.change.log2,data_limma_fdrall_withfeats)
                                
                                
                            }else{
                                
                                if(input.intensity.scale=="log2" || log2transform==TRUE){
                                    
                                    max.fold.change.log2<-maxfoldchange
                                    data_limma_fdrall_withfeats_2<-cbind(max.fold.change.log2,data_limma_fdrall_withfeats)
                                }
                                
                                
                            }
                            
                            if(logistic_reg==TRUE){
                                
                                fname4<-paste("logitreg","results_allfeatures.txt",sep="")
                            }else{
                                
                                
                                if(poisson_reg==TRUE){
                                    fname4<-paste("poissonreg","results_allfeatures.txt",sep="")
                                    
                                }else{
                                    fname4<-paste(featselmethod,"results_allfeatures.txt",sep="")
                                }
                            }
                            
                            fname4<-paste("Tables/",fname4,sep="")
                            allmetabs_res<-data_limma_fdrall_withfeats_2
                            write.table(data_limma_fdrall_withfeats_2,file=fname4,sep="\t",row.names=FALSE)
                            
                            
                            
                        }
					
					}else{
					
					if(dim(data_m_fc)[2]<kfold){
						print("Number of samples is too small to calculate cross-validation accuracy.")
					}

					}
					#feat_sigfdrthresh_cv<-c(feat_sigfdrthresh_cv,pred_acc)
					
                   
					
					print("########################################")
					print(paste("Relative standard deviation (RSD) threshold: ", log2.fold.change.thresh," %",sep=""))
					#print(paste("FDR threshold: ", fdrthresh,sep=""))
					print(paste("Number of metabolites left after RSD filtering: ", dim(data_m_fc)[1],sep=""))
					print(paste("Number of sig metabolites: ", length(goodip),sep=""))
					
					if(length(goodip)<1){
						try(dev.off(),silent=TRUE)
						next;
					}
					if(pred.eval.method=="CV"){
					feat_sigfdrthresh_cv[lf]<-pred_acc
                    feat_sigfdrthresh_permut[lf]<-permut_acc
					print(paste(kfold,"-fold CV accuracy: ", pred_acc,sep=""))
                    print(paste("Permuted ",kfold,"-fold CV accuracy: ", permut_acc,sep=""))
                    
					}else{
					if(pred.eval.method=="AUC"){
					feat_sigfdrthresh_cv[lf]<-pred_acc
                    feat_sigfdrthresh_permut[lf]<-permut_acc
					print(paste("ROC area under the curve (AUC) is : ", pred_acc,sep=""))
                    print(paste("Permuted ROC area under the curve (AUC) is : ", permut_acc,sep=""))
					}else{
						if(pred.eval.method=="BER"){
					feat_sigfdrthresh_cv[lf]<-pred_acc
                    feat_sigfdrthresh_permut[lf]<-permut_acc
					print(paste("Balanced accuracy rate is : ", pred_acc,sep=""))
                    print(paste("Permuted balanced accuracy rate is : ", permut_acc,sep=""))
					}
					}
					}
					print("######################################")

						
					t1<-table(classlabels)
		
					
					#patientcolors <- unlist(lapply(sampleclass, color.map))
					if(length(goodip)>2){
					
					goodfeats<-as.data.frame(data_m_fc_withfeats[goodip,]) #[sel.diffdrthresh==TRUE,])

					goodfeats<-unique(goodfeats)
					
					
					rnames_goodfeats<-as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
					
					if(length(which(duplicated(rnames_goodfeats)==TRUE))>0){
						goodfeats<-goodfeats[-which(duplicated(rnames_goodfeats)==TRUE),]
					rnames_goodfeats<-rnames_goodfeats[-which(duplicated(rnames_goodfeats)==TRUE)]
					}


					rownames(goodfeats)<-as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
					

					
						data_m<-as.matrix(goodfeats[,-c(1:2)])

						
						rownames(data_m)<-as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
		
					data_m<-unique(data_m)			

				    X<-t(data_m)
				    
				    #X<-replace(as.matrix(X),which(is.na(X)==TRUE),0)
				    

			

						#d1<-try(readLines(search_link),silent=TRUE)


if(heatmap.col.opt=="RdBu"){
    
    heatmap.col.opt="redblue"
}

heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
heatmap_cols<-rev(heatmap_cols)

if(heatmap.col.opt=="topo"){
    heatmap_cols<-topo.colors(256)
    heatmap_cols<-rev(heatmap_cols)
}else{
    if(heatmap.col.opt=="heat"){
        heatmap_cols<-heat.colors(256)
        heatmap_cols<-rev(heatmap_cols)
    }else{
        
        if(heatmap.col.opt=="yellowblue"){
            
            heatmap_cols<-colorRampPalette(c("yellow","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
            #heatmap_cols<-blue2yellow(256) #colorRampPalette(c("yellow","blue"))(256)
            heatmap_cols<-rev(heatmap_cols)
        }else{
            
            if(heatmap.col.opt=="redblue"){
                
                heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
                heatmap_cols<-rev(heatmap_cols)
            }else{
                
                #my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                if(heatmap.col.opt=="redyellowgreen"){
                    
                    heatmap_cols <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                    heatmap_cols<-rev(heatmap_cols)
                }else{
                    if(heatmap.col.opt=="yellowwhiteblue"){
                        
                        heatmap_cols<-colorRampPalette(c("yellow2","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                        heatmap_cols<-rev(heatmap_cols)
                    }else{
                        
                        if(heatmap.col.opt=="redwhiteblue"){
                            
                            heatmap_cols<-colorRampPalette(c("red","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                            heatmap_cols<-rev(heatmap_cols)
                        }else{
                            
                            
                            
                            heatmap_cols <- colorRampPalette(brewer.pal(10, heatmap.col.opt))(256)
                            heatmap_cols<-rev(heatmap_cols)
                            
                        }
                        
                    }
                    
                }
                
            }
            
        }
    }
    
}



						hr <- try(hclust(as.dist(1-WGCNA::cor(t(data_m),method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #metabolites
						hc <- try(hclust(as.dist(1-WGCNA::cor(data_m,method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #samples
						
						if(is(hr,"try-error") || is(hc,"try-error")){
								
                                #	print(hr)
                                #print(hc)
							print("Hierarchical clustering can not be performed. ")
						}
                        else{
						heatmap_file<-paste("heatmap_",featselmethod,".tiff",sep="")
						
						#dev.off()
						
                        #print("plotting HCA")
						
						#tiff(heatmap_file,width=plots.width,height=plots.height,res=plots.res)
						
                        heatmap_mainlabel="" #2-way HCA using all significant features"
                        
                        
                        if(output.device.type!="pdf"){
                            
                            temp_filename_1<-"Figures/HCA_All_selectedfeats.png"
                            
                            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                        }
                        
                      

						
                        ##save(list=c("data_m","hr","hc","heatmap_cols","heatmap_mainlabel","patientcolors"),file="hca_objects.Rda")
						if(znormtransform==FALSE){
							
                            if(hca_type=="two-way"){
                                
                                #   plot(c(1,1),plot=FALSE)
                                #l <- legend(0, 0, bty='n', l1,plot=FALSE, pch = pch_per_group, pt.cex = 0.6)
                                # calculate right margin width in ndc
                                #w <- grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
                                 w <- 0.1
                                par(omd=c(0, 1-w, 0, 1))
                             
                             
						h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main=heatmap_mainlabel, ColSideColors=patientcolors,labRow = FALSE, labCol = FALSE)
                        
                        print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                        
                        
                            }else{
                                #   plot(c(1,1),plot=FALSE)
                                #l <- legend(0, 0, bty='n', l1,plot=FALSE, pch = pch_per_group, pt.cex = 0.6)
                                # calculate right margin width in ndc
                                #w <- grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
                               w<-0.1
                               par(omd=c(0, 1-w, 0, 1))
                                
                                h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main=heatmap_mainlabel, ColSideColors=patientcolors,labRow = FALSE, labCol = FALSE)
                                
                                 print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                            }
						}else{
                            if(hca_type=="two-way"){
                                # plot(c(1,1),plot=FALSE)
                                #l <- legend(0, 0, bty='n', l1,plot=FALSE, pch = pch_per_group, pt.cex = 0.6)
                                # calculate right margin width in ndc
                                # w <- grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
                                w <- 0.1
                              
                                par(omd=c(0, 1-w, 0, 1))
							h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main=heatmap_mainlabel, ColSideColors=patientcolors,labRow = FALSE, labCol = FALSE)
                            
                             print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                            }else{
                                # plot(c(1,1),plot=FALSE)
                                #l <- legend(0, 0, bty='n', l1,plot=FALSE, pch = pch_per_group, pt.cex = 0.6)
                                # calculate right margin width in ndc
                                #w <- grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
                                 w <- 0.1
                                par(omd=c(0, 1-w, 0, 1))
                                h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main=heatmap_mainlabel, ColSideColors=patientcolors,labRow = FALSE, labCol = FALSE)
                                
                                 print(legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,class_labels_levels, col = col_vec,pch = rep(19,length(col_vec)), pt.cex = 0.6, title = "Class",cex=0.8))
                            }
                        }
						#dev.off()
						
						
                        
                        if(output.device.type!="pdf"){
                            
                            try(dev.off(),silent=TRUE)
                        }

						
						
						mycl_samples <- cutree(hc, h=max(hc$height)/2)
						mycl_metabs <- cutree(hr, h=max(hr$height)/2)
						
						ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],goodfeats[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])

						cnames1<-colnames(ord_data)
						cnames1[1]<-"mz_cluster_label"
						colnames(ord_data)<-cnames1
						fname1<-paste("Tables/Clustering_based_sorted_intensity_data.txt",sep="")
						write.table(ord_data,file=fname1,sep="\t",row.names=FALSE)

						fname2<-paste("Tables/Sample_clusterlabels.txt",sep="")
						
						sample_clust_num<-mycl_samples[h73$colInd]
						classlabels_temp<-as.data.frame(classlabels)
						temp1<-classlabels_temp[h73$colInd,1]
						
				
                #		print("class labels")
                #		print(classlabels)
						temp2<-classlabels
						
						#print(head(temp2))
						#temp1<-as.data.frame(temp1)
						
						#print(dim(temp1))
						temp2<-as.data.frame(temp2)
						
						match_ind<-match(temp1,temp2[,1])
						
						temp3<-temp2[match_ind,]
						
						#print(head(temp3))
						temp4<-cbind(temp1,temp3,sample_clust_num)
						
						#print(head(temp1))
						
						
						rnames1<-rownames(temp4)
						temp4<-cbind(rnames1,temp4)
						temp4<-as.data.frame(temp4)
						#temp4<-temp4[,-c(1)]
						
						
                        if(output.device.type!="pdf"){
                            
                            temp_filename_1<-"Figures/barplot_dependent_variable_ordered_by_HCA.png"
                            
                            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                        }
                        
                        

						if(analysismode=="regression"){
						
							
							names(temp3)<-as.character(temp4[,1])
							
							#tiff("Barplot_sample_cluster_ymat.tiff", width=plots.width,height=plots.height,res=plots.res, compression="lzw")
							barplot(temp3,col="brown",ylab="Y",cex.axis=0.5,cex.names=0.5,main="Dependent variable levels in samples; \n ordered based on hierarchical clustering")
							#dev.off()
							

							
						}
						
                        
                        if(output.device.type!="pdf"){
                            
                           try(dev.off(),silent=TRUE)
                        }

						#temp4<-temp4[,-c(2)]
						#colnames(temp4)<-c("SampleID","Class","HCA Cluster #")
						write.table(temp4,file=fname2,sep="\t",row.names=FALSE)
					

					
						fname3<-paste("Tables/Metabolite_clusterlabels.txt",sep="")
						
						mycl_metabs_ord<-mycl_metabs[rev(h73$rowInd)]
						mz_rt_info<-rownames(mycl_metabs_ord)
						mycl_metabs_ord<-cbind(mz_rt_info,mycl_metabs_ord)
						write.table(mycl_metabs_ord,file=fname3,sep="\t",row.names=TRUE)
					
						}
						
						
						#print("no  problem here")
					
					}
                    else{
						if(length(goodip)>1){
						print("Number of FDR significant features is too small to perform PCA and hierarchical clustering using correlation as similarity measure. Using distance measure for hierarchical clustering.")
						
						goodfeats<-as.data.frame(data_m_fc_withfeats[sel.diffdrthresh==TRUE,])


					#rownames(goodfeats)<-as.character(goodfeats[,1])
					rownames(goodfeats)<-as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
					
					#assign colors by sample class
					
					
						data_m<-as.matrix(goodfeats[,-c(1:2)])
						rownames(data_m)<-as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
						
                        
                        if(heatmap.col.opt=="RdBu"){
                            
                            heatmap.col.opt="redblue"
                        }
                        
                        heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
                        heatmap_cols<-rev(heatmap_cols)
                        
                        if(heatmap.col.opt=="topo"){
                            heatmap_cols<-topo.colors(256)
                            heatmap_cols<-rev(heatmap_cols)
                        }else{
                            if(heatmap.col.opt=="heat"){
                                heatmap_cols<-heat.colors(256)
                                heatmap_cols<-rev(heatmap_cols)
                            }else{
                                
                                if(heatmap.col.opt=="yellowblue"){
                                    
                                    heatmap_cols<-colorRampPalette(c("yellow","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                                    #heatmap_cols<-blue2yellow(256) #colorRampPalette(c("yellow","blue"))(256)
                                    heatmap_cols<-rev(heatmap_cols)
                                }else{
                                    
                                    if(heatmap.col.opt=="redblue"){
                                        
                                        heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
                                        heatmap_cols<-rev(heatmap_cols)
                                    }else{
                                        
                                        #my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                                        if(heatmap.col.opt=="redyellowgreen"){
                                            
                                            heatmap_cols <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                                            heatmap_cols<-rev(heatmap_cols)
                                        }else{
                                            if(heatmap.col.opt=="yellowwhiteblue"){
                                                
                                                heatmap_cols<-colorRampPalette(c("yellow2","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                                                heatmap_cols<-rev(heatmap_cols)
                                            }else{
                                                
                                                if(heatmap.col.opt=="redwhiteblue"){
                                                    
                                                    heatmap_cols<-colorRampPalette(c("red","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                                                    heatmap_cols<-rev(heatmap_cols)
                                                }else{
                                                    
                                                    
                                                    
                                                    heatmap_cols <- colorRampPalette(brewer.pal(10, heatmap.col.opt))(256)
                                                    heatmap_cols<-rev(heatmap_cols)
                                                    
                                                }
                                                
                                            }
                                            
                                        }
                                        
                                    }
                                    
                                }
                            }
                            
                        }
                        

						heatmap_file<-paste("heatmap_",featselmethod,".tiff",sep="")
                        hr=TRUE
                        hc=TRUE
                        heatmap_mainlabel="" #"2-way HCA using all significant features"
						
                    
                        
                        if(output.device.type!="pdf"){
                            
                            temp_filename_1<-"Figures/HCA_all_selectedfeats.png"
                            
                            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                        }
                        
                        						if(znormtransform==FALSE){
                                                    
						h73<-heatmap.2(as.matrix(data_m), col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main="", ColSideColors=patientcolors,labRow = FALSE, labCol = FALSE)
						}else{
							h73<-heatmap.2(as.matrix(data_m), col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main="", ColSideColors=patientcolors,labRow = FALSE, labCol = FALSE)
							}
						#dev.off()
					
                    
                    if(output.device.type!="pdf"){
                        
                        try(dev.off(),silent=TRUE)
                    }


						
						
						}
					
					
					}
					
					
					}
					
					else
					{
					
                    
                    #print("regression")
                    
					print("########################################")
					print(paste("RSD threshold: ", log2.fold.change.thresh,sep=""))
					#print(paste("FDR threshold: ", fdrthresh,sep=""))
					print(paste("Number of metabolites left after RSD filtering: ", dim(data_m_fc)[1],sep=""))
					print(paste("Number of sig metabolites: ", length(goodip),sep=""))
					
                    
                    
                    if(featselmethod=="lmreg"){
                        #d4<-read.table(paste(parentoutput_dir,"/Stage2/lmreg_pval_coef_stderr.txt",sep=""),sep="\t",header=TRUE,quote = "")
                    
                    d4<-read.table("lmreg_pval_coef_stderr.txt",sep="\t",header=TRUE)
                    
                    }
                    
                
                    
					if(length(goodip)>=1){
					
						
						subdata<-t(data_m_fc[goodip,])
					
                    
                    if(length(class_labels_levels)==2){
                        
                        
                        #zvec=foldchangeres
                    }else{
                        
                        zvec=NA
                        
                        if(featselmethod=="lmreg" && analysismode=="regression"){
                            
                            cnames_matrix<-colnames(d4)
                            
                           
                            
                            cnames_colindex<-grep("Estimate_",cnames_matrix)
                            
                            
                            zvec<-d4[,c(cnames_colindex[1])]
                            #zvec<-d4$Estimate_var1
                            
                            #if(length(zvec)<1){
                            #   zvec<-d4$X.Estimate_var1.
                                
                                #}
                        }
                        
                        
                    }
                    
                    
                   
                    
                   
                    
                    roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
                        if(length(x) != 1) stop("'x' must be of length 1")
                        10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
                    }
                    
                    d4<-as.data.frame(data_limma_fdrall_withfeats)
                    # d4<-as.data.frame(d1)
                    x1increment=round_any(max(d4$mz)/10,10,f=floor)
                    x2increment=round_any(max(d4$time)/10,10,f=floor)
                    
                    #manplots
                    
                    if(featselmethod=="lmreg" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat"
                    | featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="logitreg" | featselmethod=="limma2wayrepeat" | featselmethod=="wilcox" | featselmethod=="ttest" | featselmethod=="poissonreg")
                    {
                        
                        
                        
                        #print("Plotting manhattan plots")
                        
                        sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
                        
                       
                        
                        goodip<-which(sel.diffdrthresh==TRUE)
                        
                        

                        
                        classlabels<-as.data.frame(classlabels)
                        
                        
                        
                        logp<-(-1)*log((d4[,1]+(10^-20)),10)
                        ythresh<-min(logp[goodip],na.rm=TRUE)
                        
                        maintext1="Type 1 manhattan plot (-logp vs mz) \n m/z features above the dashed horizontal line meet the selection criteria"
                        maintext2="Type 2 manhattan plot (-logp vs time) \n m/z features above the dashed horizontal line meet the selection criteria"
                        
                        
                        # print("here1 A")
                        #print(zvec)
                        
                        if(is.na(zvec[1])==FALSE){
                            maintext1=paste(maintext1,"\nred: negative association "," & green: positive association ",sep="")
                            maintext2=paste(maintext2,"\nred: negative association "," & green: positive association ",sep="")
                        }
                        
                        
                        if(output.device.type!="pdf"){
                            
                            temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
                            
                            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                        }
                        
                        
                       
                        	try(get_manhattanplots(xvec=d4$mz,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="-logP",xincrement=x1increment,yincrement=1,maintext=maintext1,col_seq=c("black"),y2thresh=1.30103,colorvec=manhattanplot.col.opt),silent=TRUE)
                            
                            
                        
                            if(output.device.type!="pdf"){
                                
                                try(dev.off(),silent=TRUE)
                            }
                            
                            
                            if(output.device.type!="pdf"){
                                
                                temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                                
                                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                            }
                            
                            
                        
                            try(get_manhattanplots(xvec=d4$time,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="Retention time",ylab="-logP",xincrement=x2increment,yincrement=1,maintext=maintext2,col_seq=c("black"),y2thresh=1.30103,colorvec=manhattanplot.col.opt),silent=TRUE)
                            
                            
                            if(output.device.type!="pdf"){
                                
                               try(dev.off(),silent=TRUE)
                            }
                        
                        
                        #print("Plotting manhattan plots")
                        #get_manhattanplots(xvec=d4$mz,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="-logP",xincrement=x1increment,yincrement=1,maintext=maintext1)
                        #get_manhattanplots(xvec=d4$time,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="Retention time",ylab="-logP",xincrement=x2increment,yincrement=1,maintext=maintext2)
                        
                    }else{
                        
                        if(featselmethod=="pls" | featselmethod=="o1pls"){
                            
                            maintext1="Type 1 manhattan plot (VIP vs mz) \n m/z features above the dashed horizontal line meet the selection criteria"
                            maintext2="Type 2 manhattan plot (VIP vs time) \n m/z features above the dashed horizontal line meet the selection criteria"
                            
                            if(is.na(zvec[1])==FALSE){
                                maintext1=paste(maintext1,"\nred: negative association "," & green: positive association ",sep="")
                                maintext2=paste(maintext2,"\nred: negative association "," & green: positive association ",sep="")
                            }
                            
                            
                            if(output.device.type!="pdf"){
                                
                                temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
                                
                                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                            }
                            
                            
                           
                            
                            try(get_manhattanplots(xvec=d4$mz,yvec=data_limma_fdrall_withfeats[,1],ythresh=pls_vip_thresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="VIP",xincrement=x1increment,yincrement=0.5,maintext=maintext1,col_seq=c("black"),colorvec=manhattanplot.col.opt),silent=TRUE)
                            
                           
                            
                            if(output.device.type!="pdf"){
                                
                                try(dev.off(),silent=TRUE)
                            }
                            
                            
                            if(output.device.type!="pdf"){
                                
                                temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                                
                                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                            }
                            
                            
                            
                            try(get_manhattanplots(xvec=d4$time,yvec=data_limma_fdrall_withfeats[,1],ythresh=pls_vip_thresh,up_or_down=zvec,xlab="Retention time",ylab="VIP",xincrement=x2increment,yincrement=0.5,maintext=maintext2,col_seq=c("black"),colorvec=manhattanplot.col.opt),silent=TRUE)
                            
                            
                            
                            if(output.device.type!="pdf"){
                                
                               try(dev.off(),silent=TRUE)
                            }
                            
                        }
                        else{
                            if(featselmethod=="spls" | featselmethod=="o1spls"){
                                
                                
                                maintext1="Type 1 manhattan plot (|loading| vs mz) \n m/z features with non-zero loadings meet the selection criteria"
                                maintext2="Type 2 manhattan plot (|loading| vs time) \n m/z features with non-zero loadings meet the selection criteria"
                                if(is.na(zvec[1])==FALSE){
                                    maintext1=paste(maintext1,"\nred: negative association "," & green: positive association ",sep="")
                                    maintext2=paste(maintext2,"\nred: negative association "," & green: positive association ",sep="")
                                }
                                
                                
                                if(output.device.type!="pdf"){
                                    
                                    temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
                                    
                                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                }
                                

                                try(get_manhattanplots(xvec=d4$mz,yvec=data_limma_fdrall_withfeats[,1],ythresh=0,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="Loading",xincrement=x1increment,yincrement=0.1,maintext=maintext1,col_seq=c("black"),colorvec=manhattanplot.col.opt),silent=TRUE)
                                
                                
                                if(output.device.type!="pdf"){
                                    
                                    try(dev.off(),silent=TRUE)
                                }
                                
                                
                                if(output.device.type!="pdf"){
                                    
                                    temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                                    
                                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                }
                                
                                
                                
                                try(get_manhattanplots(xvec=d4$time,yvec=data_limma_fdrall_withfeats[,1],ythresh=0,up_or_down=zvec,xlab="Retention time",ylab="Loading",xincrement=x2increment,yincrement=0.1,maintext=maintext2,col_seq=c("black"),colorvec=manhattanplot.col.opt),silent=TRUE)
                                
                                
                                if(output.device.type!="pdf"){
                                    
                                    try(dev.off(),silent=TRUE)
                                }
                                
                            }
                        }
                        
                    }
                    
                    
					data_minval<-min(parent_data_m[,-c(1:2)],na.rm=TRUE)*0.5
					#subdata<-apply(subdata,2,function(x){naind<-which(is.na(x)==TRUE);if(length(naind)>0){ x[naind]<-median(x,na.rm=TRUE)};return(x)})
					subdata<-apply(subdata,2,function(x){naind<-which(is.na(x)==TRUE);if(length(naind)>0){ x[naind]<-data_minval};return(x)})
					
					
					#print(head(subdata))
					#print(dim(subdata))
					#print(dim(classlabels))
					
					#print(dim(classlabels))
					
					classlabels_response_mat<-as.data.frame(classlabels_response_mat)
					
					if(length(classlabels)>dim(parent_data_m)[2]){
					#classlabels<-as.data.frame(classlabels[,1])
					classlabels_response_mat<-as.data.frame(classlabels_response_mat[,1])
					}
					
					#if(dim(classlabels)[2]<2)
					{
						
					#return(list("x"=subdata,"y"=classlabels_response_mat))
					
					svm_model_reg<-try(svm(x=subdata,y=(classlabels_response_mat[,1]),type="eps",cross=kfold),silent=TRUE)
						
						if(is(svm_model_reg,"try-error")){
							print("SVM could not be performed. Skipping to the next step.")
							termA<-(-1)
							pred_acc<-termA
						}else{
							termA<-svm_model_reg$tot.MSE
							#print("termA is ")
							#print(termA)
							pred_acc<-termA
							print(paste(kfold,"-fold mean squared error: ", pred_acc,sep=""))
							
							# visualize
		     				
						}
					
					}
					
					print("######################################")
					}else{
						print("Number of selected variables is too small to perform CV.")
						
					}
					
					#print("termA is ")
					#print(termA)
                    
                    # print("dim of goodfeats")
                    goodfeats<-as.data.frame(data_m_fc_withfeats[sel.diffdrthresh==TRUE,])



                    
                    goodip<-which(sel.diffdrthresh==TRUE)
                    
                    #print(length(goodip))
                    
                    res_score<-termA
                    
                    #if(res_score<best_cv_res){
                    
                    if(length(which(sel.diffdrthresh==TRUE))>0){
                    if(res_score<best_cv_res){
                        
                        
                        best_logfc_ind<-lf
                        
                        best_feats<-goodip
                        best_cv_res<-res_score
                        best_acc<-pred_acc
                        best_limma_res<-data_limma_fdrall_withfeats[sel.diffdrthresh==TRUE,]
                        
                    }
                    }else{
                     res_score<-(9999999)
                    }
                    res_score_vec[lf]<-res_score
                    
                    
				if(length(which(sel.diffdrthresh==TRUE))>2){
					
					goodfeats<-unique(goodfeats)
					#rownames(goodfeats)<-as.character(goodfeats[,1])
					rownames(goodfeats)<-as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
					
					
						data_m<-as.matrix(goodfeats[,-c(1:2)])
		
						rownames(data_m)<-as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
					
					
				    X<-t(data_m)
               
                    pca_comp<-min(dim(X)[1],dim(X)[2])
                
					t1<-seq(1,dim(data_m)[2])
			
				    col <-col_vec[1:length(t1)]
		

						
						hr <- try(hclust(as.dist(1-cor(t(data_m),method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #metabolites
						hc <- try(hclust(as.dist(1-cor(data_m,method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #samples
						
                        
                        if(heatmap.col.opt=="RdBu"){
                            
                            heatmap.col.opt="redblue"
                        }
                        
                        heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
                        heatmap_cols<-rev(heatmap_cols)
                        
                        if(heatmap.col.opt=="topo"){
                            heatmap_cols<-topo.colors(256)
                            heatmap_cols<-rev(heatmap_cols)
                        }else{
                            if(heatmap.col.opt=="heat"){
                                heatmap_cols<-heat.colors(256)
                                heatmap_cols<-rev(heatmap_cols)
                            }else{
                                
                                if(heatmap.col.opt=="yellowblue"){
                                    
                                    heatmap_cols<-colorRampPalette(c("yellow","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                                    #heatmap_cols<-blue2yellow(256) #colorRampPalette(c("yellow","blue"))(256)
                                    heatmap_cols<-rev(heatmap_cols)
                                }else{
                                    
                                    if(heatmap.col.opt=="redblue"){
                                        
                                        heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
                                        heatmap_cols<-rev(heatmap_cols)
                                    }else{
                                        
                                        #my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                                        if(heatmap.col.opt=="redyellowgreen"){
                                            
                                            heatmap_cols <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                                            heatmap_cols<-rev(heatmap_cols)
                                        }else{
                                            if(heatmap.col.opt=="yellowwhiteblue"){
                                                
                                                heatmap_cols<-colorRampPalette(c("yellow2","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                                                heatmap_cols<-rev(heatmap_cols)
                                            }else{
                                                
                                                if(heatmap.col.opt=="redwhiteblue"){
                                                    
                                                    heatmap_cols<-colorRampPalette(c("red","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                                                    heatmap_cols<-rev(heatmap_cols)
                                                }else{
                                                    
                                                    
                                                    
                                                    heatmap_cols <- colorRampPalette(brewer.pal(10, heatmap.col.opt))(256)
                                                    heatmap_cols<-rev(heatmap_cols)
                                                    
                                                }
                                                
                                            }
                                            
                                        }
                                        
                                    }
                                    
                                }
                            }
                            
                        }
                        


						if(is(hr,"try-error") || is(hc,"try-error")){
								
							print("Hierarchical clustering can not be performed. ")
						}else{


						mycl_samples <- cutree(hc, h=max(hc$height)/2)
                        t1<-table(mycl_samples)
                        col_clust<-topo.colors(length(t1))
                        patientcolors=rep(col_clust,t1) #mycl_samples[col_clust]
						heatmap_file<-paste("heatmap_",featselmethod,"_imp_features.tiff",sep="")
						
						#tiff(heatmap_file,width=plots.width,height=plots.height,res=plots.res, compression="lzw")
                        
                        
                        if(output.device.type!="pdf"){
                            
                            temp_filename_1<-"Figures/HCA_all_selectedfeats.png"
                            
                            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                        }
                        
                        
                        
                        
						if(znormtransform==FALSE){
						h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main="Using all selected features",labRow = FALSE, labCol = FALSE)
						}else{
							h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main="Using all selected features",labRow = FALSE, labCol = FALSE)
							}
                        
                        
                        if(output.device.type!="pdf"){
                            
                            try(dev.off(),silent=TRUE)
                        }
				
						
						
						mycl_samples <- cutree(hc, h=max(hc$height)/2)
						mycl_metabs <- cutree(hr, h=max(hr$height)/2)
						
						ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],goodfeats[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])

cnames1<-colnames(ord_data)
						cnames1[1]<-"mz_cluster_label"
						colnames(ord_data)<-cnames1
						fname1<-paste("Tables/Clustering_based_sorted_intensity_data.txt",sep="")
						write.table(ord_data,file=fname1,sep="\t",row.names=FALSE)

						fname2<-paste("Tables/Sample_clusterlabels.txt",sep="")
						
						sample_clust_num<-mycl_samples[h73$colInd]
						
						
						
						classlabels<-as.data.frame(classlabels)
						
						
						temp1<-classlabels[h73$colInd,]
												
						temp3<-cbind(temp1,sample_clust_num)
						
						rnames1<-rownames(temp3)
						temp4<-cbind(rnames1,temp3)
						temp4<-as.data.frame(temp4)
						

						if(analysismode=="regression"){
							
														
							#names(temp3[,1)<-as.character(temp4[,1])
							
						
							
							temp3<-temp4[,-c(1)]
							temp3<-as.data.frame(temp3)
							temp3<-apply(temp3,2,as.numeric)
							
							
							temp_vec<-as.vector(temp3[,1])
							
							
							
							names(temp_vec)<-as.character(temp4[,1])
							
                            
                            if(output.device.type!="pdf"){
                                
                                temp_filename_1<-"Figures/Barplot_dependent_variable_ordered_by_HCA.png"
                                
                                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                            }
                            
                            
							
							#tiff("Barplot_sample_cluster_ymat.tiff", width=plots.width,height=plots.height,res=plots.res, compression="lzw")
							barplot(temp_vec,col="brown",ylab="Y",cex.axis=0.5,cex.names=0.5,main="Dependent variable levels in samples; \n ordered based on hierarchical clustering")
							#dev.off()
                            
                            
                            if(output.device.type!="pdf"){
                                
                                try(dev.off(),silent=TRUE)
                            }
							

							
						}
						
                        #	print(head(temp_vec))
						#temp4<-temp4[,-c(2)]
						write.table(temp4,file=fname2,sep="\t",row.names=FALSE)
					

					
						fname3<-paste("Metabolite_clusterlabels.txt",sep="")
						
						mycl_metabs_ord<-mycl_metabs[rev(h73$rowInd)]
						
						}
				}

						
						
						
						
					}
					
                    classlabels_orig<-classlabels_orig_parent
                    if(pairedanalysis==TRUE){
                        
                        classlabels_orig<-classlabels_orig[,-c(2)]
                        
                        
                    }else{
                        
                        if(featselmethod=="lmreg" || featselmethod=="logitreg" || featselmethod=="poissonreg"){
                            classlabels_orig<-classlabels_orig[,c(1:2)]
                            classlabels_orig<-as.data.frame(classlabels_orig)
                        }
                    }
  
                    classlabels_orig_wgcna<-classlabels_orig
                    
                    
                    
                    if(analysismode=="classification"){
                        
                        classlabels_temp<-classlabels_orig_wgcna #cbind(classlabels_sub[,1],classlabels)
                        
                        #  degree_eval_res<-try(degree_eval(X=data_m_fc_withfeats,Y=classlabels_temp,sigfeats=data_m_fc_withfeats[goodip,c(1:2)],sigfeatsind=goodip),silent=TRUE)
                        
                        # degree_eval_res<-try(diffrank_eval(X=data_m_fc_withfeats,Y=classlabels_temp,sigfeats=data_m_fc_withfeats[goodip,c(1:2)],sigfeatsind=goodip),silent=TRUE)
                        
                        degree_eval_res<-diffrank_eval(X=data_m_fc_withfeats,Y=classlabels_temp,sigfeats=data_m_fc_withfeats[goodip,c(1:2)],sigfeatsind=goodip,num_nodes=num_nodes,abs.cor.thresh=abs.cor.thresh,cor.fdrthresh=cor.fdrthresh)
        
                    }
                    
                    
                    
                    if(analysismode=="classification")
                    {
                        
                        
                            degree_rank<-rep(1,dim(data_m_fc_withfeats)[1])
                            
                            if(is(degree_eval_res,"try-error")){
                                
                                
                                degree_rank<-rep(1,dim(data_m_fc_withfeats)[1])
                                
                            }else{
                                if(degree_rank_method=="overall"){
                                    degree_rank<-rank((-1)*degree_rank[goodip])
                                }else{
                                    
                                    if(FALSE){
                                        degree_rank_mat<-degree_eval_res$all[,-c(1:4)]
                                        degree_rank_list<-new("list")
                                        
                                        for(i in 1:dim(degree_rank_mat)[2]){
                                            
                                            degree_rank_list[[i]]<-degree_rank_mat[,i]/max(degree_rank_mat[,i])
                                            
                                        }
                                        
                                        #diff_degree_measure<-abs(degree_rank_list[[1]]-degree_rank_list[[2]])
                                        
                                        diff_degree_measure<-apply(degree_rank_mat,1,function(x){res<-{};for(i in 1:length(x)){res<-c(res,(x[i]-x[-i]));};tempres<-abs(res);res_ind<-which(tempres==max(tempres,na.rm=TRUE));return(abs(res[res_ind[1]]));})
                                        
                                    }
                                    
                                    diff_degree_measure<-degree_eval_res$all$DiffRank
                                    
                                    degree_rank<-rank((-1)*diff_degree_measure)
                                }
                                
                                
                }
                            
               
               if(featselmethod=="lmreg" | featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma1way" | featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="limma1wayrepeat" | featselmethod=="limma2wayrepeat" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="wilcox" | featselmethod=="ttest" | featselmethod=="poissonreg")
                        {
                            diffexp_rank<-rank(data_limma_fdrall_withfeats[,2]) #order(data_limma_fdrall_withfeats[,2],decreasing=FALSE)
                        }else{
                            
                            
                            if(featselmethod=="rfesvm"){
                                
                                
                                diffexp_rank<-rank_vec
                                data_limma_fdrall_withfeats<-cbind(rank_vec,data_limma_fdrall_withfeats)
                                
                                
                            }else{
                            
                            
                                    if(featselmethod=="pamr"){
                                        
                                        diffexp_rank<-rank_vec
                                        data_limma_fdrall_withfeats<-cbind(rank_vec,data_limma_fdrall_withfeats)
                                        
                                        
                                    }else{
                                    
                                         if(featselmethod=="MARS"){
                                        
                                                diffexp_rank<-rank((-1)*data_limma_fdrall_withfeats[,2])
                                        
                                         }else{
                                                diffexp_rank<-rank((1)*data_limma_fdrall_withfeats[,2])
                                         }
                                    }
                            }
                        }
                        
			 if(input.intensity.scale=="raw" && log2transform==FALSE){
							
								 
								
								max.fold.change.log2<-maxfoldchange
							     data_limma_fdrall_withfeats_2<-cbind(max.fold.change.log2,degree_rank,diffexp_rank,data_limma_fdrall_withfeats)

						
					}else{

					 if(input.intensity.scale=="log2" || log2transform==TRUE){

						max.fold.change.log2<-maxfoldchange
						data_limma_fdrall_withfeats_2<-cbind(max.fold.change.log2,degree_rank,diffexp_rank,data_limma_fdrall_withfeats)
					}

								
			}
                        
                   
                        
                        if(logistic_reg==TRUE){
                            
                            fname4<-paste("logitreg","results_allfeatures.txt",sep="")
                        }else{
                            
                            
                            if(poisson_reg==TRUE){
                                fname4<-paste("poissonreg","results_allfeatures.txt",sep="")
                            }else{
                                fname4<-paste(featselmethod,"results_allfeatures.txt",sep="")
                            }
                        }
                        
                        fname4<-paste("Tables/",fname4,sep="")
			
                        allmetabs_res<-data_limma_fdrall_withfeats_2
                        
                        write.table(data_limma_fdrall_withfeats_2,file=fname4,sep="\t",row.names=FALSE)
                        
                         
                        if(length(goodip)>=1){
                            
                            data_limma_fdrall_withfeats_2<-data_limma_fdrall_withfeats_2[goodip,]
                            
                            data_limma_fdrall_withfeats_2<-as.data.frame(data_limma_fdrall_withfeats_2)
                            
                            
                       
                        
                        if(featselmethod=="lmreg" | featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma1way" | featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="limma1wayrepeat" | featselmethod=="limma2wayrepeat" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="wilcox" | featselmethod=="ttest" | featselmethod=="poissonreg")
                        {
                            diffexp_rank<-rank(data_limma_fdrall_withfeats_2[,3]) #order(data_limma_fdrall_withfeats[,2],decreasing=FALSE)
                            
                        }else{
                            
                            
                            diffexp_rank<-rank((1)*data_limma_fdrall_withfeats_2[,3])
                        }
                        
                        if(FALSE)
                        {
                        
                            if(is(degree_eval_res,"try-error")){
                                
                                
                                degree_rank<-rep(1,dim(data_limma_fdrall_withfeats_2)[1])
                                
                            }else{
                                if(degree_rank_method=="overall"){
                                    degree_rank<-rank((-1)*degree_rank[goodip])
                                }else{
                                    
                                    degree_rank_mat<-degree_eval_res$sigfeats[,-c(1:4)]
                                    degree_rank_list<-new("list")
                                    
                                    for(i in 1:dim(degree_rank_mat)[2]){
                                        
                                        degree_rank_list[[i]]<-degree_rank_mat[,i]/max(degree_rank_mat[,i])
                                        
                                    }
                                    
                                    #diff_degree_measure<-abs(degree_rank_list[[1]]-degree_rank_list[[2]])
                                    
                                    diff_degree_measure<-apply(degree_rank_mat,1,function(x){res<-{};for(i in 1:length(x)){res<-c(res,(x[i]-x[-i]));};tempres<-abs(res);res_ind<-which(tempres==max(tempres,na.rm=TRUE));return(abs(res[res_ind[1]]));})
                                    
                                    degree_rank<-rank((-1)*diff_degree_measure)
                                }
                                
                            
                            }
                        
                        }
                        
                        #degree_rank<-rep(1,nrow(data_limma_fdrall_withfeats_2))
                        
                        #data_limma_fdrall_withfeats_2$degree_rank<-degree_rank
                        data_limma_fdrall_withfeats_2$diffexp_rank<-diffexp_rank
                        if(logistic_reg==TRUE){
                            
                            fname4<-paste("logitreg","results_selectedfeatures.txt",sep="")
                            
                        }else{
                            
                            
                            if(poisson_reg==TRUE){
                                
                                fname4<-paste("poissonreg","results_selectedfeatures.txt",sep="")
                                
                            }else{
                                fname4<-paste(featselmethod,"results_selectedfeatures.txt",sep="")
                            }
                        }
                        
                        rank_list<-t(data_limma_fdrall_withfeats_2[,c(2:3)])
                        
                        #print(rank_list)
                        
                        if(length(rocfeatlist)>length(goodip)){
                            
                            rocfeatlist<-seq(1,(length(goodip)))
                            numselect<-length(goodip)
                            rocfeatlist<-rocfeatlist+1
                        }else{
                            
                            numselect<-length(rocfeatlist)
                        }
                       
	       
                        #data_limma_fdrall_withfeats_2<-cbind(maxfoldchange,data_limma_fdrall_withfeats)
                        goodfeats<-as.data.frame(data_limma_fdrall_withfeats_2)
                        
                        }
                    }else{
                        
                        
                        
                        
                      
                        
                        if(featselmethod=="lmreg" | featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma1way" | featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="limma1wayrepeat" | featselmethod=="limma2wayrepeat" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="wilcox" | featselmethod=="ttest" | featselmethod=="poissonreg")
                        {
                            diffexp_rank<-rank(data_limma_fdrall_withfeats[,1]) #order(data_limma_fdrall_withfeats[,2],decreasing=FALSE)
                        }else{
                            
                            
                            if(featselmethod=="rfesvm"){
                                
                                
                                diffexp_rank<-rank_vec
                                data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats
                                
                                
                            }else{
                            
                            
                                    if(featselmethod=="pamr"){
                                        
                                        diffexp_rank<-rank_vec
                                        data_limma_fdrall_withfeats<-cbind(rank_vec,data_limma_fdrall_withfeats)
                                        
                                        
                                    }else{
                                    
                                            if(featselmethod=="MARS"){
                                                
                                                diffexp_rank<-rank((-1)*data_limma_fdrall_withfeats[,1])
                                                
                                            }else{
                                                diffexp_rank<-rank((1)*data_limma_fdrall_withfeats[,2])
                                            }
                                    }
                             }
                            
                        }
                        
                        
                            
                            
                            # degree_rank<-rep(1,dim(data_limma_fdrall_withfeats_2)[1])
                            
                        
                        
                        data_limma_fdrall_withfeats_2<-cbind(diffexp_rank,data_limma_fdrall_withfeats)
                        # fname4<-paste(featselmethod,"_sigfeats.txt",sep="")
                        
                        goodfeats<-data_limma_fdrall_withfeats_2[goodip,] #[sel.diffdrthresh==TRUE,]
                        # write.table(goodfeats,file=fname4,sep="\t",row.names=FALSE)
                        ##
                        fname4<-paste("Tables/",featselmethod,"results_allfeatures.txt",sep="")
                    
                        allmetabs_res<-data_limma_fdrall_withfeats_2
                      
                      
                      #  write.table(data_limma_fdrall_withfeats_2[,-which(names(data_limma_fdrall_withfeats_2)%in%c("degree_rank"))],file=fname4,sep="\t",row.names=FALSE)
                    
                        write.table(data_limma_fdrall_withfeats_2,file=fname4,sep="\t",row.names=FALSE)
                    
                    
                        ###
                    }
                    
                   
                    }
                    
                    
                    
                    if(length(goodip)>1){
                    goodfeats_by_DICErank<-{}
                    
                    if(analysismode=="classification"){
                        if(featselmethod=="lmreg" | featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma1way" | featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="limma1wayrepeat" | featselmethod=="limma2wayrepeat" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="wilcox" | featselmethod=="ttest" | featselmethod=="poissonreg")
                        {
                            goodfeats<-goodfeats[order(goodfeats[,3],decreasing=FALSE),]
                            
                            if(length(goodip)>1){
                           # goodfeats_by_DICErank<-data_limma_fdrall_withfeats_2[r1$top.list,]
                            }
                        }else{
                            
                            goodfeats<-goodfeats[order(goodfeats[,3],decreasing=FALSE),]
                            
                             if(length(goodip)>1){
                            #goodfeats_by_DICErank<-data_limma_fdrall_withfeats_2[r1$top.list,]
                            
                             }
                        }
                    }else{
                        if(analysismode=="regression"){
                            
                            try(dev.off(),silent=TRUE)
                            
                            if(featselmethod=="lmreg" | featselmethod=="pls" | featselmethod=="spls" |      featselmethod=="o1pls" | featselmethod=="RF" | featselmethod=="MARS"){
                                
                                ##save(goodfeats,file="goodfeats.Rda")
                                goodfeats<-goodfeats[order(goodfeats[,1],decreasing=FALSE),]
                                
                            }else{
                                
                                
                                #goodfeats<-goodfeats[order(goodfeats[,1],decreasing=TRUE),]
                                
                            }
                        }
                        
                    }
                    
                    if(logistic_reg==TRUE){
                        
                        
                        fname4<-paste("logitreg.selected.features.final.txt",sep="")
                    }else{
                        
                        if(poisson_reg==TRUE){
                         
                        
                         fname4<-paste("poissonreg.selected.features.final.txt",sep="")
                         
                        }else{
                        fname4<-paste(featselmethod,".selected.features.final.txt",sep="")
                        }
                    }
                    
                    
                    
						
                    }
                    
                    if(logistic_reg==TRUE){
                        
                        
                        fname4<-paste("logitreg.selected.features.final.txt",sep="")
                    }else{
                        
                        
                        if(poisson_reg==TRUE){
                            
                            
                            fname4<-paste("poissonreg.selected.features.final.txt",sep="")
                            
                        }else{
                            fname4<-paste(featselmethod,".selected.features.final.txt",sep="")
                        }
                    }
                    #fname4<-paste(featselmethod,"_sigfeats_ordered_by_significance.txt",sep="")
		    
		   # fname4<-paste("Tables/",fname4,sep="")
                    
                    if(length(which(names(goodfeats)%in%c("degree_rank")))>0){
                        
                        #write.table(goodfeats[,-which(names(goodfeats)%in%c("degree_rank"))],file=fname4,sep="\t",row.names=FALSE)
                        
                        write.table(goodfeats,file=fname4,sep="\t",row.names=FALSE)
                        
                    }else{
                    write.table(goodfeats,file=fname4,sep="\t",row.names=FALSE)
                    }




class_label_A<-class_labels_levels[1]
class_label_B<-class_labels_levels[2]

goodfeats_allfields<-{}

if(length(which(sel.diffdrthresh==TRUE))>1){

goodfeats<-as.data.frame(goodfeats)
mzvec<-goodfeats$mz
timevec<-goodfeats$time

if(length(mzvec)>4){
max_per_row<-3


par_rows<-ceiling(9/max_per_row)

}else{
max_per_row<-length(mzvec)
par_rows<-1
}


goodfeats<-as.data.frame(goodfeats)

cnamesd1<-colnames(goodfeats)
time_ind<-which(cnamesd1=="time")

goodfeats_allfields<-as.data.frame(goodfeats)

    
file_ind<-1


mz_ind<-which(cnamesd1=="mz")





#if(length(class_labels_levels)<10)
if(analysismode=="classification" && nrow(goodfeats)>=1 && length(goodip)>=1)
{
    
    goodfeats_temp<-cbind(goodfeats[,mz_ind],goodfeats[,time_ind],goodfeats[,-c(1:time_ind)])
    
    cnames_temp<-colnames(goodfeats_temp)
    cnames_temp[1]<-"mz"
    cnames_temp[2]<-"time"
    colnames(goodfeats_temp)<-cnames_temp
    

if(length(class_labels_levels)==2){

#print("Generating ROC curve using top features on training set")


#try(get_roc(dataA=goodfeats_temp,classlabels=classlabels,classifier=rocclassifier,kname="radial",rocfeatlist=rocfeatlist,rocfeatincrement=rocfeatincrement,mainlabel="Training set ROC curve using top features"),silent=TRUE)

}

#print("ROC done")
best_subset<-{}
best_acc<-0

xvec<-{}
yvec<-{}
#for(i in 2:max_varsel)

if(nrow(goodfeats_temp)>30){

	max_cv_varsel<-30
}else{
	max_cv_varsel<-nrow(goodfeats_temp)
}

cv_yvec<-lapply(2:max_cv_varsel,function(i)
{
    
    subdata<-t(goodfeats_temp[1:i,-c(1:2)])
    svm_model<-try(svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
    
    if(is(svm_model,"try-error")){
        
        res1<-NA
        
    }else{
        
       
       
        res1<-svm_model$avg_acc
        
       
    }
    
    return(res1)
})

xvec<-seq(2:max_cv_varsel)
 
yvec<-unlist(cv_yvec)
 
#if(svm_model$avg_acc>best_acc){
            
  #          best_acc<-svm_model$avg_acc
    #        best_subset<-seq(1,i)
            
            
      #  }


if(pred.eval.method=="CV"){
    ylab_text=paste(pred.eval.method," accuracy (%)",sep="")
    
}else{
    if(pred.eval.method=="BER"){
    ylab_text=paste("Balanced accuracy"," (%)",sep="")
    }else{

	ylab_text=paste("AUC"," (%)",sep="")
   }
}


if(length(yvec)>0){
    
    
    if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/kfoldCV_forward_selection.png"
        
        png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
    }else{
        temp_filename_1<-"Figures/kfoldCV_forward_selection.pdf"
        
        pdf(temp_filename_1)
        
    }
    
    
   
plot(x=xvec,y=yvec,main="k-fold CV classification accuracy based on forward selection of top features",xlab="Feature index",ylab=ylab_text,type="b",col="#0072B2",cex.main=0.8)


if(output.device.type!="pdf"){
    
   try(dev.off(),silent=TRUE)
}else{
    try(dev.off(),silent=TRUE)
    
}

cv_mat<-cbind(xvec,yvec)
colnames(cv_mat)<-c("Feature Index",ylab_text)

write.table(cv_mat,file="Tables/kfold_cv_mat.txt",sep="\t")
}

if(pairedanalysis==TRUE)
{

if(featselmethod=="pls" | featselmethod=="spls"){
    classlabels_sub<-classlabels_sub[,-c(1)]
	classlabels_temp<-cbind(classlabels_sub)
}else{
     classlabels_sub<-classlabels_sub[,-c(1)]
    classlabels_temp<-cbind(classlabels_sub)
    
}

}else{
	classlabels_temp<-cbind(classlabels_sub,classlabels)
}


num_sig_feats<-nrow(goodfeats)

if(num_sig_feats<3){
    
    pca.stage2.eval=FALSE

}
	
    if(pca.stage2.eval==TRUE)
    {
        
        
        pca_comp<-min(10,dim(X)[2])
         
        #dev.off()
        
        # print("plotting")
        #pdf("sig_features_evaluation.pdf", height=2000,width=2000)
        library(pcaMethods)

        p1<-pcaMethods::pca(X,method="rnipals",center=TRUE,scale="uv",cv="q2",nPcs=pca_comp)
        
        if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/PCAdiagnostics_selectedfeats.png"
            
            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        }
        
        
        
        p2<-plot(p1,col=c("darkgrey","grey"),main="PCA diagnostics after variable selection")
        
        print(p2)
        if(output.device.type!="pdf"){
            
            try(dev.off(),silent=TRUE)
        }
        #dev.off()

    }

classlabels_orig<-classlabels_orig_parent


if(pairedanalysis==TRUE){
    
    classlabels_orig<-classlabels_orig[,-c(2)]
    
    
}else{
    
    if(featselmethod=="lmreg" || featselmethod=="logitreg" || featselmethod=="poissonreg"){
        classlabels_orig<-classlabels_orig[,c(1:2)]
        classlabels_orig<-as.data.frame(classlabels_orig)
    }
}



classlabels_orig_wgcna<-classlabels_orig


	goodfeats_temp<-cbind(goodfeats[,mz_ind],goodfeats[,time_ind],goodfeats[,-c(1:time_ind)])
	cnames_temp<-colnames(goodfeats_temp)
    cnames_temp<-c("mz","time",cnames_temp[-c(1:2)])
    colnames(goodfeats_temp)<-cnames_temp
   
    
    
    
    if(num_sig_feats>=3){
        if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/PCAplots_selectedfeats.pdf"
            
            #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
            pdf(temp_filename_1)
        }
        
       
             par(mfrow=c(1,1),family="sans",cex=cex.plots)
             try(get_pcascoredistplots(X=goodfeats_temp,Y=classlabels_orig,feature_table_file=NA,parentoutput_dir=getwd(),class_labels_file=NA,sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,col_vec=col_vec,pairedanalysis=pairedanalysis,pca.cex.val=pca.cex.val,legendlocation=legendlocation,pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,filename="selected",paireddesign=paireddesign,lineplot.col.opt=lineplot.col.opt,lineplot.lty.option=lineplot.lty.option),silent=TRUE)
       
       if(output.device.type!="pdf"){
           
           try(dev.off(),silent=TRUE)
       }
    }
       ##save(list=ls(),file="timeseries.Rda")

       #if(FALSE)
{
    
       if(log2transform==TRUE || input.intensity.scale=="log2"){
           
           if(znormtransform==TRUE){
               ylab_text_2="scale normalized"
           }else{
               if(quantile_norm==TRUE){
                   
                   ylab_text_2="quantile normalized"
               }else{
                   
                   ylab_text_2=""
               }
           }
           ylab_text=paste("log2 intensity ",ylab_text_2,sep="")
       }else{
           if(znormtransform==TRUE){
               ylab_text_2="scale normalized"
           }else{
               if(quantile_norm==TRUE){
                   
                   ylab_text_2="quantile normalized"
               }else{
                   ylab_text_2=""
               }
           }
           ylab_text=paste("Raw intensity ",ylab_text_2,sep="")
       }
       par(mfrow=c(1,1),family="sans",cex=cex.plots)

     if(pairedanalysis==TRUE)
       {
           
           if(output.device.type!="pdf"){
               
               temp_filename_1<-"Figures/Lineplots_selectedfeats.pdf"
               
               #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
               pdf(temp_filename_1)
           }
           get_lineplots(X=goodfeats_temp,Y=classlabels_orig,feature_table_file=NA,parentoutput_dir=getwd(),class_labels_file=NA,lineplot.col.opt=lineplot.col.opt,alphacol=0.3,col_vec=col_vec,pairedanalysis=pairedanalysis,pca.cex.val=pca.cex.val,legendlocation=legendlocation,pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,filename="selected",ylabel=ylab_text,error.bar=error.bar,cex.val=cex.plots,lineplot.lty.option=lineplot.lty.option) #,silent=TRUE)  #,silent=TRUE)
       
           if(output.device.type!="pdf"){
               
               try(dev.off(),silent=TRUE)
           }
          
       }
         
         
         
    
}
    
    
   
	

if(nrow(goodfeats)<1){
    
    print(paste("No features selected for ",featselmethod,sep=""))
}
#else
{

    
    write.table(goodfeats_temp,file="Tables/boxplots_file.normalized.txt",sep="\t",row.names=FALSE)
    
    goodfeats<-goodfeats[,-c(1:time_ind)]


    goodfeats_raw<-data_matrix_beforescaling_rsd[goodip,]
    write.table(goodfeats_raw,file="Tables/boxplots_file.raw.txt",sep="\t",row.names=FALSE)
    
   
     
     if(output.device.type!="pdf"){
         
         try(dev.off(),silent=TRUE)
         
         temp_filename_1<-"Figures/Boxplots.selectedfeats.normalized.pdf"
        
         pdf(temp_filename_1)
     }
       par(mfrow=c(1,1),family="sans",cex=cex.plots)
     get_boxplots(X=goodfeats_temp,Y=classlabels_orig,parentoutput_dir=output_dir,boxplot.col.opt=boxplot.col.opt,alphacol=0.3,newdevice=FALSE,cex=cex.plots,ylabel=ylab_text)
     
     
                
                
                if(output.device.type!="pdf"){
                    
                    try(dev.off(),silent=TRUE)
                }
                
                if(output.device.type!="pdf"){
                    
                    temp_filename_1<-"Figures/Boxplots.selectedfeats.raw.pdf"
                    
                    pdf(temp_filename_1)
                }
                
                
                par(mfrow=c(1,1),family="sans",cex=cex.plots)
                get_boxplots(X=goodfeats_raw,Y=classlabels_orig,parentoutput_dir=output_dir,boxplot.col.opt=boxplot.col.opt,alphacol=0.3,newdevice=FALSE,cex=cex.plots,ylabel="raw Intensity")
                
             
                
                if(output.device.type!="pdf"){
                    
                    try(dev.off(),silent=TRUE)
                }
                
                
                ##save(goodfeats_temp,file="goodfeats_temp.Rda")
                # #save(classlabels_orig,file="classlabels_orig.Rda")
                
                if(output.device.type!="pdf"){
                    
                    temp_filename_1<-"Figures/Barplots_selectedfeats.pdf"
                    
                    #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                    pdf(temp_filename_1,bg="transparent") #, height = 5.5, width = 3)
                }
                
              
                par(mfrow=c(1,1),family="sans",cex=cex.plots,pty="s")
                #try(get_barplots(feature_table_file,class_labels_file,X=goodfeats_temp,Y=classlabels_orig,parentoutput_dir=output_dir,newdevice=FALSE,ylabel=ylab_text,cex.val=cex.plots,barplot.col.opt=barplot.col.opt,error.bar=error.bar),silent=TRUE)
                
               #save(list=ls(),file="getbarplots.Rda")
               
               get_barplots(feature_table_file,class_labels_file,X=goodfeats_temp,Y=classlabels_orig,parentoutput_dir=output_dir,newdevice=FALSE,ylabel=ylab_text,cex.val=cex.plots,barplot.col.opt=barplot.col.opt,error.bar=error.bar,barplot.xaxis=barplot.xaxis)
                
                if(output.device.type!="pdf"){
                    
                    try(dev.off(),silent=TRUE)
                }
                
                if(output.device.type!="pdf"){
                    
                    temp_filename_1<-"Figures/Individual_sample_plots_selectedfeats.pdf"
                    
                    #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                    pdf(temp_filename_1)
                }
                
                #  par(mfrow=c(2,2))
                par(mfrow=c(1,1),family="sans",cex=cex.plots)
                #try(get_individualsampleplots(feature_table_file,class_labels_file,X=goodfeats_temp,Y=classlabels_orig,parentoutput_dir=output_dir,newdevice=FALSE,ylabel=ylab_text,cex.val=cex.plots,sample.col.opt=sample.col.opt),silent=TRUE)
                get_individualsampleplots(feature_table_file,class_labels_file,X=goodfeats_temp,Y=classlabels_orig,parentoutput_dir=output_dir,newdevice=FALSE,ylabel=ylab_text,cex.val=cex.plots,sample.col.opt=individualsampleplot.col.opt)
                
                
               
               if(output.device.type!="pdf"){
                    
                       try(dev.off(),silent=TRUE)
                    }
               
                
	 	if(globalclustering==TRUE){
            
             if(output.device.type!="pdf"){
                     
                     temp_filename_1<-"Figures/GlobalclusteringEM.png"
                     
                     png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                 }
                 m1<-Mclust(t(data_m_fc_withfeats[,-c(1:2)]))
                 s1<-m1$classification #summary(m1)
                 
                 EMcluster<-m1$classification
                 
                 col_vec <- colorRampPalette(brewer.pal(10, "RdBu"))(length(levels(as.factor(classlabels_orig[,2]))))
                 #col_vec<-topo.colors(length(levels(as.factor(classlabels_orig[,2])))) #patientcolors #heatmap_cols[1:length(levels(classlabels_orig[,2]))]
                 t1<-table(EMcluster,classlabels_orig[,2])
                 
                 par(mfrow=c(1,1))
                 plot(t1,col=col_vec,main="EM cluster labels\n using all features",cex.axis=1,ylab="Class",xlab="Cluster number")
                 
                 par(xpd=TRUE)
                 try(legend("bottomright",legend=levels(classlabels_orig[,2]),text.col=col_vec,pch=13,cex=0.4),silent=TRUE)
                 
                par(xpd=FALSE)
                 
                 t1<-cbind(EMcluster,classlabels_orig[,2])
                 write.table(t1,file="Tables/EM_clustering_labels_using_allfeatures.txt",sep="\t")
                 
                 if(output.device.type!="pdf"){
                     
                     try(dev.off(),silent=TRUE)
                 }
                 
                 if(output.device.type!="pdf"){
                     
                     temp_filename_1<-"Figures/GlobalclusteringHCA.png"
                     
                     png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                 }
                 
                 #if(FALSE)
                 {
                 
                 #p1<-heatmap.2(as.matrix(data_m_fc_withfeats[,-c(1:2)]),scale="row",symkey=FALSE,col=topo.colors(n=256))
                 
                 if(heatmap.col.opt=="RdBu"){
                     
                     heatmap.col.opt="redblue"
                 }
                 
                 heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
                 heatmap_cols<-rev(heatmap_cols)
                 
                 if(heatmap.col.opt=="topo"){
                     heatmap_cols<-topo.colors(256)
                     heatmap_cols<-rev(heatmap_cols)
                 }else
                 {
                     if(heatmap.col.opt=="heat"){
                         heatmap_cols<-heat.colors(256)
                         heatmap_cols<-rev(heatmap_cols)
                     }else{
                         
                         if(heatmap.col.opt=="yellowblue"){
                             
                             heatmap_cols<-colorRampPalette(c("yellow","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                             #heatmap_cols<-blue2yellow(256) #colorRampPalette(c("yellow","blue"))(256)
                             heatmap_cols<-rev(heatmap_cols)
                         }else{
                             
                             if(heatmap.col.opt=="redblue"){
                                 
                                 heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
                                 heatmap_cols<-rev(heatmap_cols)
                             }else{
                                 
                                 #my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                                 if(heatmap.col.opt=="redyellowgreen"){
                                     
                                     heatmap_cols <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
                                     heatmap_cols<-rev(heatmap_cols)
                                 }else{
                                     if(heatmap.col.opt=="yellowwhiteblue"){
                                         
                                         heatmap_cols<-colorRampPalette(c("yellow2","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                                         heatmap_cols<-rev(heatmap_cols)
                                     }else{
                                         
                                         if(heatmap.col.opt=="redwhiteblue"){
                                             
                                             heatmap_cols<-colorRampPalette(c("red","white","blue"))(256) #colorRampPalette(c("yellow","white","blue"))(256)
                                             heatmap_cols<-rev(heatmap_cols)
                                         }else{
                                             
                                             
                                             
                                             heatmap_cols <- colorRampPalette(brewer.pal(10, heatmap.col.opt))(256)
                                             heatmap_cols<-rev(heatmap_cols)
                                             
                                         }
                                         
                                     }
                                     
                                 }
                                 
                             }
                             
                         }
                     }
                     
                     
                 }
                 
                 
                 
                 #col_vec<-heatmap_cols[1:length(levels(classlabels_orig[,2]))]
                 c1<-WGCNA::cor(as.matrix(data_m_fc_withfeats[,-c(1:2)]),method="pearson",use="pairwise.complete.obs") #cor(d1[,-c(1:2)])
                 d2<-as.dist(1-c1)
                 clust1<-hclust(d2)
                 
                 hr <- try(hclust(as.dist(1-WGCNA::cor(t(data_m_fc_withfeats),method="pearson",use="pairwise.complete.obs"))),silent=TRUE) #metabolites
                 #hc <- try(hclust(as.dist(1-WGCNA::cor(data_m,method=cor.method,use="pairwise.complete.obs"))),silent=TRUE) #samples
                 
                 
                 h73<-heatmap.2(as.matrix(data_m_fc_withfeats[,-c(1:2)]), Rowv=as.dendrogram(hr), Colv=as.dendrogram(clust1),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main="Global clustering\n using all features", ColSideColors=patientcolors,labRow = FALSE, labCol = FALSE)
                 
                 # par(xpd=TRUE)
                 #legend("bottomleft",legend=levels(classlabels_orig[,2]),text.col=unique(patientcolors),pch=13,cex=0.4)
                 #par(xpd=FALSE)

                 clust_res<-cutreeDynamic(distM=as.matrix(d2),dendro=clust1)
                 
                 #mycl_samples <- cutree(clust1, h=max(clust1$height)/2)
                 
                 HCAcluster<-clust_res
                 
                 c2<-cbind(clust1$labels,HCAcluster)
                 
                 rownames(c2)<-c2[,1]
                 
                 c2<-as.data.frame(c2)

                 t1<-table(HCAcluster,classlabels_orig[,2])
                 
                 plot(t1,col=col_vec,main="HCA (Cutree Dynamic) cluster labels\n using all features",cex.axis=1,ylab="Class",xlab="Cluster number")
                 
                 par(xpd=TRUE)
                 try(legend("bottomright",legend=levels(classlabels_orig[,2]),text.col=col_vec,pch=13,cex=0.4),silent=TRUE)
                 par(xpd=FALSE)

                 t1<-cbind(HCAcluster,classlabels_orig[,2])
                write.table(t1,file="Tables/HCA_clustering_labels_using_allfeatures.txt",sep="\t")
                 }
                 
                 
                 
                 if(output.device.type!="pdf"){
                     
                     try(dev.off(),silent=TRUE)
                 }
		}

}
		#dev.off()
}
else
{
    #goodfeats_allfields<-as.data.frame(goodfeats)

    goodfeats<-goodfeats[,-c(1:time_ind)]
    
    
}

}

rm(goodfeats_temp)
				if(length(goodip)>0){
				
				try(dev.off(),silent=TRUE)
				}
			}
			else{
				try(dev.off(),silent=TRUE)
				break;
			}
	

	if(analysismode=="classification" & WGCNAmodules==TRUE){
                print("Doing wgcna")
                classlabels_temp<-classlabels_orig_wgcna #cbind(classlabels_sub[,1],classlabels)
                #print(classlabels_temp)
		data_temp<-data_matrix_beforescaling[,-c(1:2)]


        cl<-makeSOCKcluster(num_nodes)

		#clusterExport(cl,"do_rsd")
                #feat_rsds<-parApply(cl,data_temp,1,do_rsd)
		#rm(data_temp)
                #feat_rsds<-abs(feat_rsds) #round(max_rsd,2)
		#print(summary(feat_rsds))
		#if(length(which(feat_rsds>0))>0)
		{
		X<-data_m_fc_withfeats #data_matrix[which(feat_rsds>=wgcnarsdthresh),]

#	print(head(X))
#		print(dim(X))


if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/WGCNA_preservation_plot.png"
    
    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
}


		print("generating preservation plot")
                #preservationres<-try(do_wgcna(X=X,Y=classlabels_temp,sigfeats=data_m_fc_withfeats[goodip,c(1:2)]),silent=TRUE)
		pres<-try(do_wgcna(X=X,Y=classlabels_temp,sigfeats=data_m_fc_withfeats[goodip,c(1:2)]),silent=TRUE)
		

        if(output.device.type!="pdf"){
            
            try(dev.off(),silent=TRUE)
        }
        	}
	}	
		#print(lf)
	
		#print("next iteration")
		#dev.off()
	}
	
setwd(parentoutput_dir)
summary_res<-cbind(log2.fold.change.thresh_list,feat_eval,feat_sigfdrthresh,feat_sigfdrthresh_cv,feat_sigfdrthresh_permut,res_score_vec)

if(fdrmethod=="none"){
    exp_fp<-round(fdrthresh*feat_eval)
}else{
exp_fp<-round(fdrthresh*feat_sigfdrthresh)
}
rank_num<-order(summary_res[,5],decreasing=TRUE)


if(featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma2wayrepeat" | featselmethod=="lmreg" | featselmethod=="logitreg" 
| featselmethod=="lm2wayanova" | featselmethod=="lm1wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="wilcox" | featselmethod=="ttest" | featselmethod=="poissonreg")
{
	summary_res<-cbind(summary_res,exp_fp)

colnames(summary_res)<-c("RSD.thresh","Number of features left after RSD filtering","Number of features selected",paste(pred.eval.method,"-accuracy",sep=""),paste(pred.eval.method," permuted accuracy",sep=""),"Score","Expected_False_Positives")
}else{
	#exp_fp<-round(fdrthresh*feat_sigfdrthresh)
	
	#if(featselmethod=="MARS" | featselmethod=="RF" | featselmethod=="pls" | featselmethod=="o1pls" | featselmethod=="o2pls"){
		
		exp_fp<-rep(NA,dim(summary_res)[1])
	#}
	summary_res<-cbind(summary_res,exp_fp)
		colnames(summary_res)<-c("RSD.thresh","Number of features left after RSD filtering","Number of features selected",paste(pred.eval.method,"-accuracy",sep=""),paste(pred.eval.method," permuted accuracy",sep=""),"Score","Expected_False_Positives")

	}

featselmethod<-parentfeatselmethod
file_name<-paste("Results_summary_",featselmethod,".txt",sep="")
write.table(summary_res,file=file_name,sep="\t",row.names=FALSE)


if(output.device.type=="pdf"){

        try(dev.off(),silent=TRUE)
}

print("##############Level 1: processing complete###########")


if(length(best_feats)>0)
{

mz_index<-best_feats

if(analysismode=="classification"){
    
    
     
log2.fold.change.thresh=log2.fold.change.thresh_list[best_logfc_ind]
	
print(paste("Best results found at RSD threshold ", log2.fold.change.thresh,sep=""))
#print(paste(kfold,"-fold CV accuracy ", best_acc,sep=""))
	if(pred.eval.method=="CV"){
	
	print(paste(kfold,"-fold CV accuracy: ", best_acc,sep=""))
	}else{
	if(pred.eval.method=="AUC"){
	
	print(paste("Area under the curve (AUC) is : ", best_acc,sep=""))
	}
	}


	data_m<-parent_data_m
	data_m_fc<-data_m #[which(abs(mean_groups)>log2.fold.change.thresh),]

	data_m_fc_withfeats<-data_matrix[,c(1:2)]
	data_m_fc_withfeats<-cbind(data_m_fc_withfeats,data_m_fc)
	

#when using a feature table generated by apLCMS

rnames<-paste("mzid_",seq(1,dim(data_m_fc)[1]),sep="")
#print(best_limma_res[1:3,])

goodfeats<-best_limma_res[order(best_limma_res$mz),-c(1:2)]

#goodfeats<-best_limma_res[,-c(1:2)]

goodfeats_all<-goodfeats

goodfeats<-goodfeats_all
rm(goodfeats_all)
}
if(length(featselmethod)>1){
    abs.cor.thresh=NA
    globalcor=FALSE
}

if(globalcor==TRUE){
    

if(length(best_feats)>2){
if(is.na(abs.cor.thresh)==FALSE){
setwd(parentoutput_dir)
print("##############Level 2: Metabolome wide correlation network analysis of differentially expressed metabolites###########")
print(paste("Generating metabolome-wide ",cor.method," correlation network using RSD threshold ", log2.fold.change.thresh," results",sep=""))
data_m_fc_withfeats<-as.data.frame(data_m_fc_withfeats)

goodfeats<-as.data.frame(goodfeats)
#print(goodfeats[1:4,])
sigfeats_index<-which(data_m_fc_withfeats$mz%in%goodfeats$mz)
sigfeats<-sigfeats_index
if(globalcor==TRUE){

#outloc<-paste(parentoutput_dir,"/Allcornetworksigfeats","log2fcthresh",log2.fold.change.thresh,"/",sep="")
outloc<-paste(parentoutput_dir,"/Stage3","/",sep="")

dir.create(outloc)
setwd(outloc)

if(networktype=="complete"){
	mwan_fdr<-do_cor(data_matrix,subindex=sigfeats_index,targetindex=NA,outloc,networkscope="global",cor.method,abs.cor.thresh,cor.fdrthresh,max.cor.num,net_node_colors,net_legend)
	}else{
	if(networktype=="GGM"){
	mwan_fdr<-get_partial_cornet(data_matrix, sigfeats.index=sigfeats_index,targeted.index=NA,networkscope="global",cor.method,abs.cor.thresh,cor.fdrthresh,outloc=outloc,net_node_colors,net_legend)
	}else{
		print("Invalid option. Please use complete or GGM.")
	}
	}
	
print("##############Level 2: processing complete###########")
}else{
	print("##############Skipping Level 2: global correlation analysis###########")
	}

print("#########################")
	if(is.na(target.metab.file)==FALSE){
	
	print("##############Level 3: Targeted correlation network analysis of differentially expressed metabolites###########")
	setwd(parentoutput_dir)

	#outloc<-paste(parentoutput_dir,"/Targetedcornetworksigfeats","log2fcthresh",log2.fold.change.thresh,"/",sep="")
	outloc<-paste(parentoutput_dir,"/Stage4","/",sep="")

	dir.create(outloc)
	setwd(outloc)
	
	print(outloc)
	
		print(paste("Searching for metabolites matching target list",sep=""))
	
	if(is.na(target.metab.file)==FALSE){
		
		#print(target.metab.file)
	dataA<-read.table(target.metab.file,sep="\t",header=TRUE)
	dataA<-unique(dataA)
	
	#print(dim(dataA))
	
	if(dim(dataA)[2]<2){
		
		dataA<-cbind(dataA,dataA)
		colnames(dataA)<-c("mz","mzduplicate")
	}
	
	#print(target.metab.file)
	#dataA<-
	#dataA<-as.data.frame(dataA)
	#dataA<-apply(dataA,2,as.numeric)
	
	#print(head(dataA))
	
g1<-getVenn(dataA=dataA,name_a="TargetSet",name_b="ExperimentalSet",dataB=data_matrix[,c(1:2)],mz.thresh=target.mzmatch.diff,time.thresh=target.rtmatch.diff,
xMSanalyzer.outloc=outloc,alignment.tool=NA)
#names(g1)



	
	if(length(g1$common)>1){
	
	if(is.na(sigfeats)==FALSE){
		
	#	print(dim(goodfeats))
	#	print(head(goodfeats))
	com_mzs<-find.Overlapping.mzs(dataA=data_matrix,dataB=goodfeats,mz.thresh=1,time.thresh=1,alignment.tool=NA)
	
	sigfeats_index<-com_mzs$index.A #which(data_matrix$mz%in%goodfeats$mz)
	print(paste(length(unique(sigfeats_index))," variables matched the selected features list",sep=""))
	

	}else{
		sigfeats_index<-NA #seq(1,dim(data_matrix)[1])
		
		}
	
	print(paste(length(unique(g1$common$index.B))," variables matched the target list",sep=""))
	
		print(paste("Generating targeted correlation network",sep=""))
	
	targeted_mat<-dataA[g1$common$index.A,]
	data_commat<-data_matrix[g1$common$index.B,]
	overlap_mat<-merge(targeted_mat,data_commat,by="mz")
	write.table(overlap_mat,file="matching_targeted_mz_data.txt",sep="\t",row.names=FALSE)
	#print(data_matrix[1:3,1:5])
	#print(length(sigfeats_index))
	if(networktype=="complete"){
	targetedan_fdr<-do_cor(data_matrix,subindex=sigfeats_index,targetindex=g1$common$index.B,outloc,networkscope="targeted",cor.method,abs.cor.thresh,cor.fdrthresh,max.cor.num,net_node_colors,net_legend)
    
    if(FALSE){
    pdf("Cornetworkplot.pdf")
    load("metabnet.Rda")
    print(plot(net_result))
    try(dev.off(),silent=TRUE)
    }
    
	}else{
	if(networktype=="GGM"){
	targetedan_fdr<-get_partial_cornet(data_matrix, sigfeats.index=sigfeats_index,targeted.index=g1$common$index.B,networkscope="targeted",cor.method,abs.cor.thresh,cor.fdrthresh,outloc=outloc,net_node_colors)
    
    if(FALSE){
    pdf("GGMnetworkplot.pdf")
    load("metabnet.Rda")
    print(plot(net_result))
    try(dev.off(),silent=TRUE)
    }
    
	}else{
		print("Invalid option. Please use complete or GGM.")
	}
	}
	}else{
	print(paste("Targeted metabolites were not found.",sep=""))
	}
	}
	
	if(analysismode=="classification"){
		classlabels_temp<-cbind(classlabels_sub[,1],classlabels)
		#do_wgcna(X=data_matrix,Y=classlabels,sigfeats.index=sigfeats_index)
	}

	print("##############Level 3: processing complete###########")
	print("#########################")
	}

	
}
}
else{
	print(paste("Can not perform network analysis. Too few metabolites.",sep=""))
}
}

print(paste(featselmethod, " processing done.",sep=""))

}
setwd(parentoutput_dir)




#print("Note A: Please note that log2 fold-change based filtering is only applicable to two-class comparison. 
#log2fcthresh of 0 will remove only those features that have exactly sample mean intensities between the two groups.
#More features will be filtered prior to FDR as log2fcthresh increases.")

#print("Note C: Please make sure all the packages are installed. You can use the command install.packages(packagename) to install packages.")
#print("Eg: install.packages(\"mixOmics\"),install.packages(\"snow\"), install.packages(\"e1071\"), biocLite(\"limma\"), install.packages(\"gplots\").")
#print("Eg: install.packages("mixOmics""),install.packages("snow"), install.packages("e1071"), biocLite("limma"), install.packages("gplots").")
##############################
##############################
###############################



	if(length(best_feats)>0){
	
	goodfeats<-as.data.frame(goodfeats)
    #goodfeats<-data_matrix_beforescaling[which(data_matrix_beforescaling$mz%in%goodfeats$mz),]
	}else{
		goodfeats-{}
		}

	cur_date<-Sys.time()
	cur_date<-gsub(x=cur_date,pattern="-",replacement="")
	cur_date<-gsub(x=cur_date,pattern=":",replacement="")
	cur_date<-gsub(x=cur_date,pattern=" ",replacement="")
    if(saveRda==TRUE){
	fname<-paste("Analysis_",featselmethod,"_",cur_date,".Rda",sep="")
	#save(list=ls(),file=fname)
    }
			################################


	return(list("diffexp_metabs"=goodfeats_allfields,  "mw.an.fdr"=mwan_fdr,"targeted.an.fdr"=targetedan_fdr,"classlabels"=classlabels_orig,"all_metabs"=allmetabs_res))


}
