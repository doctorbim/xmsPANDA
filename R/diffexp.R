diffexp <-
function(Xmat=NA,Ymat=NA,feature_table_file,parentoutput_dir,class_labels_file,num_replicates=3,summarize.replicates=TRUE,summary.method="mean",summary.na.replacement="zeros",missing.val=0,rep.max.missing.thresh=0.3,
 all.missing.thresh=0.1,group.missing.thresh=0.7,input.intensity.scale="raw",
log2transform=TRUE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=TRUE,lowess_norm=FALSE,madscaling=FALSE,rsd.filt.list=1,
pairedanalysis=FALSE,featselmethod=c("limma","pls"),fdrthresh=0.05,fdrmethod="BH",cor.method="spearman",networktype="complete",abs.cor.thresh=0.4,cor.fdrthresh=0.05,kfold=10,
pred.eval.method="BER",globalcor=TRUE,
target.metab.file=NA,target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=100, 
numtrees=20000,analysismode="classification",net_node_colors=c("green","red"), net_legend=TRUE,
svm_kernel="radial",heatmap.col.opt="redblue",manhattanplot.col.opt=c("darkblue","red3"),boxplot.col.opt=c("grey57"),sample.col.opt="rainbow",rf_selmethod="rankbased",pls_vip_thresh=1,num_nodes=2,max_varsel=100,pls_ncomp=5,pca.stage2.eval=TRUE,scoreplot_legend=TRUE,pca.global.eval=TRUE,rocfeatlist=seq(2,11,1),rocfeatincrement=TRUE,rocclassifier="svm",foldchangethresh=2,wgcnarsdthresh=20,WGCNAmodules=TRUE,optselect=TRUE,max_comp_sel=1,saveRda=TRUE,legendlocation="topleft",pca.cex.val=4,
pca.ellipse=FALSE,ellipse.conf.level=0.95,permutations.count=1000,svm.acc.tolerance=5,limmadecideTests=TRUE,pls.vip.selection="max",globalclustering=TRUE,plots.res=600,plots.width=8,plots.height=8,plots.type="cairo",output.device.type="pdf",pvalue.thresh=0.05,individualsampleplot.col.opt=NA,pamr.threshold.select.max=FALSE,aggregation.method="RankAggreg",...)
{
    
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
    pls.permut.count=permutations.count
    
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
    
    #WORKING
    if(featselmethod=="lm2wayanovarepeat" | featselmethod=="lm1wayanovarepeat" | featselmethod=="spls1wayrepeat" | featselmethod=="limma2wayrepeat")
    {
        print("Note 3: Class labels format should be: Sample ID, Subject, Factor, Time. lm2wayanovarepeat and lm1wayanovarepeat include post-hoc comparisons.")
    }
    print("#####Starting processing now##############################################")
    print(paste("Program is running. Please check the logfile for runtime status: ",fname,sep=""))
    
    fname_params<-paste(parentoutput_dir,"/InputParameters.csv",sep="")
    #sink(fname_params)
    # save(list=ls(),file="cur.Rda")
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
    c1<-rbind(c1,"lowess_norm:")
    c2<-rbind(c2,lowess_norm)
    c1<-rbind(c1,"madscaling:")
    c2<-rbind(c2,madscaling)
    c1<-rbind(c1,"rsd.filt.list:")
    c2<-rbind(c2,rsd.filt.list)
    c1<-rbind(c1,"pairedanalysis:")
    c2<-rbind(c2,pairedanalysis)
    c1<-rbind(c1,"featselmethod:")
    c2<-rbind(c2,featselmethod)
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
    c1<-rbind(c1,"feat_weight:")
    c2<-rbind(c2,feat_weight)
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
for(i in 1:length(featselmethod))
{
    
                    outloc<-paste(parentoutput_dir,featselmethod[i],sep="/")
            diffexp.res[[i]]<-diffexp.child(Xmat,Ymat,feature_table_file,parentoutput_dir,class_labels_file,num_replicates,feat.filt.thresh,summarize.replicates,summary.method,summary.na.replacement,missing.val,rep.max.missing.thresh,
             all.missing.thresh,group.missing.thresh,input.intensity.scale,
            log2transform,medcenter,znormtransform,quantile_norm,lowess_norm,madscaling,rsd.filt.list,
            pairedanalysis,featselmethod[i],fdrthresh,fdrmethod,cor.method,networktype,abs.cor.thresh,cor.fdrthresh,kfold,pred.eval.method,feat_weight,globalcor,
            target.metab.file,target.mzmatch.diff,target.rtmatch.diff,max.cor.num,samplermindex,pcacenter,pcascale,
            numtrees,analysismode,net_node_colors,net_legend,svm_kernel,heatmap.col.opt,manhattanplot.col.opt,boxplot.col.opt,sample.col.opt,alphacol,rf_selmethod,pls_vip_thresh,num_nodes,max_varsel, pls_ncomp=pls_ncomp,pca.stage2.eval=pca.stage2.eval,scoreplot_legend=scoreplot_legend,pca.global.eval=pca.global.eval,rocfeatlist=rocfeatlist,rocfeatincrement=rocfeatincrement,rocclassifier=rocclassifier,foldchangethresh=foldchangethresh,wgcnarsdthresh=wgcnarsdthresh,WGCNAmodules=WGCNAmodules,
            optselect=optselect,max_comp_sel=max_comp_sel,saveRda=saveRda,legendlocation=legendlocation,degree_rank_method=degree_rank_method,pca.cex.val=pca.cex.val,pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,pls.permut.count=pls.permut.count,svm.acc.tolerance=svm.acc.tolerance,limmadecideTests=limmadecideTests,pls.vip.selection=pls.vip.selection,globalclustering=globalclustering,plots.res=plots.res,plots.width=plots.width,plots.height=plots.height,plots.type=plots.type,output.device.type=output.device.type,pvalue.thresh,individualsampleplot.col.opt,pamr.threshold.select.max)

        diffexp.res[[i]]$all_metabs<-diffexp.res[[i]]$all_metabs[order(diffexp.res[[i]]$all_metabs$mz),]
        mz_rt_all<-paste(diffexp.res[[i]]$all_metabs$mz,"_",diffexp.res[[i]]$all_metabs$time,sep="")
        tname<-paste("mz_rt_",i,".Rda",sep="")
        #save(mz_rt_all,file=tname)

if(i==1){
    
    #sort(rankingCriteria, index.return = TRUE)$ix
        ranked_list<-mz_rt_all[sort(diffexp.res[[i]]$all_metabs$diffexp_rank,index.return=TRUE)$ix]
}

        if(i>1){
            mz_rt_1<-paste(diffexp.res[[(i-1)]]$diffexp_metabs$mz,"_",diffexp.res[[(i-1)]]$diffexp_metabs$time,sep="")
            mz_rt_2<-paste(diffexp.res[[i]]$diffexp_metabs$mz,"_",diffexp.res[[i]]$diffexp_metabs$time,sep="")
            cnamesd1<-colnames(diffexp.res[[(i-1)]]$diffexp_metabs)
            time_ind<-which(cnamesd1=="time")
            mz_ind<-which(cnamesd1=="mz")
            common_feat_ind<-which(mz_rt_1%in%mz_rt_2)
            common_feat_ind2<-which(mz_rt_2%in%mz_rt_1)
            
            
            ranked_list<-rbind(ranked_list,mz_rt_all[sort(diffexp.res[[i]]$all_metabs$diffexp_rank,index.return=TRUE)$ix])
            
            if(length(common_feat_ind)<1){
             ermsg<-paste("No common significant features found between ",featselmethod[(i-1)]," and ",featselmethod[i],sep="")
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
            mat_A<-mat_A[,c(1:(mz_ind-1))]
            common_feats<-cbind(mat_A,mat_B[m1,]) #merge(mat_B,mat_A,by.x="mz_rt_2",by.y="mz_rt_1")
            print(paste("Number of common sig feats between ",featselmethod[(i-1)]," and ",featselmethod[i],sep=""))
            print(dim(common_feats))
            
            }
        }

}
if(consensus_analysis==TRUE){
        
        
    
    if(aggregation.method=="RankAggreg"){
    save(ranked_list,file="ranked_list.Rda")
    r1<-RankAggreg(x=ranked_list,k=max_varsel,verbose=TRUE,distance="Spearman",method="CE")
    
    #r1<-list(top.list=order(degree_rank))
    print("Aggregated rank")
    
    common_row_index<-which(mz_rt_all%in%r1$top.list)
    
    common_feats<-diffexp.res[[i]]$all_metabs[common_row_index,]
    
    cnamesd1<-colnames(common_feats)
    time_ind<-which(cnamesd1=="time")
    mz_ind<-which(cnamesd1=="mz")
    #common_feats<-common_feats[,c(1:(mz_ind-1))]
    
    print("Aggregated sig feats list")
    print(r1$top.list)
    
    }
    cnames_1<-colnames(common_feats)
    
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


        print("Dimension of consensus significant feature table")
        print(dim(common_feats))
        
        num_common_feats<-dim(common_feats)[1]
        
        if(num_common_feats<1){
            
            stop("No common features found.")
        }
        
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
        
        
        write.table(common_feats,file="Aggregated_significant_features.txt",sep="\t",row.names=FALSE)
        
        subdata<-t(Xmat[,-c(1:2)])
        classlabels<-Ymat[,2]
        classlabels<-as.data.frame(classlabels)
        svm_model<-try(svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
        #
        
        
        
        if(is(svm_model,"try-error")){
            kfold_acc<-NA
            
        }else{
            kfold_acc<-svm_model$avg_acc
        }
        
       numcores<-round(detectCores()*0.6)
        
        kfold_acc_rand<-{}
        #for(p1 in 1:100){
         cl <- makeCluster(getOption("cl.cores", numcores))
         clusterEvalQ(cl,library(e1071))
         clusterEvalQ(cl,library(pROC))
          clusterEvalQ(cl,library(ROCR))
         clusterExport(cl,"svm_cv",envir = .GlobalEnv)
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
         stopCluster(cl)
        
        kfold_acc_rand<-unlist(kfold_acc_rand)
        
        kfold_acc_rand<-mean(kfold_acc_rand,na.rm=TRUE)
        
        summary_res<-cbind(num_common_feats,kfold_acc,kfold_acc_rand)
        colnames(summary_res)<-c("Number of significant features after aggregation",paste(pred.eval.method,"-accuracy",sep=""),paste(pred.eval.method," permuted accuracy",sep=""))

        
        file_name<-paste("../Results_summary_aggregated.txt",sep="")
        write.table(summary_res,file=file_name,sep="\t",row.names=FALSE)
        
        if(output.device.type=="pdf"){
        pdf("Figures/Aggregatedresults.pdf")
        }
        
        if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/HCA_aggregated_sigfeats.png"
            
            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        }
        
        #   print("here")
        #print(head(Ymat))
        #save(Xmat,file="Xmat.Rda")
        #save(Ymat,file="Ymat.Rda")
        
        g1<-get_hca(feature_table_file=NA,parentoutput_dir=outloc,class_labels_file=NA,X=Xmat,Y=Ymat,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode=analysismode,
        sample.col.opt=sample.col.opt,plots.width=plots.width,plots.height=plots.height,plots.res=plots.res, plots.type=plots.type, alphacol=0.3, hca_type="two-way",newdevice=FALSE,input.type="intensity",mainlab="",cexRow=0.4,cexCol=0.4)
        
        if(output.device.type!="pdf"){
            
            dev.off()
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
                
                temp_filename_1<-"Figures/kfold_forward_selection_aggregated_sigfeats.png"
                
                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
            }
            
            
            plot(x=xvec,y=yvec,main=msg1,xlab="Feature index",ylab=ylab_text,type="b",col="brown")
            
            
            if(output.device.type!="pdf"){
                
                dev.off()
            }


            cv_mat<-cbind(xvec,yvec)
            colnames(cv_mat)<-c("Feature Index",ylab_text)
            
            write.table(cv_mat,file="aggregated_kfold_cv_mat.txt",sep="\t")
        }
        
        # save(Xmat,file="Xmat.Rda")
        #save(Ymat,file="Ymat.Rda")
        #save(sample.col.opt,file="sample.col.opt.Rda")
        
        #print("consensus boxplots")
        
        if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/Boxplots_aggregated_sigfeats.png"
            
            #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
            pdf(temp_filename_1)
        }
        
        
        
        get_boxplots(X=Xmat,Y=Ymat,parentoutput_dir=outloc,sample.col.opt=sample.col.opt,boxplot.col.opt=boxplot.col.opt, alphacol=0.3,newdevice=FALSE,cex=0.7)
        
       
        
        if(output.device.type!="pdf"){
            
            dev.off()
        }
        
        # print("Consensus Narplots")
         
         if(output.device.type!="pdf"){
             
             temp_filename_1<-"Figures/Barplots_aggregated_sigfeats.pdf"
             
             #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
             pdf(temp_filename_1)
         }
         
         
         
            get_barplots(feature_table_file=NA,class_labels_file=NA,X=Xmat,Y=Ymat,parentoutput_dir=outloc,newdevice=FALSE,ylabel="intensity",cex.val=0.9,barplot.col.opt=boxplot.col.opt)

if(output.device.type!="pdf"){
    
    dev.off()
}


if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/Individual_sample_plots_aggregated_sigfeats.png"
    
    #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
    pdf(temp_filename_1)
}


get_individualsampleplots(feature_table_file=NA,class_labels_file=NA,X=Xmat,Y=Ymat,parentoutput_dir=outloc,newdevice=FALSE,
ylabel="intensity",cex.val=0.9,sample.col.opt=sample.col.opt)
            
            
              if(output.device.type!="pdf"){
        
            dev.off()
              }
 
 #if(pairedanalysis==TRUE)
 {
     
     # print("Consensus timeseriesplots")
      
      if(output.device.type!="pdf"){
          
          temp_filename_1<-"Figures/Lineplots_aggregated_sigfeats.png"
          pdf(temp_filename_1)
          #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
      }
      
      
     
 get_lineplots(X=Xmat,Y=Ymat,feature_table_file=NA,parentoutput_dir=getwd(),class_labels_file=NA,sample.col.opt=sample.col.opt, alphacol=0.3,col_vec=col_vec,pairedanalysis=pairedanalysis,pca.cex.val=pca.cex.val,legendlocation=legendlocation,pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,filename="significant")  #,silent=TRUE)
 }
 
    if(output.device.type!="pdf"){
        
        dev.off()
    }
 }
 
 
 if(output.device.type=="pdf"){
        try(dev.off(),silent=TRUE)
 }
        }

suppressWarnings(sink(file=NULL))

print("###############################")
print("###############################")
print("###############################")
print("Program ended successfully.")
print(paste("All the result files are in the specified output location: ",parentoutput_dir,sep=""))
print("There will be a sub-folder for each step: pre-processing, statistical comparisons (limma/RF/MARS), and network analysis (Global and/or Targeted).")
print("Network files that can be used with Cytoscape and CIRCOS tool are provided in the output folder.")
print("An Rda file that includes all the main objects generated during the analysis for future reference or further processing is also written.")
print("Enjoy!")
print("###############################")


s1<-"Stage 1 results: Preprocessing (Normalization, transformation)"
s2<-"Stage 2 results: Feature selection & evaluation results (Manhattan plots, PCA, HCA, boxplots, table of significant features, clustering results)"
s3<-"Stage 3 results: Correlation based network analysis"
s4<-"Stage 4 results: Correlation based targeted network analysis"
s5<-"Consensus results: HCA, k-fold CV, boxplots, and barplots for aggregated significant features"
sm<-rbind(s1,s2,s3,s4,s5)
setwd(parentoutput_dir)
write.table(sm,file="Readme.txt",sep="\t",row.names=FALSE)
return(list("individual.featsel.res"=diffexp.res,"aggregated.res"=common_feats))
}else{
    
    suppressWarnings(sink(file=NULL))
    
    print("###############################")
    print("###############################")
    print("###############################")
    print("Program ended successfully. Consensus analysis could not be performed as not all features were selected by all feature selection methods.")
    print(paste("All the result files are in the specified output location: ",parentoutput_dir,sep=""))
    print("There will be a sub-folder for each step: pre-processing, statistical comparisons (limma/RF/MARS), and network analysis (Global and/or Targeted).")
    print("Network files that can be used with Cytoscape and CIRCOS tool are provided in the output folder.")
    print("An Rda file that includes all the main objects generated during the analysis for future reference or further processing is also written.")
    print("Enjoy!")
    print("###############################")
    
    
    s1<-"Stage 1 results: Preprocessing (Normalization, transformation)"
    s2<-"Stage 2 results: Feature selection & evaluation results (Manhattan plots, PCA, HCA, boxplots, table of significant features, clustering results)"
    s3<-"Stage 3 results: Correlation based network analysis"
    s4<-"Stage 4 results: Correlation based targeted network analysis"
    
    sm<-rbind(s1,s2,s3,s4)
    setwd(parentoutput_dir)
    write.table(sm,file="Readme.txt",sep="\t",row.names=FALSE)
    return(list("individual.featsel.res"=diffexp.res))
    
}


    
    }else{
        
        
        diffexp.res<-diffexp.child(Xmat,Ymat,feature_table_file,parentoutput_dir,class_labels_file,num_replicates,feat.filt.thresh,summarize.replicates,summary.method,summary.na.replacement,missing.val,rep.max.missing.thresh,
        all.missing.thresh,group.missing.thresh,input.intensity.scale,
        log2transform,medcenter,znormtransform,quantile_norm,lowess_norm,madscaling,rsd.filt.list,
        pairedanalysis,featselmethod,fdrthresh,fdrmethod,cor.method,networktype,abs.cor.thresh,cor.fdrthresh,kfold,pred.eval.method,feat_weight,globalcor,
        target.metab.file,target.mzmatch.diff,target.rtmatch.diff,max.cor.num,samplermindex,pcacenter,pcascale,
        numtrees,analysismode,net_node_colors,net_legend,svm_kernel,heatmap.col.opt,manhattanplot.col.opt,boxplot.col.opt,sample.col.opt,alphacol,rf_selmethod,pls_vip_thresh,num_nodes,max_varsel, pls_ncomp,pca.stage2.eval,scoreplot_legend,pca.global.eval,rocfeatlist,rocfeatincrement,rocclassifier,foldchangethresh,wgcnarsdthresh,WGCNAmodules,
        optselect,max_comp_sel,saveRda,legendlocation,degree_rank_method,pca.cex.val,pca.ellipse,ellipse.conf.level,pls.permut.count,svm.acc.tolerance,limmadecideTests,pls.vip.selection,globalclustering,plots.res,plots.width,plots.height,plots.type,output.device.type,pvalue.thresh,individualsampleplot.col.opt,pamr.threshold.select.max)
        
        
        suppressWarnings(sink(file=NULL))
        
        print("###############################")
        print("###############################")
        print("###############################")
        print("Program ended successfully.")
        print(paste("All the result files are in the specified output location: ",parentoutput_dir,sep=""))
        print("There will be a sub-folder for each step: pre-processing, statistical comparisons (limma/RF/MARS), and network analysis (Global and/or Targeted).")
        print("Network files that can be used with Cytoscape and CIRCOS tool are provided in the output folder.")
        print("An Rda file that includes all the main objects generated during the analysis for future reference or further processing is also written.")
        print("Enjoy!")
        print("###############################")
        
        
        s1<-"Stage 1 results: Preprocessing (Normalization, transformation)"
        s2<-"Stage 2 results: Feature selection & evaluation results (Manhattan plots, PCA, HCA, boxplots, table of significant features, clustering results)"
        s3<-"Stage 3 results: Correlation based network analysis"
        s4<-"Stage 4 results: Correlation based targeted network analysis"
        sm<-rbind(s1,s2,s3,s4)
        setwd(parentoutput_dir)
        write.table(sm,file="Readme.txt",sep="\t",row.names=FALSE)
        return(diffexp.res)
        
        
    }
}
