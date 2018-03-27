runlmreg <-
function(X,Y,fdrmethod="BH",fdrthresh=0.05,pvalue.thresh=0.05){
    
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
        
        #print(data_mat_anova)
        
        anova_res<-diffexplmreg(dataA=data_mat_anova,logistic_reg)
        
        return(anova_res)
    })
    
    save(res1,file="lmregres.Rda")
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
                dev.off()
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
