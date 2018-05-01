do_plsda <-
function(X,Y,oscmode="pls",numcomp=3,kfold=10,evalmethod="CV",keepX=15,sparseselect=FALSE,analysismode="classification",vip.thresh=1,sample.col.opt="default",sample.col.vec=c("red","green","blue","purple"),scoreplot_legend=TRUE,feat_names=NA,pairedanalysis=FALSE,optselect=FALSE,class_labels_levels_main=NA,legendlocation="bottomleft",plotindiv=TRUE,pls.vip.selection="max",output.device.type="pdf",plots.res=600,plots.width=8,plots.height=8,plots.type="cairo")
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


                                    save(linn.pls,file="linnpls.Rda")
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

                    save(linn.pls,file="linnpls.Rda")



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
            good_feats<-which(linn.vip[c1]>vip.thresh)
            
        }
        
        }else{
            
            if(opt_comp>1){
                linn.vip.mean<-apply(linn.vip,1,mean)
                
                good_feats<-which(linn.vip.mean[c1]>vip.thresh)
            }else{
                good_feats<-which(linn.vip[c1]>vip.thresh)
                
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
        #save(linn.pls2,file="linnpls2.Rda")
        #save(linn.pls3,file="linnpls3.Rda")
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
               
                write.table(r2_q2_valid_res,file="pls_r2_q2_res_usingallfeatures.txt",sep="\t",row.names=TRUE)
                    if(plotindiv==TRUE){
                barplot(r2_q2_valid_res[1:2,],beside=TRUE,main="PLS leave-one-out validation diagnostics using all features",ylab="Variation",col=c("darkgrey","lightgrey"))
                legend("topright",c("R2","Q2"),col=c("darkgrey","lightgrey"),pch=c(20))
            
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
                    barplot(r2_q2_valid_res[1:2,],beside=TRUE,main="PLS loo validation diagnostics \n using all features",ylab="Variation",col=c("darkgrey","lightgrey"))
                    legend("topright",c("R2","Q2"),col=c("darkgrey","lightgrey"),pch=c(20))
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
		
               write.table(r2_q2_valid_res,file="pls_r2_q2_res_usingselectfeats.txt",sep="\t",row.names=TRUE)
               
               if(plotindiv==TRUE){
                   barplot(r2_q2_valid_res[1:2,],beside=TRUE,main="PLS loo validation diagnostics \n using selected features",ylab="Variation",col=c("darkgrey","lightgrey"))
                   legend("topright",c("R2","Q2"),col=c("darkgrey","lightgrey"),pch=c(20))
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
       
        # save(list=ls(),file="debug.Rda")
        if(plotindiv==TRUE){
	
    
    if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/PLS_pairwise_component_plots.pdf"
        
        pdf(temp_filename_1)
        #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        }
        get_plsplots(X,plsres=linn.pls,plsvar=pls_var,samplelabels=Yclass,filename=NA,ncomp=opt_comp,center=TRUE,scale=TRUE,legendcex=0.5,outloc=getwd(),col_vec=col.stimu,sample.col.opt="default",alphacol=0.3,legendlocation="topright",class_levels=class_labels_levels)
    #,silent=TRUE)
    
    if(output.device.type!="pdf"){
        
       try(dev.off(),silent=TRUE)
    }
    
    
    
	}
        #legend = c(class_labels_levels), cex = 0.55)
    }
    
    write.table(linn.pls$variates$X,file="pls_scores.txt",sep="\t")
    write.table(linn.pls$loadings$X,file="pls_loadings.txt",sep="\t")
    #save(linn.pls,file="pls_res.Rda")
    return(list("model"=linn.pls,"vip_res"=linn.vip,"valid_res"=v1,"cv_res"=cv_res,"opt_comp"=opt_comp,"selected_variables"=good_feats,"bad_variables"=bad_variables))
    
}
