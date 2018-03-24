diffexp.child <-
function(Xmat,Ymat,feature_table_file,parentoutput_dir,class_labels_file,num_replicates,feat.filt.thresh,summarize.replicates,summary.method,
summary.na.replacement,missing.val,rep.max.missing.thresh,
 all.missing.thresh,group.missing.thresh,input.intensity.scale,
log2transform,medcenter,znormtransform,quantile_norm,lowess_norm,madscaling,rsd.filt.list,
pairedanalysis,featselmethod,fdrthresh,fdrmethod,cor.method,networktype,abs.cor.thresh,cor.fdrthresh,kfold,pred.eval.method,feat_weight,globalcor,
target.metab.file,target.mzmatch.diff,target.rtmatch.diff,max.cor.num, samplermindex,pcacenter,pcascale,
numtrees,analysismode,net_node_colors,net_legend,svm_kernel,heatmap.col.opt,manhattanplot.col.opt,boxplot.col.opt,sample.col.opt,alphacol,rf_selmethod,pls_vip_thresh,num_nodes,max_varsel,pls_ncomp,pca.stage2.eval,scoreplot_legend,pca.global.eval,rocfeatlist,rocfeatincrement,rocclassifier,foldchangethresh,wgcnarsdthresh,WGCNAmodules,optselect,max_comp_sel,saveRda,legendlocation,degree_rank_method,pca.cex.val,pca.ellipse,ellipse.conf.level,pls.permut.count,svm.acc.tolerance,limmadecideTests,pls.vip.selection,globalclustering,plots.res,plots.width,plots.height,plots.type,output.device.type,pvalue.thresh,individualsampleplot.col.opt,pamr.threshold.select.max)
{
    
   
   
		#############
        
        remove_firstrun=FALSE #TRUE or FALSE
        run_number=1
        minmaxtransform=FALSE
        pca.CV=TRUE
        max_rf_var=5000
        
        logistic_reg=FALSE
        
        goodfeats_allfields={}
        mwan_fdr={}
        targetedan_fdr={}
        data_m_fc_withfeats={}
        classlabels_orig={}
        
        if(input.intensity.scale=="log2"){
            
            log2transform=FALSE
        }

        parentfeatselmethod<-featselmethod
        
        if(featselmethod=="rf"){
            
            featselmethod="RF"
            
            rfconditional=FALSE
        }else{
            
            if(featselmethod=="rfconditional"){
                
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
						
							if(featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="rfesvm" | featselmethod=="wilcox" | featselmethod=="wilcoxrepeat" | featselmethod=="pamr"){
								print("##########Level 1: Finding discriminatory metabolites ###########")

							if(featselmethod=="logitreg"){

								featselmethod="lmreg"
								logistic_reg=TRUE
                            }else{
                                logistic_reg=FALSE
                                
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
		X<-read.table(feature_table_file,sep="\t",header=TRUE,quote = "")
		
			
            #print(X[1:3,])
		
		#X<-X[order(X$mz),]
		
        X[,1]<-round(X[,1],5)
        X[,2]<-round(X[,2],1)
        
        
		Xmat<-t(X[,-c(1:2)])
		
		Xmat<-as.data.frame(Xmat)
		
        #print(Xmat[1:3,1:5])
		
				}else{
					X<-Xmat
                    
                    X[,1]<-round(X[,1],5)
                    X[,2]<-round(X[,2],1)
                    
					Xmat<-t(X[,-c(1:2)])
		
					Xmat<-as.data.frame(Xmat)
		}
		
						print("Feature selection method:")
						print(featselmethod)			
		
        #save(Xmat,file="Xmat.Rda")
						
	if(analysismode=="regression")
	{
	
					#log2.fold.change.thresh_list<-c(0)
					
				print("Performing regression analysis")
				if(is.na(Ymat)==TRUE){
								classlabels<-read.table(class_labels_file,sep="\t",header=TRUE,quote = "")
						
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
 
                        if(is.na(all(diff(match(rnames_xmat,rnames_ymat))))==FALSE){
                                    if(all(diff(match(rnames_xmat,rnames_ymat)) > 0)==TRUE){
                                        
                                        data_matrix<-data_preprocess(Xmat,Ymat,feature_table_file,parentoutput_dir,class_labels_file=NA,num_replicates=num_replicates,feat.filt.thresh,summarize.replicates,
                                        summary.method, all.missing.thresh,group.missing.thresh=NA,
                                        log2transform,medcenter,znormtransform,quantile_norm,lowess_norm,madscaling,missing.val, samplermindex,rep.max.missing.thresh,summary.na.replacement,featselmethod)
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
                        featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="rfesvm" | featselmethod=="pamr")
                    {
                        # print("here")
                            
                        stop(paste(featselmethod," for paired analysis, featselmethod should be: limma1wayrepeat, limma2wayrepeat, lm1wayanovarepeat, lm2wayanovarepeat, spls1wayrepeat, or spls2wayrepeat",sep=""))
                        
                    }
                    
                    
                }
		
        #	print(head(Ymat))
								if(is.na(Ymat)==TRUE){
											classlabels<-read.table(class_labels_file,sep="\t",header=TRUE,quote = "")
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
									
		if(featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma2wayrepeat" | featselmethod=="limma1way" | featselmethod=="limma1wayrepeat" | featselmethod=="MARS" | featselmethod=="RF" | featselmethod=="pls" | featselmethod=="o1pls" | featselmethod=="o2pls" | featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="spls" | featselmethod=="pls1wayrepeat" | featselmethod=="spls1wayrepeat" | featselmethod=="pls2wayrepeat" | featselmethod=="spls2wayrepeat" | featselmethod=="pls2way" | featselmethod=="spls2way" | featselmethod=="o1spls" | featselmethod=="o2spls" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="rfesvm" | featselmethod=="wilcox" | featselmethod=="ttest" | featselmethod=="pamr")
		{
			#analysismode="classification"

                #if(is.na(Ymat)==TRUE)
				{
						#classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)

						if(analysismode=="classification"){
						if(featselmethod=="lmreg" | featselmethod=="logitreg")
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
							if(featselmethod=="limma" | featselmethod=="limma1way" | featselmethod=="MARS" | featselmethod=="RF" | featselmethod=="pls" | featselmethod=="o1pls" | featselmethod=="o2pls" | featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="spls" | featselmethod=="o1spls" | featselmethod=="o2spls" | featselmethod=="rfesvm" | featselmethod=="pamr")
							{
								
								
								if(featselmethod=="lmreg" | featselmethod=="logitreg")
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
								
                                #save(Xmat,file="Xmat.Rda")
                                
                                #save(classlabels,file="Xmat_classlabels.Rda")
                                
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
								
                                
                                # save(Xmat_temp,file="Xmat_temp.Rda")
                                
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
											
                                            # save(Xmat_temp,file="Xmat_temp.Rda")
                                             
                                             
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
              print(classlabels)
					#rownames(Xmat)<-as.character(classlabels[,1])
					
					#write.table(classlabels,file="organized_classlabels.txt",sep="\t",row.names=FALSE)
					
					Xmat1<-cbind(classlabels,Xmat)
                    #write.table(Xmat1,file="organized_featuretable.txt",sep="\t",row.names=TRUE)

						featselmethod="limma2way"
						pairedanalysis = TRUE
										
								}
								else{
									print(classlabels)
									stop("Only one factor specificied in the class labels file.")			
								}
										}
								
						
						
					}
                #else{
                #	classlabels<-Ymat
                    

#	}
				classlabels<-as.data.frame(classlabels)
				
				
				#colnames(classlabels)<-c("SampleID","Class")
				#f1<-table(classlabels$SampleID)
				#Ymat=classlabels
				
				#		classlabels<-as.data.frame(classlabels)
				#	classlabels_response_mat<-classlabels[,-c(1)]
				#	classlabels_response_mat<-as.data.frame(classlabels_response_mat)
			
					
				
			
			
				
				if(featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="pls2way" | featselmethod=="spls2way" | featselmethod=="wilcox"){
					
					analysismode="classification"

					classlabels<-read.table(class_labels_file,sep="\t",header=TRUE,quote = "")
					
     #               write.table(classlabels,file="original_classlabelsA.txt",sep="\t",row.names=TRUE)
                    

							#cnames[2]<-"Factor1"
														
							cnames<-colnames(classlabels)
							
							factor_inf<-classlabels[,-c(1)]
							factor_inf<-as.data.frame(factor_inf)
							
							colnames(classlabels)<-c("SampleID",paste("Factor",seq(1,dim(factor_inf)[2]),sep=""))
							
							analysismode="classification"

							Xmat_temp<-Xmat #t(Xmat)
							Xmat_temp<-cbind(classlabels,Xmat_temp)
							
                            
                            # save(Xmat_temp,file="Xmat_temp.Rda")
                             
							if(featselmethod=="lm1wayanova" | featselmethod=="wilcox"){
						
	
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
				
						if(featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="pls1wayrepeat" | featselmethod=="spls1wayrepeat" | featselmethod=="pls2wayrepeat" | featselmethod=="spls2wayrepeat" | featselmethod=="wilcoxrepeat"){
							
							analysismode="classification"

							classlabels<-read.table(class_labels_file,sep="\t",header=TRUE,quote = "")
					
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
							if(featselmethod=="lm1wayanovarepeat" | featselmethod=="pls1wayrepeat" | featselmethod=="spls1wayrepeat" | featselmethod=="wilcoxrepeat"){
								
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

#save(Ymat,file="Ymat.Rda")
#                                       save(Xmat,file="Xmat.Rda")

										if(featselmethod=="spls1wayrepeat"){
                                                                                                featselmethod="spls"
												
                                                                                        }else{
                                                                                                if(featselmethod=="pls1wayrepeat"){
                                                                                                featselmethod="pls"
                                                                                        }														
											}
                                                                                        
                                                                                        
                                                                                        if(featselmethod=="wilcoxrepeat"){
                                                                                            
                                                                                            featselmethod=="wilcox"
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
                                    
                                    if(featselmethod=="lmreg" | featselmethod=="logitreg"){
                                        
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
            
            
            #      write.table(Xmat,file="organized_featuretableA1.txt",sep="\t",row.names=TRUE)
            
            #write.table(Xmat_temp,file="organized_featuretableA2.txt",sep="\t",row.names=TRUE)
            #write.table(classlabels,file="organized_classlabelsA1.txt",sep="\t",row.names=TRUE)
            
            #write.table(Ymat,file="organized_Ymat1.txt",sep="\t",row.names=TRUE)
            
            colnames(Xmat)<-as.character(Ymat[,1])
            
			Xmat<-cbind(X[,c(1:2)],Xmat)
						
						Xmat<-as.data.frame(Xmat)
						Ymat<-as.data.frame(Ymat)
						
					#	print(Ymat[1:3,])
					#	print(Xmat[1:3,1:5])
                    
		
        #print("Before preprocessing")
        # 	print("Ymat")
        #	print(Ymat[1:3,])
        #	print("Xmat:")
        #	print(Xmat[1:3,])
            
            #   print("dim xmat")
            #print(dim(Xmat))
            
            
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
                            
                            if(is.na(all(diff(match(rnames_xmat,rnames_ymat))))==FALSE){
                                if(all(diff(match(rnames_xmat,rnames_ymat)) > 0)==TRUE){
                                    
                                    data_matrix<-data_preprocess(Xmat=Xmat,Ymat=Ymat,feature_table_file=NA,parentoutput_dir,class_labels_file=NA,num_replicates=num_replicates,feat.filt.thresh,summarize.replicates,summary.method,all.missing.thresh,group.missing.thresh,
                                    log2transform,medcenter,znormtransform,quantile_norm,lowess_norm,madscaling,missing.val, samplermindex,rep.max.missing.thresh,summary.na.replacement,featselmethod)
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

#print("here")
#	print(dim(data_matrix))
		  data_matrix_beforescaling<-data_matrix$data_matrix_prescaling
  
    data_matrix_beforescaling<-as.data.frame( data_matrix_beforescaling)
   data_matrix<-data_matrix$data_matrix_afternorm_scaling
   
   #data_matrix<-data_matrix[c(1:1000),]
   
	data_m<-data_matrix[,-c(1:2)]
	
	#print(classlabels[1:10,])
	#print(data_matrix[1:10,1:4])
	
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
			
            #	print("After preprocessing")
            #print("Ymat")
            #print(Ymat[1:3,])
            #print("Xmat:")
            #print(Xmat[1:3,1:5])
			
            
            
            #       		write.table(Ymat,file="ordered_classlabels.txt",sep="\t",row.names=TRUE)
	
    # write.table(classlabels_orig,file="original_classlabels.txt",sep="\t",row.names=TRUE)

			rnames1<-as.character(Ymat[,1])
			rnames2<-as.character(classlabels_orig[,1])
            #print("rownames")
            #print(head(rnames1))
            #print(head(rnames2))
			sorted_index<-{}
			for(i in 1:length(rnames1)){

				
				#sorted_index<-c(sorted_index,grep(rnames1[i],pattern=paste("^",rnames2,"$",sep="")))

				sorted_index<-c(sorted_index,grep(x=rnames2,pattern=paste("^",rnames1[i],"$",sep="")))

			}
			classlabels_orig<-classlabels_orig[sorted_index,]

            #write.table(classlabels_response_mat,file="original_classlabelsB.txt",sep="\t",row.names=TRUE)
			classlabelsA<-classlabels
			
			#print(classlabels)
			#print(which(duplicated(classlabels)==TRUE))
			
			if(length(which(duplicated(classlabels)==TRUE))>0){
				rownames(classlabels)<-paste("S",seq(1,dim(classlabels)[1]),sep="")
			}else{
				rownames(classlabels)<-as.character(classlabels[,1])
			}#as.character(classlabels[,1])
			#print(classlabels)
			#print(classlabels[1:10,])
            
 #           save(classlabels,file="classlabels.Rda")
  #          save(classlabels_orig,file="classlabels_orig.Rda")
   #         save(classlabels_response_mat,file="classlabels_response_mat.Rda")
            
            if(pairedanalysis==TRUE){
                
                save(subject_inf,file="subjectinf.Rda")
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

							#save(classlabels,file="classlabels_1.Rda")
							#save(class_labels_levels,file="class_labels_levels.Rda")

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
	mean_overall<-apply(data_temp,1,do_mean)

	print("mean overall")
	print(summary(mean_overall))
	bad_feat<-which(mean_overall==0)

	if(length(bad_feat)>0){

		data_matrix_beforescaling<-data_matrix_beforescaling[-bad_feat,]
		data_m<-data_m[-bad_feat,]
		data_matrix<-data_matrix[-bad_feat,]
	
	} 
	
	
	#Step 5) RSD/CV calculation
	
	if(FALSE){
		if(length(class_labels_levels)==2){
		if(log2transform==TRUE){
		
		
		data_temp<-data_matrix_beforescaling[,-c(1:2)]
		
		mean_groupA<-apply(data_temp[,1:num_samps_group[[1]]],1,do_mean)
		mean_groupB<-apply(data_temp[,(num_samps_group[[1]]+1):(num_samps_group[[1]]+num_samps_group[[2]])],1,do_mean)
		#print(summary(mean_groupA))
		#print(summary(mean_groupB))
		#print(head(data_temp))
		#mean_groups<-log2((mean_groupA/(mean_groupB)))
		
		mean_groups<-mean_groupA-mean_groupB
		
		}else
		{
			data_temp<-data_matrix_beforescaling[,-c(1:2)]

			mean_groupA<-apply(data_temp[,1:num_samps_group[[1]]],1,do_mean)
			mean_groupB<-apply(data_temp[,(num_samps_group[[1]]+1):(num_samps_group[[1]]+num_samps_group[[2]])],1,do_mean)
			mean_groups<-log2((mean_groupA)/(mean_groupB))
	
		}
		
		}else{
			print("More than 2 classes found. Skipping fold change calculation.")
			#log2.fold.change.thresh_list<-c(0)
			mean_groups<-rep(1000,dim(data_m)[1])
			
			}
	}
						
	
	}else{
		
		classlabels<-(classlabels[,-c(1)])
		}
		
        #	print("######classlabels#########")
        #print(classlabels)
	
    
    class_labels_levels_new<-levels(classlabels)
    
    if(analysismode=="classification"){
    test_classlabels<-cbind(class_labels_levels_main,class_labels_levels_new)
    }
    
    # print("test class labels")
    #print(test_classlabels)

		#data_temp<-data_matrix_beforescaling[,-c(1:2)]
		#	feat_rsds<-apply(data_temp,1,do_rsd)
			
		#	print("Summary of RSD across all features:")
		#	print(summary(feat_rsds))

	parent_data_m<-round(data_m,2)
	
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

s1<-apply(X,1,sum)

#print(length(which(s1==0)))

s1<-apply(X,2,sum)

#print(length(which(s1==0)))
#print(dim(X))

 #  print(X[1:10,1:10])
   
	#if(pca.CV==FALSE){
   
  # rm(pcaMethods)
  #try(detach("package:pcaMethods",unload=TRUE),silent=TRUE)
   library(mixOmics)
 
if(FALSE){  
  if(znormtransform==FALSE & medcenter==FALSE){
    
   # metabpcaresultlog2allmetabs5pcs<-pca(X,ncomp=3,center=pcacenter,scale=pcascale)
    
    metabpcaresultlog2allmetabs5pcs<-pca(X,ncomp=3,center=FALSE,scale=TRUE)
    
}else{
    
    if(medcenter==FALSE & znormtransform==TRUE){
        
       # metabpcaresultlog2allmetabs5pcs<-pca(X,ncomp=3,center=pcacenter,scale=FALSE)
	
	metabpcaresultlog2allmetabs5pcs<-pca(X,ncomp=3,center=FALSE,scale=TRUE)
    }else{
        
        if(medcenter==TRUE & znormtransform==TRUE){
            metabpcaresultlog2allmetabs5pcs<-pca(X,ncomp=3,center=FALSE,scale=FALSE)
	    
        }else{
        
        if(znormtransform==FALSE & medcenter==TRUE & pcacenter==FALSE){
            metabpcaresultlog2allmetabs5pcs<-pca(X,ncomp=3,center=FALSE,scale=pcascale)
        }else{
	
	
            metabpcaresultlog2allmetabs5pcs<-pca(X,ncomp=3,center=pcacenter,scale=pcascale)
        
			}
    }
    }
    
}
}
#if(pca.global.eval==TRUE)
 if(FALSE)
    {
	metabpcaresultlog2allmetabs5pcs<-pca(X,ncomp=5,center=TRUE,scale=TRUE)

#}else{
	#metabpcaresultlog2allmetabs5pcs_1<-pcaMethods::pca(X,ncomp=3,center=pcacenter)
#	}

#metabpcaresultlog2allmetabs5pcs<-pca(X,ncomp=5,center=FALSE,scale=TRUE)

#metabpcaresultnotransform10pcsallmetabs<-metabpcaresult
result<-metabpcaresultlog2allmetabs5pcs
}

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

#col_vec<-heat.colors(length(class_labels_levels)) #heat.colors(length(class_labels_levels))

if(is.na(individualsampleplot.col.opt)==TRUE){
    
    individualsampleplot.col.opt=col_vec
}

#c("mediumpurple4","mediumpurple1","blueviolet","darkblue","blue","cornflowerblue","cyan4","skyblue",
#"darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
#"red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
#"aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")

 #p2<-c("lightblue", "mistyrose", "lightcyan", "lavender")

#if(pca.global.eval==TRUE)

data_temp<-data_matrix_beforescaling[,-c(1:2)]



cl<-makeSOCKcluster(num_nodes)

clusterExport(cl,"do_rsd")
			feat_rsds<-parApply(cl,data_temp,1,do_rsd)
			
			stopCluster(cl)
			
		sum_rsd<-summary(feat_rsds,na.rm=TRUE)
		max_rsd<-max(feat_rsds,na.rm=TRUE)
		max_rsd<-round(max_rsd,2)
		
			print("Summary of RSD across all features:")
			print(sum_rsd)

if(log2.fold.change.thresh_list[length(log2.fold.change.thresh_list)]>max_rsd){
stop(paste("The maximum relative standard deviation threshold in rsd.filt.list should be below ",max_rsd,sep=""))
}

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
print(subject_inf)
}


for(lf in 1:length(log2.fold.change.thresh_list))
#for(lf in 1:length(var.thresh_list))
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
                output_dir<-paste(output_dir1,featselmethod,"signalthresh",group.missing.thresh,"RSD",log2.fold.change.thresh,"/",sep="")
                
                
            }
            
            dir.create(output_dir,showWarnings=FALSE)

			setwd(output_dir)
            
            dir.create("Figures")
            
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
            
            if(dim(data_m_fc)[2]>50){
                
                if(output.device.type!="pdf"){
                    
                    temp_filename_1<-"Figures/SampleIntensityDistribution.png"
                    
                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                
                }
                
                size_num<-min(100,dim(data_m_fc)[2])
                
                print(size_num)
                samp_index<-sample(x=1:dim(data_m_fc)[2],size=size_num)
                try(boxplot(data_m_fc[,samp_index],main="Intensity distribution across samples after preprocessing",xlab="Samples",ylab=ylab_text,col=boxplot.col.opt),silent=TRUE)
                
                if(output.device.type!="pdf"){
                    
                    dev.off()
                }
                
            }else{
                
                if(output.device.type!="pdf"){
                    
                    temp_filename_1<-"Figures/SampleIntensityDistribution.png"
                    
                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                }
                try(boxplot(data_m_fc,main="Intensity distribution across samples after preprocessing",xlab="Samples",ylab=ylab_text,col=boxplot.col.opt),silent=TRUE)
                if(output.device.type!="pdf"){
                    
                    dev.off()
                }
                
            }
           
            
            
			data_m_fc_withfeats<-cbind(data_m_fc_withfeats,data_m_fc)
			
			
			feat_eval[lf]<-0
			res_score_vec[lf]<-0
			#feat_sigfdrthresh_cv[lf]<-0
			
			#filename<-paste(fheader,log2.fold.change.thresh,".txt",sep="")
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



                                p1<-pcaMethods::pca(tempX,method="rnipals",center=TRUE,scale="uv",cv="q2",nPcs=10)
                       
                                if(output.device.type!="pdf"){

                                temp_filename_2<-"Figures/PCAdiagnostics_allfeats.png"

                                png(temp_filename_2,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                }


                                p2<-plot(p1,col=c("darkgrey","grey"),main="PCA diagnostics using all features after RSD filtering")
                                if(output.device.type!="pdf"){

                                dev.off()
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

                                        if(featselmethod=="lmreg" || featselmethod=="logitreg"){
                                            classlabels_orig<-classlabels_orig[,c(1:2)]
                                            classlabels_orig<-as.data.frame(classlabels_orig)
                                        }
                                }



                                if(output.device.type!="pdf"){

                                temp_filename_1<-"Figures/PCAscore_distribution_allfeats.pdf"

                                #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")

                                pdf(temp_filename_1)
                                }

                                #added this
                                # save(classlabels_orig,file="classlabels_orig.Rda")

                                try(get_pcascoredistplots(X=data_m_fc_withfeats,Y=classlabels_orig,feature_table_file=NA,parentoutput_dir=getwd(),class_labels_file=NA,sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,col_vec=col_vec,pairedanalysis=pairedanalysis,pca.cex.val=pca.cex.val,legendlocation=legendlocation,pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,filename="all"),silent=TRUE)

                                if(output.device.type!="pdf"){

                                    dev.off()
                                }


                                classlabels_orig<-classlabels_orig_parent
                   }else{
                       #regression
                            tempgroup<-rep("A",dim(data_m_fc)[2])
                            col_vec1<-rep("black",dim(data_m_fc)[2])
                            class_labels_levels_main1<-c("A")
                            
                            # print(sample.col.opt)
                            
                                                     get_pca(X=data_m_fc,samplelabels=tempgroup,legendlocation=legendlocation,filename="all",ncomp=3,center=TRUE,scale=TRUE,legendcex=0.5,outloc=getwd(),col_vec=col_vec1,sample.col.opt=sample.col.opt,alphacol=0.3,class_levels=NA,pca.cex.val=pca.cex.val,pca.ellipse=FALSE)
                            
                            
                            
                        }


if(featselmethod=="pamr"){
    
    #print("HERE")
    # save(data_m_fc,classlabels,file="pamdebug.Rda")
    
    pamr.res<-do_pamr(X=data_m_fc,Y=classlabels,fdrthresh=fdrthresh,nperms=pls.permut.count,pamr.threshold.select.max=pamr.threshold.select.max)
    
    goodip<-pamr.res$feature.list
    
    save(pamr.res,file="pamr.res.Rda")
    feature_rowindex<-seq(1,nrow(data_m_fc))
    
    discore<-rep(0,nrow(data_m_fc))
    
    discore[goodip]<-pamr.res$max.discore
    
    
    sel.diffdrthresh<-feature_rowindex%in%goodip
    
    rank_vec<-rank(-discore)
    
    selected_id_withmztime<-cbind(data_m_fc_withfeats[goodip,c(1:2)],pamr.res$pam_toplist)
    
    write.csv(selected_id_withmztime,file="dscores.selectedfeats.csv",row.names=FALSE)
    data_limma_fdrall_withfeats<-cbind(rank_vec,data_m_fc_withfeats)
    write.table(data_limma_fdrall_withfeats,file="pamr_ranked_feature_table.txt",sep="\t",row.names=FALSE)
    
    
}


if(featselmethod=="rfesvm"){
    
    #print(classlabels)
    
    if(length(class_labels_levels)<3){
        featureRankedList = svmrfeFeatureRanking(x=t(data_m_fc),y=classlabels,svmkernel=svm_kernel,kfold=kfold,pred.eval.method=pred.eval.method)
        
        best_subset<-featureRankedList$best_subset
        featureRankedList<-featureRankedList$featureRankedList
    }else{
        
        featureRankedList = svmrfeFeatureRankingForMulticlass(x=t(data_m_fc),y=classlabels,svmkernel=svm_kernel,kfold=kfold,pred.eval.method=pred.eval.method)
        best_subset<-featureRankedList$best_subset
        featureRankedList<-featureRankedList$featureRankedList
    }
    
    rank_vec<-seq(1,dim(data_m_fc_withfeats)[1])
    goodip<-best_subset
    dtemp1<-data_m_fc_withfeats[goodip,]
    save(dtemp1,file="svm_sel.Rda")
    save(classlabels,file="classlabels.Rda")
    
    sel.diffdrthresh<-rank_vec%in%goodip
     
     data_m_fc_withfeats<-data_m_fc_withfeats #[featureRankedList,]
    
    data_m_fc<-data_m_fc #[featureRankedList,]
    
    rank_vec<-seq(1,dim(data_m_fc_withfeats)[1])
    #rank_vec<-rank_vec[featureRankedList]
    
    rank_vec<-sort(featureRankedList,index.return=TRUE)$ix
    
    data_limma_fdrall_withfeats<-cbind(rank_vec,data_m_fc_withfeats)
    
    #save(sigfeats_order,file="rfesvm_ranked_feature_table.Rda")
    save(featureRankedList,file="rfesvm_featureRank.Rda")
    
    #write.table(data_limma_fdrall_withfeats,file="rfesvm_ranked_feature_table.txt",sep="\t",row.names=FALSE)
    
   
    
    #print(best_subset)
    #print(sel.diffdrthresh)
    #print(goodip)
    
    save(rank_vec,file="rank_vec.Rda")
    save(best_subset,file="best_subset.Rda")
    save(sel.diffdrthresh,file="sel.diffdrthresh.Rda")
    save(goodip,file="goodip.Rda")
}

						
								#print("starting this")
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
												print("Desing matrix")
												print(design)
                                                print(cont.matrix)
												#print(data_m[1:10,1:10])
											
                                            #save(list=ls(),file="limma.Rda")
												#save(f,file="f.Rda")
												#save(design,file="design.Rda")
												#save(data_m_fc,file="data_m_fc.Rda")
												#save(subject_inf,file="subject_inf.Rda")	
												#f1<-c(seq(1,num_samps_group[[1]]),seq(1,num_samps_group[[1]]))
												corfit<-duplicateCorrelation(data_m_fc,design=design,block=subject_inf,ndups=1)
												
                                                #	print(corfit$consensus)

												fit<-lmFit(data_m_fc,design,block=f1,cor=corfit$consensus)

											}else{
											
											
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
																dev.off()
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
                                                
                                                dev.off()
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
										filename<-"limma_posthoc1wayanova_results.txt"
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
																	dev.off()
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
                                               
                                               dev.off()
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
                                
                                #print(data_limma_fdrall_withfeats[1:10,1:10])
                                
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
                                    
                                    temp_filename_1<-"Figures/HCA_Factor1sigfeats.png"
                                    
                                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                }
                                                                get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X1,Y=Y1,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                                sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type="two-way",newdevice=FALSE,input.type="intensity",mainlab="Factor1")
                                
                               if(output.device.type!="pdf"){
                                    
                                    dev.off()
                                }

                                
                                #Y2=cbind(classlabels_orig[,1],classlabels_orig[,2])
                                
                                if(pairedanalysis==FALSE){
                                    Y2=classlabels_orig[,c(1,3)] #cbind(classlabels_orig[,1],classlabels_orig[,3])
                                }else{
                                    Y2=classlabels_orig[,c(1,4)] #cbind(classlabels_orig[,1],classlabels_orig[,4])
                                    
                                }
                                
                                Y2<-as.data.frame(Y2)
                                
                                if(output.device.type!="pdf"){
                                    
                                    temp_filename_1<-"Figures/HCA_Factor2sigfeats.png"
                                    
                                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                }
                               
                                                                get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X2,Y=Y2,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                                sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type="two-way",newdevice=FALSE,input.type="intensity",mainlab="Factor2")
                                
                                 if(output.device.type!="pdf"){
                                    
                                    dev.off()
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
                              
                              #save(classlabels_orig,file="classlabels_orig.Rda")
                              
                              if(output.device.type!="pdf"){
                                  
                                  temp_filename_1<-"Figures/HCA_Factor1xFacto2sigfeats.png"
                                  
                                  png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                              }
                                 get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X3,Y=Y3,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                                sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type="two-way",newdevice=FALSE,input.type="intensity",mainlab="Factor1 x Factor2")
                                
                                 if(output.device.type!="pdf"){
                                    
                                    dev.off()
                                }


                                
										filename<-"limma_2wayanova_posthocresults.txt"
										
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
								
								#tiff("pval_dist.tiff",width=plots.width,height=plots.height,res=plots.res, compression="lzw")
								#hist(d4[,1],xlab="p",main="Distribution of p-values")
								#dev.off()
								
								results2<-decideTests(fit2,method="nestedF",adjust.method=fdrmethod,p.value=fdrthresh)
									
									
									
								}
									
									
					
                        if(featselmethod=="RF")
                        {
							
							maxint<-apply(data_m_fc,1,max)
							
							print("Performing random forest analysis using the cforest function")
							
                            data_m_fc_withfeats<-as.data.frame(data_m_fc_withfeats)
                            
                            data_m_fc<-as.data.frame(data_m_fc)
                            write.table(classlabels,file="classlabels_rf.txt",sep="\t",row.names=FALSE)
                  
                            save(data_m_fc,classlabels,numtrees,analysismode,file="rfdebug.Rda")
                            
                            if(rfconditional==TRUE){
                            rfcondres1<-do_rf_conditional(X=data_m_fc,classlabels,ntrees=numtrees,analysismode) #,silent=TRUE)
                            filename<-"RFconditional_variable_importance_allfeats.txt"
                            }else{
                                
                                rfcondres1<-do_rf(X=data_m_fc,classlabels,ntrees=numtrees,analysismode)
                                filename<-"RF_variable_importance_allfeats.txt"
                            }
                            
                            
                            #print(rfcondres1)
                            varimp_res2<-rfcondres1$rf_varimp
                            
                            
                            
                            write.table(varimp_res2, file=filename,sep="\t",row.names=FALSE)
                        
                            min_varimp<-min(varimp_res2,na.rm=TRUE)
                         
                    
                            if(rf_selmethod=="absVIMthresh"){
					     
                                                 goodip<-which(varimp_res2>abs(min_varimp))
                                                 
                                                  if(length(goodip)<1){
                                                    print("No significant variables found.")
                                                }
                                                var_names<-paste(sprintf("%.3f",data_m_fc_withfeats[,1]),sprintf("%.1f",data_m_fc_withfeats[,2]),sep="_")

                                                names(varimp_res2)<-as.character(var_names)
                                                sel.diffdrthresh<-varimp_res2>abs(min_varimp)
                                                
                                                        if(length(which(sel.diffdrthresh==TRUE))<1){
                                                            print("No significant variables found.")
                                                        }
									
					     }else{
					     		
					     				
					     				if(rf_selmethod=="rankbased"){
					     				
					     					
					     					if(max_varsel>dim(data_m_fc)[1]){
					     						max_varsel=dim(data_m_fc)[1]
					     					}
					     					
					     					sorted_varimp_res<-varimp_res2[order(varimp_res2,decreasing=TRUE)[1:(max_varsel)]]
					     				
					     					goodip<-order(varimp_res2,decreasing=TRUE)[1:(max_varsel)]
					     					sel.diffdrthresh<-varimp_res2>=min(sorted_varimp_res,na.rm=TRUE)
					     					
					     				}
					     			
					     	
					     	}
					   #  print(goodip)
					    
					   
					  


						num_var_rf<-length(which(sel.diffdrthresh==TRUE))
						
						if(num_var_rf>10){
							
							num_var_rf=10
						}
						sorted_varimp_res<-varimp_res2[order(varimp_res2,decreasing=TRUE)[1:(num_var_rf)]]

						sorted_varimp_res<-sort(sorted_varimp_res)
		
						barplot_text=paste("Variable Importance measures (top ",length(sorted_varimp_res)," shown)\n",sep="")
						
				        if(output.device.type!="pdf"){
                            
                            temp_filename_1<-"Figures/RF_sigfeats_VIMbarplot.png"
                            
                            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                        }
                        

						barplot(sorted_varimp_res, xlab="Significant features", main=barplot_text,cex.axis=0.5,cex.names=0.4, ylab="raw permutation accuracy\n variable importance measure")
					
                    if(output.device.type!="pdf"){
                        
                        dev.off()
                    }




                    rank_num<-rank(-varimp_res2)
                
					data_limma_fdrall_withfeats<-cbind(varimp_res2,rank_num,data_m_fc_withfeats)
					
					cnames_tab<-colnames(data_m_fc_withfeats)
					cnames_tab<-c("Raw VIM","Rank",cnames_tab)
					
					 goodip<-which(sel.diffdrthresh==TRUE)
					    
					feat_sigfdrthresh[lf]<-length(which(sel.diffdrthresh==TRUE))
						
					   
					data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(data_limma_fdrall_withfeats$mz),]
					
					colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
					
					#data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
                    #write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
					

						}
							
						if(featselmethod=="MARS"){
							
							#print(classlabels)
						#varimp_res<-do_mars(X=data_m_fc,classlabels,ntrees=numtrees)
						marsres1<-do_mars(X=data_m_fc,classlabels, analysismode,kfold)
						
						varimp_marsres1<-marsres1$mars_varimp
						
						#print(varimp_marsres1[1:10,])
						
						mars_mznames<-rownames(varimp_marsres1)
						
						#all_names<-rownames(data_m_fc)
						all_names<-paste("mz",seq(1,dim(data_m_fc)[1]),sep="")
					
						com1<-match(all_names,mars_mznames)
						
						#g1<-grep(pattern="NA",x=mars_mznames)
						#if(length(g1)>0){
						#	varimp_marsres1<-varimp_marsres1[-g1,]
						#}
						
						#print(varimp_res2[1:10])
						#print(summary(varimp_res2))

						filename<-"MARS_variable_importance.txt"
						
						#print(com1[1:10])
						
						#sel.diffdrthresh<-com1
						sel.diffdrthresh<-varimp_marsres1[,4]>50
						
                        goodip<-which(sel.diffdrthresh==TRUE)
                        
						#print(dim(data_m_fc_withfeats))
						#print(dim(varimp_marsres1))
						#print(head(varimp_marsres1))
						#print(tail(varimp_marsres1))
						
						data_limma_fdrall_withfeats<-cbind(varimp_marsres1[,c(4,6)],data_m_fc_withfeats)
					
						cnames_tab<-colnames(data_m_fc_withfeats)
						cnames_tab<-c("GCV importance","RSS importance",cnames_tab)
						feat_sigfdrthresh[lf]<-length(which(sel.diffdrthresh==TRUE))
						
						colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
					
						#data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
                        #write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
					
					
						goodip<-which(sel.diffdrthresh==TRUE)
						
						#goodip<-sel.diffdrthresh==TRUE
						
						}
						
								if(featselmethod=="pls" | featselmethod=="o1pls" | featselmethod=="o2pls" | featselmethod=="spls" | featselmethod=="o1spls" | featselmethod=="o2spls")
								{
									
                                    #	write.table(data_m_fc,file="data_m_fc.txt",sep="\t",row.names=FALSE)
                                    #write.table(data_m_fc_withfeats,file="data_m_fc_withfeats.txt",sep="\t",row.names=FALSE)
									write.table(classlabels,file="classlabels.txt",sep="\t",row.names=FALSE)
									

									classlabels<-as.data.frame(classlabels)
									
								#	if(analysismode=="classification"){
								#		classlabels<-as.numeric(as.factor(classlabels[,1]))
								#	}
									
									#write.table(classlabels,file="classlabels_numeric.txt",sep="\t",row.names=FALSE)
									
                                    if(is.na(max_comp_sel)==TRUE){
                                            max_comp_sel=pls_ncomp
                                    }
                
							rand_pls_sel<-{} #new("list")
									if(featselmethod=="spls" | featselmethod=="o1spls" | featselmethod=="o2spls"){
								
                                #  print("here")
                                            if(featselmethod=="o1spls"){
                                    
                                                featselmethod="o1pls"
                                            }else{
                                    
                                                if(featselmethod=="o2spls"){
                                                        featselmethod="o2pls"
                                                    }
                                    
                                            }
                                
									if(pairedanalysis==TRUE){
                                        
									classlabels_temp<-cbind(classlabels_sub[,2],classlabels)

#									plsres1<-try(do_plsda(X=data_m_fc,Y=classlabels_sub,oscmode=featselmethod,numcomp=pls_ncomp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=TRUE,analysismode,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=optselect,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,pls.vip.selection=pls.vip.selection),silent=TRUE)

plsres1<-do_plsda(X=data_m_fc,Y=classlabels_sub,oscmode=featselmethod,numcomp=pls_ncomp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=TRUE,analysismode,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=optselect,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,output.device.type=output.device.type,plots.res=plots.res,plots.width=plots.width,plots.height=plots.height,plots.type=plots.type)

if (is(plsres1, "try-error")){
    print(paste("sPLS could not be performed at RSD threshold: ",log2.fold.change.thresh,sep=""))
    #break;
}

opt_comp<-plsres1$opt_comp
									#for(randindex in 1:100)

									if(is.na(pls.permut.count)==FALSE){
                                        if(pls.permut.count>0){
									rand_pls_sel<-lapply(1:pls.permut.count,function(x)
									{

										
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
                                    
                                    plsres1<-do_plsda(X=data_m_fc,Y=classlabels,oscmode=featselmethod,numcomp=pls_ncomp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=TRUE,analysismode,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=optselect,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,pls.vip.selection=pls.vip.selection,output.device.type=output.device.type,plots.res=plots.res,plots.width=plots.width,plots.height=plots.height,plots.type=plots.type)
                                    
                                    opt_comp<-plsres1$opt_comp
                                    
                                    if (is(plsres1, "try-error")){
                                        print(paste("sPLS could not be performed at RSD threshold: ",log2.fold.change.thresh,sep=""))
                                        break;
                                    }
											#for(randindex in 1:100)
									if(is.na(pls.permut.count)==FALSE){
                                        
                                        if(pls.permut.count>0){
									rand_pls_sel<-lapply(1:pls.permut.count,function(x)
                                                                        {


plsresrand<-do_plsda_rand(X=data_m_fc,Y=classlabels[sample(x=seq(1,dim(classlabels)[1]),size=dim(classlabels)[1]),],oscmode=featselmethod,numcomp=opt_comp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=TRUE,analysismode,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=FALSE,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,plotindiv=FALSE)
                                                                                
                                                                                #rand_pls_sel<-cbind(rand_pls_sel,plsresrand$vip_res[,1])
                                                                        	#return(plsresrand$vip_res[,1])		
										if (is(plsresrand, "try-error")){
                                            
									                return(rep(0,dim(data_m_fc)[1]))
                                                                                        }else{
                                                                                        return(plsresrand$vip_res[,1])
                                                                                }
									})	
									}
                                    }
									}
									pls_vip_thresh<-0
									
									if (is(plsres1, "try-error")){
										print(paste("PLS could not be performed at RSD threshold: ",log2.fold.change.thresh,sep=""))
                                        break;
									}else{	
									opt_comp<-plsres1$opt_comp
									}

									}else{
                                        #PLS
											if(pairedanalysis==TRUE){
                                                                        classlabels_temp<-cbind(classlabels_sub[,2],classlabels)
                                                                        plsres1<-do_plsda(X=data_m_fc,Y=classlabels_temp,oscmode=featselmethod,numcomp=pls_ncomp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=FALSE,analysismode=analysismode,vip.thresh=pls_vip_thresh,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=optselect,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,pls.vip.selection=pls.vip.selection,output.device.type=output.device.type,plots.res=plots.res,plots.width=plots.width,plots.height=plots.height,plots.type=plots.type)
                                                                        
                                                                        if (is(plsres1, "try-error")){
                                                                            print(paste("PLS could not be performed at RSD threshold: ",log2.fold.change.thresh,sep=""))
                                                                            break;
                                                                        }else{
                                                                            opt_comp<-plsres1$opt_comp
                                                                        }

                                                                        }else{
											plsres1<-do_plsda(X=data_m_fc,Y=classlabels,oscmode=featselmethod,numcomp=pls_ncomp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=FALSE,analysismode=analysismode,vip.thresh=pls_vip_thresh,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=optselect,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,pls.vip.selection=pls.vip.selection,output.device.type=output.device.type,plots.res=plots.res,plots.width=plots.width,plots.height=plots.height,plots.type=plots.type)
                                            
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
                                            
                                            cl <- makeCluster(getOption("cl.cores", num_nodes))
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
#t1fname<-paste("ranpls",x,".Rda",sep="")
#save(list=ls(),file=t1fname)
print(paste("PLSDA permutation number: ",x,sep=""))

plsresrand<-do_plsda_rand(X=data_m_fc,Y=classlabels[sample(x=seq(1,dim(classlabels)[1]),size=dim(classlabels)[1]),],oscmode=featselmethod,numcomp=opt_comp,kfold=kfold,evalmethod=pred.eval.method,keepX=max_varsel,sparseselect=FALSE,analysismode,sample.col.vec=col_vec,scoreplot_legend=scoreplot_legend,pairedanalysis=pairedanalysis,optselect=FALSE,class_labels_levels_main=class_labels_levels_main,legendlocation=legendlocation,plotindiv=FALSE) #,silent=TRUE)

                                                                                if (is(plsresrand, "try-error")){
                                                                                    
                                                                                    #save(plsresrand,file="plsresrand.Rda")
                                                                                    #print(plsresrand)
                                                                                    return(1)
                                                                                }else{
                                                                                    return(plsresrand$vip_res[,1])
                                                                                    
                                                                                }
                                                                                #rand_pls_sel<-cbind(rand_pls_sel,plsresrand$vip_res[,1])
										
                                                                        })
                                                        
                                                        stopCluster(cl)
                                        }
                                                        #save(rand_pls_sel,file="rand_pls_sel1.Rda")
                                                        
									}
										
									}
											opt_comp<-plsres1$opt_comp
                                    }
										
									if(length(plsres1$bad_variables)>0){
										
										data_m_fc_withfeats<-data_m_fc_withfeats[-c(plsres1$bad_variables),]
											data_m_fc<-data_m_fc[-c(plsres1$bad_variables),]
									}
																			


					
									if(is.na(pls.permut.count)==FALSE){
                                        
                                         if(pls.permut.count>0){
                                        save(rand_pls_sel,file="rand_pls_sel.Rda")
                                        
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
                                    
                                                    #save(rand_pls_sel,file="rand_pls_sel2.Rda")
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
                                                                dev.off()
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
 
 #r2_q2_valid_res<-rbind(plsres1$valid_res$R2,plsres1$valid_res$Q2,plsres1$valid_res$MSEP)
                                    
                                    #  if(length(r2_q2_valid_res)>0){
                                    #rownames(r2_q2_valid_res)<-c("R2","Q2","MSEP")
                                    
                                    
                                    #cnames_vres<-paste("PLScomp",seq(1,plsres1$opt_comp),sep="")
                                    #colnames(r2_q2_valid_res)<-cnames_vres
                                    # }
                                    #   barplot(r2_q2_valid_res[1:2,],beside=TRUE,main="PLS diagnostics",ylab="Variation",col=c("darkgrey","lightgrey"))
                                    # legend("topright",c("R2","Q2"),col=c("darkgrey","lightgrey"),pch=c(20))
									
									write.table(vip_res,file="vip_res.txt",sep="\t",row.names=FALSE)
									
									
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

										
										}else{
											
												cnames_tab<-colnames(data_m_fc_withfeats)
												cnames_tab<-c("VIP","Rank",cnames_tab)
												
						
                        if(max_comp_sel>opt_comp){
                            
                            max_comp_sel<-opt_comp
                        }
														
									
												#pls_vip<-plsres1$vip_res[,c(1:max_comp_sel)]
									
                                    #if(FALSE)
                                    {
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

                                        }
                                    
                                    #vip_res1<-plsres1$vip_res
                                    pls_vip<-vip_res1
											
                                            
                                            #  print(summary(pls_vip))
                                            #print(head(pls_vip_thresh))
                                            #print(summary(rand_pls_sel_fdr))
                                            
                                            #pls
                                                sel.diffdrthresh<-pls_vip>=pls_vip_thresh & rand_pls_sel_fdr<fdrthresh & rand_pls_sel_prob<pvalue.thresh
									
									}
									
									rank_vec<-order(pls_vip,decreasing=TRUE)
										rank_vec2<-seq(1,length(rank_vec))

				ranked_vec<-pls_vip[rank_vec]
				rank_num<-match(pls_vip,ranked_vec)


									
									data_limma_fdrall_withfeats<-cbind(pls_vip,rank_num,data_m_fc_withfeats)
					
					
						feat_sigfdrthresh[lf]<-length(which(sel.diffdrthresh==TRUE)) #length(plsres1$selected_variables) #length(which(sel.diffdrthresh==TRUE))
						
						filename<-paste(featselmethod,"_variable_importance.txt",sep="")

						colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
					
						#data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
                        #write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
					

								goodip<-which(sel.diffdrthresh==TRUE) #plsres1$selected_variables #which(sel.diffdrthresh==TRUE)	
									
								
								}
								
								#stop("Please choose limma, RF, RFcond, or MARS for featselmethod.")
								if(featselmethod=="lmreg" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat"| featselmethod=="logitreg" | featselmethod=="wilcox")
								{
								pvalues<-{}
								
								classlabels_response_mat<-as.data.frame(classlabels_response_mat)
					
					
								
								
										if(featselmethod=="lm1wayanova")
										{
										
											print("Performing one-way ANOVA analysis")
											
											#print(dim(data_m_fc))
											#print(dim(classlabels_response_mat))
											#print(dim(classlabels))
											
											#data_mat_anova<-cbind(t(data_m_fc),classlabels_response_mat)
                                            
                                            
                                            numcores<-round(detectCores()*0.6)
                                            
                                            cl <- makeCluster(getOption("cl.cores", numcores))
                                            
                                            clusterExport(cl,"diffexponewayanova",envir = .GlobalEnv)
                                            
                                             clusterExport(cl,"anova",envir = .GlobalEnv)
                                             
                                             
                                             clusterExport(cl,"TukeyHSD",envir = .GlobalEnv)
                                             
                                              clusterExport(cl,"aov",envir = .GlobalEnv)
                                            

#res1<-apply(data_m_fc,1,function(x){
                                                res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat){
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
                                        write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
											
											
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
                                            
                                            cl <- makeCluster(getOption("cl.cores", numcores))
                                            
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
                                            
                                            #pvalues<-t(pvalues)
                                            # print(dim(pvalues))
                                            #print(dim(data_m_fc_withfeats))
                                            
                                            
                                            
                                            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
                                            
                                            colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
                                            
                                            #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
                                            #  write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
                                            
                                            
                                            data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
                                            
                                        }
                                        
                                        
										if(featselmethod=="lmreg")
										{
										
											
											
											#print(dim(data_m_fc))
                                            #	print(dim(classlabels_response_mat))
                                            #      print(head(classlabels_response_mat))
                                            
                                            if(logistic_reg==TRUE){
                                                print("Performing logistic regression analysis:")
                                                
                                               
                                                classlabels_response_mat[,1]<-as.numeric((classlabels_response_mat[,1]))-1
                                                
                                                fileheader="logitreg"
                                            }else{
                                                print("Performing linear regression analysis:")
                                                fileheader="lmreg"
                                            }

#print(head(classlabels_response_mat))
   
   
   numcores<-round(detectCores()*0.5)
   
   cl <- makeCluster(getOption("cl.cores", numcores))
   
   clusterExport(cl,"diffexplmreg",envir = .GlobalEnv)
  clusterExport(cl,"lm",envir = .GlobalEnv)
  clusterExport(cl,"glm",envir = .GlobalEnv)
   clusterExport(cl,"summary",envir = .GlobalEnv)
   clusterExport(cl,"anova",envir = .GlobalEnv)

											#data_mat_anova<-cbind(t(data_m_fc),classlabels_response_mat)
													res1<-apply(data_m_fc,1,function(x,classlabels_response_mat,logistic_reg){
												
                                                
														#print(x)
														xvec<-x
														#print(length(xvec))
                                                        #print(dim(classlabels_response_mat))
                                                        
                                                        if(dim(classlabels_response_mat)[2]>1){
                                                            
                                                            for(cnum in 2:dim(classlabels_response_mat)[2]){
                                                                
                                                                # classlabels_response_mat[,cnum]<-as.numeric(classlabels_response_mat[,cnum])
                                                            
                                                            }
                                                        }
                                                        
														data_mat_anova<-cbind(xvec,classlabels_response_mat)
														
														cnames<-colnames(data_mat_anova)
														cnames[1]<-"Response"
														
														colnames(data_mat_anova)<-cnames
														
														#print(data_mat_anova)
														
														anova_res<-diffexplmreg(dataA=data_mat_anova,logistic_reg)
														
														return(anova_res)
													},classlabels_response_mat,logistic_reg)
														
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


                                        write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
                                        }
                                        
                                        filename<-paste(fileheader,"_pval_coef_stderr.txt",sep="")
                                        


										data_allinf_withfeats<-cbind(all_inf_mat,data_m_fc_withfeats)
                                        
                                        write.table(data_allinf_withfeats, file=filename,sep="\t",row.names=FALSE)

                                        #write.table(all_inf_mat,file="all_inf_mat.txt",sep="\t",row.names=FALSE)
                                        
                                        #  write.table(classlabels_response_mat,file="classlabels_response_mat.txt",sep="\t",row.names=FALSE)

										cnames_tab<-colnames(data_m_fc_withfeats)
										
                                        #cnames_tab<-c(paste("P.value_var",1:dim(classlabels_response_mat)[2],sep=""),
                                        #paste("Estimate_var",1:dim(classlabels_response_mat)[2],sep=""), paste("StdError_var",1:dim(classlabels_response_mat)[2],sep=""),
                                        #paste("statistic_var",1:dim(classlabels_response_mat)[2],sep=""),cnames_tab)
										
                                        cnames_tab<-c(cnames_1,cnames_tab)
                                        
                                           class_column_names<-colnames(classlabels_response_mat)
                                           # print(class_column_names)
                                        # cnames_tab<-c(paste("P.value_var",1:dim(classlabels_response_mat)[2],sep=""),
                                        #paste("Estimate_var",1:dim(classlabels_response_mat)[2],sep=""), paste("StdError_var",1:dim(classlabels_response_mat)[2],sep=""),
                                        #paste("statistic_var",1:dim(classlabels_response_mat)[2],sep=""),cnames_tab)
                                        
                                        
                                        #         cnames_tab<-c(paste("P.value_",class_column_names,sep=""),
                                        #paste("Estimate_",class_column_names,sep=""), #paste("StdError_",class_column_names,sep=""),
                                        #paste("statistic_",class_column_names,sep=""),cnames_tab)
                                        
                                        #print(head(cnames_tab))

										colnames(data_allinf_withfeats)<-as.character(cnames_tab)
										
										
										write.table(data_allinf_withfeats, file=filename,sep="\t",row.names=FALSE)
											
											
											
										}
											
											
											
													if(featselmethod=="lm2wayanova")
													{
										
														print("Performing two-way ANOVA analysis with Tukey post hoc comparisons")
											
                                            #print(dim(data_m_fc))
                                            #			print(dim(classlabels_response_mat))
		
        
        numcores<-round(detectCores()*0.5)
        
        cl <- makeCluster(getOption("cl.cores", numcores))
        
        
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
													
														#save(data_mat_anova,file="data_mat_anova.Rda")
	
														#diffexplmtwowayanova
														anova_res<-diffexplmtwowayanova(dataA=data_mat_anova)
														
														
														return(anova_res)
														},classlabels_response_mat)
														
														
                                                        stopCluster(cl)
                                                        #	print("done")

#save(res1,file="res1.Rda")
            
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
                                                                        
                                                                         write.table(twoanova_res,file="twoanova_with_posthoc_pvalues.txt",sep="\t",row.names=FALSE)
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
															dev.off()
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
                                        
                                        # save(data_m_fc_withfeats,file="data_m_fc_withfeats.Rda")
                                        data_limma_fdrall_withfeats<-cbind(pvalues1,fdr_adjust_pvalue1,pvalues2,fdr_adjust_pvalue2,pvalues3,fdr_adjust_pvalue3,posthoc_pval_mat,data_m_fc_withfeats)

fdr_adjust_pvalue<-cbind(fdr_adjust_pvalue1,fdr_adjust_pvalue2,fdr_adjust_pvalue3)
fdr_adjust_pvalue<-apply(fdr_adjust_pvalue,1,function(x){min(x,na.rm=TRUE)})


                                                                                colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)

                                                                                #data_limma_fdrall_withfeats<-data_limma_fdrall_withfeats[order(fdr_adjust_fpvalue),]
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
                                        Y1=cbind(classlabels_orig[,1],classlabels_response_mat[,1])
                                        Y1<-as.data.frame(Y1)



if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/HCA_Factor1sigfeats.png"
    
    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
}

get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X1,Y=Y1,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                                        sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type="two-way",newdevice=FALSE,input.type="intensity",mainlab="Factor1")
                                        
                                        if(output.device.type!="pdf"){
                                            
                                            dev.off()
                                        }
                                        }else{
                                            
                                            print("No significant features for Factor 1.")
                                        }

if(length(which(fdr_adjust_pvalue2<fdrthresh))>0){
                                        X2=data_m_fc_withfeats[which(fdr_adjust_pvalue2<fdrthresh),]
                                        
                                        
                                        Y2=cbind(classlabels_orig[,1],classlabels_response_mat[,2])
                                        Y2<-as.data.frame(Y2)
                                        
                                        if(output.device.type!="pdf"){
                                            
                                            temp_filename_1<-"Figures/HCA_Factor2sigfeats.png"
                                            
                                            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                        }
                                        
                                        
                                        get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X2,Y=Y2,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                                        sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type="two-way",newdevice=FALSE,input.type="intensity",mainlab="Factor2")
                                        

                                        if(output.device.type!="pdf"){
                                            
                                            dev.off()
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
                                            
                                            temp_filename_1<-"Figures/HCA_Factor1xFactor2sigfeats.png"
                                            
                                            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                        }
                                        
                                        get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X3,Y=Y3,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
                                        sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type="two-way",newdevice=FALSE,input.type="intensity",mainlab="Factor1 x Factor2")
                                        

                                        if(output.device.type!="pdf"){
                                            
                                            dev.off()
                                        }


}else{
    
    print("No significant features for the interaction.")
}


                                        #data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,posthoc_pval_mat,data_m_fc_withfeats)

											
                                            #data_limma_fdrall_withfeats<-cbind(pvalues1,fdr_adjust_pvalue1,data_m_fc_withfeats)
												 
                                                 #	 cnames_tab<-colnames(data_m_fc_withfeats)
                                                 #cnames_tab<-c("P.value.Factor1","adjusted.P.value.Factor1",cnames_tab)
                                                 #colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
									
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
														
                                                        
                                                        
                                                        
                                                        
                                                        numcores<-round(detectCores()*0.5)
                                                        
                                                        cl <- makeCluster(getOption("cl.cores", numcores))
                                                        
                                                        clusterExport(cl,"diffexplmonewayanovarepeat",envir = .GlobalEnv)
                                                        clusterEvalQ(cl,library(nlme))
                                                        clusterEvalQ(cl,library(multcomp))
                                                        clusterExport(cl,"lme",envir = .GlobalEnv)
                                                        clusterExport(cl,"interaction",envir = .GlobalEnv)
                                                        clusterExport(cl,"anova",envir = .GlobalEnv)
                                                        #clusterExport(cl,"classlabels_response_mat",envir = .GlobalEnv)
                                                        #clusterExport(cl,"subject_inf",envir = .GlobalEnv)
                                                        
                                                        #res1<-apply(data_m_fc,1,function(x){
                                                        
                                                        res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat,subject_inf){
                                                            
                                                            #res1<-apply(data_m_fc,1,function(x){
												
														xvec<-x
														
														colnames(classlabels_response_mat)<-paste("Factor",seq(1,dim(classlabels_response_mat)[2]),sep="")
														
														data_mat_anova<-cbind(xvec,classlabels_response_mat)
														
														cnames<-colnames(data_mat_anova)
														cnames[1]<-"Response"
														
														colnames(data_mat_anova)<-cnames
														
														#print("subject inf is")
														#print(subject_inf)
														
														
														
														anova_res<-diffexplmonewayanovarepeat(dataA=data_mat_anova,subject_inf=subject_inf)
														
                                                        

														
														return(anova_res)
                                                        
                                                        
														},classlabels_response_mat,subject_inf)
                                                        
                                                        stopCluster(cl)
                                                        
														main_pval_mat<-{}
														pvalues<-{}

#	save(res1,file="lmres1.Rda")
														
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
															dev.off()
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
                                        write.table(data_limma_fdrall_withfeats, file=filename,sep="\t",row.names=FALSE)
											
											data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)

		
									}
													
																if(featselmethod=="lm2wayanovarepeat"){
														
														print("Performing two-way ANOVA with repeated measurements analysis using nlme::lme()")
														
                                                        
                                                        numcores<-round(detectCores()*0.5)
                                                        
                                                        cl <- makeCluster(getOption("cl.cores", numcores))
                                                       
                                                        clusterExport(cl,"diffexplmtwowayanovarepeat",envir = .GlobalEnv)
                                                        clusterEvalQ(cl,library(nlme))
                                                        clusterEvalQ(cl,library(multcomp))
                                                        clusterExport(cl,"lme",envir = .GlobalEnv)
                                                        clusterExport(cl,"interaction",envir = .GlobalEnv)
                                                        clusterExport(cl,"anova",envir = .GlobalEnv)
                                                        #clusterExport(cl,"classlabels_response_mat",envir = .GlobalEnv)
                                                        #clusterExport(cl,"subject_inf",envir = .GlobalEnv)
                                                        
                                                        #res1<-apply(data_m_fc,1,function(x){

#	print(dim(data_m_fc))
#	print(dim(classlabels_response_mat))
                                                        
                                                        res1<-parApply(cl,data_m_fc,1,function(x,classlabels_response_mat,subject_inf){
                                                            
                                                            #  res1<-apply(data_m_fc,1,function(x){
												
                                                # save(classlabels_response_mat,file="classlabels_response_mat.Rda")
                                                        
                                                        #       save(subject_inf,file="subject_inf.Rda")

                                                        xvec<-x
														
                                                        #save(xvec,file="xvec.Rda")
														colnames(classlabels_response_mat)<-paste("Factor",seq(1,dim(classlabels_response_mat)[2]),sep="")
														
														data_mat_anova<-cbind(xvec,classlabels_response_mat)
														
														cnames<-colnames(data_mat_anova)
														cnames[1]<-"Response"
														
														colnames(data_mat_anova)<-cnames
														
														#print(subject_inf)
														#print(dim(data_mat_anova))
														
														
														subject_inf<-as.data.frame(subject_inf)
														#print(dim(subject_inf))
														
														anova_res<-diffexplmtwowayanovarepeat(dataA=data_mat_anova,subject_inf=subject_inf[,1])
														
														return(anova_res)
														},classlabels_response_mat,subject_inf)
                                                        
                                                        
														main_pval_mat<-{}
														
                                                        stopCluster(cl)
																		posthoc_pval_mat<-{}

#print(head(res1))
#	print("here")

															pvalues<-{}


                                                                bad_lm1feats<-{}

save(res1,file="res1.Rda")
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
																
                                                                write.table(twoanovarepeat_res,file="2wayanovarepeat_with_posthoc_pvalues.txt",sep="\t",row.names=FALSE)
                                                                

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
            dev.off()
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
Y1=cbind(classlabels_orig[,1],classlabels_response_mat[,1])
Y1<-as.data.frame(Y1)



if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/HCA_Factor1sigfeats.png"
    
    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
}


get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X1,Y=Y1,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type="two-way",newdevice=FALSE,input.type="intensity",mainlab="Factor 1")

if(output.device.type!="pdf"){
    
    dev.off()
}
}else{
    print("No significant features for Factor 1.")
    
}

if(length(which(fdr_adjust_pvalue2<fdrthresh))>0){

X2=data_m_fc_withfeats[which(fdr_adjust_pvalue2<fdrthresh),]
Y2=cbind(classlabels_orig[,1],classlabels_response_mat[,2])
Y2<-as.data.frame(Y2)


if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/HCA_Factor2sigfeats.png"
    
    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
}


get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X2,Y=Y2,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type="two-way",newdevice=FALSE,input.type="intensity",mainlab="Factor 2")

if(output.device.type!="pdf"){
    
    dev.off()
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
    
    temp_filename_1<-"Figures/HCA_Factor1xFactor2sigfeats.png"
    
    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
}

get_hca(feature_table_file=NA,parentoutput_dir=output_dir,class_labels_file=NA,X=X3,Y=Y3,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=FALSE,analysismode="classification",
sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3, hca_type="two-way",newdevice=FALSE,input.type="intensity",mainlab="Factor 1 x Factor 2")


if(output.device.type!="pdf"){
    
    dev.off()
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
					  			| featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="logitreg" | featselmethod=="limma2wayrepeat" | featselmethod=="wilcox")
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

				
					Targetvar<-classlabels[,1]
                                        dataA<-cbind(Targetvar,t(data_m_fc))
                                        dataA<-as.data.frame(dataA)
                                        #df.summary <- dataA %>% group_by(Targetvar) %>%  summarize_each(funs(mean))

                                        df.summary <- dataA %>% group_by(Targetvar) %>%  summarize_each(funs(mean))
                                        df2<-as.data.frame(df.summary[,-c(1)])

                                        save(df2,file="df2.Rda")
                                        save(dataA,file="dataA.Rda")
                                        save(Targetvar,file="Targetvar.Rda")
                                        
                                        if(log2transform==TRUE || input.intensity.scale=="log2"){

                                        foldchangeres<-apply(df2,2,function(x){res<-{};for(i in 1:length(x)){res<-c(res,(x[i]-x[-i]));};tempres<-abs(res);res_ind<-which(tempres==max(tempres,na.rm=TRUE));return(res[res_ind[1]]);})
                                        
                                        print("Using log2 fold change threshold of")
                                        print(foldchangethresh)


                                        }else{
                                            
                                           

#foldchangeres<-apply(df2[,1:5],2,function(x){res<-{};x<-abs(x);for(i in 1:length(x)){res<-c(res,(x[i]/(x[-i]+0.0001)));};tempres<-abs(log2(res));print(tempres);res_ind<-which(tempres==max(tempres,na.rm=TRUE));print(res);print(res_ind);return(res[res_ind[1]]);})

#



                                             #raw intensities
                                            if(znormtransform==FALSE)
                                            {
                                                foldchangeres<-apply(log2(df2+1),2,function(x){res<-{};for(i in 1:length(x)){res<-c(res,(x[i]-x[-i]));};tempres<-abs(res);res_ind<-which(tempres==max(tempres,na.rm=TRUE));return(res[res_ind[1]]);})
                                                
                                    

                                                #foldchangeres<-log2(foldchangeres+0.0001)
                                                
                                                foldchangethresh=foldchangethresh
                                                print("Using raw fold change threshold of")
                                                print(foldchangethresh)

                                            }else{
                                                
                                                # foldchangeres<-apply(df2,2,function(x){res<-{};x<-abs(x);for(i in 1:length(x)){res<-c(res,(x[i]/(x[-i]+0.0001)));};tempres<-abs(res);res_ind<-which(tempres==max(tempres,na.rm=TRUE));return(res[res_ind[1]]);})
                                                
                                                
                                                foldchangeres<-apply(df2,2,function(x){res<-{};for(i in 1:length(x)){res<-c(res,(x[i]-(x[-i])));};tempres<-abs(res);res_ind<-which(tempres==max(tempres,na.rm=TRUE));return(res[res_ind[1]]);})
                                                
                                               print(summary(foldchangeres))
                                               
                                               
                                                #foldchangethresh=2^foldchangethresh
                                                print("Using Z-score change threshold of")
                                                print(foldchangethresh)
                                                
                                            }
                        # save(foldchangeres,file="foldcha")

						#foldchangeres<-apply(df2,2,function(x){res<-{};x<-log2(x+1);for(i in 1:length(x)){res<-c(res,(x[i]-x[-i]));};tempres<-abs(res);res_ind<-which(tempres==max(tempres,na.rm=TRUE));return(res[res_ind[1]]);})

						
                                        }
 

	                                       

					if(length(class_labels_levels)==2){

							
							zvec=foldchangeres
					}else{
	
							zvec=NA
	
							if(featselmethod=="lmreg" && analysismode=="regression"){
                                
                            cnames_matrix<-colnames(data_limma_fdrall_withfeats)
                            cnames_colindex<-grep("Estimate_",cnames_matrix)
                            
                            print(cnames_matrix)
                            print(cnames_colindex)
                            
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

					roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}
					d4<-as.data.frame(data_limma_fdrall_withfeats)
                    
                    max_mz_val<-roundUpNice(max(d4$mz))
                    max_time_val<-roundUpNice(max(d4$time))
                    
					x1increment=round_any(max_mz_val/10,10,f=floor)
					x2increment=round_any(max_time_val/10,10,f=floor)

					 if(featselmethod=="lmreg" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat"
					  			| featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="logitreg" | featselmethod=="limma2wayrepeat" | featselmethod=="wilcox")
					  {



# print("Plotting manhattan plots")
					  
                      sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh


					  			 goodip<-which(sel.diffdrthresh==TRUE)
                       

							classlabels<-as.data.frame(classlabels)
						
						
							   
							   logp<-(-1)*log((d4[,1]+(10^-20)),10)
							   ythresh<-min(logp[goodip],na.rm=TRUE)

			maintext1="Type 1 manhattan plot (-logp vs mz) \n m/z features above the dashed horizontal line meet the selection criteria"
			maintext2="Type 2 manhattan plot (-logp vs time) \n m/z features above the dashed horizontal line meet the selection criteria"

#  print("here1")
#       print(zvec)

			if(is.na(zvec[1])==FALSE){
				 maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
							maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
			}

    yvec_val=logp
    ylabel="(-)log10p"
    yincrement=1
    y2thresh=1.30103 #(-1)*log10(0.05)
    #save(list=c("d4","yvec_val","ythresh","zvec","x1increment","yincrement","maintext1","x2increment","maintext2","ylabel","y2thresh"),file="manhattanplot_objects.Rda")
    
    
    if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
        
        png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
    }
    
    #save(list=ls(),file="m1.Rda")
    try(get_manhattanplots(xvec=d4$mz,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab=ylabel,xincrement=x1increment,yincrement=yincrement,maintext=maintext1,col_seq=c("black"),y2thresh=y2thresh,colorvec=manhattanplot.col.opt),silent=TRUE)

#get_manhattanplots(xvec=d4$mz,yvec=yvec_val,ythresh=ythresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab=ylabel,xincrement=x1increment,yincrement=yincrement,maintext=maintext1,col_seq=c("black"),y2thresh=y2thresh)
    
    if(output.device.type!="pdf"){
    
        dev.off()
    }

if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
    
    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
}

try(get_manhattanplots(xvec=d4$time,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="Retention time",ylab="-log10p",xincrement=x2increment,yincrement=1,maintext=maintext2,col_seq=c("black"),y2thresh=y2thresh,colorvec=manhattanplot.col.opt),silent=TRUE)


if(output.device.type!="pdf"){
    
    dev.off()
}






				}else{
							  
							if(featselmethod=="pls" | featselmethod=="o1pls"){

				maintext1="Type 1 manhattan plot (VIP vs mz) \n m/z features above the dashed horizontal line meet the selection criteria"
			maintext2="Type 2 manhattan plot (VIP vs time) \n m/z features above the dashed horizontal line meet the selection criteria"

			if(is.na(zvec[1])==FALSE){
				 maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
							maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
			}


            yvec_val<-data_limma_fdrall_withfeats[,1]
            ythresh=pls_vip_thresh

            ylabel="VIP"
            yincrement=0.5
            y2thresh=NA
            #  save(list=c("d4","yvec_val","ythresh","zvec","x1increment","yincrement","maintext1","x2increment","maintext2","ylabel","y2thresh"),file="manhattanplot_objects.Rda")
            


if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
    
    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
}

try(get_manhattanplots(xvec=d4$mz,yvec=data_limma_fdrall_withfeats[,1],ythresh=pls_vip_thresh,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="VIP",xincrement=x1increment,yincrement=0.5,maintext=maintext1,col_seq=c("black"),colorvec=manhattanplot.col.opt),silent=TRUE)


if(output.device.type!="pdf"){
    
    dev.off()
}

                            
                            
                            if(output.device.type!="pdf"){
                                
                                temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                                
                                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                            }

                            
                            
                            try(get_manhattanplots(xvec=d4$time,yvec=data_limma_fdrall_withfeats[,1],ythresh=pls_vip_thresh,up_or_down=zvec,xlab="Retention time",ylab="VIP",xincrement=x2increment,yincrement=0.5,maintext=maintext2,col_seq=c("black"),colorvec=manhattanplot.col.opt),silent=TRUE)
							
                            if(output.device.type!="pdf"){
                                
                                dev.off()
                            }
                            
                            }
							else{
										if(featselmethod=="spls" | featselmethod=="o1spls"){
											

maintext1="Type 1 manhattan plot (|loading| vs mz) \n m/z features with non-zero loadings meet the selection criteria"
			maintext2="Type 2 manhattan plot (|loading| vs time) \n m/z features with non-zero loadings meet the selection criteria"
	if(is.na(zvec[1])==FALSE){
				 maintext1=paste(maintext1,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
							maintext2=paste(maintext2,"\n",manhattanplot.col.opt[2],": lower in class ",class_labels_levels_main[1]," & ",manhattanplot.col.opt[1],": higher in class ",class_labels_levels_main[1],sep="")
			}					
    yvec_val<-data_limma_fdrall_withfeats[,1]
    ythresh=0
ylabel="Loading (absolute)"
yincrement=0.1
y2thresh=NA
#save(list=c("d4","yvec_val","ythresh","zvec","x1increment","yincrement","maintext1","x2increment","maintext2","ylabel","y2thresh"),file="manhattanplot_objects.Rda")

if(output.device.type!="pdf"){
    
    temp_filename_1<-"Figures/ManhattanPlot_Type1.png"
    
    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
}

						try(get_manhattanplots(xvec=d4$mz,yvec=data_limma_fdrall_withfeats[,1],ythresh=0,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="Loading (absolute)",xincrement=x1increment,yincrement=0.1,maintext=maintext1,col_seq=c("black"),colorvec=manhattanplot.col.opt),silent=TRUE)
							
                            if(output.device.type!="pdf"){
                                
                                dev.off()
                            }
                            
                            if(output.device.type!="pdf"){
                                
                                temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                                
                                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                            }
                            
                            
                            try(get_manhattanplots(xvec=d4$time,yvec=data_limma_fdrall_withfeats[,1],ythresh=0,up_or_down=zvec,xlab="Retention time",ylab="Loading (absolute)",xincrement=x2increment,yincrement=0.1,maintext=maintext2,col_seq=c("black"),colorvec=manhattanplot.col.opt),silent=TRUE)
										
                                        
                                        if(output.device.type!="pdf"){
                                            
                                            dev.off()
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
                    
                    #print(subdata)
                    # print("good ip is")
                    #print(goodip)
                    #subdata<-apply(subdata,2,function(x){naind<-which(is.na(x)==TRUE);if(length(naind)>0){ x[naind]<-data_minval};return(x)})
#					svm_model<-try(svm(x=subdata,y=classlabels,type="C",cross=kfold),silent=TRUE)
					

#print("svm_cv")
            #  print(dim(subdata))
            #       print(length(classlabels))
            #svm_model<-svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95)
                    svm_model<-try(svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
					
					classlabels<-as.data.frame(classlabels)
                    
                    #print("subdata")
                    #print(subdata[1:3,])
                    #print("classlabels")
                    #print(classlabels)
                    #print(dim(classlabels))




					if(is(svm_model,"try-error")){
						print("SVM could not be performed. Please try lowering the kfold or set kfold=total number of samples for Leave-one-out CV. Skipping to the next step.")
						termA<-(-1)
						pred_acc<-termA
                        permut_acc<-(-1)
					}else{
					#pred_acc<-svm_model$tot.accuracy
					
					#if(pred.eval.method=="AUC"){
					#pred_acc<-multiclass.roc(classlabels,as.numeric(svm_model$fitted))
					#pred_acc_orig<-pred_acc$auc[1]
					#pred_acc<-pred_acc_orig*100
					
					#}
						
						pred_acc<-svm_model$avg_acc
					
                            print("Accuracy is:")
                            print(pred_acc)
                            
                            print("Calculating permuted CV accuracy")
                            
                            permut_acc<-{}
                            #permut_acc<-lapply(1:100,function(j){
                            numcores<-round(detectCores()*0.5)
                            cl <- makeCluster(getOption("cl.cores", numcores))
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
					
					#if(featselmethod=="limma")
					
					if(featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma2wayrepeat" | featselmethod=="lmreg" | featselmethod=="logitreg"
| featselmethod=="lm2wayanova" | featselmethod=="lm1wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="wilcox")
					{
                        if(fdrmethod=="none"){
                            exp_fp<-(dim(data_m_fc)[1]*fdrthresh)+1
                        }else{
						exp_fp<-(feat_sigfdrthresh[lf]*fdrthresh)+1
                        }
                    }
					
                    
                    #  print(termA)
                    #print(permu_acc)
					termB<-(dim(parent_data_m)[1]*dim(parent_data_m)[1])/(dim(data_m_fc)[1]*dim(data_m_fc)[1]*100)
					
					
					res_score<-(100*(termA-permut_acc))-(feat_weight*termB*exp_fp)
					res_score<-round(res_score,2)
                    #	print("res score is")
                    #print(res_score)
                    #print(best_cv_res)
						if(lf==0)
						{
						best_logfc_ind<-lf
						
						best_feats<-goodip
						best_cv_res<-res_score
						best_acc<-pred_acc
						best_limma_res<-data_limma_fdrall_withfeats[goodip,] #[sel.diffdrthresh==TRUE,]
						
						}else{
						
					if(res_score>best_cv_res){
    #metabpcaresultnotransform10pcsallmetabs<-metabpcaresult
						
						
						
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
                                fname4<-paste(featselmethod,"results_allfeatures.txt",sep="")
                            }
                            
                            #allmetabs_res<-data_limma_fdrall_withfeats_2
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
						try(dev.off(),silent=FALSE)
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
					
				

						#sampleclass<-{}
						#patientcolors<-{}
						
						t1<-table(classlabels)
						#	print(t1)
				#	col <- rep(col_vec[1:length(t1)], t1)
				#col_all=topo.colors(256)

					


							
							#print("here")
							#print(sampleclass)
							#print(patientcolors)
					
					#patientcolors <- unlist(lapply(sampleclass, color.map))
					if(length(goodip)>2){
					
					goodfeats<-as.data.frame(data_m_fc_withfeats[goodip,]) #[sel.diffdrthresh==TRUE,])

					goodfeats<-unique(goodfeats)
					
					#print("here")
					#print(goodfeats[1:10,1:5])
					#rownames(goodfeats)<-as.character(goodfeats[,1])
					rnames_goodfeats<-as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
					
					if(length(which(duplicated(rnames_goodfeats)==TRUE))>0){
						goodfeats<-goodfeats[-which(duplicated(rnames_goodfeats)==TRUE),]
					rnames_goodfeats<-rnames_goodfeats[-which(duplicated(rnames_goodfeats)==TRUE)]
					}


					rownames(goodfeats)<-as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
					
					#print("here 2")
					#assign colors by sample class
					#print(dim(goodfeats))
					
						data_m<-as.matrix(goodfeats[,-c(1:2)])
						#data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE);if(length(naind)>0){ x[naind]<-median(x,na.rm=TRUE)};return(x)})
						#data_m<-t(data_m)
						
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
						
                        heatmap_mainlabel="2-way HCA using all significant features"
                        
                        
                        if(output.device.type!="pdf"){
                            
                            temp_filename_1<-"Figures/HCA_All_sigfeats.png"
                            
                            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                        }
                        
                      

						
                        #save(list=c("data_m","hr","hc","heatmap_cols","heatmap_mainlabel","patientcolors"),file="hca_objects.Rda")
						if(znormtransform==FALSE){
							
						h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main=heatmap_mainlabel, ColSideColors=patientcolors)
						}else{
							h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main=heatmap_mainlabel, ColSideColors=patientcolors)
							}
						#dev.off()
						
						
                        
                        if(output.device.type!="pdf"){
                            
                            dev.off()
                        }

						
						
						mycl_samples <- cutree(hc, h=max(hc$height)/2)
						mycl_metabs <- cutree(hr, h=max(hr$height)/2)
						
						ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],goodfeats[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])

						cnames1<-colnames(ord_data)
						cnames1[1]<-"mz_cluster_label"
						colnames(ord_data)<-cnames1
						fname1<-paste("Clustering_based_sorted_intensity_data.txt",sep="")
						write.table(ord_data,file=fname1,sep="\t",row.names=FALSE)

						fname2<-paste("Sample_clusterlabels.txt",sep="")
						
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
                            
                            dev.off()
                        }

						#temp4<-temp4[,-c(2)]
						#colnames(temp4)<-c("SampleID","Class","HCA Cluster #")
						write.table(temp4,file=fname2,sep="\t",row.names=FALSE)
					

					
						fname3<-paste("Metabolite_clusterlabels.txt",sep="")
						
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
						
                        #save(list=c("data_m","hr","hc","heatmap_cols","heatmap_mainlabel","patientcolors"),file="hca_objects.Rda")
                        
						#tiff(heatmap_file,width=plots.width,height=plots.height,res=plots.res, compression="lzw")
                        
                        
                        if(output.device.type!="pdf"){
                            
                            temp_filename_1<-"Figures/HCA_all_sigfeats.png"
                            
                            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                        }
                        
                        						if(znormtransform==FALSE){
						h73<-heatmap.2(as.matrix(data_m), col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main="2-way HCA using all significant features", ColSideColors=patientcolors)
						}else{
							h73<-heatmap.2(as.matrix(data_m), col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main="2-way HCA using all significant features", ColSideColors=patientcolors)
							}
						#dev.off()
					
                    
                    if(output.device.type!="pdf"){
                        
                        dev.off()
                    }

						names(h73)
						
						
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
                    
                    d4<-read.table("lmreg_pval_coef_stderr.txt",sep="\t",header=TRUE,quote = "")
                    
                    }
                    
                    # print(dim(d4))
                    #print(head(d4))
                    
					if(length(goodip)>=1){
						
						
					#	goodip<-sel.diffdrthresh==TRUE
						
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
                    | featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="logitreg" | featselmethod=="limma2wayrepeat" | featselmethod=="wilcox")
                    {
                        
                        
                        
                        #print("Plotting manhattan plots")
                        
                        sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
                        
                       
                        
                        goodip<-which(sel.diffdrthresh==TRUE)
                        
                        print("good ip is ")
                        
                        print(goodip)
                        
                        print(summary(fdr_adjust_pvalue))
                        print(final.pvalues)


                        
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
                                
                                dev.off()
                            }
                            
                            
                            if(output.device.type!="pdf"){
                                
                                temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                                
                                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                            }
                            
                            
                        
                            try(get_manhattanplots(xvec=d4$time,yvec=logp,ythresh=ythresh,up_or_down=zvec,xlab="Retention time",ylab="-logP",xincrement=x2increment,yincrement=1,maintext=maintext2,col_seq=c("black"),y2thresh=1.30103,colorvec=manhattanplot.col.opt),silent=TRUE)
                            
                            
                            if(output.device.type!="pdf"){
                                
                                dev.off()
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
                                
                                dev.off()
                            }
                            
                            
                            if(output.device.type!="pdf"){
                                
                                temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                                
                                png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                            }
                            
                            
                            
                            try(get_manhattanplots(xvec=d4$time,yvec=data_limma_fdrall_withfeats[,1],ythresh=pls_vip_thresh,up_or_down=zvec,xlab="Retention time",ylab="VIP",xincrement=x2increment,yincrement=0.5,maintext=maintext2,col_seq=c("black"),colorvec=manhattanplot.col.opt),silent=TRUE)
                            
                            
                            
                            if(output.device.type!="pdf"){
                                
                                dev.off()
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
                                    
                                    dev.off()
                                }
                                
                                
                                if(output.device.type!="pdf"){
                                    
                                    temp_filename_1<-"Figures/ManhattanPlot_Type2.png"
                                    
                                    png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                                }
                                
                                
                                
                                try(get_manhattanplots(xvec=d4$time,yvec=data_limma_fdrall_withfeats[,1],ythresh=0,up_or_down=zvec,xlab="Retention time",ylab="Loading",xincrement=x2increment,yincrement=0.1,maintext=maintext2,col_seq=c("black"),colorvec=manhattanplot.col.opt),silent=TRUE)
                                
                                
                                if(output.device.type!="pdf"){
                                    
                                    dev.off()
                                }
                                
                            }
                        }
                        
                    }
                    
                    #write.table(subdata,file="test.txt",sep="\t")
                    #subdata<-na.omit(subdata)
					#write.table(subdata,file="test2.txt",sep="\t")
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
						print("Number of significant metabolites is too small to perform CV.")
						
					}
					
					#print("termA is ")
					#print(termA)
                    
                    # print("dim of goodfeats")
                    goodfeats<-as.data.frame(data_m_fc_withfeats[sel.diffdrthresh==TRUE,])


#print(dim(goodfeats))
                    
                    goodip<-which(sel.diffdrthresh==TRUE)
                    
                    #print(length(goodip))
                    
				if(length(which(sel.diffdrthresh==TRUE))>2){
					 #plot(, as.numeric(classlabels))
				    # points(x, log(x), col = 2)
				     #points(x, new, col = 4)
					
					#res_score<-10000000-termA	
					
					
					res_score<-termA
					
					#if(res_score<best_cv_res){
						
						if(res_score<best_cv_res){
						
						{
						best_logfc_ind<-lf
						
						best_feats<-goodip
						best_cv_res<-res_score
						best_acc<-pred_acc
						best_limma_res<-data_limma_fdrall_withfeats[sel.diffdrthresh==TRUE,]
						}
					}
					res_score_vec[lf]<-res_score
					

					
					
						

					goodfeats<-unique(goodfeats)
					#rownames(goodfeats)<-as.character(goodfeats[,1])
					rownames(goodfeats)<-as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
					
						#assign colors by sample class
					
					
						data_m<-as.matrix(goodfeats[,-c(1:2)])
						#data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE);if(length(naind)>0){ x[naind]<-median(x,na.rm=TRUE)};return(x)})
						#data_m<-t(data_m)
						
						rownames(data_m)<-as.character(paste(goodfeats[,1],goodfeats[,2],sep="_"))
					
					
				    X<-t(data_m)
				   
                   #print("HCA analysis")
                    
               
                pca_comp<-min(dim(X)[1],dim(X)[2])
                
					t1<-seq(1,dim(data_m)[2])
					
					#col_vec=topo.colors(dim(data_m)[2])
				    
				    ## 1) raw data
				    #tiff("PC1PC2sigmetabs_colorbyclasses.tiff",width=plots.width,height=plots.height,res=plots.res)
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
                            
                            temp_filename_1<-"Figures/HCA_all_sigfeats.png"
                            
                            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                        }
                        
                        
                        
                        
						if(znormtransform==FALSE){
						h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main="Using all significant features")
						}else{
							h73<-heatmap.2(data_m, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main="Using all significant features")
							}
                        
                        
                        if(output.device.type!="pdf"){
                            
                            dev.off()
                        }
					#	dev.off()
						
						
						mycl_samples <- cutree(hc, h=max(hc$height)/2)
						mycl_metabs <- cutree(hr, h=max(hr$height)/2)
						
						ord_data<-cbind(mycl_metabs[rev(h73$rowInd)],goodfeats[rev(h73$rowInd),c(1:2)],data_m[rev(h73$rowInd),h73$colInd])

cnames1<-colnames(ord_data)
						cnames1[1]<-"mz_cluster_label"
						colnames(ord_data)<-cnames1
						fname1<-paste("Clustering_based_sorted_intensity_data.txt",sep="")
						write.table(ord_data,file=fname1,sep="\t",row.names=FALSE)

						fname2<-paste("Sample_clusterlabels.txt",sep="")
						
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
                                
                                dev.off()
                            }
							

							
						}
						
                        #	print(head(temp_vec))
						#temp4<-temp4[,-c(2)]
						write.table(temp4,file=fname2,sep="\t",row.names=FALSE)
					

					
						fname3<-paste("Metabolite_clusterlabels.txt",sep="")
						
						mycl_metabs_ord<-mycl_metabs[rev(h73$rowInd)]
						#mz_rt_info<-rownames(mycl_metabs_ord)
						#mycl_metabs_ord<-cbind(mz_rt_info,mycl_metabs_ord)
						#write.table(mycl_metabs_ord,file=fname3,sep="\t",row.names=TRUE)

					
						}
				}

						
						
						
						
					}
					
                    classlabels_orig<-classlabels_orig_parent
                    if(pairedanalysis==TRUE){
                        
                        classlabels_orig<-classlabels_orig[,-c(2)]
                        
                        
                    }else{
                        
                        if(featselmethod=="lmreg" || featselmethod=="logitreg"){
                            classlabels_orig<-classlabels_orig[,c(1:2)]
                            classlabels_orig<-as.data.frame(classlabels_orig)
                        }
                    }
                    
                    # print("dim of classlabels orig xyplots B")
                    #print(dim(classlabels_orig))
                    #print(head(classlabels_orig))
                    
                    #  print("HERE")
                    

                    #print(head(classlabels_orig_parent))
                    classlabels_orig_wgcna<-classlabels_orig
                    
                    
                    
                    if(analysismode=="classification"){
                        
                        classlabels_temp<-classlabels_orig_wgcna #cbind(classlabels_sub[,1],classlabels)
                        
                        
                        
                        
                        #   print("generating DICE plot")
                        
                        #pres<-do_wgcna(X=data_m_fc_withfeats,Y=classlabels_temp,sigfeatsind=goodip)
                        degree_eval_res<-try(degree_eval(X=data_m_fc_withfeats,Y=classlabels_temp,sigfeats=data_m_fc_withfeats[goodip,c(1:2)],sigfeatsind=goodip),silent=TRUE)
                        
                        
                    }
                    
                    
                    
                    if(analysismode=="classification"){
                        
                        
                        if(is(degree_eval_res,"try-error")){
                            
                            
                            degree_rank<-rep(1,dim(data_m_fc_withfeats)[1])
                        }else{
                            
                        if(degree_rank_method=="overall"){
                            
                            
                            
                            degree_rank<-rank((-1)*degree_eval_res$all[,3]) #,decreasing=TRUE)
                            
                        }else{
                            
                            degree_rank_mat<-degree_eval_res$all[,-c(1:4)]
                            degree_rank_list<-new("list")
                            
                            for(i in 1:dim(degree_rank_mat)[2]){
                                
                                degree_rank_list[[i]]<-degree_rank_mat[,i]/max(degree_rank_mat[,i])
                                degree_rank_mat[,i]<-100*(degree_rank_mat[,i]/max(degree_rank_mat[,i]))
                            }
                            
                            #diff_degree_measure<-abs(degree_rank_list[[1]]-degree_rank_list[[2]])
                            
                            diffexpres<-apply(degree_rank_mat,1,function(x){res<-{};for(i in 1:length(x)){res<-c(res,(x[i]-x[-i]));};tempres<-abs(res);res_ind<-which(tempres==max(tempres,na.rm=TRUE));return(abs(res[res_ind[1]]));})
                            
                            
                            
                            degree_rank<-rank((-1)*diffexpres) #diff_degree_measure)
                            
                            # degree_rank2<-order((-1)*diffexpres)
                            
                            #dtemp<-cbind(degree_eval_res$all,degree_rank_mat,diffexpres,degree_rank2,degree_rank)
                            # write.table(dtemp,file="degree_debug.txt",sep="\t",row.names=FALSE)
                            
                            
                        }
                        }
                        
                        
                        
                        # print(data_limma_fdrall_withfeats[1:10,1:10])
                        if(featselmethod=="lmreg" | featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma1way" | featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="limma1wayrepeat" | featselmethod=="limma2wayrepeat" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="wilcox")
                        {
                            diffexp_rank<-rank(data_limma_fdrall_withfeats[,2]) #order(data_limma_fdrall_withfeats[,2],decreasing=FALSE)
                        }else{
                            
                            #goodfeats<-goodfeats[order(goodfeats[,2],decreasing=TRUE),]
                            
                            #sigfeats_rank<-order(data_limma_fdrall_withfeats[,2],decreasing=TRUE)
                            if(featselmethod=="rfesvm"){
                                
                                
                                
                                data_limma_fdrall_withfeats<-cbind(rank_vec,data_limma_fdrall_withfeats)
                                
                                
                            }
                            
                            
                            if(featselmethod=="pamr"){
                                
                                data_limma_fdrall_withfeats<-cbind(rank_vec,data_limma_fdrall_withfeats)
                                
                                
                            }
                            
                             if(featselmethod=="MARS"){
                            
                                    diffexp_rank<-rank((-1)*data_limma_fdrall_withfeats[,2])
                            
                             }else{
                            diffexp_rank<-rank((1)*data_limma_fdrall_withfeats[,2])
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
                            fname4<-paste(featselmethod,"results_allfeatures.txt",sep="")
                        }
                        
                        allmetabs_res<-data_limma_fdrall_withfeats_2
                        
                        write.table(data_limma_fdrall_withfeats_2,file=fname4,sep="\t",row.names=FALSE)
                        
                        
                        
                         
                         
                        if(length(goodip)>1){
                        data_limma_fdrall_withfeats_2<-data_limma_fdrall_withfeats_2[goodip,]
                        
                        data_limma_fdrall_withfeats_2<-as.data.frame(data_limma_fdrall_withfeats_2)
                        
                        
                       
                        
                        if(featselmethod=="lmreg" | featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma1way" | featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="limma1wayrepeat" | featselmethod=="limma2wayrepeat" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="wilcox")
                        {
                            diffexp_rank<-rank(data_limma_fdrall_withfeats_2[,3]) #order(data_limma_fdrall_withfeats[,2],decreasing=FALSE)
                            
                        }else{
                            
                            #goodfeats<-goodfeats[order(goodfeats[,2],decreasing=TRUE),]
                            
                            #sigfeats_rank<-order(data_limma_fdrall_withfeats[,2],decreasing=TRUE)
                            
                            diffexp_rank<-rank((1)*data_limma_fdrall_withfeats_2[,3])
                        }
                        
                        
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
                            
                            diff_degree_measure<-abs(degree_rank_list[[1]]-degree_rank_list[[2]])
                            
                            
                            degree_rank<-rank((-1)*diff_degree_measure)
                        }
                        
                        
                        }
                        
                        
                        # print("Degree")
                        #print(dim(data_limma_fdrall_withfeats_2))
                        
                        #print(length(goodip))
                        # save(data_limma_fdrall_withfeats_2,file="data_limma_fdrall_withfeats_2.Rda")
                        
                        #save(degree_rank_list,file="degree_rank_list.Rda")
                        
                        # print(length(data_limma_fdrall_withfeats_2$degree_rank))
                        
                        #print(length(degree_rank))
                        
                        data_limma_fdrall_withfeats_2$degree_rank<-degree_rank
                        data_limma_fdrall_withfeats_2$diffexp_rank<-diffexp_rank
                        if(logistic_reg==TRUE){
                            
                            fname4<-paste("logitreg","results_significantfeatures.txt",sep="")
                            
                        }else{
                            fname4<-paste(featselmethod,"results_significantfeatures.txt",sep="")
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
                       
			if(FALSE)
            {
                        r1<-RankAggreg(x=rank_list,k=numselect,verbose=FALSE,distance="Kendall")
                        
                        #r1<-list(top.list=order(degree_rank))
                        print("Aggregated rank")
                        
                        print(r1$top.list)
                        
                        save(r1,file="aggregrank.Rda")
            }

#write.table(data_limma_fdrall_withfeats_2,file=fname4,sep="\t",row.names=FALSE)
                        
                        #data_limma_fdrall_withfeats_2<-cbind(maxfoldchange,data_limma_fdrall_withfeats)
                        goodfeats<-as.data.frame(data_limma_fdrall_withfeats_2)
                        
                        }
                    }else{
                        
                        
                        data_limma_fdrall_withfeats_2<-data_limma_fdrall_withfeats
                        fname4<-paste(featselmethod,"_sigfeats.txt",sep="")
                        
                        goodfeats<-data_limma_fdrall_withfeats[goodip,] #[sel.diffdrthresh==TRUE,]
                        # write.table(goodfeats,file=fname4,sep="\t",row.names=FALSE)
                        
                    }
                    
                    #best_feats<-goodip
                    
                    #if(length(which(sel.diffdrthresh==TRUE))>1){
                    #	goodfeats<-goodfeats[,-c(1:2)]
                    #}
                    }
                    
                    
                    
                    if(length(goodip)>1){
                    goodfeats_by_DICErank<-{}
                    
                    if(analysismode=="classification"){
                        if(featselmethod=="lmreg" | featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma1way" | featselmethod=="lmreg" | featselmethod=="logitreg" | featselmethod=="limma1wayrepeat" | featselmethod=="limma2wayrepeat" | featselmethod=="lm1wayanova" | featselmethod=="lm2wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="wilcox")
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
                            
                            dev.off()
                            
                            if(featselmethod=="lmreg"){
                                
                                goodfeats<-goodfeats[order(goodfeats[,1],decreasing=FALSE),]
                                
                            }else{
                                goodfeats<-goodfeats[order(goodfeats[,1],decreasing=TRUE),]
                                
                            }
                        }
                        
                    }
                    
                    fname4<-paste(featselmethod,"_sigfeats_ordered_by_significance.txt",sep="")
                    
						
                    }
                    fname4<-paste(featselmethod,"_sigfeats_ordered_by_significance.txt",sep="")
                    write.table(goodfeats,file=fname4,sep="\t",row.names=FALSE)





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
if(analysismode=="classification" && nrow(goodfeats)>2 && length(goodip)>1)
{
    
    goodfeats_temp<-cbind(goodfeats[,mz_ind],goodfeats[,time_ind],goodfeats[,-c(1:time_ind)])
    
    cnames_temp<-colnames(goodfeats_temp)
    cnames_temp[1]<-"mz"
    cnames_temp[2]<-"time"
    colnames(goodfeats_temp)<-cnames_temp
    
    
    # save(goodfeats_temp,file="goodfeats_by_DiffExprank.Rda")
    #save(classlabels,file="classlabels.Rda")


if(length(class_labels_levels)==2){

print("Generating ROC curve using top features on training set")


try(get_roc(dataA=goodfeats_temp,classlabels=classlabels,classifier=rocclassifier,kname="radial",rocfeatlist=rocfeatlist,rocfeatincrement=rocfeatincrement,mainlabel="Training set ROC curve using top features"),silent=TRUE)

 if(length(goodip)>1){

if(FALSE){
goodfeats_by_DICErank<-data_limma_fdrall_withfeats_2[r1$top.list,]

#save(goodfeats_by_DICErank,file="goodfeats_by_DICErank.Rda")

#goodfeats_by_DICErank<-na.omit(goodfeats_by_DICErank)

print("Generating ROC curve using diff exp and diff centrality ranking")
goodfeats_temp<-cbind(goodfeats_by_DICErank[,mz_ind],goodfeats_by_DICErank[,time_ind],goodfeats_by_DICErank[,-c(1:time_ind)])

cnames_temp<-colnames(goodfeats_temp)
cnames_temp[1]<-"mz"
cnames_temp[2]<-"time"
colnames(goodfeats_temp)<-cnames_temp

get_roc(dataA=goodfeats_temp,classlabels=classlabels,classifier=rocclassifier,kname="radial",rocfeatlist=rocfeatlist,rocfeatincrement=rocfeatincrement,mainlabel="ROC curve using DICE ranking")
 }

	}
}

#print("ROC done")
best_subset<-{}
best_acc<-0

xvec<-{}
yvec<-{}
for(i in 2:max_varsel){
    
    subdata<-t(goodfeats_temp[1:i,-c(1:2)])
    svm_model<-try(svm_cv(v=kfold,x=subdata,y=classlabels,kname=svm_kernel,errortype=pred.eval.method,conflevel=95),silent=TRUE)
    
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
    
    if(output.device.type!="pdf"){
        
        temp_filename_1<-"Figures/kfoldCV_forward_selection.png"
        
        png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
    }
    
    
   
plot(x=xvec,y=yvec,main="k-fold CV classification accuracy based on forward selection of top features",xlab="Feature index",ylab=ylab_text,type="b",col="brown")


if(output.device.type!="pdf"){
    
    dev.off()
}

cv_mat<-cbind(xvec,yvec)
colnames(cv_mat)<-c("Feature Index",ylab_text)

write.table(cv_mat,file="kfold_cv_mat.txt",sep="\t")
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

#print("Stage 2 PCA eval")
        # print(dim(X))
        
        p1<-pcaMethods::pca(X,method="rnipals",center=TRUE,scale="uv",cv="q2",nPcs=pca_comp)
        
        if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/PCAdiagnostics_sigfeats.png"
            
            png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        }
        
        
        
        p2<-plot(p1,col=c("darkgrey","grey"),main="PCA diagnostics after variable selection")
        
        print(p2)
        if(output.device.type!="pdf"){
            
            dev.off()
        }
        #dev.off()
        
        
        
        try(detach("package:pcaMethods",unload=TRUE),silent=TRUE)
        library(mixOmics)
        pca_comp<-min(dim(X)[1],dim(X)[2])
        metabpcaresultlog2allmetabs5pcs<-mixOmics::pca(X,ncomp=pca_comp,center=TRUE,scale=TRUE)
				    #metabpcaresultnotransform10pcsallmetabs<-metabpcaresult
                    result<-metabpcaresultlog2allmetabs5pcs

#        pcavar1<-result$sdev/sum(result$sdev)
        
        s1<-summary(result)
        
        pcavar1<-s1$importance[2,]
        names(pcavar1)<-paste("PC",seq(1,length(pcavar1)),sep="")
        pcavar1<-pcavar1[1:12]*100
        
        if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/PCVariance_sigfeats.png"
            
            #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
        }
        
        
        
        
        mp <- barplot(pcavar1, beside = TRUE,
        col =pca_col_vec,
        legend = colnames(pcavar1), ylim = c(0,100),
        main = "% Variance explained by each PC", font.main = 4,
        cex.names = 0.7)
        text(mp, round(pcavar1,2), labels = round(pcavar1,2), pos = 3)
        
         print(mp)
        
        
        if(output.device.type!="pdf"){
            
            #  dev.off()
        }
        try(detach("package:pcaMethods",unload=TRUE),silent=TRUE)
        #tiff("R2Q2_usingsigfeats.tiff",width=plots.width*1.5,height=plots.height*1.5,res=plots.res, compression="lzw")
        
        r1<-pcavar1
        r1<-round(r1,2)
        
        #col_vec<-col_vec[sample(1:length(col_vec),length(col_vec))]
        
        
        
        col<- col_vec[1:length(t1)]
        cex <- rep(0.4,length(t1))
        pch <- rep(15,length(t1))
        
        ## The first two parameters are the x and y coordinates of the legend on the graph
        ## The third one is the text of the legend
        #legend("bottomleft", l1, col = col,pch = pch, pt.cex = cex, title = "class #", cex=0.8)
        l1<-levels(as.factor(class_labels_levels))
        
        if(scoreplot_legend==TRUE){
            #	print(legend("bottomleft", l1, col = col,pch = pch, pt.cex = cex, title = "class #", cex=0.6))
        }
        #dev.off()
        
        
    }

classlabels_orig<-classlabels_orig_parent


#print("dim of classlabels_temp")
#print(dim(classlabels_temp))

#print("dim of classlabels orig xyplots A")
#print(dim(classlabels_orig))
#print(head(classlabels_orig))
if(pairedanalysis==TRUE){
    
    classlabels_orig<-classlabels_orig[,-c(2)]
    
    
}else{
    
    if(featselmethod=="lmreg" || featselmethod=="logitreg"){
        classlabels_orig<-classlabels_orig[,c(1:2)]
        classlabels_orig<-as.data.frame(classlabels_orig)
    }
}

#print("dim of classlabels orig xyplots B")
#print(dim(classlabels_orig))
#print(head(classlabels_orig))

classlabels_orig_wgcna<-classlabels_orig

#print("dim of classlabels sub")
#print(dim(classlabels_sub))
#print(head(classlabels_sub))


#classlabels_orig<-classlabels_orig[seq(1,dim(classlabels_orig)[1],num_replicates),]


    print("here2")
	goodfeats_temp<-cbind(goodfeats[,mz_ind],goodfeats[,time_ind],goodfeats[,-c(1:time_ind)])
	
    print("dim good feats")
       print(dim(goodfeats_temp))
    
    print(goodfeats_temp[1:3,1:6])
    
 
    
    
    
    if(num_sig_feats>3){
        if(output.device.type!="pdf"){
            
            temp_filename_1<-"Figures/PCAplots_sigfeats.pdf"
            
            #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
            pdf(temp_filename_1)
        }
        
        #  save(list=ls(),file="pcascoresigplots.Rda")
              try(get_pcascoredistplots(X=goodfeats_temp,Y=classlabels_orig,feature_table_file=NA,parentoutput_dir=getwd(),class_labels_file=NA,sample.col.opt=sample.col.opt,plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,col_vec=col_vec,pairedanalysis=pairedanalysis,pca.cex.val=pca.cex.val,legendlocation=legendlocation,pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,filename="significant"),silent=TRUE)
       
       if(output.device.type!="pdf"){
           
           try(dev.off(),silent=TRUE)
       }
       
       #save(list=ls(),file="timeseries.Rda")
       
       if(output.device.type!="pdf"){
           
           temp_filename_1<-"Figures/Lineplots_sigfeats.pdf"
           
           #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
           pdf(temp_filename_1)
       }
       
       #if(pairedanalysis==TRUE)
         {
            
            # save(list=ls(),file="lineplots.Rda")
           
           try(get_lineplots(X=goodfeats_temp,Y=classlabels_orig,feature_table_file=NA,parentoutput_dir=getwd(),class_labels_file=NA,sample.col.opt=sample.col.opt,alphacol=0.3,col_vec=col_vec,pairedanalysis=pairedanalysis,pca.cex.val=pca.cex.val,legendlocation=legendlocation,pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,filename="significant"),silent=TRUE)  #,silent=TRUE)
         }
         
         if(output.device.type!="pdf"){
             
             try(dev.off(),silent=TRUE)
         }
    }
    
    
    
   
	

if(nrow(goodfeats)<1){
    
    print(paste("No significant features found for ",featselmethod,sep=""))
}
#else
{

goodfeats<-goodfeats[,-c(1:time_ind)]


     par_rows=2
     
     max_per_row=2
     
     write.table(goodfeats,file="boxplots_file.txt",sep="\t",row.names=FALSE)
     
     if(output.device.type!="pdf"){
         
         temp_filename_1<-"Figures/Boxplots_sigfeats.pdf"
         
         #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
     
        pdf(temp_filename_1)
     }
     
     
     
           par(mfrow=c(par_rows,max_per_row))
                
                for(m in 1:dim(goodfeats)[1])
                {
                    if(m%%9==0){
                        #dev.off()
                        file_ind<-file_ind+1
                        boxplots_fname<-paste("boxplots_file",file_ind,".tiff",sep="")
                        #tiff(boxplots_fname, width=plots.width,height=plots.height,res=plots.res, compression="lzw")
                        #tiff(boxplots_fname, width=2000,height=3000,res=plots.res, compression="lzw")
                        #pdf(boxplots_fname)
                        par(mfrow=c(par_rows,max_per_row))
                    }
                    x1<-as.vector(t(goodfeats[m,c(1:num_samps_group[[1]])]))
                    x1<-replace_outliers(x1)
                    x2<-as.vector(t(goodfeats[m,c((num_samps_group[[1]]+1):(num_samps_group[[1]]+num_samps_group[[2]]))]))
                    x2<-replace_outliers(x2)
                    round_mzval<-sprintf("%.4f",mzvec[m])
                    round_timeval<-sprintf("%.1f",timevec[m])
                    
                    mzname<-paste("mz_time: ",round_mzval,"_",round_timeval,sep="")
                    if(length(class_labels_levels)>2){
                        
                        if(length(class_labels_levels)==(-3)){
                            x3<-as.vector(t(goodfeats[m,c((num_samps_group[[1]]+num_samps_group[[2]]+1):(num_samps_group[[1]]+num_samps_group[[2]]+num_samps_group[[3]]))]))
                            x3<-replace_outliers(x3)
                            class_label_C<-class_labels_levels[3]
                            boxplot(x1,x2, x3, ylab="Intensity",main=mzname,xaxt="n",col=boxplot.col.opt)
                            axis(side=1,at=seq(1),labels=class_label_A, col=col_vec[1])
                            axis(side=1,at=2,labels=class_label_B,col=col_vec[2])
                            axis(side=1,at=3,labels=class_label_C,col=col_vec[3])
                        }else{
                            if(length(class_labels_levels)==(-4)){
                                
                                t1<-table(sampleclass)
                                num_samps_group[[1]]<-t1[1]
                                num_samps_group[[2]]<-t1[2]
                                num_samps_group[[3]]<-t1[3]
                                num_samps_group[[4]]<-t1[4]
                                x3<-as.vector(t(goodfeats[m,c((num_samps_group[[1]]+num_samps_group[[2]]+1):(num_samps_group[[1]]+num_samps_group[[2]]+num_samps_group[[3]]))]))
                                class_label_C<-class_labels_levels[3]
                                x4<-as.vector(t(goodfeats[m,c((num_samps_group[[1]]+num_samps_group[[2]]+num_samps_group[[3]]+1):(num_samps_group[[1]]+num_samps_group[[2]]+num_samps_group[[3]]+num_samps_group[[4]]))]))
                                class_label_D<-class_labels_levels[4]
                                x4<-replace_outliers(x4)
                                boxplot(x1,x2, x3, x4,ylab="Intensity",main=mzname,xaxt="n",col=boxplot.col.opt)
                                axis(side=1,at=seq(1),labels=class_label_A, col=col_vec[1])
                                axis(side=1,at=2,labels=class_label_B,col=col_vec[2])
                                axis(side=1,at=3,labels=class_label_C,col=col_vec[3])
                                axis(side=1,at=4,labels=class_label_D,col=col_vec[4])
                                
                                
                            }else{
                                if(length(class_labels_levels)>=1){
                                    
                                    t1<-table(sampleclass)
                                    cur_d<-{}
                                    for(c in 1:length(class_labels_levels)){
                                        
                                        
                                        num_samps_group[[1]]<-t1[1]
                                        
                                        cvec<-as.vector(t(goodfeats[m,c(groupwiseindex[[c]])]))
                                        
                                        
                                        cvec<-replace_outliers(cvec)
                                        cur_d<-cbind(cur_d,cvec)
                                    }
                                    #x3<-as.vector(t(goodfeats[m,c((num_samps_group[[1]]+num_samps_group[[2]]+1):(num_samps_group[[1]]+num_samps_group[[2]]+num_samps_group[[3]]))]))
                                    #class_label_C<-class_labels_levels[3]
                                    #x4<-as.vector(t(goodfeats[m,c((num_samps_group[[1]]+num_samps_group[[2]]+num_samps_group[[3]]+1):(num_samps_group[[1]]+num_samps_group[[2]]+num_samps_group[[3]]+num_samps_group[[4]]))]))
                                    #class_label_D<-class_labels_levels[4]
                                    
                                    cur_d<-as.data.frame(cur_d)
                                    
                                    #class_labels_boxplot<-levels(classlabelsA) #paste(seq(1,length(class_labels_levels)),sep="")
                                   
				 
                                    #print(class_labels_boxplot)
                                    colnames(cur_d)<-NULL #paste(seq(1,length(class_labels_levels)),sep="") #as.character(class_labels_levels)
                                    cur_d<-round(cur_d,2)
                                   # print("i")
				   # print(m)
				   # print(dim(cur_d))
                                    #print(cur_d[1:4,])
                                   # print(length(cur_d))
				#	print(cur_d)
                                    
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
					
				    boxplot(cur_d,ylab=ylab_text,main=mzname,xaxt="n",col=boxplot.col.opt)
                                    
                                    for(i in 1:length(class_labels_levels)){
                                        axis(side=1,at=c(i),labels=class_labels_levels[i], col=col_vec[i],cex.axis=0.75)
                                        
                                        #text(, round(pcavar1,2), labels = round(pcavar1,2), pos = 3)
                                    }
                                    
                                }
                                
                            }
                        }
                    }else{
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

                                   # boxplot(cur_d,ylab=ylab_text,main=mzname,xaxt="n")
                        boxplot(x1,x2, ylab=ylab_text,main=mzname,xaxt="n",col=boxplot.col.opt)
                        axis(side=1,at=seq(1),labels=class_label_A, col="red")
                        axis(side=1,at=2,labels=class_label_B,col="green")
                    }
                }
                
                
                
                if(output.device.type!="pdf"){
                    
                    dev.off()
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
                
                #save(goodfeats_temp,file="goodfeats_temp.Rda")
                # save(classlabels_orig,file="classlabels_orig.Rda")
                
                if(output.device.type!="pdf"){
                    
                    temp_filename_1<-"Figures/Barplots_sigfeats.pdf"
                    
                    #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                    pdf(temp_filename_1)
                }
                
                
                try(get_barplots(feature_table_file,class_labels_file,X=goodfeats_temp,Y=classlabels_orig,parentoutput_dir=output_dir,newdevice=FALSE,ylabel=ylab_text,cex.val=0.9,barplot.col.opt=boxplot.col.opt),silent=TRUE)
                
                
                if(output.device.type!="pdf"){
                    
                    dev.off()
                }
                
                if(output.device.type!="pdf"){
                    
                    temp_filename_1<-"Figures/Individual_sample_plots_sigfeats.pdf"
                    
                    #png(temp_filename_1,width=plots.width,height=plots.height,res=plots.res,type=plots.type,units="in")
                    pdf(temp_filename_1)
                }
                
               
              
              #save(list=ls(),file="individualsampleplots.Rda")
           
           try(get_individualsampleplots(feature_table_file,class_labels_file,X=goodfeats_temp,Y=classlabels_orig,parentoutput_dir=output_dir,newdevice=FALSE,ylabel=ylab_text,cex.val=0.9,sample.col.opt=sample.col.opt),silent=TRUE)
                #,silent=TRUE)
                
                
               
               if(output.device.type!="pdf"){
                    
                       dev.off()
                    }
                #  print("Classlabels orig")
                # print(classlabels_orig)
                 
                 
                
                 
                 
                
                
	 	if(globalclustering==TRUE){
            
             if(output.device.type!="pdf"){
                     
                     temp_filename_1<-"Figures/Globalclustering.png"
                     
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
                 write.table(t1,file="EM_clustering_labels_using_allfeatures.txt",sep="\t")
                 
                 
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
                 
                 
                 h73<-heatmap.2(as.matrix(data_m_fc_withfeats[,-c(1:2)]), Rowv=as.dendrogram(hr), Colv=as.dendrogram(clust1),  col=heatmap_cols, scale="row",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1,xlab="",ylab="", main="Using all features", ColSideColors=patientcolors)
                 
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
                write.table(t1,file="HCA_clustering_labels_using_allfeatures.txt",sep="\t")
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
            
            dev.off()
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



#if(featselmethod=="limma")
if(featselmethod=="limma" | featselmethod=="limma2way" | featselmethod=="limma2wayrepeat" | featselmethod=="lmreg" | featselmethod=="logitreg" 
| featselmethod=="lm2wayanova" | featselmethod=="lm1wayanova" | featselmethod=="lm1wayanovarepeat" | featselmethod=="lm2wayanovarepeat" | featselmethod=="wilcox")
{
	summary_res<-cbind(summary_res,exp_fp)
#print(summary_res)
colnames(summary_res)<-c("RSD.thresh","Number of features left after RSD filtering","Number of FDR significant features",paste(pred.eval.method,"-accuracy",sep=""),paste(pred.eval.method," permuted accuracy",sep=""),"Score","Expected_False_Positives")
}else{
	#exp_fp<-round(fdrthresh*feat_sigfdrthresh)
	
	#if(featselmethod=="MARS" | featselmethod=="RF" | featselmethod=="pls" | featselmethod=="o1pls" | featselmethod=="o2pls"){
		
		exp_fp<-rep(NA,dim(summary_res)[1])
	#}
	summary_res<-cbind(summary_res,exp_fp)
		colnames(summary_res)<-c("RSD.thresh","Number of features left after RSD filtering","Number of significant features",paste(pred.eval.method,"-accuracy",sep=""),paste(pred.eval.method," permuted accuracy",sep=""),"Score","Expected_False_Positives")

	}

featselmethod<-parentfeatselmethod
file_name<-paste("Results_summary_",featselmethod,".txt",sep="")
write.table(summary_res,file=file_name,sep="\t",row.names=FALSE)


	


print("##############Level 1: processing complete###########")

#print("best feats is")
#print(best_feats)

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

print("network analysis")
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
	dataA<-read.table(target.metab.file,sep="\t",header=TRUE,quote = "")
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
	print(paste(length(unique(sigfeats_index))," metabolites matched the significant features list",sep=""))
	

	}else{
		sigfeats_index<-NA #seq(1,dim(data_matrix)[1])
		
		}
	
	print(paste(length(unique(g1$common$index.B))," metabolites matched the target list",sep=""))
	
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
    dev.off()
    }
    
	}else{
	if(networktype=="GGM"){
	targetedan_fdr<-get_partial_cornet(data_matrix, sigfeats.index=sigfeats_index,targeted.index=g1$common$index.B,networkscope="targeted",cor.method,abs.cor.thresh,cor.fdrthresh,outloc=outloc,net_node_colors)
    
    if(FALSE){
    pdf("GGMnetworkplot.pdf")
    load("metabnet.Rda")
    print(plot(net_result))
    dev.off()
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
	#print(head(goodfeats))	
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
	save(list=ls(),file=fname)
    }
			################################


	return(list("diffexp_metabs"=goodfeats_allfields,  "mw.an.fdr"=mwan_fdr,"targeted.an.fdr"=targetedan_fdr,"classlabels"=classlabels_orig,"all_metabs"=allmetabs_res))


}
