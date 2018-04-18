diffexplmonewayanovarepeat <-
function(dataA,subject_inf,analysismode="classification"){


dataA<-as.data.frame(dataA)


if(analysismode=="classification"){
dataA$Factor1<-as.factor(dataA$Factor1)
}


Subject<-subject_inf
dataA<-cbind(dataA,subject_inf)


#save(list=ls(),file="res.Rda")
res <- try(lme(as.numeric(Response) ~ Factor1, random = ~ 1 | subject_inf/Factor1, data=dataA,control=lmeControl(opt="optim")),silent=TRUE)

#res <- lme(as.numeric(Response) ~ Factor1, random = ~ 1 | subject_inf/Factor1, data=dataA,control=lmeControl(opt="optim"))

#print(res)

if(is(res,"try-error")){

    #return(res)
    
    return(list("mainpvalues"=NA,"posthoc"=NA))
}else{
anova_res<-anova(res)


num_rows<-dim(anova_res)
pvalues_factors<-data.frame(t(anova_res["p-value"][-c(1),]))

if(analysismode=="classification"){
posthoc<-summary(glht(res,linfct=mcp(Factor1="Tukey")))

posthoc_pvalues<-posthoc$test$pvalues
names(posthoc_pvalues)<-names(posthoc$test$tstat)


names(pvalues_factors)<-rownames(anova_res)[-c(1)]



return(list("mainpvalues"=pvalues_factors,"posthoc"=posthoc_pvalues))
}else{
	return(list("mainpvalues"=pvalues_factors))
	}
}

}
