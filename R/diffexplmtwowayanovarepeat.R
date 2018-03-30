diffexplmtwowayanovarepeat <-
function(dataA,subject_inf){

dataA<-as.data.frame(dataA)

dataA$Factor1<-as.factor(dataA$Factor1)
dataA$Factor2<-as.factor(dataA$Factor2)

subject_inf<-as.vector(subject_inf)

Subject<-subject_inf
dataA<-cbind(dataA,subject_inf)

res <- try(lme(as.numeric(Response) ~ Factor1 + Factor2 + Factor1 * Factor2, random = ~ 1 | subject_inf/Factor2, data=dataA,control=lmeControl(opt="optim")),silent=TRUE)

#res <- lme(as.numeric(Response) ~ Factor1 + Factor2 + Factor1 * Factor2, random = ~ 1 | subject_inf/Factor2, data=dataA,control=lmeControl(opt="optim"))



if(is(res,"try-error")){
    
    #return(res)
    
    return(list("mainpvalues"=NA,"posthoc"=NA))
}else{
    
    anova_res<-anova(res)


#print(anova_res)

dataA$SHD<-interaction(dataA$Factor1,dataA$Factor2)
mod2<-lme(Response~-1+SHD, data=dataA, random=~1|subject_inf/Factor2,control=lmeControl(opt="optim"))
posthoc<-summary(glht(mod2,linfct=mcp(SHD="Tukey")))

#summary(posthoc)

posthoc_pvalues<-posthoc$test$pvalues
names(posthoc_pvalues)<-names(posthoc$test$tstat)

num_rows<-dim(anova_res)
pvalues_factors<-data.frame(t(anova_res["p-value"][-c(1),]))

names(pvalues_factors)<-rownames(anova_res)[-c(1)]

return(list("mainpvalues"=pvalues_factors,"posthoc"=posthoc_pvalues))
}

}
