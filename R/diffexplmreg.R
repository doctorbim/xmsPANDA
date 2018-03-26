diffexplmreg <-
function(dataA,logistic_reg=FALSE){
	
dataA<-as.data.frame(dataA)

save(dataA,file="lmreg_func.Rda")


if(logistic_reg==TRUE){
    cnames1<-colnames(dataA)
    cnames1[2]<-"Class"
    colnames(dataA)<-cnames1
	
    labels_1<-levels(as.factor(dataA$Class))

dataA$Class<-replace(dataA$Class,which(dataA$Class==labels_1[1]),0)
 dataA$Class<-replace(dataA$Class,which(dataA$Class==labels_1[2]),1)



	a1 <- glm(dataA$Class ~ .,family=binomial(logit),data=dataA)
}else{
a1 <- lm(dataA$Response ~ .,data=dataA) # aov(dataA$Response ~ .,data=dataA) # + chocolate$Factor1*chocolate$Factor2)
}
s1<-summary(a1)

if(FALSE){
anova_res<-anova(a1)
num_rows<-dim(anova_res)
pvalues_factors<-data.frame(t(anova_res["Pr(>F)"][-c(num_rows),]))
}

if(logistic_reg==FALSE){
    
    r2<-s1$adj.r.squared
}else{
    r2<-NA
}
s1<-s1$coefficients

s1<-s1[-c(1),]

if(dim(dataA)[2]<3){ # && dim(dataA)[1]<3){
    #s1<-as.data.frame(s1)
s1<-t(s1)

}


    confint_lower<-s1[,1]-(1.96*s1[,2])
    confint_upper<-s1[,1]+(1.96*s1[,2])
    

return(list("mainpvalues"=s1[,4],"estimates"=s1[,1],"statistic"=s1[,3],"stderr"=s1[,2],"r2"=r2,"confint"=c(confint_lower,confint_upper)))


}
