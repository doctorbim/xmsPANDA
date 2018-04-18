do_multilevelpls <-
function(dataA,stimu.time,repeat.stimu2,keepX=20)
{

t500<-tune.multilevel(t(data_m), cond = stimu.time,
sample = repeat.stimu2, method = 'splsda',test.keepX = c(500, 500,500))



res.2level_top10 <- multilevel(t(data_m), cond = stimu.time,
sample = repeat.stimu2, ncomp = 3,
keepX = c(10, 10,10), tab.prob.gene = NULL, method = 'splsda')

res.2level <- multilevel(t(data_m), cond = stimu.time,
sample = repeat.stimu2, ncomp = 3,
keepX = c(200, 200,200), tab.prob.gene = NULL, method = 'splsda')

goodcomp1<-which(res.2level$loadings$X[,1]!=0)
goodcomp2<-which(res.2level$loadings$X[,2]!=0)
goodcomp3<-which(res.2level$loadings$X[,3]!=0)
good_metabs<-data_matrix[c(goodcomp1,goodcomp2,goodcomp3),]
good_metabs<-good_metabs[order(good_metabs$mz),]

#res.2level_rand<-new("list")
good_metabs_rand<-new("list")
rand_mz<-new("list")
rand_mz_length<-{}
rand_mz_list<-{}

for(i in 1:1000){

stimu.time_rand<-stimu.time[sample(x=1:dim(stimu.time)[1],size=dim(stimu.time)[1]),]
repeat.stimu2_rand<-repeat.stimu2[sample(x=1:dim(stimu.time)[1],size=dim(stimu.time)[1])]

res.2level_rand <- multilevel(t(data_m), cond = stimu.time_rand,
sample = repeat.stimu2_rand, ncomp =3,
keepX = c(200, 200, 200), tab.prob.gene = NULL, method = 'splsda')


goodcomp1<-which(res.2level_rand$loadings$X[,1]!=0)
goodcomp2<-which(res.2level_rand$loadings$X[,2]!=0)
goodcomp3<-which(res.2level_rand$loadings$X[,3]!=0)
good_metabs_rand_cur<-data_matrix[c(goodcomp1,goodcomp2,goodcomp3),1:2]
good_metabs_rand_cur<-good_metabs_rand_cur[order(good_metabs_rand_cur$mz),1:2]

good_metabs_rand[[i]]<-good_metabs_rand_cur

rand_mz[[i]]<-which(good_metabs$mz%in%good_metabs_rand[[i]]$mz)


rand_mz_length<-c(rand_mz_length,length(which(good_metabs_rand[[i]]$mz%in%good_metabs$mz)))
rand_mz_list<-c(rand_mz_list,good_metabs_rand[[i]]$mz[which(good_metabs_rand[[i]]$mz%in%good_metabs$mz)])

}

#save(list=ls(),file="1000permut_anal_msplsda.Rda")

t1<-table(rand_mz_list)
 t2<-t1/1000
 
 t3<-as.vector(t2)
 # length(which(t3<0.2)) =414
 #length(which(t3<0.05)) =71
 
 good_rand_eval_list<-names(t2[which(t2<0.05)])

good_metabs_filt_by_permut<-good_metabs[which(good_metabs$mz%in%as.numeric(good_rand_eval_list)),]


X1<-res.2level$variates$X %*% t(res.2level$loadings$X)



X2<-t(X1)
s1<-apply(X3,1,sum)
X3<-X2[-which(s1==0),]
S1<-cor(t(X3),use="pairwise.complete.obs")
A1<-eigen(S1)
(A1$values[1])/sum(A1$values)

#33.84, 33.67, 32.48

# color for plotIndiv
col.stimu = as.numeric(classlabels[,4])

color_vec<-c("green","purple")
col.stimu<-color_vec[col.stimu]

# pch for plots
#pch.time = rep(20, 48)
#pch.time[time == 't2'] = 4

pch.time = c(rep(c(3,5,7,9,12,13,2,17,21),7),c(3,7,9,12,13,2,17,21),rep(c(3,5,7,9,12,13,2,17,21),4))
#
#c(rep(15,18),rep(16,18),rep(3,18),rep(17,17),rep(18,18),rep(20,18))

cex <- rep(0.5,length(col.stimu))
#pdf("splsda.pdf")

tiff("pca_plot_alltop200_pc2and3_600dpi.tiff",width=3000,height=3000,res=300)
plotIndiv(res.2level, col = col.stimu, pch = pch.time, ind.names = FALSE,pt.cex=cex,cex=0.6,comp=c(2,3))
#legend('bottomright', col = c(rep("black",9),"green","purple"),
#legend = c(levels(as.factor(classlabels[,3])), levels(as.factor(classlabels[,4]))), pch = c(3,5,7,9,12,13,2,17,21,15,15), cex = 0.55)

dev.off()

d1<-pca(t(data_m),center=TRUE,scale=TRUE,ncomp=10)
tiff("pca1and2_plot_allmz_600dpi.tiff",width=1500,height=1500,res=300)
plotIndiv(d1, col = col.stimu, pch = pch.time, ind.names = FALSE,pt.cex=cex,cex=0.6,comp=c(1,2))
#legend('bottomright', col = c(rep("black",9),"green","purple"),
#legend = c(levels(as.factor(classlabels[,3])), levels(as.factor(classlabels[,4]))), pch = c(3,5,7,9,12,13,2,17,21,15,15), cex = 0.35)

dev.off()

tiff("pca2and3_plot_allmz_600dpi.tiff",width=1500,height=1500,res=300)
plotIndiv(d1, col = col.stimu, pch = pch.time, ind.names = FALSE,pt.cex=cex,cex=0.6,comp=c(2,3))
#legend('bottomright', col = c(rep("black",9),"green","purple"),
#legend = c(levels(as.factor(classlabels[,3])), levels(as.factor(classlabels[,4]))), pch = c(3,5,7,9,12,13,2,17,21,15,15), cex = 0.35)

dev.off()

tiff("pca1and5_plot_allmz_600dpi.tiff",width=1500,height=1500,res=300)
plotIndiv(d1, col = col.stimu, pch = pch.time, ind.names = FALSE,pt.cex=cex,cex=0.6,comp=c(1,5))
#legend('bottomright', col = c(rep("black",9),"green","purple"),
#legend = c(levels(as.factor(classlabels[,3])), levels(as.factor(classlabels[,4]))), pch = c(3,5,7,9,12,13,2,17,21,15,15), cex = 0.35)

dev.off()


#legend('topright', col = 'black', legend = levels(as.factor(classlabels[,3])),
#pch = unique(pch.time), cex = 0.8)



write.table(good_metabs,file="splsda_twofactor_top200metabs_pc1and2and3_raw.txt",sep="\t",row.names=TRUE)

}
