library(xmsPANDA)


demetabs_res<-diffexp(feature_table_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/feature_table_time_series.txt",
parentoutput_dir="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/testlive_xmsPANDAv1.0.3.1limma2wayrepeatclass/",
class_labels_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/cress_time_and_group_classlabels_lm2wayanovarepeat.txt",
num_replicates = 1,
feat.filt.thresh =NA, summarize.replicates =FALSE, summary.method="mean",summary.na.replacement="zeros",rep.max.missing.thresh=0.5,
all.missing.thresh=0.1,
group.missing.thresh=0.7,
log2transform = TRUE, medcenter=FALSE, znormtransform = FALSE,
quantile_norm = TRUE, lowess_norm = FALSE, madscaling = FALSE,
rsd.filt.list = c(1,5,10,15), pairedanalysis = TRUE, featselmethod="lm2wayanovarepeat",
fdrthresh = 0.05, fdrmethod="BH",cor.method="spearman", abs.cor.thresh = 0.5, cor.fdrthresh=0.2,
kfold=10,feat_weight=1,globalcor=FALSE,target.metab.file=NA,
target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=200,missing.val=0,networktype="complete",
samplermindex=NA,numtrees=1000,analysismode="classification",net_node_colors=c("green","red"),
net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3,pca.stage2.eval=TRUE,scoreplot_legend=TRUE)

