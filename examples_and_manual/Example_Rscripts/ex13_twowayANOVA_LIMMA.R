

library(xmsPANDA)

demetabs_res<-diffexp(feature_table_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/lcms_table.txt",
parentoutput_dir="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/testlive_xmsPANDAv0.0.7limma2wayanovaclass/",
class_labels_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/classlabels_gender.txt",
num_replicates = 1,
    feat.filt.thresh =NA, summarize.replicates =FALSE, summary.method="mean",summary.na.replacement="zeros",rep.max.missing.thresh=0.5,
    all.missing.thresh=0.25,
    group.missing.thresh=NA, 
    log2transform = FALSE, medcenter=FALSE, znormtransform = FALSE, 
    quantile_norm = FALSE, lowess_norm = FALSE, madscaling = FALSE, 
    rsd.filt.list = c(-1), pairedanalysis = FALSE, featselmethod="limma2way", 
    fdrthresh = 0.2, fdrmethod="BH",cor.method="spearman", abs.cor.thresh = 0.5, cor.fdrthresh=0.2,
    kfold=10,feat_weight=1,globalcor=FALSE,target.metab.file=NA,
    target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=200,missing.val=0,networktype="complete",
    samplermindex=NA,numtrees=1000,analysismode="classification",net_node_colors=c("green","red"), 
    net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3) 
    
