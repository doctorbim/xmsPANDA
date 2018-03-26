library(xmsPANDA)


demetabs_res<-diffexp(feature_table_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/lcms_table.txt",
parentoutput_dir="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/testlive_xmsPANDAv0.0.6_v35pdfsfdr0.5o2pls/",
class_labels_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/classlabels.txt",
num_replicates = 1,
    feat.filt.thresh =NA, summarize.replicates =FALSE, summary.method="mean",summary.na.replacement="zeros",rep.max.missing.thresh=0.3,group.missing.thresh=0.7, 
    log2transform = TRUE, medcenter=FALSE, znormtransform = FALSE, 
    quantile_norm = TRUE, lowess_norm = FALSE, madscaling = FALSE, 
    rsd.filt.list = seq(0, 30, 5), pairedanalysis = FALSE, featselmethod="o2pls", 
    fdrthresh = 0.2, fdrmethod="BH",cor.method="spearman", abs.cor.thresh = 0.5, cor.fdrthresh=0.2,
    kfold=10,feat_weight=1,globalcor=TRUE,target.metab.file=NA,
    target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=200,plots.width=2000,plots.height=2000,plots.res=300, missing.val=0,networktype="complete",
    samplermindex=NA,numtrees=1000,analysismode="classification",net_node_colors=c("green","red"), net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3) 
    
names(demetabs_res)


