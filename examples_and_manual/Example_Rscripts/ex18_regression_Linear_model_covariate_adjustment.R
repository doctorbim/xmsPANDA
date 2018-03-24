
#Ex7: Regression using linear model with or without covariates
library(xmsPANDA)


demetabs_res<-diffexp(feature_table_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/lcms_table.txt",
parentoutput_dir="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/testlive_xmsPANDAv0.0.7lmreg/",
class_labels_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/proteinAB_levels.txt",
num_replicates = 1,
    feat.filt.thresh =NA, summarize.replicates =FALSE, summary.method="mean",summary.na.replacement="zeros",rep.max.missing.thresh=0.5,
    all.missing.thresh=0.25,
    group.missing.thresh=NA, 
    log2transform = FALSE, medcenter=FALSE, znormtransform = FALSE, 
    quantile_norm = FALSE, lowess_norm = FALSE, madscaling = FALSE, 
    rsd.filt.list = seq(0, 30, 5), pairedanalysis = FALSE, featselmethod="lmreg", 
    fdrthresh = 0.2, fdrmethod="BH",cor.method="spearman", abs.cor.thresh = 0.4, cor.fdrthresh=0.2,
    kfold=10,feat_weight=1,globalcor=FALSE,target.metab.file=NA,
    target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=200,missing.val=0,networktype="complete",
    samplermindex=NA,numtrees=1000,analysismode="regression",net_node_colors=c("green","red"), 
    net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3) 
    

    
names(demetabs_res)
