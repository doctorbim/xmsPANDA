
#Ex3: Regression using Random Forest
library(xmsPANDA)


demetabs_res<-diffexp(feature_table_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/lcms_table.txt",
parentoutput_dir="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/RF_reg_Ex/",
class_labels_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/proteinA_levels.txt",
num_replicates = 1,
    feat.filt.thresh =NA, summarize.replicates =FALSE, summary.method="mean",summary.na.replacement="zeros",rep.max.missing.thresh=0.3,group.missing.thresh=0.7, 
    log2transform = FALSE, medcenter=FALSE, znormtransform = FALSE, 
    quantile_norm = FALSE, lowess_norm = FALSE, madscaling = FALSE, 
    rsd.filt.list = seq(0, 20, 5), pairedanalysis = TRUE, featselmethod="RF", 
    fdrthresh = 0.2, fdrmethod="BH",cor.method="spearman", abs.cor.thresh = 0.5, cor.fdrthresh=0.2,
    kfold=10,feat_weight=1,globalcor=FALSE,target.metab.file=NA,
    target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=200,plots.width=2000,plots.height=2000,plots.res=300, missing.val=0,networktype="complete",
    samplermindex=NA,numtrees=1000,analysismode="regression",net_node_colors=c("green","red"), net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3) 
    

    
names(demetabs_res)
