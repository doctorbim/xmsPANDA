

library(xmsPANDA)

avg_data<-data_preprocess(feature_table_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/lcms_table.txt",
parentoutput_dir="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/preprocessing_results/",
class_labels_file=NA,
num_replicates=1,
feat.filt.thresh=NA,summarize.replicates=TRUE,summary.method="mean",
 all.missing.thresh=NA,group.missing.thresh=NA,
log2transform=FALSE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=FALSE,
lowess_norm=FALSE,madscaling=FALSE,missing.val=0,samplermindex=NA, 
rep.max.missing.thresh=0.5,
summary.na.replacement="zeros")