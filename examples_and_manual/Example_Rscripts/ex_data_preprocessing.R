

library(xmsPANDA)

avg_data<-data_preprocess(feature_table_file="/Users/Documents/feature_table.txt",
parentoutput_dir="/Users/Documents/preprocessing_results/",
class_labels_file=NA,
num_replicates=3,
feat.filt.thresh=NA,summarize.replicates=TRUE,summary.method="mean",
 all.missing.thresh=0.5,group.missing.thresh=0.8,
log2transform=FALSE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=FALSE,
lowess_norm=FALSE,madscaling=FALSE,missing.val=0,samplermindex=NA, 
rep.max.missing.thresh=0.3,
summary.na.replacement="zeros")
