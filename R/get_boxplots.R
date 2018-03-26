get_boxplots <-
function(X=NA,Y=NA,feature_table_file,parentoutput_dir,class_labels_file,boxplot.col.opt="grey57",sample.col.opt="rainbow",alphacol=0.3,newdevice=TRUE,cex=0.8,replace.by.NA=FALSE,pairedanalysis=FALSE,filename="")
{

get_boxplots_child(X=X,Y=Y,feature_table_file=feature_table_file,parentoutput_dir=parentoutput_dir,class_labels_file=class_labels_file,boxplot.col.opt,sample.col.opt=sample.col.opt,alphacol=alphacol,newdevice=newdevice,cex=cex,replace.by.NA=replace.by.NA,pairedanalysis=pairedanalysis,filename=filename)

}
