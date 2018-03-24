get_pcascoredistplots <-
function(X=NA,Y=NA,feature_table_file,parentoutput_dir,class_labels_file,sample.col.opt="rainbow",plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3,col_vec,pairedanalysis=FALSE,pca.cex.val=3,legendlocation="topright",pca.ellipse=TRUE,ellipse.conf.level=0.5,filename="all")
{

#print("Score plots")
get_pcascoredistplots_child(X=X,Y=Y,feature_table_file=feature_table_file,parentoutput_dir=parentoutput_dir,class_labels_file=class_labels_file,sample.col.opt=sample.col.opt,plots.width=plots.width,plots.height=plots.height,plots.res=plots.res, alphacol=alphacol,col_vec=col_vec,pairedanalysis=pairedanalysis,pca.cex.val=pca.cex.val,legendlocation=legendlocation,pca.ellipse=pca.ellipse,ellipse.conf.level=ellipse.conf.level,filename=filename)

}
