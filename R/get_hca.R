get_hca <-
function(feature_table_file=NA,parentoutput_dir,class_labels_file=NA,X=NA,Y=NA,heatmap.col.opt="RdBu",cor.method="spearman",is.data.znorm=FALSE,analysismode="classification",
sample.col.opt="rainbow",plots.width=8,plots.height=8,plots.res=600, plots.type="cairo", alphacol=0.3, hca_type="two-way",newdevice=FALSE,input.type="intensity",mainlab="",cexRow=1, cexCol=1)
{
    
    h73<-get_hca_child(X=X,Y=Y,feature_table_file=feature_table_file,parentoutput_dir=parentoutput_dir,class_labels_file=class_labels_file,heatmap.col.opt=heatmap.col.opt,cor.method=cor.method,is.data.znorm=is.data.znorm,analysismode=analysismode,
    sample.col.opt=sample.col.opt,plots.width=plots.width,plots.height=plots.height,plots.res=plots.res, plots.type=plots.type, alphacol=alphacol, hca_type=hca_type,newdevice=newdevice,input.type=input.type,mainlab=mainlab,cexRow=cexRow, cexCol=cexCol)
    return(h73)
}
