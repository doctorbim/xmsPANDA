
#load xmsPANDA v1.0.7.2
library(xmsPANDA)


feature_table_file<-"/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/xmsPaNDA/xmsPANDA_v1.0.7.5/examples_and_manual/Example_feature_table_and_classlabels/feature_table_one_or_two_factor_analysis.txt"
class_labels_file<-"/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/xmsPaNDA/xmsPANDA_v1.0.7.5/examples_and_manual/Example_feature_table_and_classlabels/classlabels_regression.txt"
outloc<-"/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/xmsPaNDA/xmsPANDA_v1.0.7.5/examples_and_manual/Example_feature_table_and_classlabels/reg/"

#start
demetabs_res<-diffexp(
        #1) arguments for input files
        feature_table_file=feature_table_file,
        parentoutput_dir=outloc,
        class_labels_file=class_labels_file,
        input.intensity.scale="raw",

        ##2) data preprocessing order: 1) summarization, 2) filtering by missing values, 3) imputation; 4) transformation and normalization
        num_replicates = 1,
        summarize.replicates =TRUE, summary.method="median",summary.na.replacement="halffeaturemin",
        rep.max.missing.thresh=0.5,
        all.missing.thresh=0.5, group.missing.thresh=0.8, missing.val=0,
        log2transform = FALSE, medcenter=FALSE, znormtransform = FALSE,
        quantile_norm = FALSE, lowess_norm = FALSE, madscaling = FALSE,
        rsd.filt.list = c(1),

        ##3) arguments for feature seletion:
        pairedanalysis = FALSE, featselmethod=c("lmreg","pls"),
        fdrthresh = 0.1, fdrmethod="BH",
        kfold=5,networktype="complete",
        samplermindex=NA,numtrees=5000,analysismode="regression",pls_vip_thresh = 2, num_nodes = 3,
        max_varsel = 10, pls_ncomp = 5,pred.eval.method="BER",rocfeatlist=seq(2,10,1),
        rocfeatincrement=TRUE,
        rocclassifier="svm",foldchangethresh=0,
        optselect=TRUE,max_comp_sel=5,saveRda=FALSE,pls.permut.count=1000,
        pca.ellipse=TRUE,ellipse.conf.level=0.95,svm.acc.tolerance=5,pamr.threshold.select.max=FALSE,
        aggregation.method="RankAggreg",

        #4) arguments for WGCNA and global clustering analysis (HCA and EM clustering)
        wgcnarsdthresh=30,WGCNAmodules=FALSE,globalclustering=FALSE,

        #5) arguments for correlation and network analysis using the selected features
        cor.method="spearman", abs.cor.thresh = 0.4, cor.fdrthresh=0.2,
        globalcor=FALSE,target.metab.file=NA,
        target.mzmatch.diff=10,target.rtmatch.diff=NA,max.cor.num=NA,

        #6) arguments for graphical options
        pca.cex.val=4,legendlocation="bottomleft",
        net_node_colors=c("green","red"),
        net_legend=FALSE,heatmap.col.opt="RdBu",sample.col.opt="rainbow",alphacol=0.3
)
sink(file=NULL)
#end
#####################################################


####################################################################
#Options for featselmethod:
#"limma": for one-way ANOVA using LIMMA (mode=classification)
#"limma2way": for two-way ANOVA using LIMMA (mode=classification)
#"limma1wayrepeat": for one-way ANOVA repeated measures using LIMMA (mode=classification)
#"limma2wayrepeat": for two-way ANOVA repeated measures using LIMMA (mode=classification)
#"lm1wayanova": for one-way ANOVA using linear model (mode=classification)
#"lm2wayanova": for two-way ANOVA using linear model (mode=classification)
#"lm1wayanovarepeat": for one-way ANOVA repeated measures using linear model (mode=classification)
#"lm2wayanovarepeat": for two-way ANOVA repeated measures using linear model (mode=classification)
#"lmreg": variable selection based on p-values calculated using a linear regression model; 
#allows adjustment for covariates (mode= regression or classification)
#"logitreg": variable selection based on p-values calculated using a logistic regression model; 
# allows adjustment for covariates (mode= classification)
#"rfesvm": uses recursive feature elimination SVM algorithm for variable selection; 
#(mode=classification)
#"wilcox": uses Wilcoxon tests for variable selection; 
#(mode=classification)
#"RF": for random forest based feature selection (mode= regression or classification)
#"RFconditional": for conditional random forest based feature selection (mode= regression or classification)
#"pamr": for prediction analysis for microarrays algorithm based on the nearest shrunked centroid method (mode=classification)
#"MARS": for multiple adaptive regression splines (MARS) based feature selection
#(mode= regression or classification)
#"pls": for partial least squares (PLS) based feature selection
#(mode= regression or classification)
#"spls": for sparse partial least squares (PLS) based feature selection
#(mode= regression or classification)
#"spls1wayanovarepeat": for one-way ANOVA repeated measures using linear model (mode=classification)
#"lm2wayanovarepeat": for two-way ANOVA repeated measures using linear model (mode=classification)
#"o1pls": for orthogonal partial least squares (OPLS) based feature selection
#(mode= regression or classification)
####################################################################
