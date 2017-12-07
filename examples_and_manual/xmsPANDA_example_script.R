
#load xmsPANDA v1.0.6
library(xmsPANDA)


feature_table_file<-"C:\\Users\\KUPPAL2\\Documents\\Mzmine_smokers_nonsmokers_PANDA.txt"
class_labels_file<-"C:\\Users\\KUPPAL2\\Documents\\classlabels.txt"
outloc<-"C:\\Users\\KUPPAL2\\Documents\\testpandav1.0.6\\"

demetabs_res<-diffexp(feature_table_file=feature_table_file,
parentoutput_dir=outloc,
class_labels_file=class_labels_file,
num_replicates = 3,
    feat.filt.thresh =NA, summarize.replicates =TRUE, summary.method="median",summary.na.replacement="zeros",
    rep.max.missing.thresh=0.5,
    all.missing.thresh=NA, group.missing.thresh=NA, input.intensity.scale="raw", 
    log2transform = FALSE, medcenter=FALSE, znormtransform = FALSE, 
    quantile_norm = FALSE, lowess_norm = FALSE, madscaling = FALSE, 
    rsd.filt.list = c(0),  pairedanalysis = FALSE, featselmethod = "limma", 
    fdrthresh = 0.05, fdrmethod = "none", cor.method = "spearman", 
    networktype = "complete", abs.cor.thresh = 0.4, cor.fdrthresh = 0.05, 
    kfold = 10, pred.eval.method = "AUC", feat_weight = 1, globalcor =TRUE, 
    target.metab.file = NA, target.mzmatch.diff = 10, target.rtmatch.diff = NA, 
    max.cor.num = 100, samplermindex = NA, pcacenter =TRUE, 
    pcascale = TRUE, numtrees = 20000, analysismode = "classification", 
    net_node_colors = c("green", "red"), net_legend = TRUE, svm_kernel = "radial", 
    heatmap.col.opt = "RdBu", sample.col.opt = "rainbow", alphacol = 0.3, 
    rf_selmethod = "rawVIMsig", pls_vip_thresh = 1, num_nodes = 2, 
    max_varsel = 100, pls_ncomp = 5, pca.stage2.eval = FALSE, 
    scoreplot_legend = TRUE, pca.global.eval = FALSE, rocfeatlist = seq(2, 10, 1), 
    rocfeatincrement = TRUE, rocclassifier = "svm", 
    foldchangethresh = 0, wgcnarsdthresh = 20, WGCNAmodules = FALSE, 
    optselect = TRUE, max_comp_sel = 1, saveRda = TRUE, legendlocation = "topleft", 
    degree_rank_method = "diffK", pca.cex.val = 4, pca.ellipse = FALSE, 
    ellipse.conf.level = 0.95, pls.permut.count = NA, svm.acc.tolerance = 5, 
    limmadecideTests = TRUE) 


sink(file=NULL)
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
#"MARS": for multiple adaptive regression splines (MARS) based feature selection
#(mode= regression or classification)
#"pls": for partial least squares (PLS) based feature selection
#(mode= regression or classification)
#"spls": for sparse partial least squares (PLS) based feature selection
#(mode= regression or classification)
#"o1pls": for orthogonal partial least squares (OPLS) based feature selection
#(mode= regression or classification)
####################################################################