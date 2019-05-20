\name{diffexp}
\alias{diffexp}
\title{diffexp}
\description{
This function performs biomarker discovery and generates a correlation network based on the 
metabolome-wide (and targeted) correlation analysis of the differentially expressed features. 
The "featselmethod" allows users to select the method for selecting discriminatory features.
The function evaluates the k-fold cross-validation accuracy using Support Vector Machine,
performs hierarchical clustering analysis, PCA analysis (R2/Q2 diagnostics), and generates 
boxplots for discriminatory features identified at each relative standard deviation (coefficient
 of variation) threshold across all samples (if one feature selection method is used). 
An optimization score that minimizes the number of false positives and increases the 
classification accuracy is used to select the best set of features. The best set is 
then used for correlation (complete or partial) based metabolome-wide network analysis. 
Additionally, users have the option to provide a list of mzs corresponding to chemicals
of interest such as (phenylalanine, choline, etc). The function uses the getVenn 
function in xMSanalyzer to find the mzs matching the
target list based on a user defined mass search threshold (+/- ppm). 


}
\usage{
diffexp(Xmat = NA, Ymat = NA, feature_table_file, parentoutput_dir, class_labels_file, 
num_replicates = 3, summarize.replicates = TRUE, summary.method = "mean", 
summary.na.replacement = "zeros", missing.val = 0, rep.max.missing.thresh = 0.3, 
all.missing.thresh = 0.5, group.missing.thresh = 0.8, input.intensity.scale = "raw", 
log2transform = TRUE, medcenter = FALSE, znormtransform = FALSE, quantile_norm = TRUE, 
lowess_norm = FALSE, madscaling = FALSE, rsd.filt.list = seq(0, 75, 5), pairedanalysis = FALSE,
 featselmethod = "limma", pvalue.thresh = 0.05, fdrthresh = 0.05, fdrmethod = "BH", 
 cor.method = "spearman", networktype = "complete", abs.cor.thresh = 0.4, cor.fdrthresh = 0.05,
  kfold = 10, pred.eval.method = "BER", globalcor = TRUE, target.metab.file = NA, 
  target.mzmatch.diff = 10, target.rtmatch.diff = NA, max.cor.num = 100, 
  numtrees = 20000, analysismode = "classification", net_node_colors = c("green", "red"),
   net_legend = TRUE, svm_kernel = "radial", heatmap.col.opt = "redblue", 
   manhattanplot.col.opt = c("darkblue", "red3"), boxplot.col.opt = c("grey57"), 
   sample.col.opt = "rainbow", rf_selmethod = "absVIMthresh", pls_vip_thresh = 1,
    num_nodes = 2, max_varsel = 100, pls_ncomp = 5, pca.stage2.eval = TRUE, 
    scoreplot_legend = TRUE, pca.global.eval = TRUE, rocfeatlist = seq(2, 11, 1), 
    rocfeatincrement = TRUE, rocclassifier = "svm", foldchangethresh = 2, wgcnarsdthresh = 20, 
    WGCNAmodules = TRUE, optselect = TRUE, max_comp_sel = 1, saveRda = TRUE, 
    legendlocation = "topleft", pca.cex.val = 4, pca.ellipse = FALSE, ellipse.conf.level = 0.95,
     permutations.count = 1000, svm.acc.tolerance = 5, limmadecideTests = TRUE, 
     pls.vip.selection = "max", globalclustering = TRUE, plots.res = 600, plots.width = 8, 
     plots.height = 8, plots.type = "cairo", output.device.type = "pdf", ,
      individualsampleplot.col.opt = NA, pamr.threshold.select.max = FALSE, 
      aggregation.method = "RankAggreg",aggregation.max.iter=1000,mars.gcv.thresh=10, ...)
}
\arguments{
 \item{Xmat}{
 R object for feature table. If this is given, then feature table can be set to NA.
 }
  \item{Ymat}{
  R object for response/class labels matrix. If this is given, then class can be set to NA.
  }
  \item{feature_table_file}{
Feature table that includes the mz, retention time, and measured intensity in each sample 
for each analyte. The first 2 columns should be the mz and time. The remaining columns
should correspond to the samples in the class labels file with each column including the intensity profile
of a sample.
Full path required. Eg: C:/My Documents/test.txt
The feature table should be in a tab-delimited format. An example of the input file is provided under the
"example" folder.
}
  \item{parentoutput_dir}{
Provide full path of the folder where you want the results to be written.
Eg: C:/My Documents/ProjectA/results/
}
  \item{class_labels_file}{
File with class labels information for each sample. Samples should be in the same order
as in the feature table. Please use the same format as in the example folder. If you
want to adjust for covariates in "lmreg" option, then you can add additional columns,
one per covariate. Categorical variables should be strings (eg: "male", "female").
Please see "classlabels_gender.txt" file as an example.
}
  \item{num_replicates}{
Number of technical replicates
}
  \item{feat.filt.thresh}{
Percent Intensity Difference or Coefficient of variation threshold; feature filtering
Use NA to skip this step. 
}
  \item{summarize.replicates}{
Do the technical replicates per sample need to be averaged or median summarized? 
}
 \item{summary.method}{
Method for summarizing the replicates. Options: "mean" or "median"
}
 \item{summary.na.replacement}{
How should the missing values be represented? 
Options: "zeros", "halffeaturemin", "halfsamplemin","halfdatamin", "none"
"zeros": replaces missing values by 0
"halfsamplemin": replaces missing value by one-half of the lowest signal intensity in the
corresponding sample
"halfdatamin": replaces missing value by one-half of the lowest signal intensity in the
complete dataset
"halffeaturemin": replaces missing value by one-half of the lowest signal intensity for the
current feature
"none": keeps missing values as NAs

Users are recommended to perform imputation prior to performing biomarker discovery.
}
\item{missing.val}{
How are the missing values represented in the input data? Options: "0" or "NA"
}
 \item{rep.max.missing.thresh}{
What propotion of replicates are allowed to have missing values during the averaging or 
median summarization step of each biological sample? If the number of replicates with
missing values is greater than the defined threshold, then the summarized value is 
represented by the "missing.val" parameter. If the number of replicates with missing values
is less than or equal to the defined threshold, then the summarized value is equal to the 
mean or the median of the non-missing values. Default: 0.5
}
  \item{all.missing.thresh}{
What propotion of total number of samples should have an intensity?
Default: 0.5
}
 \item{input.intensity.scale}{
Are the intensities in the input feature table at raw scale or log2 scale?
eg: "raw" or "log2"
Default: "raw"
}
  \item{group.missing.thresh}{
What propotion of samples in either of the two groups should have an intensity?
If at least x% of the samples in either group have a signal, then the feature is retained
for further analysis. Default: 0.7
}
  \item{log2transform}{
Data transformation: Please refer to http://www.biomedcentral.com/1471-2164/7/142
Try different combinations; such as log2transform=TRUE, znormtransfrom=FALSE
or log2transform=FALSE, znormtransfrom=TRUE
}
  \item{medcenter}{
Median centering of metabolites
}
  \item{znormtransform}{
Auto scaling; each metabolite will have a mean of 0 and unit variance
}
  \item{quantile_norm}{
Performs quantile normalization. Normalization options: Please set only one of the options to be TRUE
}
  \item{lowess_norm}{
Performs lowess normalization. Normalization options: Please set only one of the options to be TRUE
}
  \item{madscaling}{
Performs median adjusted scale normalization. Normalization options: Please set only one of the options to be TRUE
}
  \item{rsd.filt.list}{
This parameter allows to perform feature filtering based on overall variance (across all samples)
prior to performing hypothesis testing. Eg: seq(0,30,5). 
}
  \item{pairedanalysis}{
 Is this a paired-study design? TRUE or FALSE
 If samples are paired, then the feature table and the class labels file should be organized so that the paired samples 
 are arranged in the same order in each group. For example, the first sample in group A and the first sample in 
 group B should be paired.  
}
\item{featselmethod}{
Options:
"limma": for one-way ANOVA using LIMMA (mode=classification)
"limma2way": for two-way ANOVA using LIMMA (mode=classification)
"limma1wayrepeat": for one-way ANOVA repeated measures using LIMMA (mode=classification)
"limma2wayrepeat": for two-way ANOVA repeated measures using LIMMA (mode=classification)
"lm1wayanova": for one-way ANOVA using linear model (mode=classification)
"lm2wayanova": for two-way ANOVA using linear model (mode=classification)
"lm1wayanovarepeat": for one-way ANOVA repeated measures using linear model (mode=classification)
"lm2wayanovarepeat": for two-way ANOVA repeated measures using linear model (mode=classification)
"lmreg": variable selection based on p-values calculated using a linear regression model; 
allows adjustment for covariates (mode= regression or classification)
"logitreg": variable selection based on p-values calculated using a logistic regression model; 
 allows adjustment for covariates (mode= classification)
"rfesvm": uses recursive feature elimination SVM algorithm for variable selection; 
(mode=classification)
"RF": for random forest based feature selection (mode= regression or classification)
"RFconditional": for conditional random forest based feature selection (mode= regression or classification)
"pamr": for prediction analysis for microarrays algorithm based on the nearest shrunken centroid
method (mode= classification)
"MARS": for multiple adaptive regression splines (MARS) based feature selection
(mode= regression or classification)
"pls": for partial least squares (PLS) based feature selection
(mode= regression or classification)
"spls": for sparse partial least squares (PLS) based feature selection
(mode= regression or classification)
"spls1wayrepeat": for sparse partial least squares (PLS) based feature selection
for one-way repeated measures (mode= regression or classification)
"spls2wayrepeat": for sparse partial least squares (PLS) based feature selection
for two-way repeated measures (mode= regression or classification)
"o1pls": for orthogonal partial least squares (OPLS) based feature selection
(mode= regression or classification)
}
  \item{fdrthresh}{
False discovery rate threshold. Eg: 0.05
}
 \item{fdrmethod}{
Options: "BH", "ST", "Strimmer", "none"
"BH": Benjamini-Hochberg (1995) (Default: more conservative than "ST" and "Strimmer")
"ST": Storey & Tibshirani (Storey 2001, PNAS) algorithm implemented in the qvalue package
"Strimmer": (Strimmer 2008, Bioinformatics) algorithm implemented in the fdrtool package
"none": No FDR correction will be performed. fdrthresh will be treated as raw p-value 
cutoff
}
\item{cor.method}{
Correlation method. Options: "pearson" or "spearman". Default: "spearman"
}
 \item{networktype}{
Options: "complete" or "GGM"
"complete": performs network analysis using ordinary Pearson or Spearman correlation
statistic
"GGM": generates network based on partial correlation analysis using the 
GeneNet package
}
  \item{abs.cor.thresh}{
Absolute Pearson correlation coefficient for network analysis. Default: 0.4
}
  \item{cor.fdrthresh}{
False discovery rate threshold for correlation analysis. Default: 0.05
}
\item{kfold}{
	Number of subsets in which the data should be divided 
	for cross-validation. If kfold=10, then the data set will
	be divided into 10 subsets of size (N/10), where N is 
	the total number of samples. 9 subsets are used for
	training and the remaining 1 is used for testing. This
	process is repeated 10 times and the CV-accuracy
	would be the mean of the classification accuracy
	from the 10 iterations. The same will be true for any
	other value of k. 
	WARNING: The kfold value should be
	less than or equal to the total number of samples.
}
 \item{pred.eval.method}{
	Criteria for evaluating the performance of the
	model. CV: Overall Cross-validation classification accuracy, 
	balanced error rate (BER): (sum of accuracy in each class)/(number of classes)
	area under the curve (AUC)
	Eg: "CV", "BER", or "AUC". Default: "BER"
 }
\item{globalcor}{
Do you want to perform correlation analysis after biomarker discovery?
Options: "TRUE" or "FALSE"
}
  \item{target.metab.file}{
File that includes the mz and/or retention time of the targeted metabolites.
See example.
}
  \item{target.mzmatch.diff}{
+/- ppm mass tolerance for searching the target m/z in the current
	  feature table
}
  \item{target.rtmatch.diff}{
+/- retention time tolerance for searching the target m/z in the current
	  feature table
}
  \item{max.cor.num}{
Maximum number of correlated metabolites to be included in the network
figure. Default: 100
}
  \item{pcacenter}{
Data centering for PCA. Options: "TRUE" or "FALSE". Default=TRUE
}
  \item{pcascale}{
Data scaling for PCA. Options: "TRUE" or "FALSE". Default=TRUE
}
  \item{samplermindex}{
Column index of any additional or irrelevant columns to be deleted.
Options: "NA" or list of column numbers. eg: c(1,3,4) Default=NA
}
  \item{numtrees}{
Number of trees to be used for random forest method. Default=500
}
  \item{analysismode}{
"classification" for group-wise comparison (case vs control) or
"regression" for continuous response variables. Default: "classification"
}
  \item{net_node_colors}{
 Colors of nodes in the correlation networks. Eg: c("pink", "skyblue"), 
or ("red","green") 
}
  \item{net_legend}{
Should the network be displayed for the correlation network? eg:
TRUE or FALSE
}
  \item{max_var}{
Max number of variables to be used for sPLS, rfesvm, and Random Forest? 
eg:150
}
  \item{svm_kernel}{
SVM kernel eg: "radial" or "linear"
}
  \item{heatmap.col.opt}{
Color scheme for HCA hetmap eg: "RdBu", "topo", "heat", or "terrain"
}
  \item{sample.col.opt}{
Color scheme for PCA and heatmap sample axis 
eg: "rainbow", "heat" or "topo"
}
  \item{rf_selmethod}{
Random forest VIP based selection method. 
If rankbased option is selected, variables are ranked 
based on the Variable Importance Measure. 
Only the top "max_varsel" variables are selected. If
absVIMthresh is selected, then all features with VIM greater 
than the absolute value of the lowest VIM are selected.
eg: "absVIMthresh" or "rankbased"
}
  \item{pls_vip_thresh}{
Threshold for VIP score from PLS/O1PLS.
eg: 1
}
  \item{max_varsel}{
Maximum number of variables to keep if "rankbased"
RF or spls is used. eg: 100
}
\item{pls_ncomp}{
Maximum number of components to be considered 
during the PLS optimal number of components 
selection step.
eg: 2
}
\item{pca.stage2.eval}{
Should PCA diagnostics be performed in stage 2?
eg: TRUE or FALSE
}
\item{scoreplot_legend}{
Should legends be included in score plots?
eg: TRUE or FALSE
}
\item{pca.global.eval}{
Should global PCA evaluation be performed? Default:TRUE
eg: TRUE or FALSE
}
\item{rocfeatlist}{
Vector indicating number of features to be used for ROC
evaluation: eg: c(2,4,6) will generate ROC for top 2, top 4,
and top 6 feautres. Default: seq(2,10,1)
}
\item{rocclassifier}{
Classifier to be used for ROC evaluation. Options: "svm" or
"logitreg". Default: "svm"
}
\item{foldchangethresh}{
Secondary feature selection criteria based on fold change
threshold. This is performed after
statistical significance or importance evaluation.
}
\item{wgcnarsdthresh}{
Relative standard deviation or coefficient of variation
(across all samples) based filtering threshold before performing
WGCNA module preservation analysis. Default: 20
}
\item{WGCNAmodules}{
Perform WGCNA module preservation analysis. TRUE or FALSE
Default: TRUE
}
\item{optselect}{
Determine optimal number of PLS components. Default: TRUE
}
\item{max_comp_sel}{
Number of PLS components to use for VIP or sparse loading
selection (sPLS). Default=1
}
\item{saveRda}{
Should the results be saved in a binary R object. Default: TRUE
}
\item{legendlocation}{
Legend location for PLS or PCA plots
}
 \item{pca.cex.val}{
Size of points on PCA plots. eg: 4
}
\item{pca.ellipse}{
Should ellipse be plotted on PCA plots? 
eg: TRUE or FALSE
}
\item{ellipse.conf.level}{
Confidence interval for PCA ellipses
eg: 0.95
}
\item{pls.permut.count}{
Number of permutations for calculating 
p-values for PLS or sPLS models.
eg: 1000
}
\item{svm.acc.tolerance}{
Stopping criteria for forward feature selection
using "rfeSVM" method. If the difference between
best accuracy and current accuracy based on the newly
added feature drops below the tolerance level, the
forward selection process is terminated.
eg: 5
}
\item{pamr.threshold.select.max}{
If two or more thresholds for shrinking the d statistic
in the PAM algorithm (Tibshirani et al. Statistical
Science 2003) have equal accuracy,
should the maximum value (lowest number of features) be used? 
Default: FALSE
}
\item{aggregation.method}{
Method for combining the results from mutliple feature selection methods
Options: Consensus: will only keep features that are selected in all methods
RankAggreg: will use the cross entropy algorithm with Spearman footrule
distance as the distance measure (RankAggreg; Pihur et al. BMC Bioinformatics 2009)
}
\item{aggregation.max.iter}{
	Maximum number of iterations used in the rank aggregation step.
	Default: 1000
}
\item{mars.gcv.thresh}{
	Minimum generalized cross-validation value (range: 0 to 100) for
	a feature to be selected. Default: 10
}
}
\details{
This function performs data transformation, normalization, FDR analysis using LIMMA, 
variable selection using random forests, evaluates the predictive accuracy of
the FDR significant features using k-fold cross-validation with a Support Vector
Machine classifier, performs two-way hierarchical clustering analysis, and principal
component analysis. An optimizaiton scheme is used to select the best set of features
from different log2 fold change filtering thresholds. Finally, metabolome-wide and
targeted correlation based network analysis of the FDR significant features is performed.
}
\value{
  \item{diffexp_metabs}{Best set of discriminatory features.}
  \item{all_metabs}{Results for all features.}
  \item{mw.an.fdr}{Metabolome-wide significant correlation network of 
  differentially expressed metabolites.}
  \item{targeted.an.fdr}{Correlation network of  differentially expressed
  metabolites with targeted metabolites.}
  
  Following files are generated in the parent output location:
Manhattan plots: showing metabolome wide p-values; 
Heatmap from Two-way hierarchical clustering analysis;
Pairwise score plots from Principal Component Analysis;
PCA score distribution plots;
ROC plots;
List of differentially expressed metabolites;
Boxplots of differentially expressed metabolites;
Correlation network figure and matrix;
Pairwise correlation matrix CIRCOS format ready to be uploaded to:
  http://mkweb.bcgsc.ca/tableviewer/visualize/
Or uploaded to Cytoscape gml format

}
\author{Karan Uppal <kuppal2@emory.edu>}
