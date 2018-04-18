
library(xmsPANDA)


##training set files###
feature_table_file="/Users/karanuppal/Mzmine_smokers_nonsmokers_PANDA.txt"
class_labels_file="/Users/karanuppal/classlabels.txt"
###################

##test set files### Same files are used in this example
testfeature_table_file="/Users/karanuppal/Mzmine_smokers_nonsmokers_PANDA.txt"
testclass_labels_file="/Users/karanuppal/classlabels.txt"
###############

##Read files####
featuretable<-read.table(feature_table_file,sep="\t",header=TRUE)
testfeaturetable<-read.table(testfeature_table_file,sep="\t",header=TRUE)

testclasslabels<-read.table(testclass_labels_file,sep="\t",header=TRUE)
classlabels<-read.table(class_labels_file,sep="\t",header=TRUE)

################

classlabels<-classlabels[,2] #column with class information
testclasslabels<-testclasslabels[,2] #column with class information

kfold=dim(featuretable[,-c(1:2)])[2] #LOOCV recommended for only small studies N<30; or can be set to 5 or 10 for k-fold CV

errortype="BAR"; #other options: "AUC", "total"

res<-get_classification.accuracy(kfold=kfold,featuretable=featuretable,classlabels=classlabels,classifier="naivebayes",testfeaturetable=testfeaturetable,testclasslabels=testclasslabels,errortype="BAR",kernelname="linear")

print(res)
