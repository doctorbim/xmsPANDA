
library(xmsPANDA)

feature_table_file="/Users/karanuppal/Mzmine_smokers_nonsmokers_PANDA.txt"
parentoutput_dir="/Users/karanuppal/test_boxplots/"
class_labels_file="/Users/karanuppal/classlabels.txt"

plots.width=2000
plots.height=2000
plots.res=300

sample.col.opt="rainbow"

	
get_boxplots(feature_table_file,parentoutput_dir,class_labels_file,sample.col.opt="rainbow",plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3)


