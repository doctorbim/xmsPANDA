
library(xmsPANDA)

feature_table_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/lcms_table.txt"
parentoutput_dir="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/testlive_boxplots3/"
class_labels_file="/Users/karanuppal/Documents/Emory/JonesLab/Projects/DifferentialExpression/metaboanalyst/classlabels_gender.txt"

plots.width=2000
plots.height=2000
plots.res=300

sample.col.opt="rainbow"

	
get_boxplots(feature_table_file,parentoutput_dir,class_labels_file,sample.col.opt="rainbow",plots.width=2000,plots.height=2000,plots.res=300, alphacol=0.3)


