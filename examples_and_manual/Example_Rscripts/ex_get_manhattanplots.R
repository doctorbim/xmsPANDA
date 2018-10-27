
library(xmsPANDA)

#################################
#file location
lmreg_file<-"/Users/Projects/sample_input_file_get_manhattanplots.txt"

#read file
d1<-read.table(lmreg_file,sep="\t",header=TRUE)

xvec1<-d1[,1] #get m/z information for type 1 Manhattan plot: column A
xvec2<-d1[,2] #get time information for type 2 Manhattan plot: column B
yvec<-(-1)*log(d1[,3]+0.001,10) #Calculate (-1)*log10(pvalue): column C
zvec<-d1[,4]  #or NA; directionality; up or down; this could also be log2 fold change: optional column D
ythresh=1.301 #y-axis threshold for significance; e.g. raw p-value or FDR
#y2thresh=NA #optional second y-axis threshold for secondary criteria for significance; e.g. FDR

options(warn=-1)
pdf("Manhattanplot_type1.pdf")

maintext1="Type 1 manhattan plot (-logp vs mz) \n m/z features above the dashed horizontal line meet the selection criteria"
get_manhattanplots(xvec=xvec1,yvec=yvec,ythresh=ythresh,y2thresh=NA,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="-logP",xincrement=150,yincrement=0.5,colorvec=c("darkblue","red3"),col_seq=c("black"),maintext=maintext1,background.points.col="gray40")
dev.off()

pdf("Manhattanplot_type2.pdf")
maintext2="Type 2 manhattan plot (-logp vs time) \n m/z features above the dashed horizontal line meet the selection criteria"
get_manhattanplots(xvec=xvec2,yvec=yvec,ythresh=ythresh,y2thresh=2.5,up_or_down=zvec,xlab="mass-to-charge (m/z)",ylab="-logP",xincrement=150,yincrement=0.5,colorvec=c("darkblue","red3"),col_seq=c("gray2"),maintext=maintext2,background.points.col="gray40")
dev.off()


#################################
