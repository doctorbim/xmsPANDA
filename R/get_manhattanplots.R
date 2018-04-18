get_manhattanplots <-
function(xvec,yvec,up_or_down,maintext="",ythresh=0.05,ylab,xlab,colorvec=c("darkblue","red3"),col_seq=c("brown","chocolate3","orange3","coral","pink","skyblue","blue","darkblue","purple","violet"),xincrement=100,yincrement=1,y2thresh=NA,pchvec=c(21,21)){
    
    d4<-xvec
    min_val<-min(0,d4)
    max_val<-max(d4)
    
   
    
    windowsize=xincrement
    
    d4<-as.vector(d4)
    pvalues<-as.vector(yvec)
    
    
    
    logp<-as.vector(yvec) #(-1)*log((pvalues+0.0001),10)
    
    if(is.na(up_or_down)==TRUE){
        up_or_down<-rep(1,length(yvec))
    }
    
    plot(d4,logp,xaxt="n",ylab=ylab,xlab=xlab,xaxt="n",yaxt="n",cex=0.2,cex.main=1,main=maintext)
    axis(1, at=seq(min_val , max_val, by=xincrement) , las=2)
    axis(2, at=seq(0 , (max(logp)+2), by=yincrement) , las=2)
    #col_seq<-c("brown","red","orange3","coral","pink","skyblue","blue","darkblue","purple","violet")
    
    if(length(col_seq)>1){
        
    s1<-seq(windowsize,max_val,windowsize)
    points(d4[which(d4>=0 & d4<=windowsize)],logp[which(d4>=0 & d4<=windowsize)],col=col_seq[1],cex=0.2)
    for(i in 1:(length(s1)-1))
    {
        points(d4[which(d4>s1[i] & d4<=s1[i+1])],logp[which(d4>s1[i] & d4<=s1[i+1])],col=col_seq[i+1],cex=0.2)
    }
    }else{
        
        points(d4[which(d4>=0 & d4<=windowsize)],logp[which(d4>=0 & d4<=windowsize)],col="black",cex=0.2)
    }
    
    if(is.na(y2thresh)==TRUE){
        
        goodip<-which(yvec>ythresh)
    }else{
        
        goodip<-which(yvec>y2thresh)
        
    }
    
    
    
    for(i in goodip){
        if(up_or_down[i]>0){
            points(d4[i],logp[i],col=colorvec[1],cex=0.7,pch=pchvec[1],bg=colorvec[1]); points(d4[i],logp[i],col=colorvec[1],cex=0.2,bg=colorvec[1])
        }else{
            
            points(d4[i],logp[i],col=colorvec[2],cex=0.7,pch=pchvec[2],bg=colorvec[2]); points(d4[i],logp[i],col=colorvec[2],cex=0.2,bg=colorvec[2])
        }
    }
    
    if(length(goodip)>0){
								hfdrfdrthresh<-logp[which(logp==min(logp[which(yvec>ythresh)],na.rm=TRUE))]
                                
                                abline(h=hfdrfdrthresh,col="gray8",lty=2,lwd=2)
                                if(is.na(y2thresh)==FALSE){
                                    abline(h=y2thresh,col="gray40",lty=2,lwd=2)
                                }
    }
    #legend("topleft", "(x,y)", pch = 1, title = "",inset = .05)
 
 #legend("topleft",c("Higher in GroupA","Lower in GroupB"),col=c(colorvec[1],colorvec[2]),pch=c(24,25),cex=0.5)
}
