get_volcanoplots <-
function(xvec,yvec,up_or_down,maintext="",ythresh=0.05,y2thresh=NA,ylab,xlab,colorvec=c("darkblue","red3"),col_seq=c("brown","chocolate3","orange3","coral","pink","skyblue","blue","darkblue","purple","violet"),xincrement=1,yincrement=1,xthresh=1,pchvec=c(21,21)){
    
    d4<-xvec
    min_val<-round(min(d4,na.rm=TRUE)+0.5)
    max_val<-round(max(d4,na.rm=TRUE)+0.5)
    
    windowsize=xincrement
    
    d4<-as.vector(d4)
    
    logp<-as.vector(yvec)
    
    if(is.na(up_or_down)==TRUE){
        up_or_down<-rep(1,length(yvec))
    }
    
    plot(d4,logp,xaxt="n",ylab=ylab,xlab=xlab,xaxt="n",yaxt="n",cex=0.2,cex.main=1,main=maintext)
    axis(1, at=seq(min_val , max_val, by=xincrement) , las=2)
    axis(2, at=seq(0 , (max(logp)+2), by=yincrement) , las=2)

    points(d4[which(d4>=0 & d4<=windowsize)],logp[which(d4>=0 & d4<=windowsize)],col="black",cex=0.2)
    
    
    
    goodip<-which(yvec>ythresh & abs(xvec)>xthresh)
   
    
    
    for(i in goodip){
        if(up_or_down[i]>0){
            points(d4[i],logp[i],col=colorvec[1],cex=0.7,pch=pchvec[1],bg=colorvec[1]); points(d4[i],logp[i],col=colorvec[1],cex=0.2,bg=colorvec[1])
        }else{
            
            points(d4[i],logp[i],col=colorvec[2],cex=0.7,pch=pchvec[2],bg=colorvec[2]); points(d4[i],logp[i],col=colorvec[2],cex=0.2,bg=colorvec[2])
        }
    }
    
    if(length(goodip)>0){
        
        
        abline(v=(-1)*xthresh,col="gray8",lty=2,lwd=2)
        abline(v=xthresh,col="gray8",lty=2,lwd=2)
        
        abline(h=ythresh,col="gray8",lty=2,lwd=2)
        
        if(is.na(y2thresh)==FALSE){
            abline(h=y2thresh,col="gray8",lty=2,lwd=2)
            
        }
        
    }
 
 
}