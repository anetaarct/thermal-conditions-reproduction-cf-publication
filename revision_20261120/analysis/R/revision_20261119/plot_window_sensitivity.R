args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args))args[1] else '.'
d<-read.csv(file.path(root,'window_sensitivity_20261119/reference_contrasts.csv'))
stopifnot(nrow(d)==24)
window_order<-c('original','apr01_may07','apr01_apr30','apr15_may07')
labels<-c('26 Mar-7 May (43 d)','1 Apr-7 May (37 d)','1-30 Apr (30 d)','15 Apr-7 May (23 d)')
pdf(file.path(root,'figures/Fig_window_sensitivity_20261119.pdf'),width=10,height=11,useDingbats=FALSE)
par(mfrow=c(4,2),mar=c(3.2,8.3,2.4,1),mgp=c(1.9,.6,0),cex=.8)
panel<-0
for(trait in c('clutch','fledging'))for(metric in c('mean','consistency'))for(model in c('baseline','within_between')){
  panel<-panel+1
  x<-d[d$trait==trait & d$model==model & d$predictor==if(model=='baseline')'temp_mean_c' else 'temp_within',]
  x<-x[match(window_order,x$window),]
  if(metric=='mean'){
    multiplier<-if(trait=='fledging')100 else 1
    estimate<-x$mean_change*multiplier;lo<-x$lower*multiplier;hi<-x$upper*multiplier;null<-0
    xlab<-if(trait=='clutch')'Mean clutch-size change (eggs / degree C)' else 'Full-mixture success change (percentage points / degree C)'
  }else{
    estimate<-x$consistency_ratio;lo<-x$ratio_lower;hi<-x$ratio_upper;null<-1
    xlab<-if(trait=='clutch')'Conditional clutch-size SD ratio / degree C' else 'Component precision ratio / degree C'
  }
  plot(estimate,4:1,xlim=range(c(lo,hi,null)),ylim=c(.6,4.4),yaxt='n',ylab='',xlab=xlab,pch=19,col='#2171b5',bty='l')
  abline(v=null,lty=2,col='gray60');segments(lo,4:1,hi,4:1,col='#2171b5',lwd=1.4)
  axis(2,at=4:1,labels=labels,las=1,tick=FALSE,cex.axis=.85)
  title(main=paste0(LETTERS[panel],'. ',if(model=='baseline')'Pooled temperature' else 'Within-female temperature'),adj=0,cex.main=1)
}
dev.off()
