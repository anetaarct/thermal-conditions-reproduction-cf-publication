# Draw regenerated predictions without overwriting publication figure files.
out<-'figures/reproduced';dir.create(out,recursive=TRUE,showWarnings=FALSE)
draw_panels<-function(d,groups,xcolumn,stem) {
 panels<-split(d,interaction(d[groups],drop=TRUE))
 pdf(file.path(out,paste0(stem,'.pdf')),width=10,height=max(6,ceiling(length(panels)/2)*3.5))
 on.exit(dev.off())
 par(mfrow=c(ceiling(length(panels)/2),2),mar=c(4,4,3,1))
 for(name in names(panels)) {
  z<-panels[[name]];x<-suppressWarnings(as.numeric(z[[xcolumn]]))
  discrete<-anyNA(x)
  if(discrete) x<-seq_len(nrow(z))
  ylabel<-if('component' %in% names(z)) as.character(z$component[1]) else 'Expected outcome'
  plot(range(x),range(z$lower,z$upper),type='n',xlab=xcolumn,ylab=ylabel,main=name,bty='l',xaxt=if(discrete)'n' else 's')
  if(discrete) {axis(1,at=x,labels=z[[xcolumn]]);arrows(x,z$lower,x,z$upper,angle=90,code=3,length=.05)} else {
   polygon(c(x,rev(x)),c(z$lower,rev(z$upper)),col=adjustcolor('#0072B2',.15),border=NA)
  }
  lines(x,z$estimate,col='#0072B2',lwd=2);points(x,z$estimate,col='#0072B2',pch=16,cex=.4)
 }
}
d<-read.csv('tables/reproduced/baseline_predictions/curves.csv')
draw_panels(d,c('trait','predictor','component'),'x','baseline_temperature_date_age')
d<-read.csv('tables/within_between_20261026_reproduced/prediction_curves.csv')
draw_panels(d,c('trait','predictor','component'),'x','within_between_temperature')
d<-read.csv('tables/reproduced/joint_calendar_year/curves.csv')
draw_panels(d,c('trait','model','component'),'year','temperature_adjusted_calendar_year')
temporal<-list()
for(trait in c('clutch','fledging')) for(component in if(trait=='clutch') c('mean','sd') else c('success','phi')) {
 d<-read.csv(paste0('tables/reproduced/temporal/',trait,'_year_',component,'.csv'))
 temporal[[paste(trait,component)]]<-data.frame(trait,component,x=d$x,estimate=d$estimate,lower=d$CI_low,upper=d$CI_high)
}
draw_panels(do.call(rbind,temporal),c('trait','component'),'x','calendar_year_without_temperature')
