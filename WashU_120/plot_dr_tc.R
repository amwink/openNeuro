# source("plot_dr_tc.R")

library(psych)
library(ggplot2)

# initialise
rm(list=ls());
graphics.off();

dr1_results=array(0,c(100,10,120));

# load files
for (n in 0:114) {

    fname=sprintf("dual_regression_sex_binage115_twosided/dr_stage1_subject00%03d.txt",n);
    tmp=matrix(scan(fname),nrow=10);
    dr1_results[,,n+1]=t(tmp);
    
    }
# dr1_results[,,120]

# show group mean timecourses
#for (n in 1:10) {
#    plot(apply(dr1_results[,n,],1,mean,trim=0.3), type="l")
#    }
#

# mean time courses
mean_ic_tc=apply(dr1_results,c(1,2),mean,trim=0.1);

# Smith RSN names
colnames(mean_ic_tc)=c("medial visual","occippital visual","lateral visual","default mode","cerebellum","sensorimotor","auditory","executive control","frontoparietal R","frontoparietal L")

numic=dim(dr1_results)[2];
types=1:(numic+2);
types=types[-8];
types=types[-7];
par(xpd = T, mar = par()$mar + c(0,0,0,9));
matplot(mean_ic_tc, type="l", lty=, col=types, xlab="window position (TR)", ylab="ECM (a.u.)");
legend(par("usr")[2]+5,par("usr")[4], legend=colnames(mean_ic_tc), lty=types, col=types);

no_cer=T;
# remove cerebellum
if (no_cer) {
    dr1_results=dr1_results[,-5,];
    mean_ic_tc=apply(dr1_results,c(1,2),mean,trim=0.1);
    colnames(mean_ic_tc)=c("medial visual","occippital visual","lateral visual","default mode","sensorimotor","auditory","executive control","frontoparietal R","frontoparietal L");
}

par(xpd = T, mar = par()$mar + c(0,0,0,9));
matplot(mean_ic_tc, type="l", lty=types, col=types, xlab="window position (TR)", ylab="ECM (a.u.)");
legend(par("usr")[2]+5,par("usr")[4], legend=colnames(mean_ic_tc), lty=types, col=types);

design=read.table("design_dr.txt",header=F);
colnames(design)=c("filename","male/female","young/old");

maledr1=dr1_results[,,design$"male/female"==-1];
femaledr1=dr1_results[,,design$"male/female"==1];
youngdr1=dr1_results[,,design$"young/old"==-1];
olddr1=dr1_results[,,design$"young/old"==1];

youngmaledr1=dr1_results[,,design$"young/old"==-1 & design$"male/female"==-1];
youngfemaledr1=dr1_results[,,design$"young/old"==-1 & design$"male/female"==1];
oldfemaledr1=dr1_results[,,design$"young/old"==1 & design$"male/female"==1];
oldmaledr1=dr1_results[,,design$"young/old"==1 & design$"male/female"==-1];

# see https://www.personality-project.org/r/html/winsor.html
wtrim=.1;

mean_ym_tc=as.data.frame(apply(youngmaledr1,c(1,2),winsor.mean,trim=wtrim));colnames(mean_ym_tc)=colnames(mean_ic_tc);
mean_yf_tc=as.data.frame(apply(youngfemaledr1,c(1,2),winsor.mean,trim=wtrim));colnames(mean_yf_tc)=colnames(mean_ic_tc);
mean_of_tc=as.data.frame(apply(oldfemaledr1,c(1,2),winsor.mean,trim=wtrim));colnames(mean_of_tc)=colnames(mean_ic_tc);
mean_om_tc=as.data.frame(apply(oldmaledr1,c(1,2),winsor.mean,trim=wtrim));colnames(mean_om_tc)=colnames(mean_ic_tc);

sd_ym_tc=as.data.frame(apply(youngmaledr1,c(1,2),winsor.sd,trim=wtrim));colnames(sd_ym_tc)=colnames(mean_ic_tc);
sd_yf_tc=as.data.frame(apply(youngfemaledr1,c(1,2),winsor.sd,trim=wtrim));colnames(sd_yf_tc)=colnames(mean_ic_tc);
sd_of_tc=as.data.frame(apply(oldfemaledr1,c(1,2),winsor.sd,trim=wtrim));colnames(sd_of_tc)=colnames(mean_ic_tc);
sd_om_tc=as.data.frame(apply(oldmaledr1,c(1,2),winsor.sd,trim=wtrim));colnames(sd_om_tc)=colnames(mean_ic_tc);

par(xpd = T, mar = par()$mar + c(0,0,0,9));
matplot(mean_ym_tc, type="l", lty=types, col=types, xlab="window position (TR)", ylab="ECM (a.u.)");
legend(par("usr")[2]+5,par("usr")[4], legend=colnames(mean_ic_tc), lty=types, col=types);

# see https://stackoverflow.com/questions/29743503
# p<-ggplot(data=data, aes(x=interval, y=OR, colour=Drug)) + geom_point() + geom_line()
# p<-p+geom_ribbon(aes(ymin=data$lower, ymax=data$upper), linetype=2, alpha=0.1)
# see https://stackoverflow.com/questions/19921842
sdfac=5;

# young males
plot.new();
pl=ggplot() +
    ggtitle ("young males", subtitle="mean EC time seris within RSN") +
    labs(x = "start position (TR)", y = "EC map (a.u.)" ) +
    #labs(x = "start position (TR)", y = "EC map (a.u.)", caption = "https://github.com/amwink/bias/tree/master/matlab/fastECM") +
    #theme( plot.caption = element_text(hjust=1.85, vjust=0, face = "italic", size = 6, angle=0) ) + 
    geom_line(aes(x=1:100,y=mean_ym_tc[,1],color=colnames(mean_ic_tc)[1])) + geom_line() +
    geom_ribbon(aes(x=1:100,ymin=mean_ym_tc[,1]-sd_ym_tc[,1]/sdfac,ymax=mean_ym_tc[,1]+sd_ym_tc[,1]/sdfac), linetype=2, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_ym_tc[,2],color=colnames(mean_ic_tc)[2])) +
    geom_ribbon(aes(x=1:100,ymin=mean_ym_tc[,2]-sd_ym_tc[,2]/sdfac,ymax=mean_ym_tc[,2]+sd_ym_tc[,2]/sdfac), linetype=2, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_ym_tc[,3],color=colnames(mean_ic_tc)[3])) +
    geom_ribbon(aes(x=1:100,ymin=mean_ym_tc[,3]-sd_ym_tc[,3]/sdfac,ymax=mean_ym_tc[,3]+sd_ym_tc[,3]/sdfac), linetype=3, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_ym_tc[,4],color=colnames(mean_ic_tc)[4])) +
    geom_ribbon(aes(x=1:100,ymin=mean_ym_tc[,4]-sd_ym_tc[,4]/sdfac,ymax=mean_ym_tc[,4]+sd_ym_tc[,4]/sdfac), linetype=4, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_ym_tc[,5],color=colnames(mean_ic_tc)[5])) +
    geom_ribbon(aes(x=1:100,ymin=mean_ym_tc[,5]-sd_ym_tc[,5]/sdfac,ymax=mean_ym_tc[,5]+sd_ym_tc[,5]/sdfac), linetype=5, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_ym_tc[,6],color=colnames(mean_ic_tc)[6])) +
    geom_ribbon(aes(x=1:100,ymin=mean_ym_tc[,6]-sd_ym_tc[,6]/sdfac,ymax=mean_ym_tc[,6]+sd_ym_tc[,6]/sdfac), linetype=6, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_ym_tc[,7],color=colnames(mean_ic_tc)[7])) +
    geom_ribbon(aes(x=1:100,ymin=mean_ym_tc[,7]-sd_ym_tc[,7]/sdfac,ymax=mean_ym_tc[,7]+sd_ym_tc[,7]/sdfac), linetype=7, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_ym_tc[,8],color=colnames(mean_ic_tc)[8])) +
    geom_ribbon(aes(x=1:100,ymin=mean_ym_tc[,8]-sd_ym_tc[,8]/sdfac,ymax=mean_ym_tc[,8]+sd_ym_tc[,8]/sdfac), linetype=8, alpha=0.1) +
    geom_line(aes(x=1:100,y=mean_ym_tc[,9],color=colnames(mean_ic_tc)[9])) +
    geom_ribbon(aes(x=1:100,ymin=mean_ym_tc[,9]-sd_ym_tc[,9]/sdfac,ymax=mean_ym_tc[,9]+sd_ym_tc[,9]/sdfac), linetype=9, alpha=0.1) +
    ylim(-2.5,7.5) + labs(colour = "resting networks");
pl;
ggsave("mean_ym_tc.pdf")

# young females
plot.new();
pl=ggplot() +
    ggtitle ("young females", subtitle="mean EC time seris within RSN") +
    labs(x = "start position (TR)", y = "EC map (a.u.)" ) +
    #labs(x = "start position (TR)", y = "EC map (a.u.)", caption = "https://github.com/amwink/bias/tree/master/matlab/fastECM") +
    #theme( plot.caption = element_text(hjust=1.85, vjust=0, face = "italic", size = 6, angle=0) ) + 
    geom_line(aes(x=1:100,y=mean_yf_tc[,1],color=colnames(mean_ic_tc)[1])) + geom_line() +
    geom_ribbon(aes(x=1:100,ymin=mean_yf_tc[,1]-sd_yf_tc[,1]/sdfac,ymax=mean_yf_tc[,1]+sd_yf_tc[,1]/sdfac), linetype=2, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_yf_tc[,2],color=colnames(mean_ic_tc)[2])) +
    geom_ribbon(aes(x=1:100,ymin=mean_yf_tc[,2]-sd_yf_tc[,2]/sdfac,ymax=mean_yf_tc[,2]+sd_yf_tc[,2]/sdfac), linetype=2, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_yf_tc[,3],color=colnames(mean_ic_tc)[3])) +
    geom_ribbon(aes(x=1:100,ymin=mean_yf_tc[,3]-sd_yf_tc[,3]/sdfac,ymax=mean_yf_tc[,3]+sd_yf_tc[,3]/sdfac), linetype=3, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_yf_tc[,4],color=colnames(mean_ic_tc)[4])) +
    geom_ribbon(aes(x=1:100,ymin=mean_yf_tc[,4]-sd_yf_tc[,4]/sdfac,ymax=mean_yf_tc[,4]+sd_yf_tc[,4]/sdfac), linetype=4, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_yf_tc[,5],color=colnames(mean_ic_tc)[5])) +
    geom_ribbon(aes(x=1:100,ymin=mean_yf_tc[,5]-sd_yf_tc[,5]/sdfac,ymax=mean_yf_tc[,5]+sd_yf_tc[,5]/sdfac), linetype=5, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_yf_tc[,6],color=colnames(mean_ic_tc)[6])) +
    geom_ribbon(aes(x=1:100,ymin=mean_yf_tc[,6]-sd_yf_tc[,6]/sdfac,ymax=mean_yf_tc[,6]+sd_yf_tc[,6]/sdfac), linetype=6, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_yf_tc[,7],color=colnames(mean_ic_tc)[7])) +
    geom_ribbon(aes(x=1:100,ymin=mean_yf_tc[,7]-sd_yf_tc[,7]/sdfac,ymax=mean_yf_tc[,7]+sd_yf_tc[,7]/sdfac), linetype=7, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_yf_tc[,8],color=colnames(mean_ic_tc)[8])) +
    geom_ribbon(aes(x=1:100,ymin=mean_yf_tc[,8]-sd_yf_tc[,8]/sdfac,ymax=mean_yf_tc[,8]+sd_yf_tc[,8]/sdfac), linetype=8, alpha=0.1) +
    geom_line(aes(x=1:100,y=mean_yf_tc[,9],color=colnames(mean_ic_tc)[9])) +
    geom_ribbon(aes(x=1:100,ymin=mean_yf_tc[,9]-sd_yf_tc[,9]/sdfac,ymax=mean_yf_tc[,9]+sd_yf_tc[,9]/sdfac), linetype=9, alpha=0.1) +
    ylim(-2.5,7.5) + labs(colour = "resting networks");
pl;
ggsave("mean_yf_tc.pdf")

# old females
plot.new();
pl=ggplot() +
    ggtitle ("old females", subtitle="mean EC time seris within RSN") +
    labs(x = "start position (TR)", y = "EC map (a.u.)" ) +
    #labs(x = "start position (TR)", y = "EC map (a.u.)", caption = "https://github.com/amwink/bias/tree/master/matlab/fastECM") +
    #theme( plot.caption = element_text(hjust=1.85, vjust=0, face = "italic", size = 6, angle=0) ) + 
    geom_line(aes(x=1:100,y=mean_of_tc[,1],color=colnames(mean_ic_tc)[1])) + geom_line() +
    geom_ribbon(aes(x=1:100,ymin=mean_of_tc[,1]-sd_of_tc[,1]/sdfac,ymax=mean_of_tc[,1]+sd_of_tc[,1]/sdfac), linetype=2, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_of_tc[,2],color=colnames(mean_ic_tc)[2])) +
    geom_ribbon(aes(x=1:100,ymin=mean_of_tc[,2]-sd_of_tc[,2]/sdfac,ymax=mean_of_tc[,2]+sd_of_tc[,2]/sdfac), linetype=2, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_of_tc[,3],color=colnames(mean_ic_tc)[3])) +
    geom_ribbon(aes(x=1:100,ymin=mean_of_tc[,3]-sd_of_tc[,3]/sdfac,ymax=mean_of_tc[,3]+sd_of_tc[,3]/sdfac), linetype=3, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_of_tc[,4],color=colnames(mean_ic_tc)[4])) +
    geom_ribbon(aes(x=1:100,ymin=mean_of_tc[,4]-sd_of_tc[,4]/sdfac,ymax=mean_of_tc[,4]+sd_of_tc[,4]/sdfac), linetype=4, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_of_tc[,5],color=colnames(mean_ic_tc)[5])) +
    geom_ribbon(aes(x=1:100,ymin=mean_of_tc[,5]-sd_of_tc[,5]/sdfac,ymax=mean_of_tc[,5]+sd_of_tc[,5]/sdfac), linetype=5, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_of_tc[,6],color=colnames(mean_ic_tc)[6])) +
    geom_ribbon(aes(x=1:100,ymin=mean_of_tc[,6]-sd_of_tc[,6]/sdfac,ymax=mean_of_tc[,6]+sd_of_tc[,6]/sdfac), linetype=6, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_of_tc[,7],color=colnames(mean_ic_tc)[7])) +
    geom_ribbon(aes(x=1:100,ymin=mean_of_tc[,7]-sd_of_tc[,7]/sdfac,ymax=mean_of_tc[,7]+sd_of_tc[,7]/sdfac), linetype=7, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_of_tc[,8],color=colnames(mean_ic_tc)[8])) +
    geom_ribbon(aes(x=1:100,ymin=mean_of_tc[,8]-sd_of_tc[,8]/sdfac,ymax=mean_of_tc[,8]+sd_of_tc[,8]/sdfac), linetype=8, alpha=0.1) +
    geom_line(aes(x=1:100,y=mean_of_tc[,9],color=colnames(mean_ic_tc)[9])) +
    geom_ribbon(aes(x=1:100,ymin=mean_of_tc[,9]-sd_of_tc[,9]/sdfac,ymax=mean_of_tc[,9]+sd_of_tc[,9]/sdfac), linetype=9, alpha=0.1) +
    ylim(-2.5,7.5) + labs(colour = "resting networks")
pl;
ggsave("mean_of_tc.pdf")

# old males
plot.new();
pl=ggplot() +
    ggtitle ("old males", subtitle="mean EC time seris within RSN") +
    labs(x = "start position (TR)", y = "EC map (a.u.)" ) +
    #labs(x = "start position (TR)", y = "EC map (a.u.)", caption = "https://github.com/amwink/bias/tree/master/matlab/fastECM") +
    #theme( plot.caption = element_text(hjust=1.85, vjust=0, face = "italic", size = 6, angle=0) ) + 
    geom_line(aes(x=1:100,y=mean_om_tc[,1],color=colnames(mean_ic_tc)[1])) + geom_line() +
    geom_ribbon(aes(x=1:100,ymin=mean_om_tc[,1]-sd_om_tc[,1]/sdfac,ymax=mean_om_tc[,1]+sd_om_tc[,1]/sdfac), linetype=2, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_om_tc[,2],color=colnames(mean_ic_tc)[2])) +
    geom_ribbon(aes(x=1:100,ymin=mean_om_tc[,2]-sd_om_tc[,2]/sdfac,ymax=mean_om_tc[,2]+sd_om_tc[,2]/sdfac), linetype=2, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_om_tc[,3],color=colnames(mean_ic_tc)[3])) +
    geom_ribbon(aes(x=1:100,ymin=mean_om_tc[,3]-sd_om_tc[,3]/sdfac,ymax=mean_om_tc[,3]+sd_om_tc[,3]/sdfac), linetype=3, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_om_tc[,4],color=colnames(mean_ic_tc)[4])) +
    geom_ribbon(aes(x=1:100,ymin=mean_om_tc[,4]-sd_om_tc[,4]/sdfac,ymax=mean_om_tc[,4]+sd_om_tc[,4]/sdfac), linetype=4, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_om_tc[,5],color=colnames(mean_ic_tc)[5])) +
    geom_ribbon(aes(x=1:100,ymin=mean_om_tc[,5]-sd_om_tc[,5]/sdfac,ymax=mean_om_tc[,5]+sd_om_tc[,5]/sdfac), linetype=5, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_om_tc[,6],color=colnames(mean_ic_tc)[6])) +
    geom_ribbon(aes(x=1:100,ymin=mean_om_tc[,6]-sd_om_tc[,6]/sdfac,ymax=mean_om_tc[,6]+sd_om_tc[,6]/sdfac), linetype=6, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_om_tc[,7],color=colnames(mean_ic_tc)[7])) +
    geom_ribbon(aes(x=1:100,ymin=mean_om_tc[,7]-sd_om_tc[,7]/sdfac,ymax=mean_om_tc[,7]+sd_om_tc[,7]/sdfac), linetype=7, alpha=0.1) + 
    geom_line(aes(x=1:100,y=mean_om_tc[,8],color=colnames(mean_ic_tc)[8])) +
    geom_ribbon(aes(x=1:100,ymin=mean_om_tc[,8]-sd_om_tc[,8]/sdfac,ymax=mean_om_tc[,8]+sd_om_tc[,8]/sdfac), linetype=8, alpha=0.1) +
    geom_line(aes(x=1:100,y=mean_om_tc[,9],color=colnames(mean_ic_tc)[9])) +
    geom_ribbon(aes(x=1:100,ymin=mean_om_tc[,9]-sd_om_tc[,9]/sdfac,ymax=mean_om_tc[,9]+sd_om_tc[,9]/sdfac), linetype=9, alpha=0.1) +
    ylim(-2.5,7.5) + labs(colour = "resting networks")
pl;
ggsave("mean_om_tc.pdf")

# now run in the shell:
# $ doall.sh gv mean*pdf
# $ for f in mean*pdf; do cm="pdftoppm -r 600 -png -singlefile $f ${f//.pdf/}"; echo $cm;$cm;done
# $ convert +append mean_ym_tc.png mean_yf_tc.png young.png
# $ convert +append mean_om_tc.png mean_of_tc.png old.png
# $ convert +append young.png old.png young_old.png



