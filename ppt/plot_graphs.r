xcm <- read.csv('xcmVt.csv',header=F)
ycm <- read.csv('ycmVt.csv',header=F)


png("XexpDataR.png")

plot(xcm$V1,xcm$V2,lwd=1.5,pch=15,cex=1.5,xlab=sprintf('Time(s)'),xlim=c(0,0.6),
     ylim=c(0.1,0.35),ylab=sprintf("x/L"),axes=T,ann=T)
lines(xcm$V1,xcm$V2,lwd=1.5)
legend(0.001,0.35,legend=sprintf("Experiment (Mean Value)"),pch=15,lty=1)

box()

dev.off()

png("YexpDataR.png")

plot(ycm$V1,ycm$V2,lwd=1.5,pch=15,cex=1.5,xlab=sprintf('Time(s)'),xlim=c(0,0.6),
     ylim=c(0.02,0.12),axes=T,ann=T ,ylab=sprintf("y/L")) 
lines(ycm$V1,ycm$V2,lwd=1.5)
legend(0.3,0.12,legend=sprintf("Experiment (Mean Value)"),pch=15,lty=1)

dev.off()
