library("plotrix")
png('normal_interaction.png')
plot.new()
plot.window(xlim=c(-7,7),ylim=c(-9,5))
draw.circle(-2,0,3,lty=1,lwd=1.5)
draw.circle(3,0,4,lty=1,lwd=1.5)
lines(c(-2,3),c(0,0),lty=2)
lines(c(-2,0),c(0,0),lty=2)

arrows(-2,0,-1.5,0,length = 0.05, angle = 30,lwd=1.75,lty=1)

points(-1,0,pch=16,cex=1.5)
text(-2,0.35,expression('n'['ij']),cex=1.35)
points(1,0,pch=16,cex=1.5)

lines(c(-1,-4),c(0,-6))
lines(c(1,4),c(0,-6))
points(-4,-6,pch=16,cex=1.5)
points(4,-6,pch=16,cex=1.5)

lines(c(-2,-2),c(-7,-5))
lines(c(2,2),c(-7,-5))

lines(c(2,1),c(-5,-5))
lines(c(-2,-1),c(-5,-5))
lines(c(-1,-0.6,0,0.6,1),c(-5,-4.5,-5.25,-4.5,-5))
lines(c(-.5,-.5),c(-6.5,-7.5))
lines(c(-0.5,0.2),c(-6.5,-6.5))
lines(c(-0.5,0.2),c(-7.5,-7.5))
points(0,-7,pch="|",lwd=1.5)
lines(c(0,2),c(-7,-7))
lines(c(-0.5,-2),c(-7,-7))
lines(c(4,2),c(-6,-6))
lines(c(-4,-2),c(-6,-6))

text(-3,-2,'i',cex=1.75)
text(4,-2,'j',cex=1.75)

text(0,-4,expression('k'['n']),cex=1.75)
text(0,-8,expression(gamma['n']),cex=1.75)

dev.off()

png("tangential_interaction.png")

plot.new()
plot.window(xlim=c(-1.5,6),ylim=c(-5,5))
draw.arc(-2,0,3,-1.3,1.3,lwd=1.5)
draw.arc(3,0,4,-2,-4.25,lwd=1.5)
points(-0.76,2,pch=16,cex=1.5)
points(0.2,-3,pch=16,cex=1.5)

lines(c(-0.76,3),c(2,4))
lines(c(0.2,2.85),c(-3,-4))

lines(c(2.85,2.85),c(-2,-4))
lines(c(2.95,2.95),c(-1,-3))

lines(c(2.95,2.875),c(-1,-0.85))

lines(c(2.875,2.875),c(-0.85,0.5))
lines(c(3,3),c(4,2))

lines(c(2.875,3.15,2.65,3.15,2.65,3.15,2.65,3.15,3),c(0.5,0.75,0.85,1.05,1.25,1.45,1.65,1.85,2))

text(-0.75,-3.5,'i',cex=1.75)
text(0.75,-4.25,'j',cex=1.75)

text(3.5,1,expression('k'['t']),cex=1.75)
text(2.5,-2.5,expression(gamma['n']),cex=1.75)

dev.off()
