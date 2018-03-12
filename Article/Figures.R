# Figures article
# 12 March 2018

#### 2 herbivores one plant ####
# parametres
S1=5;
a=0.2;
qv1r1=1.5;
qv1r2=1.8;

# herbivores 
m1=0.5;
qh1r1=3;
qh1r2=6;
g11=1.5;

m2=0.5;
qh2r1=6;
qh2r2=2.5;
g21=1.5;

# ZNGI
penteR1=m1/g11*qh1r1
penteR2=m1/g11*qh1r2
xseq=seq(penteR1,20,length.out=20)
yseq=seq(penteR2,20,length.out=20)
R1=rep(penteR1,length(xseq))
R2=rep(penteR2,length(yseq))

penteR1b=m2/g21*qh2r1
penteR2b=m2/g21*qh2r2
xseqb=seq(penteR1b,20,length.out=20)
yseqb=seq(penteR2b,20,length.out=20)
R1b=rep(penteR1b,length(xseqb))
R2b=rep(penteR2b,length(yseqb))

# vector
pentvect=qv1r2/qv1r1
supR1=S1*qv1r1
supR2=S1*qv1r2
ordo=supR2-pentvect*supR1
x0=supR1#11
y0=pentvect*x0+ordo
x1=0.75
y1=pentvect*x1+ordo
cross=(penteR2-ordo)/pentvect

# plot
setwd("~/Desktop/Article_Tanguy/Article")

pdf("TwoHerbivoresOnePlant.pdf",width=12,height=10)
par(mar=c(5,5,4,2))
plot(R2~xseq,type='l',xlim=c(0,10),ylim=c(0,10),xlab="",ylab="",lwd=3,axes=F,xaxs='i',yaxs='i')
axis(1,lwd.ticks=0,at=c(0,5,10),label=c("","",""),cex.axis=2)
mtext(side = 1, text = expression(paste("R"[1])), line = 3,cex=4)
#axis(1,at=c(0,5,10),labels=c("",expression(paste("R"[1])),""),cex.axis=2,lwd.ticks=0)
axis(2,las=1,at=c(0,5,10),labels=c("",expression(paste("R"[2])),""),cex.axis=4,lwd.ticks=0)
lines(yseq~R1,lwd=3)
lines(R2b~xseqb,lty=2,lwd=3.5)
lines(yseqb~R1b,lty=2,lwd=3.5)
arrows(x0,y0,x1,y1,lwd=4,col='grey70')
points(supR2~supR1,lwd=3)
points(penteR2~cross,lwd=3)
text(0.5,7,"I",cex=3.5)
text(1.5,7,"II",cex=3.5)
text(8,7,"III",cex=3.5)
text(8,1.5,"IV",cex=3.5)
text(8,0.4,"V",cex=3.5)
text(8.6,supR2,"Supply",cex=3.5)
text(cross,2.5,"E",cex=3.5)
box()
dev.off()

#### coexistence ####
# plants
a1 = 0.01;
a2 = 0.01;
Qv1r1 = 0.2;
Qv1r2 = 0.08;
Qv2r1 = 0.1;
Qv2r2 = 0.25;

# herbivore 1
g11 = 0.15;
g12 = 0.11; # 0.1
Qh1r1 = 0.3;
Qh1r2 = 0.2;
m1 = 0.1;

# herbivore 2
g21 = 0.08;
g22 = 0.15;
Qh2r1 = 0.2;
Qh2r2 = 0.3;
m2 = 0.1;

# lim
limx=1.2
limy=1.2

# equilibrium point
V1barre=(m1*g22*Qh1r1*Qv2r2-g12*m2*Qh2r2*Qv2r1)/(g11*g22*Qv1r1*Qv2r2-g12*g21*Qv2r1*Qv1r2)
V2barre=(m2*Qh2r2-g21*V1barre*Qv1r2)/(g22*Qv2r2)
R1equ=V1barre*Qv1r1+V2barre*Qv2r1
R2equ=V1barre*Qv1r2+V2barre*Qv2r2

# feasibility cone
V1bound=Qv1r2/Qv1r1
V2bound=Qv2r2/Qv2r1

# ZNGI herbivore 1
R1stx=Qh1r1*m1/g11;
R1sty=Qh1r1*m1/g12;
R2stx=Qh1r2*m1/g11;
R2sty=Qh1r2*m1/g12;

zngixR1y1=seq(R2stx,limy,length.out=20)
zngixR1x1=rep(R1stx,length(zngixR1y1))
zngixR2x1=seq(R1stx,limy,length.out=20)
zngixR2y1=rep(R2stx,length(zngixR2x1))

zngiyR1y1=seq(R2sty,limy,length.out=20)
zngiyR1x1=rep(R1sty,length(zngiyR1y1))
zngiyR2x1=seq(R1sty,limy,length.out=20)
zngiyR2y1=rep(R2sty,length(zngiyR2x1))

# ZNGI herbivore 2
R1stx=Qh2r1*m2/g21;
R1sty=Qh2r1*m2/g22;
R2stx=Qh2r2*m2/g21;
R2sty=Qh2r2*m2/g22;

zngixR1y2=seq(R2stx,limy,length.out=20)
zngixR1x2=rep(R1stx,length(zngixR1y2))
zngixR2x2=seq(R1stx,limy,length.out=20)
zngixR2y2=rep(R2stx,length(zngixR2x2))

zngiyR1y2=seq(R2sty,limy,length.out=20)
zngiyR1x2=rep(R1sty,length(zngiyR1y2))
zngiyR2x2=seq(R1sty,limy,length.out=20)
zngiyR2y2=rep(R2sty,length(zngiyR2x2))

# polygon vertices
predx1=limy/V2bound
predy1=limx*V1bound

xpolygon=c(0,predx1,limx,limx)
ypolygon=c(0,limy,limy,predy1)

# vectors
pente1=(g11*V1barre*Qv1r2+g12*V2barre*Qv2r2)/(g11*V1barre*Qv1r1+g12*V2barre*Qv2r1)
pente2=(g21*V1barre*Qv1r2+g22*V2barre*Qv2r2)/(g21*V1barre*Qv1r1+g22*V2barre*Qv2r1)

# Vxequ=(m2*Qh22*gy1*Qvy1-gy2*Qvy2*m1*Qh11)/(gx2*Qvx2*gy1*Qvy1-gy2*Qvy2*gx1*Qvx1)
# Vyequ=(m1*Qh11-gx1*Vxequ*Qvx1)/(gy1*Qvy1)
# R1equ=Vxequ*Qvx1+Vyequ*Qvy1
# R2equ=Vxequ*Qvx2+Vyequ*Qvy2
# 
# S1equ1=seq(1,10)
# S1equ2=seq(1,10)
# 
# S2equ1=(S1equ1-ax)*(gy1*Vyequ)/(gx1*Vxequ)+ay*Vyequ
# S2equ2=(S1equ1-ax)*(gy2*Vyequ)/(gx2*Vxequ)+ay*Vyequ
# 
# ch1=(S1equ1*Qvx2+S2equ1*Qvy2)/(S1equ1*Qvx1+S2equ1*Qvy1)
# ch2=(S1equ2*Qvx2+S2equ2*Qvy2)/(S1equ2*Qvx1+S2equ2*Qvy1)
# 
# pente1=mean(ch1)
# pente2=mean(ch2)
# 
# ordo1=R2equ-pente1*R1equ
# ordo2=R2equ-pente2*R1equ

# plot
setwd("~/Desktop/Article_Tanguy/Article")

pdf("Coexistence2Herbivores.pdf",width=12,height=10)
par(mar=c(5,5,4,2))
plot(0,0,type='n',xlim=c(0,limx),ylim=c(0,limy),xlab='',ylab='',xaxs='i',yaxs='i',lwd=2,axes=F)
axis(1,lwd.ticks=0,at=c(0,0.65,2),label=c("","",""),cex.axis=2)
mtext(side = 1, text = expression(paste("R"[1])), line = 3,cex=4)
axis(2,lwd.ticks=0,at=c(0,0.65,2),label=c("",expression(paste("R"[2])),""),las=1,cex.axis=4)
polygon(xpolygon,ypolygon,col='grey87',lty=0)
lines(zngixR1y1~zngixR1x1,lwd=3.5)
lines(zngixR2y1~zngixR2x1,lwd=3.5)
lines(zngiyR1y1~zngiyR1x1,lwd=3.5)
lines(zngiyR2y1~zngiyR2x1,lwd=3.5)
lines(zngixR1y2~zngixR1x2,lty=2,lwd=3.5)
lines(zngixR2y2~zngixR2x2,lty=2,lwd=3.5)
lines(zngiyR1y2~zngiyR1x2,lty=2,lwd=3.5)
lines(zngiyR2y2~zngiyR2x2,lty=2,lwd=3.5)

ordo1=R2equ-pente1*R1equ
ordo2=R2equ-pente2*R1equ
x0=R1equ-0.125
y0=pente1*x0+ordo1
x1=2
y1=x1*pente1+ordo1
arrows(x0,y0,x1,y1,code=1,lwd=4,col='grey45',length=0.2)
x0=R1equ-0.125#35
y0=pente2*x0+ordo2
x1=x0+0.01
y1=x1*pente2+ordo2
arrows(x0,y0,x1,y1,code=1,lwd=4,col='grey45',length=0.2)
x2=2
y2=x2*pente2+ordo2
lines(c(y1,y2)~c(x1,x2),lty=2,lwd=4,col='grey45')

text(0.07,1,"I",cex=3.5)
text(0.17,1,"II",cex=3.5)
text(0.5,1,"III",cex=3.5)
text(0.9,1,"IV",cex=3.5)
text(0.9,0.6,"V",cex=3.5)
text(0.9,0.16,"VI",cex=2.5)
text(0.9,0.08,'VII',cex=3.5)
box()

dev.off()

# x0=R1equ-0.125
# y0=pente1*x0+ordo1
# x1=2
# y1=x1*pente1+ordo1
# arrows(x0,y0,x1,y1,code=1,lwd=2.5,col='grey45',length=0.2)
# x0=R1equ-0.135
# y0=pente2*x0+ordo2
# x1=x0+0.01
# y1=x1*pente2+ordo2
# arrows(x0,y0,x1,y1,code=1,lwd=2.5,col='grey45',length=0.2)
# x2=2
# y2=x2*pente2+ordo2
# lines(c(y1,y2)~c(x1,x2),lty=2,lwd=2.5,col='grey45')


#### exclusion ####

# plants
a1 = 0.01;
a2 = 0.01;
Qv1r1 = 0.2;
Qv1r2 = 0.08;
Qv2r1 = 0.1;
Qv2r2 = 0.25;

# herbivore 1
g11 = 0.11;
g12 = 0.15; # 0.1
Qh1r1 = 0.3;
Qh1r2 = 0.2;
m1 = 0.1;

# herbivore 2
g21 = 0.15;
g22 = 0.08;
Qh2r1 = 0.2;
Qh2r2 = 0.3;
m2 = 0.1;

# lim
limx=1.2
limy=1.2

# equilibrium point
V1barre=(m1*g22*Qh1r1*Qv2r2-g12*m2*Qh2r2*Qv2r1)/(g11*g22*Qv1r1*Qv2r2-g12*g21*Qv2r1*Qv1r2)
V2barre=(m2*Qh2r2-g21*V1barre*Qv1r2)/(g22*Qv2r2)
R1equ=V1barre*Qv1r1+V2barre*Qv2r1
R2equ=V1barre*Qv1r2+V2barre*Qv2r2

# feasibility cone
V1bound=Qv1r2/Qv1r1
V2bound=Qv2r2/Qv2r1

# ZNGI herbivore 1
R1stx=Qh1r1*m1/g11;
R1sty=Qh1r1*m1/g12;
R2stx=Qh1r2*m1/g11;
R2sty=Qh1r2*m1/g12;

zngixR1y1=seq(R2stx,limy,length.out=20)
zngixR1x1=rep(R1stx,length(zngixR1y1))
zngixR2x1=seq(R1stx,limy,length.out=20)
zngixR2y1=rep(R2stx,length(zngixR2x1))

zngiyR1y1=seq(R2sty,limy,length.out=20)
zngiyR1x1=rep(R1sty,length(zngiyR1y1))
zngiyR2x1=seq(R1sty,limy,length.out=20)
zngiyR2y1=rep(R2sty,length(zngiyR2x1))

# ZNGI herbivore 2
R1stx=Qh2r1*m2/g21;
R1sty=Qh2r1*m2/g22;
R2stx=Qh2r2*m2/g21;
R2sty=Qh2r2*m2/g22;

zngixR1y2=seq(R2stx,limy,length.out=20)
zngixR1x2=rep(R1stx,length(zngixR1y2))
zngixR2x2=seq(R1stx,limy,length.out=20)
zngixR2y2=rep(R2stx,length(zngixR2x2))

zngiyR1y2=seq(R2sty,limy,length.out=20)
zngiyR1x2=rep(R1sty,length(zngiyR1y2))
zngiyR2x2=seq(R1sty,limy,length.out=20)
zngiyR2y2=rep(R2sty,length(zngiyR2x2))

# polygon vertices
predx1=limy/V2bound
predy1=limx*V1bound

xpolygon=c(0,predx1,limx,limx)
ypolygon=c(0,limy,limy,predy1)

# vectors
pente1=(g11*V1barre*Qv1r2+g12*V2barre*Qv2r2)/(g11*V1barre*Qv1r1+g12*V2barre*Qv2r1)
pente2=(g21*V1barre*Qv1r2+g22*V2barre*Qv2r2)/(g21*V1barre*Qv1r1+g22*V2barre*Qv2r1)

# plot
setwd("~/Desktop/Article_Tanguy/Article")

pdf("Exclusion2Herbivores.pdf",width=12,height=10)
par(mar=c(5,5,4,2))
plot(0,0,type='n',xlim=c(0,limx),ylim=c(0,limy),xlab='',ylab='',xaxs='i',yaxs='i',lwd=2,axes=F)
axis(1,lwd.ticks=0,at=c(0,0.65,2),label=c("","",""),cex.axis=2)
mtext(side = 1, text = expression(paste("R"[1])), line = 3,cex=4)
#axis(1,lwd.ticks=0,at=c(0,0.65,2),label=c("",expression(paste("R"[1])),""),cex.axis=2)
axis(2,lwd.ticks=0,at=c(0,0.65,2),label=c("",expression(paste("R"[2])),""),las=1,cex.axis=4)
polygon(xpolygon,ypolygon,col='grey87',lty=0)
lines(zngixR1y1~zngixR1x1,lwd=3.5)
lines(zngixR2y1~zngixR2x1,lwd=3.5)
lines(zngiyR1y1~zngiyR1x1,lwd=3.5)
lines(zngiyR2y1~zngiyR2x1,lwd=3.5)
lines(zngixR1y2~zngixR1x2,lty=2,lwd=3.5)
lines(zngixR2y2~zngixR2x2,lty=2,lwd=3.5)
lines(zngiyR1y2~zngiyR1x2,lty=2,lwd=3.5)
lines(zngiyR2y2~zngiyR2x2,lty=2,lwd=3.5)

ordo1=R2equ-pente1*R1equ
ordo2=R2equ-pente2*R1equ
x0=R1equ-0.125
y0=pente1*x0+ordo1
x1=2
y1=x1*pente1+ordo1
arrows(x0,y0,x1,y1,code=1,lwd=4,col='grey45',length=0.2)
x0=R1equ-0.135
y0=pente2*x0+ordo2
x1=x0+0.01
y1=x1*pente2+ordo2
arrows(x0,y0,x1,y1,code=1,lwd=4,col='grey45',length=0.2)
x2=2
y2=x2*pente2+ordo2
lines(c(y1,y2)~c(x1,x2),lty=2,lwd=4,col='grey45')

text(0.07,1,"I",cex=3.5)
text(0.17,1,"II",cex=3.5)
text(0.5,1,"III",cex=3.5)
text(0.75,1,"IV",cex=3.5)
text(0.75,0.6,"V",cex=3.5)
text(0.75,0.16,"VI",cex=2.5)
text(0.75,0.08,'VII',cex=3.5)
box()

dev.off()

