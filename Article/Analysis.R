# Competition herbivores 
# 9 March 2018
library(deSolve)
# parameters
# plant 1
a1 = 0.01;
Qv1r1 = 0.2;
Qv1r2 = 0.08;

# plant 2
a2 = 0.01;
Qv2r1 = 0.1;
Qv2r2 = 0.25;

# herbivore 1
g11 = 0.15;
g12 = 0.11; # 0.1
Qh1r1 = 0.3;
Qh1r2 = 0.2;
m1 = 0.1;

# herbivore 2
g21 = 0.09;
g22 = 0.15;
Qh2r1 = 0.2;
Qh2r2 = 0.3;
m2 = 0.1;

# supply
S1 = 5;
S2 = 5;
supplyR1=S1*Qv1r1+S2*Qv2r1
supplyR2=S1*Qv1r2+S2*Qv2r2

# initial conditions
Vx0 = 10;
Vy0 = 10;
H10 = 0.1;
H20 = 0.1;

time=10000;
dt=1;
tvoulu=seq(0,time,by=dt);

limx=1.5
limy=1.5

# equilibrium point
V1barre=(m1*g22*Qh1r1*Qv2r2-g12*m2*Qh2r2*Qv2r1)/(g11*g22*Qv1r1*Qv2r2-g12*g21*Qv2r1*Qv1r2)
V2barre=(m2*Qh2r2-g21*V1barre*Qv1r2)/(g22*Qv2r2)

# feasibility cone
V1bound=Qv1r2/Qv1r1
V2bound=Qv2r2/Qv2r1

# ZNGI herbivore 1
R1st1=Qh1r1*m1/g11;
R1st2=Qh1r1*m1/g12;
R2st1=Qh1r2*m1/g11;
R2st2=Qh1r2*m1/g12;

zngixR1y1=seq(R2st1,limy,length.out=20)
zngixR1x1=rep(R1st1,length(zngixR1y1))
zngixR2x1=seq(R1st1,limy,length.out=20)
zngixR2y1=rep(R2st1,length(zngixR2x1))

zngiyR1y1=seq(R2st2,limy,length.out=20)
zngiyR1x1=rep(R1st2,length(zngiyR1y1))
zngiyR2x1=seq(R1st2,limy,length.out=20)
zngiyR2y1=rep(R2st2,length(zngiyR2x1))

# ZNGI herbivore 2
R1st1=Qh2r1*m2/g21;
R1st2=Qh2r1*m2/g22;
R2st1=Qh2r2*m2/g21;
R2st2=Qh2r2*m2/g22;

zngixR1y2=seq(R2st1,limy,length.out=20)
zngixR1x2=rep(R1st1,length(zngixR1y2))
zngixR2x2=seq(R1st1,limy,length.out=20)
zngixR2y2=rep(R2st1,length(zngixR2x2))

zngiyR1y2=seq(R2st2,limy,length.out=20)
zngiyR1x2=rep(R1st2,length(zngiyR1y2))
zngiyR2x2=seq(R1st2,limy,length.out=20)
zngiyR2y2=rep(R2st2,length(zngiyR2x2))

# polygon vertices
predx1=limy/V2bound
predy1=limx*V1bound

xpolygon=c(0,predx1,limx,limx)
ypolygon=c(0,limy,limy,predy1)

#plot
plot(0,0,type='l',xlim=c(0,limx),ylim=c(0,limy),xlab='',ylab='',xaxs='i',yaxs='i',lwd=2,axes=F)
axis(1,lwd.ticks=0,at=c(0,0.65,2),label=c("","",""),cex.axis=2)
mtext(side = 1, text = expression(paste("R"[1])), line = 3,cex=2)
#axis(1,lwd.ticks=0,at=c(0,0.65,2),label=c("",expression(paste("R"[1])),""),cex.axis=2)
axis(2,lwd.ticks=0,at=c(0,0.65,2),label=c("",expression(paste("R"[2])),""),las=1,cex.axis=2)
polygon(xpolygon,ypolygon,col='grey87',lty=0)
lines(zngixR1y1~zngixR1x1,lwd=2)
lines(zngixR2y1~zngixR2x1,lwd=2)
lines(zngiyR1y1~zngiyR1x1,lwd=2)
lines(zngiyR2y1~zngiyR2x1,lwd=2)
lines(zngixR1y2~zngixR1x2,lty=2,lwd=2)
lines(zngixR2y2~zngixR2x2,lty=2,lwd=2)
lines(zngiyR1y2~zngiyR1x2,lty=2,lwd=2)
lines(zngiyR2y2~zngiyR2x2,lty=2,lwd=2)
box()

#abline(0,V1bound,col='grey70',lwd=2)
#abline(0,V2bound,col='grey70',lwd=2)


#### real ZNGI ####
S1vect=seq(1,50)
S2vect=seq(1,50)
# ZNGI R1 H1
R1isoH1vect=matrix(0,nrow=2500,ncol=2)
R1isoH1vect=as.data.frame(R1isoH1vect)
index=1
for (i in 1:50){
  S1=S1vect[i]
  for (j in 1:50){
    S2=S2vect[j]
    apolinom = g12*g11*m1*Qh1r1
    bpolinom = g11*g12*Qv2r1*S2+a2*g11*m1*Qh1r1-g12*g11*S1*Qv1r1+g12*a1*m1*Qh1r1
    cpolinom = a1*g12*Qv2r1*S2-a2*S1*g11*Qv1r1+a2*a1*m1*Qh1r1
    
    delta1= bpolinom^2 - 4*apolinom*cpolinom
    sol1=(-bpolinom+sqrt(delta1))/(2*apolinom)
    sol2=(-bpolinom-sqrt(delta1))/(2*apolinom)
    
    H1barre=sol1
    V2barre=(S1*g11*Qv1r1-a1*m1*Qh1r1-g11*m1*Qh1r1*H1barre)/(a1*g12*Qv2r1+g11*g12*Qv2r1*H1barre)
    V1barre=(m1*Qh1r1-g12*Qv2r1*V2barre)/(g11*Qv1r1)
    
    R1isoH1=V1barre*Qv1r1+V2barre*Qv2r1
    R1isoH1vect[index,1]=R1isoH1
    R1isoH1vect[index,2]=V1barre*Qv1r2+V2barre*Qv2r2
    index=index+1
  }
}

# ZNGI R2 H1
R2isoH1vect=matrix(0,nrow=2500,ncol=2)
R2isoH1vect=as.data.frame(R2isoH1vect)

index=1
for (i in 1:50){
  S1=S1vect[i]
  for (j in 1:50){
    S2=S2vect[j]
    apolinom = g12*g11*m1*Qh1r2
    bpolinom = g11*g12*Qv2r2*S2+a2*g11*m1*Qh1r2-g12*g11*S1*Qv1r2+g12*a1*m1*Qh1r2
    cpolinom = a1*g12*Qv2r2*S2-a2*S1*g11*Qv1r2+a2*a1*m1*Qh1r2
    
    delta1= bpolinom^2 - 4*apolinom*cpolinom
    sol1=(-bpolinom+sqrt(delta1))/(2*apolinom)
    sol2=(-bpolinom-sqrt(delta1))/(2*apolinom)
    
    H1barre=sol1
    V2barre=(S1*g11*Qv1r2-a1*m1*Qh1r2-g11*m1*Qh1r2*H1barre)/(a1*g12*Qv2r2+g11*g12*Qv2r2*H1barre)
    V1barre=(m1*Qh1r2-g12*Qv2r2*V2barre)/(g11*Qv1r2)
    
    R2isoH1=V1barre*Qv1r2+V2barre*Qv2r2
    R2isoH1vect[index,2]=R2isoH1
    R2isoH1vect[index,1]=V1barre*Qv1r1+V2barre*Qv2r1
    index=index+1
  }
}

# plot(R1isoH1vect[,2]~R1isoH1vect[,1],type='l',ylim=c(-1,10),xlim=c(-1,10))
# lines(R2isoH1vect[,2]~R2isoH1vect[,1],col='red')
mod1=lm(R1isoH1vect[,2]~R1isoH1vect[,1])
mod2=lm(R2isoH1vect[,2]~R2isoH1vect[,1])
plot(0,0,type='n',xlim=c(0,4),ylim=c(0,4))
abline(mod1$coef)
abline(mod2$coef)


#### test of vectors ####
# equilibrium point
V1barre=(m1*g22*Qh1r1*Qv2r2-g12*m2*Qh2r2*Qv2r1)/(g11*g22*Qv1r1*Qv2r2-g12*g21*Qv2r1*Qv1r2)
V2barre=(m2*Qh2r2-g21*V1barre*Qv1r2)/(g22*Qv2r2)

# EDO
equadiff=function(t,x,parms){
  res=rep(0,length(x))
  S1=parms[1]
  S2=parms[2]
  Min1=(g11*Qv1r1*x[3]+g12*Qv2r1*x[4])/Qh1r1
  Min2=(g11*Qv1r2*x[3]+g12*Qv2r2*x[4])/Qh1r2
  gh1=min(Min1,Min2)
  res[1]=(gh1-m1)*x[1]
  Min1=(g21*Qv1r1*x[3]+g22*Qv2r1*x[4])/Qh2r1
  Min2=(g21*Qv1r2*x[3]+g22*Qv2r2*x[4])/Qh2r2
  gh2=min(Min1,Min2)
  res[2]=(gh2-m2)*x[2]
  res[3]=S1-a1*x[3]-g11*x[3]*x[1]-g21*x[3]*x[2]
  res[4]=S2-a2*x[4]-g12*x[4]*x[1]-g22*x[3]*x[2]
  return(list(res))
}

suptest=seq(0.6,10,length.out=15)
leng=length(suptest)*length(suptest)
resu=matrix(nrow=leng,ncol=8)
resu=as.data.frame(resu)
names(resu)=c("H1","H2","V1","V2",'R1','R2','SR1','SR2')
count=1
for (i in 1:length(suptest)){
  for (j in 1:length(suptest)){
    S1=suptest[i]
    S2=suptest[j]
    parms=c(S1,S2)
    result=lsoda(c(H10,H20,Vx0,Vy0),tvoulu,equadiff,parms=parms,rtol=1e-12)
    endval=dim(result)[1]-1
    resu[count,1]=result[endval,2]
    resu[count,2]=result[endval,3]
    resu[count,3]=result[endval,4]
    resu[count,4]=result[endval,5]
    resu[count,5]=result[endval,4]*Qv1r1+result[endval,5]*Qv2r1
    resu[count,6]=result[endval,4]*Qv1r2+result[endval,5]*Qv2r2
    resu[count,7]=S1*Qv1r1+S2*Qv2r1
    resu[count,8]=S1*Qv1r2+S2*Qv2r2
    count=count+1
    print(count)
  }
}

thres=1 #0.01
resu2=na.omit(resu)
for (i in 1:dim(resu2)[1]){
  if (resu2[i,1]>thres && resu2[i,2]>thres){
    points(resu2[i,8]~resu2[i,7],col='green')
  }else{
    if (resu2[i,1]<thres && resu2[i,2]>thres){
      points(resu2[i,8]~resu2[i,7],col='blue')
    }else{
      if (resu2[i,1]>thres && resu2[i,2]<thres){
        points(resu2[i,8]~resu2[i,7],col='red')
      }
    }
  }
}

Vxequ=(m2*Qh2r2*g12*Qv2r1-g22*Qv2r2*m1*Qh1r1)/(g21*Qv1r2*g12*Qv2r1-g22*Qv2r2*g11*Qv1r1)
Vyequ=(m1*Qh1r1-g11*Vxequ*Qv1r1)/(g12*Qv2r1)
R1equ=Vxequ*Qv1r1+Vyequ*Qv2r1
R2equ=Vxequ*Qv1r2+Vyequ*Qv2r2

S1equ1=seq(1,10)
S1equ2=seq(1,10)

S2equ1=(S1equ1-a1)*(g12*Vyequ)/(g11*Vxequ)+a2*Vyequ
S2equ2=(S1equ1-a1)*(g22*Vyequ)/(g21*Vxequ)+a2*Vyequ

ch1=(S1equ1*Qv1r2+S2equ1*Qv2r2)/(S1equ1*Qv1r1+S2equ1*Qv2r1)
ch2=(S1equ2*Qv1r2+S2equ2*Qv2r2)/(S1equ2*Qv1r1+S2equ2*Qv2r1)

pente1=mean(ch1)
pente2=mean(ch2)

ordo1=R2equ-pente1*R1equ
ordo2=R2equ-pente2*R1equ
x0=R1equ-0.125
y0=pente1*x0+ordo1
x1=2
y1=x1*pente1+ordo1
arrows(x0,y0,x1,y1,code=1,lwd=2.5,col='grey45',length=0.2)
x0=R1equ-0.135
y0=pente2*x0+ordo2
x1=x0+0.01
y1=x1*pente2+ordo2
arrows(x0,y0,x1,y1,code=1,lwd=2.5,col='grey45',length=0.2)
x2=2
y2=x2*pente2+ordo2
lines(c(y1,y2)~c(x1,x2),lty=2,lwd=2.5,col='grey45')

pente3=(g11*V1barre*Qv1r2+g12*V2barre*Qv2r2)/(g11*V1barre*Qv1r1+g12*V2barre*Qv2r1)
pente4=(g21*V1barre*Qv1r2+g22*V2barre*Qv2r2)/(g21*V1barre*Qv1r1+g22*V2barre*Qv2r1)

ordo1=R2equ-pente3*R1equ
ordo2=R2equ-pente4*R1equ
x0=R1equ-0.125
y0=pente3*x0+ordo1
x1=2
y1=x1*pente3+ordo1
arrows(x0,y0,x1,y1,code=1,lwd=2.5,col='red',length=0.2)
x0=R1equ-0.135
y0=pente4*x0+ordo2
x1=x0+0.01
y1=x1*pente4+ordo2
arrows(x0,y0,x1,y1,code=1,lwd=2.5,col='red',length=0.2)
x2=2
y2=x2*pente4+ordo2
lines(c(y1,y2)~c(x1,x2),lty=2,lwd=2.5,col='red')

