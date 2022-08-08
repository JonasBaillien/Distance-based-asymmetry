### Examples for the proposed asymmetry measures in univariate and multivariate
### setting. Plots for these are also included:
### - univariate gamma distribution
### - normal mixture for multimodal distribution
### - beta distribution
### - rule of thumb counterexample for skewness 
### - non-rooted distribution
### - univariate skew-normal distribution
### - bivariate skew-normal distribution
### - univariate QBA-distribution
### - bivariate QBA (linear combination)
### - non-parametric version for AIS data
### - parametric version for the bivariate examples in the LC paper
### - ToS distribution of Abe & Fujisama



### gamma distribution
alpha=2 # >1
beta=2 # >0

m=(alpha-1)/beta
sm=dgamma(x = m,shape = alpha,rate = beta)

n=10000
x=seq(-10,15,length.out = n)
y=dgamma(x = x,shape = alpha,rate = beta)
ymirror=dgamma(x = 2*m-x,shape = alpha,rate = beta)
gamma=(y-ymirror)/sm
MS=x[which.max(gamma)]-m
gammastar=max(gamma)*sign(MS)

x11()
par(mfrow=c(2,1))
plot(x,y,type="l",lwd=2,xlab="x",ylab="f(x)")
lines(x,ymirror,col=2,lwd=2)
plot(x-m,gamma,type="l",col=1,lwd=2,xlab="s",ylab=expression(gamma[X](s)))

skew=2/sqrt(alpha)



### normal mixture
library(LaplacesDemon)
mu=c(0,5)
sigma=c(1,1.0001)
p=c(0.5,0.5)

n=10000
x=seq(-10,15,length.out = n)
y=dnormm(x = x,p = p,mu = mu,sigma = sigma)

m=x[which.max(y)]
sm=max(y)

ymirror=dnormm(x = 2*m-x,p = p,mu = mu,sigma = sigma)
gamma=(y-ymirror)/sm
MS=x[which.max(gamma)]-m
gammastar=max(gamma)*sign(MS)

x11()
par(mfrow=c(2,1))
plot(x,y,type="l",lwd=2,xlab="x",ylab="f(x)")
lines(x,ymirror,col=2,lwd=2)
plot(x-m,gamma,type="l",col=1,lwd=2,xlab="s",ylab=expression(gamma[X](s)))





### beta distribution
alpha=4 # >1
beta=2 # >1

m=(alpha-1)/(alpha+beta-2)
sm=dbeta(x = m,shape1 =  alpha,shape2 =  beta)

n=100000
x=seq(-1,2,length.out = n)
y=dbeta(x = x,shape1 =  alpha,shape2 =  beta)
y2=dbeta(x = x-(alpha-beta)/(alpha+beta-2),shape1 = beta,shape2 = alpha)
ymirror=dbeta(x = 2*m-x,shape1 =  alpha,shape2 =  beta)
gamma=(y-ymirror)/sm
MS=x[which.max(gamma)]-m
gammastar=max(gamma)*sign(MS)


x11()
par(mfrow=c(2,1))
plot(x,y,type="l",lwd=2,xlab="x",ylab="f(x)",xlim=c(-1,2))
lines(x,ymirror,col=2,lwd=2)
plot(x-m,gamma,type="l",col=1,lwd=2,xlab="x",ylab=expression(gamma[X](s)))

skew=2*(beta-alpha)*sqrt(alpha+beta+1)/((alpha+beta+2)*sqrt(alpha*beta))

# alpha=22
# beta=21
# x=seq((alpha-beta)/(alpha+beta-2),1,length.out=10000)
# g=x^(alpha-2)*(1-x)^(beta-2)*((alpha-1)*(1-x)-(beta-1)*x)-(x+(beta-alpha)/(alpha+beta-2))^(beta-2)*(1-x-(beta-alpha)/(alpha+beta-2))^(alpha-2)*((beta-1)*(1-x-(beta-alpha)/(alpha+beta-2))-(alpha-1)*(x+(beta-alpha)/(alpha+beta-2)))
# x11()
# plot(x,g,type='l',lwd=1)
# abline(h=0,col=2)



### problematic distribution where mean, mode, median rule does not hold
# p=0.4926946
p=0.5
n=10000
x=seq(-5,5,length.out = n)
f=function(x,p){
  fx=(1-p)*(1+(1-p)/(2*p)*x)*1*(-2*p/(1-p)<=x)*1*(x<0) +  (1-p)*exp(-x)*1*(x>=0)
  return(fx)
}
y=f(x,p)
x11()
plot(x,y,type='l')

m=0
sm=f(m,p)

ymirror=f(-x,p)

gamma=(y-ymirror)/sm
MS=x[which.max(gamma)]-m
gammastar=max(gamma)*sign(MS)


x11(height=5,width=7)
par(mfrow=c(2,1),mai = c(0.8, 1, 0.1, 0.1))
plot(x,y,type="l",lwd=2,xlab="x",ylab=expression(f[X](x)),cex.axis=1.3,cex.lab=1.3)#,main=expression(beta~"=0.9"))
lines(x,ymirror,col=2,lwd=3,lty=2)
plot(x-m,gamma,type="l",col=1,lwd=2,xlab="s",ylab=expression(gamma[X](s)),cex.lab=1.3,cex.axis=1.3)
# abline(v=-log(2*p/(1-p)))

rootf=function(p){return(exp(-2*p/(1-p))-1-(1-p)/(2*p)*(log((1-p)/(2*p))-1))}
pr=uniroot(f = rootf,interval = c(0.3,0.5))$root
p=seq(0.01,0.99,by=0.001)
gammastarf=exp(-2*p/(1-p))*1*(p<pr)-(1+(1-p)/(2*p)*(log((1-p)/(2*p))-1))*1*(p>=pr)
x11(height=5,width=7)
par(mai = c(0.8, 1, 0.1, 0.1))
plot(p[1:483],gammastarf[1:483],type="l",lwd=2,xlab=expression(beta),ylab=expression(paste(gamma,"*(X)")),cex.lab=1.3,ylim=c(-1,1),xlim=c(p[1],p[981]),cex.axis=1.3)
lines(p[484:981],gammastarf[484:981],lwd=2)




### problematic distribution where density does not go to zero
n=10000
x=seq(-0.5,1.2,length.out = n)
f=function(x){
  fx=(0.9+3/5*x)*1*(x<1/3)*1*(x>0) +  (1.2-3/10*x)*1*(x>=1/3)*(x<=1)
  return(fx)
}
y=f(x)
x11()
plot(x,y,type='l')

m=1/3
sm=f(m)

ymirror=f(2/3-x)

gamma=(y-ymirror)/sm
MS=x[which.max(gamma)]-m
gammastar=max(gamma)*sign(MS)


x11()
par(mfrow=c(2,1))
plot(x,y,type="l",lwd=2,xlab="x",ylab=expression(f[X](x)))
lines(x,ymirror,col=2,lwd=2)
plot(x-m,gamma,type="l",col=1,lwd=2,xlab="s",ylab=expression(gamma[X](s)))




### univariate skew-normal
library(sn)
xi=0
omega=2
alpha=3
dp=c(xi,omega,alpha)

m=modeSECdistr(dp = dp,family = "SN")
sm=dsn(x = m,dp = dp)

n=10000
x=seq(-5,5,length=n)
y=dsn(x = x,dp = dp)
ymirror=dsn(x = 2*m-x,dp = dp)

x11()
par(mfrow=c(2,1),mai = c(0.8, 1, 0.1, 0.1))
plot(x,y,type="l",ylim=c(-0.02,0.4),lwd=2,ylab=expression(f[X](x)),cex.lab=1.3,xlab="x",cex.axis=1.3)
lines(x,ymirror,col=2,lwd=2,lty=2)
shape::Arrows(x0 = m,y0 = 0,x1 = m,y1 = sm,col=4,code = 3,arr.type="T",lwd=2)
shape::Arrows(x0 = x[5000],y0 = ymirror[5000],x1 = x[5000],y1 = y[5000],col=4,code = 3,arr.type="T",lwd=2)
shape::Arrows(x0 = m,y0 = sm/2,x1 = x[5000],y1 = sm/2,col=4,code = 3,arr.type="T",lwd=2)
text(x = m,y = sm/4,labels = expression(f[X](M[X])),pos = 4,col = 4,cex = 1.3)
text(x = x[5000],y = 0.27,labels = expression(f[X](M[X]+s)-f[X](M[X]-s)),pos = 2,col = 4,cex = 1.3)
text(x = (x[5000]+m)/2,y = y[5000]/2,labels = expression(s),pos = 3,col = 4,cex = 1.3)

s=seq(-5,5,length=n)
sl=dsn(x = m-s,dp = dp)
sr=dsn(x = m+s,dp = dp)
MS=s[which.max(sr-sl)]
VMS=max(sr-sl)/sm


plot(s,(sr-sl)/sm,type="l",xlab="s",ylab=expression(gamma[X](s)),lwd=2,cex.lab=1.3,cex.axis=1.3)
shape::Arrows(x0 = MS,y0 = 0,x1 = MS,y1 = VMS,col=4,code = 3,arr.type="T",lwd=2)
points(x = MS,y = VMS,col=2,pch=19)
text(x = MS,y = 0,labels = expression(M[S]),pos=1,cex=1.3)
text(x = MS,y = VMS/2,col=4,labels = expression(Gamma(X)),pos=4,cex=1.3)


r=1000
alphas=seq(-10,10,length.out = r)
MS=rep(NA,r)
VMS=rep(NA,r)
for(i in 1:r){
  xi=0
  omega=2
  alpha=alphas[i]
  dp=c(xi,omega,alpha)
  
  m=modeSECdistr(dp = dp,family = "SN")
  sm=dsn(x = m,dp = dp)
  s=seq(-10,10,length=n)
  sl=dsn(x = m-s,dp = dp)
  sr=dsn(x = m+s,dp = dp)
  MS[i]=s[which.max(sr-sl)]
  VMS[i]=max(sr-sl)/sm*sign(m)
}
x11()
par(mfrow=c(2,1),mai = c(1, 1, 0.1, 0.1))
plot(alphas,MS,type="l",lwd=2,xlab=expression(beta),ylab=expression(M[S]),cex.lab=1.2)
plot(alphas,VMS,type="l",lwd=2,xlab=expression(beta),ylab=expression(paste(gamma,"*(X)")),cex.lab=1.2)




### bivariate skew-normal
library(sn)
xi=c(1.25,-2.6)
Omega=matrix(c(3.5,-0.9,-0.9,5.8),nrow=2,ncol=2)
alpha=c(6,-6)
nu=4
dp=list("xi"=xi,"Omega"=Omega,"alpha"=alpha)#,"nu"=nu)

m=modeSECdistr(dp = dp,family = "SN")
sm=dmsn(x = m,dp = dp)

n=200
x1=seq(-5,10,length.out = n)
x2=seq(-10,5,length.out = n)
xgrid=expand.grid(x1,x2)
y=matrix(dmsn(x = xgrid,dp = dp),nrow=n,ncol=n)
ymirror=matrix(dmsn(x = sweep(x = -xgrid,MARGIN = 2,STATS = 2*m,FUN = "+"),dp = dp),nrow=n,ncol=n)

x11(height = 6,width=6)
contour(x1,x2,y,lwd=2,xlim=c(-5,10),ylim=c(-10,5),xlab=expression(x[1]),ylab=expression(x[2]),cex.lab=1.5)
contour(x1,x2,ymirror,col=4,add=T,lwd=2)
points(m[1],m[2],col=2,pch=19,cex=2)

x11(height = 6,width=6)
gamma=(y-ymirror)/sm
contour(x1-m[1],x2-m[2],gamma,lwd=2,cex=2,xlim=c(-5,5),ylim=c(-5,5),xlab=expression(s[1]),ylab=expression(s[2]),cex.lab=1.5)


s1=seq(-5,5,length.out = n)
s2=s1
s=expand.grid(s1,s2)
gamma2=(matrix(dmsn(x = sweep(x = s,MARGIN = 2,STATS = m,FUN = "+"),dp = dp),nrow=n,ncol=n)-
         matrix(dmsn(x = sweep(x = -s,MARGIN = 2,STATS = m,FUN = "+"),dp = dp),nrow=n,ncol=n))/sm
Gamma=max(gamma2)
indM=which(gamma2 == max(gamma2), arr.ind = TRUE)
MS=c(s1[indM[1]],s2[indM[2]])
gammastar=Gamma*MS/sqrt(sum(MS^2))
x11()
contour(s1,s2,gamma2,xlim=c(-5,10),ylim=c(-10,5))


library(ggplot2)
library(plotly)
library(reshape2)
library(gridExtra)
x11()
surface <- plot_ly(showscale = FALSE) %>% add_surface(x=x1,y=x2,z = y,opacity=0.99,color=I("red")) %>% add_surface(x=x1,y=x2,z = ymirror,opacity=0.8) %>% layout(scene = list(aspectmode='manual',aspectratio = list(x=1, y=1, z=0.5)))
surface

surface2 <- plot_ly(z=gamma,x=s1,y=s2,type="surface") %>% layout(scene = list(aspectmode='manual',aspectratio = list(x=1, y=1, z=0.5)))
surface2

# univariate QBA
library(QBAsyDist) 

dens=function(x,alpha,mu,phi,nu,basefunc){
  switch(basefunc,"normal"=dAND(y = x,mu = mu,alpha = alpha,phi = phi)
         ,"laplace"=dALaD(y = x,mu = mu,alpha = alpha,phi = phi)
         ,"logistic"=dALoD(y = x,mu = mu,alpha = alpha,phi = phi)
         ,"t"=dATD(y = x,mu = mu,nu = nu,alpha = alpha,phi = phi))
}

basefunc="normal"
alpha=0.25
mu=0
nu=4
# phi=switch(basefunc,"normal"=alpha^2*(1-alpha)^2/((1-2*alpha)^2*(1-2/pi)+alpha*(1-alpha))
#                    ,"laplace"=alpha^2*(1-alpha)^2/((1-2*alpha)^2*+alpha*(1-alpha)*2)
#                    ,"logistic"=alpha^2*(1-alpha)^2/((1-2*alpha)^2*(pi^2/3-4*ln(2)^2)+alpha*(1-alpha)*pi^2/3)
#                    ,"t"=alpha^2*(1-alpha)^2/((1-2*alpha)^2*(nu/(nu-2)-nu/pi*(gamma(nu/2-1/2)/gamma(nu/2))^2)+alpha*(1-alpha)*nu/(nu-2)))
phi=1



m=mu
sm=dens(m,alpha,mu,phi,nu,basefunc)

n=10000
x=seq(-5,5,length=n)
y=dens(x,alpha,mu,phi,nu,basefunc)
ymirror=dens(2*m-x,alpha,mu,phi,nu,basefunc)

x11()
plot(x,y,type="l")
lines(x,ymirror,col=2)
lines(x,(y-ymirror)/sm,type="l",col=5)

nn=1001
alphaseq=seq(0.1,0.9,length.out = nn)
difm=rep(NA,nn)
mas=rep(NA,nn)
for(i in 1:nn){
  alpha=alphaseq[i]
  phi=switch(basefunc,"normal"=alpha^2*(1-alpha)^2/((1-2*alpha)^2*(1-2/pi)+alpha*(1-alpha))
             ,"laplace"=alpha^2*(1-alpha)^2/((1-2*alpha)^2*+alpha*(1-alpha)*2)
             ,"logistic"=alpha^2*(1-alpha)^2/((1-2*alpha)^2*(pi^2/3-4*ln(2)^2)+alpha*(1-alpha)*pi^2/3)
             ,"t"=alpha^2*(1-alpha)^2/((1-2*alpha)^2*(nu/(nu-2)-nu/pi*(gamma(nu/2-1/2)/gamma(nu/2))^2)+alpha*(1-alpha)*nu/(nu-2)))
  
  y=dens(x,alpha,mu,phi,nu,basefunc)
  ymirror=dens(2*m-x,alpha,mu,phi,nu,basefunc)
  mas[i]=max(y-ymirror)/max(y)
  difm[i]=m-x[which.max(y-ymirror)]
}
plot(alphaseq,mas*difm/abs(difm),type="l")


### exact value of Gamma(X), MS and gamma*(x)
source('~/PhD/code/multivariate skewness measures/mode gamma.R')
MS=switch(basefunc,"normal"=snorm(alpha=alpha,phi=phi)
                  ,"laplace"=slap(alpha=alpha,phi=phi)
                  ,"logistic"=slog(alpha=alpha,phi=phi)
                  ,"t"=st(alpha=alpha,phi=phi,nu=nu))
GAMMA=switch(basefunc ,"normal"=Gnorm(alpha=alpha)
                      ,"laplace"=Glap(alpha=alpha)
                      ,"logistic"=Glog(alpha=alpha,s=MS,phi=phi,mu=mu)
                      ,"t"=Gt(alpha=alpha,nu=nu))
gammastar=sign(MS)*GAMMA



### bivariate QBA
source("densityplot2D.R")
alpha=c(0.25,0.75)
mu=c(1.25,-2.6)
A=matrix(c(2,-1.5,0.5,1.5),nrow=2)
tpars=NA
basefunc=c("normal","logistic")
n=2000

# set up the grid on which to plot the density
x1=seq(-30,30,length.out=n) # grid in 1 direction
x2=seq(-30,30,length.out=n)
xgrid=as.matrix(expand.grid(x1,x2)) # square grid

# compute and plot the density
f=matrix(densityf(x = xgrid,basefunc = basefunc,alpha = alpha,mu = mu,A = A,tpars = tpars),nrow=n,ncol=n)
fmirror=matrix(densityf(x = sweep(x = -xgrid,MARGIN = 2,STATS = 2*mu,FUN = "+"),basefunc = basefunc,alpha = alpha,mu = mu,A = A,tpars = tpars),nrow=n,ncol=n)

x11(height=6,width=6)
contour(x1,x2,f,lwd=2,xlab=expression(x[1]),ylab=expression(x[2]),cex.lab=1.5)
contour(x1,x2,fmirror,col=4,add=T,lwd=2)
points(mu[1],mu[2],col=2,pch=19,cex=2)

x11(height=6,width=6)
gamma2=(f-fmirror)/max(f)
contour(x1-mu[1],x2-mu[2],gamma2,lwd=2,cex=2,xlab=expression(s[1]),ylab=expression(s[2]),cex.lab=1.5)



s1=seq(-30,30,length.out=n)
s2=seq(-30,30,length.out=n)
s=as.matrix(expand.grid(s1,s2))
fm=densityf(x = matrix(mu,nrow=2,ncol=2,byrow=T),basefunc = basefunc,alpha = alpha,mu = mu,A = A,tpars = tpars)[1]
gamma=(matrix(densityf(x = sweep(x = s,MARGIN = 2,STATS = mu,FUN = "+"),basefunc = basefunc,alpha = alpha,mu = mu,A = A,tpars = tpars),nrow=n,ncol=n)-
         matrix(densityf(x = sweep(x = -s,MARGIN = 2,STATS = mu,FUN = "+"),basefunc = basefunc,alpha = alpha,mu = mu,A = A,tpars = tpars),nrow=n,ncol=n))/fm

Gamma=max(gamma)
indM=which(gamma==max(gamma),arr.ind = T)
SM=c(s1[indM[1]],s2[indM[2]])
gammastar=Gamma*SM/sqrt(sum(SM^2))

x11()
contour(s1,s2,gamma,levels = seq(min(gamma),max(gamma),length.out = 500))
points(SM[1],SM[2],pch=18,col=2)

# -0.6754127 -0.1564060 0.5
# -0.6754102 -0.1369216 1
# -0.6754100 -0.1608588 1.1
# -0.6754101 -0.1650469 1.2
# -0.6754099 -0.1506902 1.5
# -0.6754096 -0.1666820 2.5



# slope zero contour for 2 normals
C=-2*((1-2*alpha[1])*B[1,1]*B[2,1]+(1-2*alpha[2])*B[1,2]*B[2,2]-sqrt(-(1-2*alpha[1])*(1-2*alpha[2]))*abs(det(B)))/(2*B[1,1]^2*(1-2*alpha[1])+2*B[1,2]^2*(1-2*alpha[2]))
abline(a=-5,b = 1/C,col=3)

func=function(s1,s2,B,alpha){
  s=matrix(c(s1,s2),nrow=2,ncol=1)
  return((1-2*alpha[1])*(t(s)%*%B[,1])^2 + (1-2*alpha[2])*(t(s)%*%B[,2])^2)
}
p1=uniroot(f = func,lower = 0,upper = 10,s2=1,B=B,alpha=alpha)
p2=uniroot(f = func,lower = 0,upper = 10,s2=2,B=B,alpha=alpha)


### data example ###
library(DAAG)
library(ks)
source("~/PhD/code/multivariate skewness measures/bivariate asymmetry measure.R")
dat=as.matrix(ais[,c(6,9)])
a=ddAsymmetry2(X = dat,gridpoints = c(500,500),mingrid = c(10,20),maxgrid = c(40,120),plot.contour = T,H.fact = 2)
X = dat;gridpoints = c(1000,700);mingrid = c(0,0);maxgrid = c(50,120);plot.contour = T;H.fact = 3


n=1000
d=ncol(dat)

test=kde(x = dat,gridsize = c(n,n),xmin = c(0,0),positive = T,xmax = c(50,120),compute.cont = T)
H=test$H
densest=kde(x = dat,H=15*H,gridsize = c(n,n),xmin = c(0,0),positive = T,xmax = c(50,120),compute.cont = T)
densest$eval.points=lapply(densest$eval.points,round,8)

x11()
plot(densest,cont=c(1,5,10,20,30,40,50,60,70,80,90,95,99))

ms=max(densest$estimate)
indM=which(densest$estimate == ms, arr.ind = TRUE)
M=rep(NA,d)
for(i in 1:d){
  M[i]=(densest$eval.points[[i]])[indM[i]]
}

x=rev(2*M[1]-densest$eval.points[[1]])
y=rev(2*M[2]-densest$eval.points[[2]])
gridM=list(x,y)


### bivariate examples from the linear combinations paper

## example 1
source("densityplot2D.R")
alpha=c(0.25,0.65)
mu=c(20,20)
A=matrix(c(12,-5,4,8),nrow=2)
tpars=NA
basefunc=c("normal","logistic")
n=4000

# set up the grid on which to plot the density
x1=seq(-100,150,length.out=n) # grid in 1 direction
x2=seq(-100,150,length.out=n)
xgrid=as.matrix(expand.grid(x1,x2)) # square grid

# compute and plot the density
f=matrix(densityf(x = xgrid,basefunc = basefunc,alpha = alpha,mu = mu,A = A,tpars = tpars),nrow=n,ncol=n)
fmirror=matrix(densityf(x = sweep(x = -xgrid,MARGIN = 2,STATS = 2*mu,FUN = "+"),basefunc = basefunc,alpha = alpha,mu = mu,A = A,tpars = tpars),nrow=n,ncol=n)

x11()
contour(x1,x2,f,lwd=2)
contour(x1,x2,fmirror,col=3,add=T,lwd=2)
points(mu[1],mu[2],col=2,pch=19,cex=2)

x11()
gamma2=(f-fmirror)/max(f)
contour(x1-mu[1],x2-mu[2],gamma2,lwd=2,cex=2)

MS=which((f-fmirror)==max(f-fmirror),arr.ind = T)
gammastar=max(f-fmirror)/max(f)*(c(x1[MS[1]],x2[MS[2]])-mu)/sqrt(sum((c(x1[MS[1]],x2[MS[2]])-mu)^2))


## example 2
alpha=c(0.25,0.65)
mu=c(20,20)
A=matrix(c(7,0,-6,3),nrow=2)
tpars=NA
basefunc=c("normal","logistic")
n=4000

# set up the grid on which to plot the density
x1=seq(-100,150,length.out=n) # grid in 1 direction
x2=seq(-100,150,length.out=n)
xgrid=as.matrix(expand.grid(x1,x2)) # square grid

# compute and plot the density
f=matrix(densityf(x = xgrid,basefunc = basefunc,alpha = alpha,mu = mu,A = A,tpars = tpars),nrow=n,ncol=n)
fmirror=matrix(densityf(x = sweep(x = -xgrid,MARGIN = 2,STATS = 2*mu,FUN = "+"),basefunc = basefunc,alpha = alpha,mu = mu,A = A,tpars = tpars),nrow=n,ncol=n)

x11()
contour(x1,x2,f,lwd=2)
contour(x1,x2,fmirror,col=3,add=T,lwd=2)
points(mu[1],mu[2],col=2,pch=19,cex=2)

x11()
gamma2=(f-fmirror)/max(f)
contour(x1-mu[1],x2-mu[2],gamma2,lwd=2,cex=2)

MS=which((f-fmirror)==max(f-fmirror),arr.ind = T)
gammastar=max(f-fmirror)/max(f)*(c(x1[MS[1]],x2[MS[2]])-mu)/sqrt(sum((c(x1[MS[1]],x2[MS[2]])-mu)^2))


## example 3
alpha=c(0.25,0.65)
mu=c(20,20)
A=matrix(c(12,0,0,8),nrow=2)
tpars=c(NA,5)
basefunc=c("normal","t")
n=4000

# set up the grid on which to plot the density
x1=seq(-100,150,length.out=n) # grid in 1 direction
x2=seq(-100,150,length.out=n)
xgrid=as.matrix(expand.grid(x1,x2)) # square grid

# compute and plot the density
f=matrix(densityf(x = xgrid,basefunc = basefunc,alpha = alpha,mu = mu,A = A,tpars = tpars),nrow=n,ncol=n)
fmirror=matrix(densityf(x = sweep(x = -xgrid,MARGIN = 2,STATS = 2*mu,FUN = "+"),basefunc = basefunc,alpha = alpha,mu = mu,A = A,tpars = tpars),nrow=n,ncol=n)

x11()
contour(x1,x2,f,lwd=2)
contour(x1,x2,fmirror,col=3,add=T,lwd=2)
points(mu[1],mu[2],col=2,pch=19,cex=2)

x11()
gamma2=(f-fmirror)/max(f)
contour(x1-mu[1],x2-mu[2],gamma2,lwd=2,cex=2)

MS=which((f-fmirror)==max(f-fmirror),arr.ind = T)
gammastar=max(f-fmirror)/max(f)*(c(x1[MS[1]],x2[MS[2]])-mu)/sqrt(sum((c(x1[MS[1]],x2[MS[2]])-mu)^2))


## simulation 1
alpha=c(0.35,0.7)
mu=c(0,0)
A=matrix(c(4,-3,1,4),nrow=2)
tpars=c(NA,6)
basefunc=c("normal","t")
n=4000

# set up the grid on which to plot the density
x1=seq(-50,50,length.out=n) # grid in 1 direction
x2=seq(-50,50,length.out=n)
xgrid=as.matrix(expand.grid(x1,x2)) # square grid

x1=seq(-20,40,length.out=n) # grid in 1 direction
x2=seq(-40,20,length.out=n)
xgrid=as.matrix(expand.grid(x1,x2)) # square grid


# compute and plot the density
f=matrix(densityf(x = xgrid,basefunc = basefunc,alpha = alpha,mu = mu,A = A,tpars = tpars),nrow=n,ncol=n)
fmirror=matrix(densityf(x = sweep(x = -xgrid,MARGIN = 2,STATS = 2*mu,FUN = "+"),basefunc = basefunc,alpha = alpha,mu = mu,A = A,tpars = tpars),nrow=n,ncol=n)

x11()
contour(x1,x2,f,lwd=2)
contour(x1,x2,fmirror,col=3,add=T,lwd=2)
points(mu[1],mu[2],col=2,pch=19,cex=2)

x11()
gamma2=(f-fmirror)/max(f)
contour(x1-mu[1],x2-mu[2],gamma2,lwd=2,cex=2)

MS=which((f-fmirror)==max(f-fmirror),arr.ind = T)
gammastar=max(f-fmirror)/max(f)*(c(x1[MS[1]],x2[MS[2]])-mu)/sqrt(sum((c(x1[MS[1]],x2[MS[2]])-mu)^2))




### ToS distribution of Abe & Fujisama ###
##########################################
rx <- function(x,lambda){
  if(lambda==0){
    return(x)
  } else {
    al <- 1-exp(-lambda^2)
    r <- (lambda*x + al-al*sqrt((lambda*x+al)^2+1-al^2))/(lambda*(1-al^2))
    return(r)
  }
}
S <- function(x,lambda){
  return(sinh(lambda*asinh(x)))
}
fx <- function(x,lambda){
  s <- S(x,lambda)
  f <- lambda/(2*pi*(1+x^2))*(1+s^2)^(1/2)*exp(-s^2/2)
  return(f)
}  
n <- 2500
x <- seq(-2,2,length.out=n)
lambda <- 2
## normal f
# sigma <- 1
# f <- dnorm(x = rx(x = x,lambda = lambda),mean = 0,sd = sigma)
# fmirror <- dnorm(x = rx(x = -x,lambda = lambda),mean = 0,sd = sigma)
# fm <- dnorm(0)

## sinh-arcsinh f
f <- fx(x = rx(x = x,lambda = lambda),lambda = lambda)
fmirror <- fx(x = rx(x = -x,lambda = lambda),lambda = lambda)
fm <- fx(x = 0,lambda = lambda)
gamma <- (f-fmirror)/fm

x11()
par(mfrow=c(2,1),mai = c(1, 1, 0.1, 0.1))
plot(x,f,type="l",lwd=2,xlab=expression(x),ylab=expression(f(x)),cex.lab=1.3,cex.axis=1.3)
lines(x,fmirror,lwd=2,lty=2,col=2)

plot(x,gamma,type="l",lwd=2,xlab=expression(x),ylab=expression(paste(gamma,"(x)")),cex.lab=1.3,cex.axis=1.3)

MS <- x[which.max(gamma)]
Gamma <- max(gamma)
gammastar <- sign(MS)*Gamma

DBA_AF <- function(x,lambda,dist="SA",delta=1){
  ## function to calculate the skewness function of the ToS proposed in F&A 2015
  ## f(r(x))
  ## x is a n-vector containing the points in which to evaluate the SF
  ## lambda is the scalar skewing parameter
  ## dist is a character specifying the symmetric reference density, 
  #   SA for sinh-arcsinh, L for logistic, N for normal
  ## delta is the skewing parameter of the SA (positive)
  
  # function r(x)
  rx <- function(x,lambda){
    if(lambda==0){
      return(x)
    } else {
      al <- 1-exp(-lambda^2)
      r <- (lambda*x + al-al*sqrt((lambda*x+al)^2+1-al^2))/(lambda*(1-al^2))
      return(r)
    }
  }
  
  # density and mirrored density
  if(dist=='N'){
    f <- dnorm(x = rx(x = x,lambda = lambda),mean = 0,sd = 1)
    fmirror <- dnorm(x = rx(x = -x,lambda = lambda),mean =  0,sd = 1)
    fm <- dnorm(x = 0,mean =0,sd = 1)
  } else if(dist=='L'){
    f <- dlogis(x = rx(x = x,lambda = lambda),location = 0,scale = 1)
    fmirror <- dlogis(x = rx(x = -x,lambda = lambda),location = 0,scale = 1)
    fm <- dlogis(x = 0,location = 0,scale = 1)
  } else {
    # sinh-arcsinh function as f
    S <- function(x,delta){
      return(sinh(delta*asinh(x)))
    }
    
    fx <- function(x,delta){
      s <- S(x,delta)
      f <- delta/(2*pi*(1+x^2))*(1+s^2)^(1/2)*exp(-s^2/2)
      return(f)
    } 
    f <- fx(x = rx(x = x,lambda = lambda),delta = delta)
    fmirror <- fx(x = rx(x = -x,lambda = lambda),delta = delta)
    fm <- fx(x = 0,delta = delta)
  }

  
  # gamma as a function of s
  gamma <- (f-fmirror)/fm
  
  # summarizing measure gamma star
  MS <- x[which.max(gamma)]
  Gamma <- max(gamma)
  gammastar <- sign(MS)*Gamma
  return(list('f'=f,'fmirror'=fmirror,'gamma'=gamma,'MS'=MS,'Gamma'=Gamma,'gammastar'=gammastar))
}

n <- 2500
x <- seq(-5,5,length.out=n)
lambda1 <- 0.75

out1 <- DBA_AF(x = x,lambda = lambda1,dist = "L")

x11(width = 6,height = 3)
par(mai = c(0.8, 1, 0.3, 0.1))
plot(x,out1$f,type="l",lwd=2,lty=1,cex.lab=1.3,cex.axis=1.3,ylim=c(0,0.3),ylab=expression(f[X](x)),xlab="x")
lines(x,out1$fmirror,lwd=2,lty=2,col=2)

x11(width = 6,height = 3)
par(mai = c(0.8, 1, 0.3, 0.1))
plot(x,out1$gamma,type="l",lwd=2,lty=1,cex.lab=1.3,cex.axis=1.3,ylim=c(-1,1),ylab=expression(gamma(s)),xlab="s")


lambda2 <- -2

out2 <- DBA_AF(x = x,lambda = lambda2,dist = "L")

x11(width = 6,height = 3)
par(mai = c(0.8, 1, 0.3, 0.1))
plot(x,out2$f,type="l",lwd=2,lty=1,cex.lab=1.3,cex.axis=1.3,ylim=c(0,0.3),ylab=expression(f[X](x)),xlab="x")
lines(x,out2$fmirror,lwd=2,lty=2,col=2)

x11(width = 6,height = 3)
par(mai = c(0.8, 1, 0.3, 0.1))
plot(x,out2$gamma,type="l",lwd=2,lty=1,cex.lab=1.3,cex.axis=1.3,ylim=c(-1,1),ylab=expression(gamma(s)),xlab="s")

### linear combination of ToS
A <- matrix(c(2,-1.5,0.5,1.5))
mu <- c(1.25,-2.6)
lambda <- c(0.75,-2)
distr <- c("N","L")
n <- 2000

x1 <- seq(-20,20,length.out=n)
x2 <- seq(-20,20,length.out=n)
xgrid <- as.matrix(expand.grid(x1,x2))

f <- function(x,distr,lambda,mu,A){
  
  B <- solve(matrix(A,nrow=2,ncol=2))
  Y <- sweep(x,2,mu)%*%B
  
  d1 <- DBA_AF(x = Y[,1],lambda = lambda[1],dist = distr[1])$f
  d2 <- DBA_AF(x = Y[,2],lambda = lambda[2],dist = distr[2])$f
  dens <- abs(det(B))*d1*d2
  return(dens)
}

fmat <- matrix(f(x = xgrid,distr = distr,lambda = lambda,mu = mu,A = A),nrow=n,ncol=n)
fmirror <- matrix(f(x = sweep(x=-xgrid,MARGIN = 2,STATS = 2*mu,FUN = '+'),distr = distr,lambda = lambda,mu = mu,A = A),nrow = n,ncol = n)
x11(width = 10,height = 10)
par(mai = c(0.8, 1, 0.1, 0.1))
contour(x1,x2,fmat,lwd=2,lty=1,cex.lab=1.3,cex.axis=1.3,ylab=expression(x[2]),xlab=expression(x[1]),xlim=c(-15,15),ylim=c(-15,15),labcex = 1)
contour(x1,x2,fmirror,lwd=2,lty=1,col=4,labcex = 1,add = T)
points(mu[1],mu[2],col=2,pch=18,cex=2)


fm <- f(x = matrix(mu,nrow=1,ncol=2),distr = distr,lambda = lambda,mu = mu,A = A)
gamma <- (fmat - fmirror)/fm
x11(width = 10,height = 10)
par(mai = c(0.8, 1, 0.1, 0.1))
contour(x1-mu[1],x2-mu[2],gamma,lwd=2,lty=1,cex.lab=1.3,cex.axis=1.3,ylab=expression(s[2]),xlab=expression(s[1]),xlim = c(-10,10),ylim = c(-10,10),labcex = 1)
