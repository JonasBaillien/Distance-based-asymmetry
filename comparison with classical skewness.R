###############################################
###  Comparison with classical skewness     ###
###############################################



### Gamma distribution ###
##########################

optfungamma=function(s,alpha,beta){
  return(dgamma(x = (alpha-1)/beta-s,shape = alpha,rate = beta)-dgamma(x = (alpha-1)/beta+s,shape = alpha,rate = beta))
}
skewnessgamma=function(alpha,beta){
  if(alpha<=1){
    gammastar=1
  } else {
      MS=optim(par = 0.35,fn = optfungamma,method = "L-BFGS-B",alpha=alpha,beta=beta)$par
    Gamma=(dgamma(x = (alpha-1)/beta+MS,shape = alpha,rate = beta)-dgamma(x = (alpha-1)/beta-MS,shape = alpha,rate = beta))/(dgamma(x = (alpha-1)/beta,shape = alpha,rate = beta))
    gammastar=sign(MS)*Gamma
  }

  return(gammastar)
}


n=10000
alpha=seq(0,10,length.out = n)
beta=5


set.seed(87)
mardiagamma=2/sqrt(alpha)
dbgamma=unlist(lapply(X = alpha,FUN = skewnessgamma,beta=beta))

x11()
plot(alpha,mardiagamma,type="l",xlab=expression(alpha),ylab="Classical skewness",lwd=2,cex.lab=1.5,cex=1.5,cex.axis=1.5)
x11()
plot(alpha,dbgamma,type="l",xlab=expression(alpha),ylab="Proposed skewness",lwd=2,cex.lab=1.5,cex=1.5,cex.axis=1.5)


n=1000
alpha=seq(0,10,length.out = n)
beta=seq(0.1,10, by=0.1)
dbg=matrix(NA,nrow=length(beta),ncol=n)
for(i in 1:length(beta)){
  dbg[i,]=unlist(lapply(X = alpha,FUN = skewnessgamma,beta=beta[i]))
}
library(plot3D)
x11()
contour2D(z = dbg,y = alpha,x = beta,xlab=expression(beta),ylab=expression(alpha),lwd=2,cex.lab=1.5,cex=1.5,cex.axis=1.5,col=gray.colors(n = 9,start = 0.8,end = 0),labcex=1.2)


### beta distribution ###
#########################


optfunbeta=function(s,alpha,beta){
  return(-dbeta(x = (alpha-1)/(alpha+beta-2)+s,shape1 = alpha,shape2 = beta)+dbeta(x = (alpha-1)/(alpha+beta-2)-s,shape1 = alpha,shape2 = beta))
}
skewnessbeta=function(alpha,beta){
  if(alpha<=1 & beta>1){
    gammastar=1
  } else if(beta<=1 & alpha > 1){
    gammastar=-1
  } else if(alpha==1 & beta==1){
    gammastar=0
  } else if(alpha<1 & beta<1){
    warning("bimodal distribution detected")
    gammastar=NA
  } else {
    if(alpha<beta){
      MS=optim(par = 0.1,fn = optfunbeta,method = "L-BFGS-B",alpha=alpha,beta=beta)$par
    } else {
      MS=optim(par = -0.1,fn = optfunbeta,method = "L-BFGS-B",alpha=alpha,beta=beta)$par
    }
    GAMMA=(dbeta(x = (alpha-1)/(alpha+beta-2)+MS,shape1 = alpha,shape2 = beta)-dbeta(x = (alpha-1)/(alpha+beta-2)-MS,shape1 = alpha,shape2 = beta))/dbeta(x = (alpha-1)/(alpha+beta-2),shape1 = alpha,shape2 = beta)
    gammastar=sign(MS)*GAMMA
  }
  return(gammastar)
}

CSKbeta=function(alpha,beta){
  return(2*(beta-alpha)*sqrt(alpha+beta+1)/((alpha+beta+2)*sqrt(alpha*beta)))
}

n=10000
beta=4
alpha=seq(0.01,10,length.out=n)

mardiabeta=CSKbeta(alpha,beta)
dbbeta=unlist(lapply(alpha,skewnessbeta,beta=beta))

x11()
plot(alpha,mardiabeta,type="l",xlab=expression(alpha),ylab="Classical skewness",lwd=2,cex.lab=1.5,cex=1.5,cex.axis=1.5)
# abline(h=0,v=beta,col=2,lwd=0.5)
x11()
plot(alpha,dbbeta,type="l",xlab=expression(alpha),ylab="Proposed skewness",lwd=2,cex.lab=1.5,cex=1.5,cex.axis=1.5)
#abline(h=0,v=beta,col=2,lwd=0.5)


n=1000
alpha=seq(1.1,10,length.out = n)
beta=seq(1.1,10, length.out = n)
dbb=matrix(NA,nrow=n,ncol=n)
mb=matrix(NA,nrow=n,ncol=n)
for(i in 1:n){
  dbb[i,]=unlist(lapply(X = alpha,FUN = skewnessbeta,beta=beta[i]))
  mb[i,]=unlist(lapply(X = alpha,FUN = CSKbeta, beta=beta[i]))
}
library(plot3D)
x11()
contour2D(z = mb,y = alpha,x = beta,xlab=expression(beta),ylab=expression(alpha),lwd=3,cex.lab=1.5,cex=2,cex.axis=1.5,col=gray.colors(n = 9,start = 0.8,end = 0),labcex=1.2,levels=c(-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8))
x11()
contour2D(z = dbb,y = alpha,x = beta,xlab=expression(beta),ylab=expression(alpha),lwd=3,cex.lab=1.5,cex=2,cex.axis=1.5,col=gray.colors(n = 9,start = 0.8,end = 0),labcex=1.2,levels=c(-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8))



### QBA normal distribution ###
###############################

source('~/PhD/code/multivariate skewness measures/QBA mode gamma.R')
library(QBAsyDist)
n=10000
alpha=seq(0.01,0.99,length.out = n)

mardiaAND=skewAND(alpha)
dbAND=unlist(lapply(X = alpha,FUN = Gnorm))*sign(0.5-alpha)

x11()
plot(alpha,mardiaAND,type="l",xlab=expression(alpha),ylab="Classical skewness",lwd=2,cex.lab=1.5,cex=1.5,cex.axis=1.5)
x11()
plot(alpha,dbAND,type="l",xlab=expression(alpha),ylab="Proposed skewness",lwd=2,cex.lab=1.5,cex=1.5,cex.axis=1.5)




### problematic distribution where mean, mode, median rule does not hold ###
############################################################################

n=10000
x=seq(-4,4,length.out = n)
f=function(x,p){
  fx=(1-p)*(1+(1-p)/(2*p)*x)*1*(-2*p/(1-p)<=x)*1*(x<0) +  (1-p)*exp(-x)*1*(x>=0)
  return(fx)
}
d1=function(p){return(2*p^2/(3*(p-1)))}
d2=function(p){return(2*p^3/(3*(1-p)^2))}
d3=function(p){return(4*p^4/(5*(p-1)^3))}
sf=function(p){
  s=( d3(p)+6-6*p -3*(d1(p)+1-p)*(d2(p)+2-2*p-(d1(p)+1-p)^2)-(d1(p)+1-p)^3 )/((d2(p)+2-2*p-(d1(p)+1-p)^2)^{3/2})
  return(s)
}  


p=c(0.25,0.375,0.5,0.625,0.75)
y1=f(x,p[1])
y2=f(x,p[2])
y3=f(x,p[3])
y4=f(x,p[4])
y5=f(x,p[5])
x11(width=6,height=5)
par(mai = c(0.8, 1, 0.1, 0.1))
plot(x,y1,type='l',col=1,lwd=2,xlab="x",ylab=expression(f[X](x)),cex.lab=1.3,cex.axis=1.3)
lines(x,y2,col=2,lwd=3,lty=2)
lines(x,y3,col=gray.colors(1, start = 0.8, end = 0.8),lwd=3,lty=1)
lines(x,y4,col=4,lwd=3,lty=4)
lines(x,y5,col=5,lwd=3,lty=5)
legend("topright",legend = c(expression(beta~"= 0.25"),expression(beta~"= 0.375"),expression(beta~"= 0.5"),expression(beta~"= 0.625"),expression(beta~"= 0.75")),lty=c(1,2,1,4,5),lwd = 3,col = c(1,2,gray.colors(1, start = 0.8, end = 0.8),4,5),cex=1.3)

# p=0.4926946
p=0.5
m=0
sm=f(m,p)
y=f(x,p)
ymirror=f(-x,p)

ps=seq(0,1,length.out=10000)
sfp=sf(ps)
x11(width=6,height=5)
plot(ps,sfp,xlab=expression(beta),ylab="Classical skewness",type="l",lwd=2,cex.axis=1.3,cex.lab=1.3)
# abline(v=0.755)
# abline(h=0)





gamma=(y-ymirror)/sm
MS=x[which.max(gamma)]-m
gammastar=max(gamma)*sign(MS)

rootf=function(p){return(exp(-2*p/(1-p))-1-(1-p)/(2*p)*(log((1-p)/(2*p))-1))}
pr=uniroot(f = rootf,interval = c(0.3,0.5))$root
p=seq(0.01,0.99,by=0.001)
gammastarf=exp(-2*p/(1-p))*1*(p<pr)-(1+(1-p)/(2*p)*(log((1-p)/(2*p))-1))*1*(p>=pr)

x11()
par(mai=c(0.8,1,0.1,0.1))
plot(p[1:483],gammastarf[1:483],type="l",lwd=2,xlab="p",ylab=expression(paste(gamma,"*(X)")),cex.lab=1.2,ylim=c(-1,1),xlim=c(p[1],p[981]),cex.axis=1.3,cex.lab=1.3)
lines(p[484:981],gammastarf[484:981],lwd=2)