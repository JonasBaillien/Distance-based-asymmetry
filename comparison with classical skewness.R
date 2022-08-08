###############################################
###  Comparison with classical skewness     ###
###############################################



### Gamma distribution ###
##########################

# function evaluate gamma(s)
optfungamma <- function(s,alpha,beta){
  # s: n-vector of points for evaluation
  # alpha, beta: parameters of gamma distribution
  
  return(dgamma(x = (alpha-1)/beta-s,shape = alpha,rate = beta)-dgamma(x = (alpha-1)/beta+s,shape = alpha,rate = beta))
}

# function to determine gamma star
skewnessgamma <- function(alpha,beta){
  # alpha, beta: parameters of gamma distribution
  
  if(alpha<=1){
    gammastar <- 1
  } else {
      MS <- optim(par = 0.35,fn = optfungamma,method = "L-BFGS-B",alpha=alpha,beta=beta)$par
    Gamma <- (dgamma(x = (alpha-1)/beta+MS,shape = alpha,rate = beta)-dgamma(x = (alpha-1)/beta-MS,shape = alpha,rate = beta))/(dgamma(x = (alpha-1)/beta,shape = alpha,rate = beta))
    gammastar <- sign(MS)*Gamma
  }

  return(gammastar)
}

### vary only alpha

# number of points
n <- 10000
# vector with values for parameter alpha
alpha <- seq(0,10,length.out = n)
# parameter beta
beta <- 5

# seed
set.seed(87)

# Pearson skewness
skewgamma <- 2/sqrt(alpha)
# our summarizing measure of asymmetry
dbgamma <- unlist(lapply(X = alpha,FUN = skewnessgamma,beta=beta))

# plotting
x11()
plot(alpha,skewgamma,type="l",xlab=expression(alpha),ylab="Classical skewness",lwd=2,cex.lab=1.5,cex=1.5,cex.axis=1.5)
x11()
plot(alpha,dbgamma,type="l",xlab=expression(alpha),ylab="Proposed skewness",lwd=2,cex.lab=1.5,cex=1.5,cex.axis=1.5)


### vary both alpha and beta
# number of points for evaluation
n <- 1000
# vector of values for alpha
alpha <- seq(0,10,length.out = n)
# vector of values for beta
beta <- seq(0.1,10, by=0.1)

# our summarizing asymmetry measure (Pearson skewness only depends on alpha)
dbg <- matrix(NA,nrow=length(beta),ncol=n)
for(i in 1:length(beta)){
  dbg[i,] <- unlist(lapply(X = alpha,FUN = skewnessgamma,beta=beta[i]))
}
library(plot3D)
x11()
contour2D(z = dbg,y = alpha,x = beta,xlab=expression(beta),ylab=expression(alpha),lwd=2,cex.lab=1.5,cex=1.5,cex.axis=1.5,col=gray.colors(n = 9,start = 0.8,end = 0),labcex=1.2)
# also independent of beta



### beta distribution ###
#########################

# function to evaluate gamma(s)
optfunbeta <- function(s,alpha,beta){
  # s: n-vector of points for evaluation
  # alpha/beta: parameters of beta distribution
  
  return(-dbeta(x = (alpha-1)/(alpha+beta-2)+s,shape1 = alpha,shape2 = beta)+dbeta(x = (alpha-1)/(alpha+beta-2)-s,shape1 = alpha,shape2 = beta))
}

# function to find gammastar
skewnessbeta <- function(alpha,beta){
  # alpha/beta: parameters of beta distribution
  
  if(alpha<=1 & beta>1){
    gammastar <- 1
  } else if(beta<=1 & alpha > 1){
    gammastar <- -1
  } else if(alpha==1 & beta==1){
    gammastar <- 0
  } else if(alpha<1 & beta<1){
    warning("bimodal distribution detected")
    gammastar=NA
  } else {
    if(alpha<beta){
      MS <- optim(par = 0.1,fn = optfunbeta,method = "L-BFGS-B",alpha=alpha,beta=beta)$par
    } else {
      MS <- optim(par = -0.1,fn = optfunbeta,method = "L-BFGS-B",alpha=alpha,beta=beta)$par
    }
    GAMMA <- (dbeta(x  <-  (alpha-1)/(alpha+beta-2)+MS,shape1 = alpha,shape2 = beta)-dbeta(x = (alpha-1)/(alpha+beta-2)-MS,shape1 = alpha,shape2 = beta))/dbeta(x = (alpha-1)/(alpha+beta-2),shape1 = alpha,shape2 = beta)
    gammastar <- sign(MS)*GAMMA
  }
  return(gammastar)
}

# Pearson skewness of beta distribution
CSKbeta <- function(alpha,beta){
  # alpha/beta: parameters of beta distribution
  
  return(2*(beta-alpha)*sqrt(alpha+beta+1)/((alpha+beta+2)*sqrt(alpha*beta)))
}

### Varying only alpha
# number of points
n <- 10000
# beta parameter
beta <- 4
# alpha parameter
alpha <- seq(0.01,10,length.out=n)

# Pearson skewness
skewbeta <- CSKbeta(alpha,beta)
# our summarizing asymmetry measure
dbbeta <- unlist(lapply(alpha,skewnessbeta,beta=beta))

x11()
plot(alpha,skewbeta,type="l",xlab=expression(alpha),ylab="Classical skewness",lwd=2,cex.lab=1.5,cex=1.5,cex.axis=1.5)
# abline(h=0,v=beta,col=2,lwd=0.5)
x11()
plot(alpha,dbbeta,type="l",xlab=expression(alpha),ylab="Proposed skewness",lwd=2,cex.lab=1.5,cex=1.5,cex.axis=1.5)
#abline(h=0,v=beta,col=2,lwd=0.5)

# varying both alpha and beta
# number of points
n <- 1000
# values for alpha
alpha <- seq(1.1,10,length.out = n)
# values for beta
beta <- seq(1.1,10, length.out = n)

# calculating both asymmetry measures
dbb <- matrix(NA,nrow=n,ncol=n)
sb <- matrix(NA,nrow=n,ncol=n)
for(i in 1:n){
  dbb[i,] <- unlist(lapply(X = alpha,FUN = skewnessbeta,beta=beta[i]))
  sb[i,] <- unlist(lapply(X = alpha,FUN = CSKbeta, beta=beta[i]))
}

# plotting 
library(plot3D)
x11()
contour2D(z = sb,y = alpha,x = beta,xlab=expression(beta),ylab=expression(alpha),lwd=3,cex.lab=1.5,cex=2,cex.axis=1.5,col=gray.colors(n = 9,start = 0.8,end = 0),labcex=1.2,levels=c(-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8))
x11()
contour2D(z = dbb,y = alpha,x = beta,xlab=expression(beta),ylab=expression(alpha),lwd=3,cex.lab=1.5,cex=2,cex.axis=1.5,col=gray.colors(n = 9,start = 0.8,end = 0),labcex=1.2,levels=c(-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8))



### QBA normal distribution ###
###############################

# load in required functions and library
source('~/PhD/code/multivariate skewness measures/QBA mode gamma.R')
library(QBAsyDist)

# values of alpha
n <- 10000
alpha <- seq(0.01,0.99,length.out = n)

# Pearson skewness
skewnessAND <- skewAND(alpha)
# our summarizing measure of asymmetry
dbAND <- unlist(lapply(X = alpha,FUN = Gnorm))*sign(0.5-alpha)

# plotting 
x11()
plot(alpha,mardiaAND,type="l",xlab=expression(alpha),ylab="Classical skewness",lwd=2,cex.lab=1.5,cex=1.5,cex.axis=1.5)
x11()
plot(alpha,dbAND,type="l",xlab=expression(alpha),ylab="Proposed skewness",lwd=2,cex.lab=1.5,cex=1.5,cex.axis=1.5)




### problematic distribution where mean, mode, median rule does not hold ###
############################################################################

# number of points
n <- 10000

# range of x-values of the distribution (domain)
x <- seq(-4,4,length.out = n)
# density function
f <- function(x,p){
  fx <- (1-p)*(1+(1-p)/(2*p)*x)*1*(-2*p/(1-p)<=x)*1*(x<0) +  (1-p)*exp(-x)*1*(x>=0)
  return(fx)
}

# moments of density function
d1 <- function(p){return(2*p^2/(3*(p-1)))}
d2 <- function(p){return(2*p^3/(3*(1-p)^2))}
d3 <- function(p){return(4*p^4/(5*(p-1)^3))}
# Pearson skewness of distribution
sf <- function(p){
  s <- ( d3(p)+6-6*p -3*(d1(p)+1-p)*(d2(p)+2-2*p-(d1(p)+1-p)^2)-(d1(p)+1-p)^3 )/((d2(p)+2-2*p-(d1(p)+1-p)^2)^{3/2})
  return(s)
}  

# plotting density for different parameter values
p <- c(0.25,0.375,0.5,0.625,0.75)
y1 <- f(x,p[1])
y2 <- f(x,p[2])
y3 <- f(x,p[3])
y4 <- f(x,p[4])
y5 <- f(x,p[5])
x11(width=6,height=5)
par(mai = c(0.8, 1, 0.1, 0.1))
plot(x,y1,type='l',col=1,lwd=2,xlab="x",ylab=expression(f[X](x)),cex.lab=1.3,cex.axis=1.3)
lines(x,y2,col=2,lwd=3,lty=2)
lines(x,y3,col=gray.colors(1, start = 0.8, end = 0.8),lwd=3,lty=1)
lines(x,y4,col=4,lwd=3,lty=4)
lines(x,y5,col=5,lwd=3,lty=5)
legend("topright",legend = c(expression(beta~"= 0.25"),expression(beta~"= 0.375"),expression(beta~"= 0.5"),expression(beta~"= 0.625"),expression(beta~"= 0.75")),lty=c(1,2,1,4,5),lwd = 3,col = c(1,2,gray.colors(1, start = 0.8, end = 0.8),4,5),cex=1.3)

# specific setting 

# find for which value of the parameter Pearson skewness is zero
rootf <- function(p){return(exp(-2*p/(1-p))-1-(1-p)/(2*p)*(log((1-p)/(2*p))-1))}
pr <- uniroot(f = rootf,interval = c(0.3,0.5))$root

# parameter
p <- 0.4926946 #(this is pr)
# mode
m <- 0
# function value in mode
sm <- f(m,p)
# density
y <- f(x,p)
# reflected density
ymirror <- f(-x,p)

# gamma(s) function
gamma <- (y-ymirror)/sm
# mode of gamma(s)
MS <- x[which.max(gamma)]-m
# summarizing measure
gammastar <- max(gamma)*sign(MS)

# summarizing measure for sequence of parameter values
p <- seq(0.01,0.99,by=0.001)
gammastarf <- exp(-2*p/(1-p))*1*(p<pr)-(1+(1-p)/(2*p)*(log((1-p)/(2*p))-1))*1*(p>=pr)

# plotting
x11()
par(mai=c(0.8,1,0.1,0.1))
plot(p[1:483],gammastarf[1:483],type="l",lwd=2,xlab="p",ylab=expression(paste(gamma,"*(X)")),cex.lab=1.2,ylim=c(-1,1),xlim=c(p[1],p[981]),cex.axis=1.3,cex.lab=1.3)
lines(p[484:981],gammastarf[484:981],lwd=2)

# range of parameters
ps <- seq(0,1,length.out=10000)
# Pearson skewness for said parameters
sfp <- sf(ps)
# plot
x11(width=6,height=5)
plot(ps,sfp,xlab=expression(beta),ylab="Classical skewness",type="l",lwd=2,cex.axis=1.3,cex.lab=1.3)
# abline(v=0.755)
# abline(h=0)
