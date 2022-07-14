# mode of gamma for linear combinations with diagonal matrix A (all A_ii>0 and A_ij=0)
Modegamma=function(alpha,A,basefuncs,nu){
  d=length(alpha)
  MS=rep(NA,d)
  for(i in 1:d){
    MS[i]=switch(basefuncs[i],"Laplace"=slap(alpha=alpha[i],phi=A[i,i]),
                              "normal"=snorm(alpha=alpha[i],phi=A[i,i]),
                              "logistic"=slog(alpha=alpha[i],phi=A[i,i]),
                              "t"=st(alpha=alpha[i],phi=A[i,i],nu=nu[i])
                 )
  }
  return(MS)
}

### mode of gamma for univariate QBA distributions
# laplace
slap=function(alpha,phi){
  # skewness parameter alpha, scale parameter phi
  if(alpha>1/2){
    Ms=phi/(2*alpha-1)*log((1-alpha)/alpha)
  } else if(alpha<1/2) {
    Ms=phi/(2*alpha-1)*log(alpha/(1-alpha))
  } else {
    Ms=0
  }
  return(Ms)
}

Glap=function(alpha){
  if(alpha<1/2){
    G=(alpha/(1-alpha))^(-(1-alpha)/(2*alpha-1))*((1-2*alpha)/alpha)
  } else if(alpha>1/2){
    G=((1-alpha)/alpha)^(-alpha/(1-2*alpha))*((2*alpha-1)/(1-alpha))
  } else {
    G=0
  }
  return(G)
}

# normal
snorm=function(alpha,phi){
  # skewness parameter alpha, scale parameter phi
  if(alpha>1/2){
    Ms=-sqrt(-4*phi^2/(2*alpha-1)*log((1-alpha)/alpha))
  } else if(alpha<1/2) {
    Ms=sqrt(-4*phi^2/(2*alpha-1)*log((1-alpha)/alpha))
  } else {
    Ms=0
  }
  return(Ms)
}

Gnorm=function(alpha){
  if(alpha<1/2){
    G=((1-alpha)/alpha)^(2*alpha^2/(2*alpha-1))*(1-(alpha/(1-alpha))^2)
  } else if(alpha>1/2){
    G=(alpha/(1-alpha))^(2*(1-alpha)^2/(1-2*alpha))*(1-((1-alpha)/alpha)^2)
  } else {
    G=0
  }
  return(G)
}


# logistic
slog=function(alpha,phi){
  # skewness parameter alpha, scale parameter phi
  if(alpha>1/2){
    func=function(s,alpha,phi){
      exp((2*alpha-1)*s/phi)*(exp(alpha*s/phi)-1)/(exp((1-alpha)*s/phi)-1)*((1+exp((1-alpha)*s/phi))/(1+exp(alpha*s/phi)))^3-(1-alpha)/alpha
    }
    Ms=uniroot(f = func,interval = c(-100*phi,-0.00004*phi),alpha=alpha,phi=phi)$root
  } else if(alpha<1/2){
    func=function(s,alpha,phi){
      exp(-(2*alpha-1)*s/phi)*(exp(-alpha*s/phi)-1)/(exp(-(1-alpha)*s/phi)-1)*((1+exp(-(1-alpha)*s/phi))/(1+exp(-alpha*s/phi)))^3-(1-alpha)/alpha
    }
    Ms=uniroot(f = func,interval = c(0.00004*phi,100*phi),alpha=alpha,phi=phi)$root
  } else {
    Ms=0
  }
  return(Ms)
}

Glog=function(alpha,s,phi,mu){
  G=(QBAsyDist::dALoD(y = mu+s,mu = mu,phi = phi,alpha = alpha)-QBAsyDist::dALoD(y = mu-s,mu = mu,phi = phi,alpha = alpha))/QBAsyDist::dALoD(y = mu,mu = mu,phi = phi,alpha = alpha)
  return(G)
}

# Student's t
st=function(alpha,phi,nu){
  # skewness parameter alpha, scale parameter phi
  if(alpha>1/2){
    Ms=-sqrt(phi^2*nu*(((1-alpha)/alpha)^(-4/(nu+3))-1)/(alpha^2-(1-alpha)^2*((1-alpha)/alpha)^(-4/(nu+3))))
  } else if(alpha<1/2) {
    Ms=sqrt(phi^2*nu*(((1-alpha)/alpha)^(-4/(nu+3))-1)/(alpha^2-(1-alpha)^2*((1-alpha)/alpha)^(-4/(nu+3))))
  } else {
    Ms=0
  }
  return(Ms)
}

Gt=function(alpha,nu){
  if(alpha<1/2){
    G=(((1-alpha)/alpha)^((2*nu+2)/(nu+3))-1)*((2*alpha-1)/(alpha^2-((1-alpha)/alpha)^(-4/(nu+3))*(1-alpha)^2))^(-(nu+1)/2)
  } else if(alpha>1/2){
    G=((alpha/(1-alpha))^((2*nu+2)/(nu+3))-1)*((1-2*alpha)/((1-alpha)^2-(alpha/(1-alpha))^(-4/(nu+3))*alpha^2))^(-(nu+1)/2)
  } else {
    G=0
  }
  return(G)
}




### multivariate QBA distribution with A=I ###
##############################################

library(QBAsyDist)
weight=function(basefunc,alpha,phi,s){
  out=switch(basefunc,"normal"=dAND(y = s,alpha = alpha,phi = 1,mu = 0)
             ,"laplace"=dALaD(y = s,alpha = alpha,phi = 1,mu = 0)
             ,"logistic"=dALoD(y = s,alpha = alpha,phi = 1,mu = 0)
             ,"t"=dATD(y = s,alpha = alpha,phi = 1,mu = 0,nu = nu))
  return(out)
}

### mode of gamma for weighted univariate QBA
# laplace
swlap=function(alpha,wplus,wmin){
  # skewness parameter alpha, scale parameter phi, weight for s term wplus, weight for -s term wmin
  if(alpha>1/2){
    Ms=1/(2*alpha-1)*log((1-alpha)*wmin/(alpha*wplus))
  } else if(alpha<1/2) {
    Ms=1/(2*alpha-1)*log(alpha*wmin/((1-alpha)*wplus))
  } else {
    Ms=0
  }
  return(Ms)
}

# normal
swnorm=function(alpha,wplus,wmin){
  # skewness parameter alpha, scale parameter phi, weight for s term wplus, weight for -s term wmin
  if(alpha>1/2){
    Ms=-sqrt(-2/(2*alpha-1)*log((1-alpha)^2*wmin/(alpha^2*wplus)))
  } else if(alpha<1/2) {
    Ms=sqrt(-2/(2*alpha-1)*log((1-alpha)^2*wmin/(alpha^2*wplus)))
  } else {
    Ms=0
  }
  return(Ms)
}


# logistic
swlog=function(alpha,wplus,wmin){
  # skewness parameter alpha, scale parameter phi, weight for s term wplus, weight for -s term wmin
  if(alpha>1/2){
    func=function(s,alpha){
      exp((2*alpha-1)*s)*(exp(alpha*s)-1)/(exp((1-alpha)*s)-1)*((1+exp((1-alpha)*s))/(1+exp(alpha*s)))^3-(1-alpha)*wmin/(alpha*wplus)
    }
    Ms=uniroot(f = func,interval = c(-100,-0.00001),alpha=alpha)$root
  } else if(alpha<1/2){
    func=function(s,alpha){
      exp(-(2*alpha-1)*s)*(exp(-alpha*s)-1)/(exp(-(1-alpha)*s)-1)*((1+exp(-(1-alpha)*s))/(1+exp(-alpha*s)))^3-(1-alpha)*wmin/(alpha*wplus)
    }
    Ms=uniroot(f = func,interval = c(0.00001,100),alpha=alpha)$root
  } else {
    Ms=0
  }
  return(Ms)
}



# Student's t
swt=function(alpha,nu,wplus,wmin){
  # skewness parameter alpha, scale parameter phi, weight for s term wplus, weight for -s term wmin
  if(alpha>1/2){
    Ms=-sqrt(nu*(((1-alpha)/alpha)^(-4/(nu+3))*(wplus/wmin)^(2/(nu+3))-1)/(alpha^2-(1-alpha)^2*((1-alpha)/alpha)^(-4/(nu+3))(wplus/wmin)^(2/(nu+3))))
  } else if(alpha<1/2) {
    Ms=sqrt(phi^2*nu*(((1-alpha)/alpha)^(-4/(nu+3))(wplus/wmin)^(2/(nu+3))-1)/(alpha^2-(1-alpha)^2*((1-alpha)/alpha)^(-4/(nu+3))(wplus/wmin)^(2/(nu+3))))
  } else {
    Ms=0
  }
  return(Ms)
}

optfun=function(s,alpha,basefuncs,nu){
  d=length(alpha)
  tp=1
  tm=1
  for(j in 1:d){
    tp=tp*switch(basefuncs[j],"normal"=dAND(y = s[j],alpha = alpha[j],phi = 1,mu = 0)
               ,"laplace"=dALaD(y = s[j],alpha = alpha[j],phi = 1,mu = 0)
               ,"logistic"=dALoD(y = s[j],alpha = alpha[j],phi = 1,mu = 0)
               ,"t"=dATD(y = s[j],alpha = alpha[j],phi = 1,mu = 0,nu = nu[j]))
    tm=tm*switch(basefuncs[j],"normal"=dAND(y = -s[j],alpha = alpha[j],phi = 1,mu = 0)
                  ,"laplace"=dALaD(y = -s[j],alpha = alpha[j],phi = 1,mu = 0)
                  ,"logistic"=dALoD(y = -s[j],alpha = alpha[j],phi = 1,mu = 0)
                  ,"t"=dATD(y = -s[j],alpha = alpha[j],phi = 1,mu = 0,nu = nu[j]))
  }
  return(tm-tp)
}



# starting values are optimal ones for univariate situation
s_init=function(alpha,basefuncs,nu){
  d=length(alpha)
  s=rep(NA,d)
  for(j in 1:d){
    s[j]=switch(basefuncs[j],"normal"=swnorm(alpha=alpha[j],wplus = 1,wmin = 1)
                  ,"laplace"=swlap(alpha=alpha[j],wplus = 1,wmin = 1)
                  ,"logistic"=swlog(alpha=alpha[j],wplus = 1,wmin = 1)
                  ,"t"=swt(alpha=alpha[j],nu=nu[j],wplus = 1,wmin = 1))
  }
  return(s)
}
basefuncs=c("logistic","laplace")
alpha=c(0.25,0.25)
sstart=s_init(alpha = alpha,basefuncs = basefuncs,nu = rep(NA,length(alpha)))

MS2=optim(par = sstart,fn = optfun,alpha=alpha,basefuncs=basefuncs,nu=rep(NA,2))$par

# 
# # univariate QBA
# library(QBAsyDist) 
# 
# dens=function(x,alpha,mu,phi,nu,basefunc){
#   switch(basefunc,"normal"=dAND(y = x,mu = mu,alpha = alpha,phi = phi)
#          ,"laplace"=dALaD(y = x,mu = mu,alpha = alpha,phi = phi)
#          ,"logistic"=dALoD(y = x,mu = mu,alpha = alpha,phi = phi)
#          ,"t"=dATD(y = x,mu = mu,nu = nu,alpha = alpha,phi = phi))
# }
# 
# basefunc="t"
# alpha=0.76
# mu=2
# nu=4
# phi=0.2
# sm=dens(x = mu,alpha = alpha,mu = mu,phi = phi,nu = nu,basefunc = basefunc)
# 
# n=10000
# x=seq(-5,5,length.out = n)
# y=dens(x = x,alpha = alpha,mu = mu,phi = phi,nu = nu,basefunc = basefunc)
# ymirror=dens(x = 2*mu-x,alpha = alpha,mu = mu,phi = phi,nu = nu,basefunc = basefunc)
# gamma=(y-ymirror)/sm
# Gamma=max(gamma)
# MS=x[which.max(gamma)]-mu
# gammastar=max(gamma)*sign(MS)
# 
# x11()
# plot(x-mu,gamma,type="l",col=1,lwd=2,xlab="s",ylab=expression(gamma[X](s)))
# 
# s=st(alpha = alpha,phi=phi,nu=nu)
# g=Gt(alpha = alpha,nu=nu)




 
# deriv=function(basefunc,s,alpha){
#   if(basefunc=="normal"){
#     if(s<0){
#       out=-(1-alpha)^2*s/sqrt(2*pi)*exp(-1/2*(1-alpha)^2*s^2)
#     } else if(s>0){
#       out=-alpha^2*s/sqrt(2*pi)*exp(-1/2*alpha^2*s^2)
#     } else {
#       out=0
#     }
#   } else if(basefunc=="laplace"){
#     if(s<0){
#       out=(1-alpha)/2*exp(-abs(-(1-alpha)*s))
#     } else if(s>0){
#       out=-alpha/2*exp(-abs(alpha*s))
#     } else {
#       out=0
#     }
#   } else if(basefunc=="logistic"){
#     if(s<0){
#       out=-(1-alpha)*exp((1-alpha)*s)*(exp((1-alpha)*s)-1)/(1+exp((1-alpha)*s))^3
#     } else if(s>0){
#       out=alpha*exp(-alpha*s)*(exp(-alpha*s)-1)/(1+exp(-alpha*s))^3
#     } else {
#       out=0
#     }
#   } else if(basefunc=="t"){
#     if(s<0){
#       out=-(nu+1)*gamma((nu+1)/2)/(sqrt(pi*nu)*gamma(nu/2))*(1-alpha)^2*s/nu*(1+(1-alpha)^2*s^2/nu)^(-(nu+3)/2)
#     } else if(s>0){
#       out=-(nu+1)*gamma((nu+1)/2)/(sqrt(pi*nu)*gamma(nu/2))*alpha^2*s/nu*(1+alpha^2*s^2/nu)^(-(nu+3)/2)
#     } else {
#       out=0
#     }
#   }
#   return(out)
# }
# 
# rootfun=function(s,basefunc,alpha,wplus,wmin){
#   return(deriv(basefunc = basefunc,s = s,alpha = alpha)*wplus+deriv(basefunc = basefunc,s = -s,alpha = alpha)*wmin)
# }
