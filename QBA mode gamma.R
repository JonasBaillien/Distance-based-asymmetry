##################################################
### mode of gamma for linear combinations with ###
### identity matrix A (all A_ii=1 and A_ij=0)  ###
### but with individual scaling phi            ###
##################################################

Modegamma <- function(alpha,A,basefuncs,nu){
  # alpha: d-vector of skewing parameters in (0,1)
  # phi: d-vector containing individual scale parameters
  # basefuncs: d-character vector containing reference densities
  #            options: "t", "laplace", "normal" and "logistic"
  # nu: d-vector with degrees of freedom for possible Stud. t-reference
  
  d=length(alpha)
  MS=rep(NA,d)
  for(i in 1:d){
    MS[i]=switch(basefuncs[i],"Laplace"=slap(alpha=alpha[i],phi=phi[i]),
                              "normal"=snorm(alpha=alpha[i],phi=phi[i]),
                              "logistic"=slog(alpha=alpha[i],phi=phi[i]),
                              "t"=st(alpha=alpha[i],phi[i],nu=nu[i])
                 )
  }
  return(MS)
}


######################################################
### mode of gamma for univariate QBA distributions ###
######################################################

### laplace

# Mode gamma
slap <- function(alpha,phi){
  # alpha: skewing parameter in (0,1)
  # phi: scaling parameter >0

  if(alpha>1/2){
    Ms <- phi/(2*alpha-1)*log((1-alpha)/alpha)
  } else if(alpha<1/2) {
    Ms <- phi/(2*alpha-1)*log(alpha/(1-alpha))
  } else {
    Ms <- 0
  }
  return(Ms)
}

# Gamma
Glap <- function(alpha){
  # alpha: skewing parameter in (0,1)

  if(alpha<1/2){
    G <- (alpha/(1-alpha))^(-(1-alpha)/(2*alpha-1))*((1-2*alpha)/alpha)
  } else if(alpha>1/2){
    G <- ((1-alpha)/alpha)^(-alpha/(1-2*alpha))*((2*alpha-1)/(1-alpha))
  } else {
    G <- 0
  }
  return(G)
}

### normal
# Mode gamma
snorm <- function(alpha,phi){
  # alpha: skewing parameter in (0,1)
  # phi: scaling parameter >0
  
  if(alpha>1/2){
    Ms <- -sqrt(-4*phi^2/(2*alpha-1)*log((1-alpha)/alpha))
  } else if(alpha<1/2) {
    Ms <- sqrt(-4*phi^2/(2*alpha-1)*log((1-alpha)/alpha))
  } else {
    Ms <- 0
  }
  return(Ms)
}

# Gamma
Gnorm <- function(alpha){
  # alpha: skewing parameter in (0,1)
  
  if(alpha<1/2){
    G <- ((1-alpha)/alpha)^(2*alpha^2/(2*alpha-1))*(1-(alpha/(1-alpha))^2)
  } else if(alpha>1/2){
    G <- (alpha/(1-alpha))^(2*(1-alpha)^2/(1-2*alpha))*(1-((1-alpha)/alpha)^2)
  } else {
    G <- 0
  }
  return(G)
}


### logistic
# Mode gamma
slog <- function(alpha,phi){
  # alpha: skewing parameter in (0,1)
  # phi: scaling parameter >0
  
  if(alpha>1/2){
    func <- function(s,alpha,phi){
      exp((2*alpha-1)*s/phi)*(exp(alpha*s/phi)-1)/(exp((1-alpha)*s/phi)-1)*((1+exp((1-alpha)*s/phi))/(1+exp(alpha*s/phi)))^3-(1-alpha)/alpha
    }
    Ms <- uniroot(f = func,interval = c(-100*phi,-0.00004*phi),alpha=alpha,phi=phi)$root
  } else if(alpha<1/2){
    func <- function(s,alpha,phi){
      exp(-(2*alpha-1)*s/phi)*(exp(-alpha*s/phi)-1)/(exp(-(1-alpha)*s/phi)-1)*((1+exp(-(1-alpha)*s/phi))/(1+exp(-alpha*s/phi)))^3-(1-alpha)/alpha
    }
    Ms <- uniroot(f = func,interval = c(0.00004*phi,100*phi),alpha=alpha,phi=phi)$root
  } else {
    Ms <- 0
  }
  return(Ms)
}

# Gamma
Glog <- function(alpha,s,phi,mu){
  # alpha: skewing parameter in (0,1)
  # s: mode of gamma(s)
  # phi: scaling parameter >0
  # mu: location parameter
  
  G <- (QBAsyDist::dALoD(y = mu+s,mu = mu,phi = phi,alpha = alpha)-QBAsyDist::dALoD(y = mu-s,mu = mu,phi = phi,alpha = alpha))/QBAsyDist::dALoD(y = mu,mu = mu,phi = phi,alpha = alpha)
  return(G)
}

### Student's t
# Mode gamma
st <- function(alpha,phi,nu){
  # alpha: skewing parameter in (0,1)
  # phi: scaling parameter >0
  # nu: degrees of freedom
  
  if(alpha>1/2){
    Ms <- -sqrt(phi^2*nu*(((1-alpha)/alpha)^(-4/(nu+3))-1)/(alpha^2-(1-alpha)^2*((1-alpha)/alpha)^(-4/(nu+3))))
  } else if(alpha<1/2) {
    Ms <- sqrt(phi^2*nu*(((1-alpha)/alpha)^(-4/(nu+3))-1)/(alpha^2-(1-alpha)^2*((1-alpha)/alpha)^(-4/(nu+3))))
  } else {
    Ms <- 0
  }
  return(Ms)
}
# Gamma
Gt <- function(alpha,nu){
  # alpha: skewing parameter in (0,1)
  # nu: degrees of freedom
  
  if(alpha<1/2){
    G <- (((1-alpha)/alpha)^((2*nu+2)/(nu+3))-1)*((2*alpha-1)/(alpha^2-((1-alpha)/alpha)^(-4/(nu+3))*(1-alpha)^2))^(-(nu+1)/2)
  } else if(alpha>1/2){
    G <- ((alpha/(1-alpha))^((2*nu+2)/(nu+3))-1)*((1-2*alpha)/((1-alpha)^2-(alpha/(1-alpha))^(-4/(nu+3))*alpha^2))^(-(nu+1)/2)
  } else {
    G <- 0
  }
  return(G)
}


