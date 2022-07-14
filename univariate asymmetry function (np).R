### functions for calculation:
### -gamma(s)
### -M_X
### -Gamma(X)
### -M_S
### -gamma*(X)
### by means of KDE in a univariate setting
### this solely comprises of non-parametric estimation of the asymmetry measure

### required packages
library(MASS)



NPasymmetry=function(X,h=NULL,factor=2){
  # takes as input 
  # - random sample "X" of size n (numeric vector)
  # - bandwidth "h" (numeric), when NULL, defaults to nrd0 bandwidth
  # - an inflation factor "factor" (numeric) of the bandwidth to ensure unimodality
  
  # returns 
  # - the mode of the empirical density
  # - value of the density at its mode
  # - asymmetry function gamma(s) over a grid x 
  # - mode of gamma(s) (on data axis)
  # - value of gamma(s) at its mode
  # - summarizing asymmetry measure
  # when one wants the asymmetry function over a grid
  
  # setting bandwidth
  if(is.null(h)){
    h=factor*bw.nrd0(X)
  }
  
  # KDE of the supplied sample
  fit=density(x = X,bw = h)
  
  # f(M_X)
  maxf=max(fit$y)
  
  # M_X
  indMX=which.max(fit$y)
  MX=fit$x[indMX]
  
  # creating a grid over which to evaluate the asymmetry function
  l=length(fit$x)
  dx=fit$x[2]-fit$x[1]
  if(indMX>l/2){
    x=seq(fit$x[1],fit$x[1]+(2*indMX-2)*dx,by=dx)
    s=x-MX
    y=c(fit$y,rep(0,2*indMX-l-1))
    ymirror=rev(y)
  } else {
    x=seq(fit$x[l]-(2*(l-indMX))*dx,fit$x[l],by=dx)
    s=x-MX
    y=c(rep(0,l-2*indMX+1),fit$y)
    ymirror=rev(y)
  }
  
  # asymmetry function
  gamma=(y-ymirror)/maxf
  # value of asymmetry function at its mode
  Gamma=max(gamma)
  # mode of asymmetry function
  MS=s[which.max(gamma)]
  # summarizing asymmetry measure
  gammastar=sign(MS)*Gamma
  
  return(list("xgrid"=x,"kde"=y,"MX"=MX,"fMX"=maxf,"s"=s,"gamma"=gamma,
              "MS"=MS,"Gamma"=Gamma,"gammastar"=gammastar))
  
}
