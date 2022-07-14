### function to calculate and plot the density-distance based asymmetry function ###
### for bivariate data. Accepts as input a data matrix containing the            ###
### observations as rows and returns the asymmetry measure and if desired also   ###
### a contourplot                                                                ###
####################################################################################
####################################################################################
####################################################################################

# depends on ks package
library(ks)

### help function used in the main function ###
###############################################

# rotate matrix A 90° clockwise
rotR=function(A){return( t(apply(A,2,rev)) )}

# rotate matrix A 90° counterclockwise
rotL=function(A){return( apply(t(A),2,rev) )}

# transpose matrix A over antidiagonal (not used)
mirrorAD=function(A){return( rotL(t(rotR(A))) )}

# extend 2 partially overlapping grids to fully overlap eachother
newgrid=function(grid1,grid2,indM){
  # grid1 is list containing grid point locations (d vectors) of original data
  # grid2 is list containing grid point locations (d vectors) of mirrored data
  # indM is the index of the mode of the data
  n=unlist(lapply(grid1,length))
  d=length(grid1)
  
  # list containing the joint grid
  newgrid=list()
  # indicator whether original grid is extended to the right or left
  indD=rep(NA,d)
  
  # extending
  for(i in 1:2){
    x=grid1[[i]]
    y=grid2[[i]]
    if(indM[i]>n[i]/2){
      newgrid[[i]]=c(x,y[(which(y>max(x))[1]):n[i]])
      indD[i]="R"
    } else {
      newgrid[[i]]=c(y[1:(which(y>min(x))[1]-2)],x)
      indD[i]="L"
    }
  }
  
  return(list("newgrid"=newgrid,"indD"=indD))
}


### main function ###
#####################

ddAsymmetry2=function(X,gridpoints=c(1000,1000),mingrid=NULL,maxgrid=NULL,plot.contour=F,H.fact=2){
  # X: a nx2 matrix containing the observations
  # gridpoints: numeric vector of length 2 containing the number of grid points to consider in the 2 dimensions
  # mingrid: numeric vector of length 2 containing the lower bounds for the grid
  # maxgrid: numeric vector of length 2 containing the upper bounds for the grid
  # plot.contour: flag if the contourplot of the kde and the asymmetry measure need to be plotted
  # H.fact: numeric inflation factor of the plug-in bandwidth matrix for the kde
  
  
  # if no grid limits are supplied, they are taken to be the range of the data +-2 times the standard deviation
  if(is.null(mingrid) | is.null(maxgrid)){
    varX=var(X)
    minX=apply(X,2,min)
    maxX=apply(X,2,max)
  }
  if(is.null(mingrid)){
    mingrid=minX-2*sqrt(diag(varX))
  }
  if(is.null(maxgrid)){
    maxgrid=maxX+2*sqrt(diag(varX))
  }
  
  
  # plug-in bandwidth
  H=Hpi(X)
  # kernel density estimate
  densest=kde(x = X,H=H.fact*H,gridsize = gridpoints,xmin = mingrid,xmax = maxgrid,compute.cont = plot.contour)
  
  # value at the mode of the kde
  md=max(densest$estimate)
  # index of the mode
  indM=which(densest$estimate == md, arr.ind = TRUE)
  # location (on the data scale) of the mode
  M=rep(NA,2)
  for(i in 1:2){
    M[i]=(densest$eval.points[[i]])[indM[i]]
  }
  indD=c()
  
  step=c(densest$eval.points[[1]][2]-densest$eval.points[[1]][1],densest$eval.points[[2]][2]-densest$eval.points[[2]][1])
  eval.points=list()
  dens=densest$estimate
  for(i in 1:2){
    if(indM[i]>(gridpoints[i]+1)/2){
      eval.points[[i]]=c(densest$eval.points[[i]][1:indM[i]],seq(M[i]+step[i],M[i]+(indM[i]-1)*step[i],by=step[i]))
      indD[i]="R"
    } else if(indM[i]==(gridpoints[i]+1)/2){
      eval.points[[i]]=densest$eval.points[[i]]
      indD[i]="C"
    } else {
      eval.points[[i]]=c(seq(M[i]-(gridpoints[i]-indM[i])*step[i],M[i]-step[i],by=step[i]),densest$eval.points[[i]][indM[i]:gridpoints[i]])
      indD[i]="L"
    }
  }
  
  n=unlist(lapply(eval.points,length))
  Odens=matrix(0,nrow=n[1],ncol=n[2])
  if(indD[1]=="R"){
    if(indD[2]=="R"){
      Odens[1:gridpoints[1],1:gridpoints[2]]=densest$estimate
    } else {
      Odens[1:gridpoints[1],(n[2]-gridpoints[2]+1):n[2]]=densest$estimate
    }
  } else {
    if(indD[2]=="R"){
      Odens[(n[1]-gridpoints[1]+1):n[1],1:gridpoints[2]]=densest$estimate
    } else {
      Odens[(n[1]-gridpoints[1]+1):n[1],(n[2]-gridpoints[2]+1):n[2]]=densest$estimate
    }
  }
  
  Mdens=rotR(rotR(Odens))
  
  # calculate asymmetry measure on the extended grid
  gamma=(Odens-Mdens)/md
  
  # value and index of maximum of gamma
  valG=max(gamma)
  indG=which(gamma == valG, arr.ind = TRUE)
  # location (on the data scale) of the mode
  MD=rep(NA,2)
  for(i in 1:2){
    MD[i]=(eval.points[[i]])[indG[i]]
  }
  sum_measure=valG*(MD-M)/sqrt(sum((M-MD)^2))
  
  # plot contour of kde and asymmetry measure
  if(plot.contour==T){
    x11()
    contour(eval.points[[1]],eval.points[[2]],Odens,main="Kernel density estimate")
    x11()
    contour(eval.points[[1]],eval.points[[2]],Odens,main="Kernel density estimate")
    contour(eval.points[[1]],eval.points[[2]],Mdens,col=2,add=T)
    indM2=which(Odens==max(Odens),arr.ind = T)
    points(eval.points[[1]][indM2[1]],eval.points[[2]][indM2[2]],col=3,pch=19)
    x11()
    contour(eval.points[[1]]-M[1],eval.points[[2]]-M[2],gamma,main="Asymmetry function",xlab=expression(d[1]),ylab=expression(d[2]))
    arrows(x0 = 0,y0 = 0,x1 = sum_measure[1],y1 = sum_measure[2],col = 2,lwd = 1.5)
  }
  
  return(list("kde"=densest,"gamma"=gamma,"grid"=eval.points,"mode"=M,"extend.direction"=indD,"summarizing measure"=sum_measure))
  

}