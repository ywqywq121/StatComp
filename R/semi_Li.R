#' @title Semiparametric partially linear model with Li method
#' @description Use R to do  Semiparametric partially linear model with Li method
#' @param y the response vector
#' @param X the data of parametric part
#' @param Z the data of non-parametric part
#' @param type the type of kernal function
#' @param h the bandwidth of kernal function
#' @param verbose the type of kernal function
#' @return a list of the regression parameters
#' @examples
#' \dontrun{
#' semi_Robinson(y,X,Z,h=1)
#' } 
#' @export

semi_Li=function(y,X,Z,type=1,h,verbose=FALSE){
  X=as.matrix(X)
  Z=as.matrix(Z)
  if(type==1){
    K=function(x,h){
      return((1-abs(x/h))*(x<h)*(x>-h)/h)
    }
    print("Use Triangle kernal function")
  }
  if(type==2){
    K=function(x,h){
      return(3/4*(1-(x/h)^2)*(x<h)*(x>-h)/h)
    }
    
    print("Use Epanechnikow kernal function")
  }
  if(type==3){
    K=function(x,h){
      return(dnorm(x/h)/h)
    }
    
    print("Use Gaussian kernal function")
  }
  if(type>3){
    print("The kernal function is wrongly set, then use defult Triangle kernal function ")
  }
  
  
  n=nrow(X)
  p=ncol(X)
  q=ncol(Z)
  
  if(p>n){
    print("The demision of X is wrong")
  }
  if(q>n){
    print("The demision of Z is wrong")
  }
  
  locus1=1:q
  locus2=1:p
  locus3=1:n
  
  y_save=numeric(n)
  X_save=matrix(nrow=n,ncol=p)
  Z_f=numeric(n)
  
  for(i in 1:n){
    
    Z_save=sweep(Z,2,Z[i,])
    KZ=matrix(nrow=n,ncol=q)
    for(j in 1:q){
      KZ[,j]=vapply(locus3,function(i)K(Z_save[i,j],h),numeric(1))
    }
    Kh= vapply(locus3,function(j)prod(KZ[j,]),numeric(1))
    Z_f[i]=mean(Kh)
    X_median=apply(X , 2 , `*` , Kh)/Z_f[i]
    X_save[i,]=vapply(locus2,function(i)mean(X_median[,i]),numeric(1))
    y_save[i]=t(y)%*%Kh/(n*Z_f[i])
    
  }
  X_tlide=X-X_save
  y_tlide=y-y_save
  
  X_linear=apply(X_tlide , 2 , `*` , Z_f)
  y_linear=y_tlide*Z_f
  linear=lm(y_linear~X_linear)
  beta_tlide=coef(linear)
  P_value=summary(linear)$coefficients
  
  if(verbose==FALSE){
    print("The asymptotic distribution of beta will not be computed for its complicacy")
    return(list(beta=beta_tlide,X_tlide=X_tlide,y_tlide=y_tlide,P_linear=P_value))
  }else{
    Theta=matrix(nrow = p,ncol=p)
    Phi=matrix(nrow = p,ncol=p)
    u=numeric(n)
    for(j in 1:n){
      Theta=Theta+X_tlide[i,]%*%t(X_tlide[i,])*(Z_f[i])^2
      u[i]=y_slide-t(X_slide)%*%beta_tlide
      Phi=Phi+X_tlide[i,]%*%t(X_tlide[i,])*(Z_f[i])^2*u[i]^2
    }
    Sigma=solve(Theta)%*%Phi%*%solve(Theta)
    print("The asymptotic distribution of beta is multivariate normal distribution with the mean vectors being beta_0, the covariance matrix being Sigma, which is in the list")
    return(list(beta=beta_tlide,X_tlide=X_tlide,y_tlide=y_tlide,P_linear=P_value,Sigma=Siama))
  }
}