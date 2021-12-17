## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,warning = FALSE, message = FALSE)

## -----------------------------------------------------------------------------
x=seq(-10,10,0.01)
plot(x,exp(((-1/2)*x^2))/sqrt(2*pi),xlim=c(-11,11), ylim=c(0,1), main="标准正态图 ",  
xlab="x", ylab="y")

## -----------------------------------------------------------------------------
x <- c(1,2,3,4)
y <-c(2,3,4,5)
x_html=knitr::kable((rbind(x,y)),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)

## -----------------------------------------------------------------------------
set.seed(123)
Ray=function(sigma,n){
  U=runif(n)
  X=sqrt(-2*sigma^2*log(1-U))
  
  hist(X,probability = TRUE,col="pink",main=sigma)
  y=seq(0,10,0.1)
  lines(y,(y/sigma)*exp(-y^2/(2*sigma^2)))
  
  }
par(mfrow=c(1,1))
for(i in seq(0.5,1.4,0.1)){Ray(i,10000)}


## -----------------------------------------------------------------------------
set.seed(1234)
mix=function(p1){
n <- 1e4
X1 <- rnorm(n)
X2 <- rnorm(n,mean=3)
r <- sample(c(0,1),n,replace=TRUE,prob=c(1-p1,p1))
Z <- r*X1+(1-r)*X2
hist(Z,probability = TRUE,col="Green",main=p1)
}
mix(0.75)
par(mfrow=c(1,1))
for(i in seq(0.1,0.9,0.05)){mix(i)}

## -----------------------------------------------------------------------------
set.seed(12345)

PGprocess=function(n,lamda,shape,size){
X=numeric(n)
for(t in 1:n){
  N=rpois(1,lambda=lamda*t)
  y=rgamma(N,shape,size)
  X[t]=sum(y)}
  return(X)
}
x1000=PGprocess(1000,1,1,1)
plot(1:1000,x1000,xlab="t",ylab="X(t)",main="compound Poisson(1)–Gamma（1,1） process",type="l")

## -----------------------------------------------------------------------------
PG10=function(lamda,shape,size){ 
  N=rpois(10000,lambda=lamda*10)
  x10=numeric(10000)
  for(i in 1:10000){y=rgamma(N[i],shape,size)
  x10[i]=sum(y)}
  
  est.mean=mean(x10)
  est.var=var(x10)
  real.mean=10*lamda*shape/size
  real.var=10*lamda*((shape/size)^2+shape/(size^2))
  res=array(c(est.mean,est.var,real.mean,real.var),dim=c(2,2))
  dimnames(res)[[1]]=c("est.mean","est.var")
  dimnames(res)[[2]]=c("real.mean","real.var")
  print(res)}
for(i in 1:3){PG10(i,1,1)} ##固定Gamma(1,1),possion参数从1取到3.
for(i in 1:3){PG10(i,1/2,2)} ##固定Gamma(1/2,2),possion参数从1取到3.
for(i in 1:3){PG10(i,2,2)} ##固定Gamma(2,2),possion参数从1取到3.
for(i in 1:3){PG10(i,3,3)} ##固定Gamma(3,3),possion参数从1取到3.

## -----------------------------------------------------------------------------
set.seed(1234)
CDF=function(x,n,alpha,beta){
  z=runif(n,0,x)
  y=30*z^(alpha-1)*(1-z)^(beta-1)*x
  cdf.x=mean(y)
  return(cdf.x)
}
cdf=numeric(99)
for(i in seq(0.01,0.99,0.01)){cdf[100*i]=CDF(i,10000,3,3)}
plot(seq(0.01,0.99,0.01),pbeta(seq(0.01,0.99,0.01),3,3),main="the CDF from 0.01 to 0.99",xlab="x",ylab="F(X)",type="l",col="red")
lines(seq(0.01,0.99,0.01),cdf,col="blue")
legend("bottomright", c("pbeta", "MC estimate"), pch = "", col=c("red", "blue"), lwd = 1)
MC=numeric(9)
for(i in seq(0.1,0.9,0.1)){MC[10*i]=CDF(i,10000,3,3)}
X=seq(0.1,0.9,0.1)
pbeta=pbeta(seq(0.1,0.9,0.1),3,3)
library(kableExtra)
x_html=knitr::kable((rbind(X,pbeta,MC)),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)#set the size of words in the table




## -----------------------------------------------------------------------------
set.seed(123)
Ray = function(sig=sigma,R = 10000, antithetic = TRUE) {
  u = runif(R)
  if (!antithetic){
  v=runif(R)
  x1=sig*sqrt(-2*log(1-u))
  x2=sig*sqrt(-2*log(1-v))
  
  ray=(x1+x2)/2
  }  
  else{
  v = 1 - u
  x1=sig*sqrt(-2*log(1-u))
  x2=sig*sqrt(-2*log(1-v))
  ray=(x1+x2)/2
  } 
 
  return(ray)
}

sd.norm=numeric(10)
for(i in 1:10){ray.norm=Ray(sig=i,R = 10000,antithetic = FALSE)
sd.norm[i]=var(ray.norm)}
sd.anti=numeric(10)
for(i in 1:10){ray.anti=Ray(sig=i,R = 10000,antithetic = TRUE)
sd.anti[i]=var(ray.anti)}
library(kableExtra)
sigma=1:10
percent.reduction=(sd.norm-sd.anti)/sd.norm
x_html=knitr::kable((rbind(sigma,sd.norm,sd.anti,percent.reduction)),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)#set the size of words in the table

## -----------------------------------------------------------------------------
    x <- seq(0, 1, .01)
    w <- 2

    g <- x^2*exp(-x^2/2) / sqrt(2*pi)
    f0 <- rep(1,length(x))
    f1 <- exp(-x)
    f2 <- (1 / pi) / (1 + x^2)
    f3 <- exp(-x) / (1 - exp(-1))
    f4 <- 4 / ((1 + x^2) * pi)
    gs <- c(expression(g(x)==e^{-x}/(1+x^2)),
            expression(f[0](x)==1),
            expression(f[1](x)==e^{-x}),
            expression(f[2](x)==1/pi/(1+x^2)),
            expression(f[3](x)==e^{-x}/(1-e^{-1})),
            expression(f[4](x)==4/((1+x^2)*pi)))
    #for color change lty to col
    par(mfrow=c(1,1))
    #figure (a)
    plot(x, g, type = "l", ylab = "",
         ylim = c(0,2), lwd = w,col=1,main='(A)')
    lines(x, f0, lty = 2, lwd = w,col=2)
    lines(x, f1, lty = 3, lwd = w,col=3)
    lines(x, f2, lty = 4, lwd = w,col=4)
    lines(x, f3, lty = 5, lwd = w,col=5)
    lines(x, f4, lty = 6, lwd = w,col=6)
    legend("topright", legend = gs,
           lty = 1:6, lwd = w, inset = 0.02,col=1:6)

    #figure (b)
    plot(x, g/f0, type = "l", ylab = "",
        ylim = c(0,3.2), lwd = w, lty = 2,col=2,main='(B)')
    lines(x, g/f1, lty = 3, lwd = w,col=3)
    lines(x, g/f2, lty = 4, lwd = w,col=4)
    lines(x, g/f3, lty = 5, lwd = w,col=5)
    lines(x, g/f4, lty = 6, lwd = w,col=6)
    legend("topright", legend = gs[-1],
           lty = 2:6, lwd = w, inset = 0.02,col=2:6)


## -----------------------------------------------------------------------------
set.seed(123)
  m <- 10000
  est <- sd <- numeric(5)
  g <- function(x) {
  x^2*exp(-x^2/2) / sqrt(2*pi) * (x > 0) * (x < 1)
  }
  x <- runif(m) #using f0
  fg <- g(x)
  est[1] <- mean(fg)
  sd[1] <- sd(fg)
  x <- rexp(m, 1) #using f1
  fg <- g(x) / exp(-x)
  est[2] <- mean(fg)
  sd[2] <- sd(fg)
  x <- rcauchy(m) #using f2
  i <- c(which(x > 1), which(x < 0))
  x[i] <- 2 #to catch overflow errors in g(x)
  fg <- g(x) / dcauchy(x)
  est[3] <- mean(fg)
  sd[3] <- sd(fg)
  u <- runif(m) #f3, inverse transform method
  x <- - log(1 - u * (1 - exp(-1)))
  fg <- g(x) / (exp(-x) / (1 - exp(-1)))
  est[4] <- mean(fg)
  sd[4] <- sd(fg)
  u <- runif(m) #f4, inverse transform method
  x <- tan(pi * u / 4)
  fg <- g(x) / (4 / ((1 + x^2) * pi))
  est[5] <- mean(fg)
  sd[5] <- sd(fg)
  est.real=0.5-est
x_html=knitr::kable((cbind(est.real,sd)),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)



## -----------------------------------------------------------------------------
set.seed(123)
  m <- 10000
  est <- sd <- numeric(5)
  g <- function(x) {
  x^2*exp(-x^2/2) / sqrt(2*pi) * (x > 0) * (x < 1)
  }
  x <- runif(m) #using f0
  fg <- g(x)
  est[1] <- mean(fg)
  sd[1] <- sd(fg)
  x <- rexp(m, 1) #using f1
  fg <- g(x) / exp(-x)
  est[2] <- mean(fg)
  sd[2] <- sd(fg)
  x <- rcauchy(m) #using f2
  i <- c(which(x > 1), which(x < 0))
  x[i] <- 2 #to catch overflow errors in g(x)
  fg <- g(x) / dcauchy(x)
  est[3] <- mean(fg)
  sd[3] <- sd(fg)
  u <- runif(m) #f3, inverse transform method
  x <- - log(1 - u * (1 - exp(-1)))
  fg <- g(x) / (exp(-x) / (1 - exp(-1)))
  est[4] <- mean(fg)
  sd[4] <- sd(fg)
  u <- runif(m) #f4, inverse transform method
  x <- tan(pi * u / 4)
  fg <- g(x) / (4 / ((1 + x^2) * pi))
  est[5] <- mean(fg)
  sd[5] <- sd(fg)
  est.real=0.5-est
x_html=knitr::kable((cbind(est.real,sd)),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)



## -----------------------------------------------------------------------------
set.seed(123)
np=0
for(i in 1:100000){
  x=rchisq(20,2)
  if(abs(mean(x)-2)<sqrt(var(x)/length(x))*qt(0.975,(length(x)-1))){np=np+1}
}
p=np/100000
p

## -----------------------------------------------------------------------------
set.seed(123)
np=0
for(i in 1:100000){
  x=rchisq(20,2)
  if(4<(length(x)-1)*var(x)/qchisq(0.95,(length(x)-1))){np=np+1}
}
p=np/100000
p

## -----------------------------------------------------------------------------
set.seed(12)
n=100
m=10000
alpha=0.05
mu0=1
p1=numeric(m)
p2=numeric(m)
p3=numeric(m)
for(i in 1:m){
  x1=rchisq(n,1)
  x2=runif(n,0,2)
  x3=rexp(n)
  ttest1 <- t.test(x1, alternative = "two.sided", mu = mu0,conf.level = 1-alpha)
  ttest2 <- t.test(x2, alternative = "two.sided", mu = mu0,conf.level = 1-alpha)
  ttest3 <- t.test(x3, alternative = "two.sided", mu = mu0,conf.level = 1-alpha)
  p1[i]=ttest1$p.value
  p2[i]=ttest2$p.value
  p3[i]=ttest3$p.value
}
Ierror.chisq=mean(p1<alpha)
Ierror.unif=mean(p2<alpha)
Ierror.exp=mean(p3<alpha)
x_html=knitr::kable((cbind(Ierror.chisq,Ierror.unif,Ierror.exp)),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)


## -----------------------------------------------------------------------------
set.seed(123)
critical=function(d){
  return(qchisq(0.95,df=d*(d+1)*(d+2)/6))
}
cvalue=c(critical(1),critical(2),critical(3))
statis<- function(x) {
  #computes the sample skewness .
  statistic=0
  sigma=cov(x,x)
  n=dim(x)[1]
  d=dim(x)[2]
  mat=(scale(x, center = TRUE, scale = TRUE))%*%t(scale(x, center = TRUE, scale = TRUE))
  statistic=sum(mat^3)
  return( statistic/n/6 )
}
n=c(10,20,30,50,100,500)
m=10000
Ierror=matrix(nrow=length(n),ncol=3)
for(i in 1:length(n)){
  for(d in 1:3){
    result=numeric(m)
    for(j in 1:m){
      x=array(rnorm(n[i]*d),dim=c(n[i],d))
      result[j]=as.integer(abs(statis(x)) >=cvalue[d] )
    }
    Ierror[i,d]=mean(result)
  }
}
dimnames(Ierror)[[1]]=c("size=10","size=20","size=30","size=50","size=100","size=500")
dimnames(Ierror)[[2]]=c("d=1","d=2","d=3")
print(Ierror)



## -----------------------------------------------------------------------------
set.seed(12)
alpha <- .1
n <- 30
m <- 2500
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
pwr1 <- numeric(N)
pwr2 <- numeric(N)
pwr3<- numeric(N)
critical=function(d){
  return(qchisq(0.95,df=d*(d+1)*(d+2)/6))
}
cvalue=c(critical(1),critical(2),critical(3))
statis<- function(x) {
  #computes the sample skewness .
  statistic=0
  sigma=cov(x,x)
  n=dim(x)[1]
  mat=(scale(x, center = TRUE, scale = TRUE))%*%t(scale(x, center = TRUE, scale = TRUE))
  statistic=sum(mat^3)
  return( statistic/n/6 )
}
for (j in 1:N) { #for each epsilon
 e <- epsilon[j]
 sktests1 <- numeric(m)
  sktests2 <- numeric(m)
   sktests3 <- numeric(m)
 for (i in 1:m) { #for each replicate
 sig1 <- sample(c(1, 10), replace = TRUE,size = n*1, prob = c(1-e, e))
 x=array(rnorm(30*1,0,sig1),dim=c(30,1))
 sktests1[i] <- as.integer(statis(x) >= critical(1))
  sig2 <- sample(c(1, 10), replace = TRUE,size = n*2, prob = c(1-e, e))
 x=array(rnorm(30*2,0,sig2),dim=c(30,2))
 sktests2[i] <- as.integer(statis(x) >= critical(2))
  sig3 <- sample(c(1, 10), replace = TRUE,size = n*3, prob = c(1-e, e))
 x=array(rnorm(30*3,0,sig3),dim=c(30,3))
 sktests3[i] <- as.integer(statis(x) >= critical(3))
}
pwr1[j] <- mean(sktests1)
pwr2[j] <- mean(sktests2)
pwr3[j] <- mean(sktests3)
}
#plot po
plot(epsilon, pwr1, type = "b",
xlab = bquote(epsilon), ylim = c(0,1),main="Power of demension 1")
abline(h = .05, lty = 3)
plot(epsilon, pwr2, type = "b",
xlab = bquote(epsilon), ylim = c(0,1),main="Power of demension 2")
abline(h = .05, lty = 3)
plot(epsilon, pwr3, type = "b",
xlab = bquote(epsilon), ylim = c(0,1),main="Power of demension 3")
abline(h = .05, lty = 3)




## -----------------------------------------------------------------------------
set.seed(123)
library(boot)
library(bootstrap)
lambda_hat <- eigen(cov(scor))$values
theta_hat <- lambda_hat[1] / sum(lambda_hat)
B <-5000 
n <- nrow(scor) 
func <- function(dat, index){
  x <- dat[index,]
  lambda <- eigen(cov(x))$values
  theta <- lambda[1] / sum(lambda)
  return(theta)
}
bootstrap_result <- boot(
  data = cbind(scor$mec, scor$vec, scor$alg, scor$ana, scor$sta),
  statistic = func, R = B)
theta_b <- bootstrap_result$t
bias_boot <- mean(theta_b) - theta_hat
se_boot <- sd(theta_b)
x_html=knitr::kable((rbind(theta_hat,bias_boot,se_boot)),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)

## -----------------------------------------------------------------------------
theta_j <- rep(0, n)
for (i in 1:n) {
  x <- scor [-i,]
  lambda <- eigen(cov(x))$values
  theta_j[i] <- lambda[1] / sum(lambda)
}
bias_jack <- (n - 1) * (mean(theta_j) - theta_hat)
se_jack <- (n - 1) * sqrt(var(theta_j) / n)
x_html=knitr::kable((rbind(bias_jack,se_jack)),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)

## -----------------------------------------------------------------------------
set.seed(123)
library(boot)
library(bootstrap)
lambda_hat <- eigen(cov(scor))$values
theta_hat <- lambda_hat[1] / sum(lambda_hat)
B <-5000 
n <- nrow(scor) 
func <- function(dat, index){
  x <- dat[index,]
  lambda <- eigen(cov(x))$values
  theta <- lambda[1] / sum(lambda)
  return(theta)
}
bootstrap_result <- boot(
  data = cbind(scor$mec, scor$vec, scor$alg, scor$ana, scor$sta),
  statistic = func, R = B)
ci <- boot.ci(bootstrap_result,type=c("norm","basic","perc","bca"))
ci

## -----------------------------------------------------------------------------
set.seed(1234)
library(boot)
n=100
B=1000

func=function(x,i){
  mu2=var(x[i,])
  mu3=mean((x[i,]-mean(x[i,]))^3)
  skew=mu3/mu2^(3/2)
  return(skew)
}
cvoer.rat.norm=function(n=100,B=1000,m=5000,a){
  cr_norm1 = numeric(m)
  left_norm1 = numeric(m)
  right_norm1 = numeric(m)
  cr_perc1 = numeric(m)
  left_perc1 = numeric(m)
  right_perc1 = numeric(m)
  cr_basic1 = numeric(m) 
  left_basic1 = numeric(m)
  right_basic1 = numeric(m)
  for(i in 1:m){
    x = as.matrix(rnorm(n))
    ci = boot.ci(boot(data = x, statistic = func, R =B), conf = a, type = c("norm", "basic", "perc"))
    cr_norm1[i] = ci$norm[2] <= 0 && ci$norm[3] >= 0
    left_norm1[i] = ci$norm[2]>0
     right_norm1[i]=ci$norm[3]<0
    cr_perc1[i] = ci$perc[4] <= 0 && ci$perc[5] >= 0
       left_perc1[i] = ci$perc[4]>0
     right_perc1[i]=ci$perc[5]<0
    cr_basic1[i] = ci$basic[4] <= 0 && ci$basic[5] >= 0
    left_basic1[i] = ci$basic[4]>0
    right_basic1[i]=ci$basic[5]<0
  }
  norm.norm=mean(cr_norm1)
  norm.norm.left=mean(left_norm1)
  norm.norm.right=mean(right_norm1)
  norm.perc=mean(cr_perc1)
  norm.perc.left=mean(left_perc1)
  norm.perc.right=mean(right_perc1)
  norm.basic=mean(cr_basic1)
  norm.basic.left=mean(left_basic1)
  norm.basic.right=mean(right_basic1)
return(c(norm.norm,norm.norm.left,norm.norm.right,norm.perc,norm.perc.left,norm.perc.right,norm.basic, norm.basic.left, norm.basic.right))
   }
cvoer.rat.chisq=function(n=100,B=1000,m=5000,a){
 cr_norm2 = numeric(m)
 left_norm2 = numeric(m)
  right_norm2 = numeric(m)
  cr_perc2 = numeric(m)
  left_perc2 = numeric(m)
  right_perc2 = numeric(m)
  cr_basic2 = numeric(m) 
  left_basic2 = numeric(m)
  right_basic2 = numeric(m)
  for(i in 1:m){
    Y = as.matrix(rchisq(n,df=5))
    ci = boot.ci(boot(data = Y, statistic = func, R =B), conf =a, type = c("norm", "basic", "perc"))
    cr_norm2[i] = ci$norm[2] <= 2*sqrt(2/5) && ci$norm[3] >= 2*sqrt(2/5)
    left_norm2[i] = ci$norm[2]>2*sqrt(2/5)
     right_norm2[i]=ci$norm[3]<2*sqrt(2/5)
    cr_perc2[i] = ci$perc[4] <= 2*sqrt(2/5) && ci$perc[5] >= 2*sqrt(2/5)
    left_perc2[i] = ci$perc[4]>2*sqrt(2/5)
     right_perc2[i]=ci$perc[5]<2*sqrt(2/5)
    cr_basic2[i] = ci$basic[4] <= 2*sqrt(2/5) && ci$basic[5] >= 2*sqrt(2/5)
    left_basic2[i] = ci$basic[4]>2*sqrt(2/5)
     right_basic2[i]=ci$basic[5]<2*sqrt(2/5)
  }
  chisq.norm=mean(cr_norm2)
  chisq.norm.left=mean(left_norm2)
  chisq.norm.right=mean(right_norm2)
  chisq.perc=mean(cr_perc2)
  chisq.perc.left=mean(left_perc2)
  chisq.perc.right=mean(right_perc2)
  chisq.basic=mean(cr_basic2)
  chisq.basic.left=mean(left_basic2)
  chisq.basic.right=mean(right_basic2)
return(c(chisq.norm,chisq.norm.left,chisq.norm.right,chisq.perc,chisq.perc.left,chisq.perc.right,chisq.basic, chisq.basic.left, chisq.basic.right))

}

compare = function(n=100,B=1000,m=5000,a){
  # compare the results and construct a table
  r1 = cvoer.rat.norm(n=n, B=B, m=m, a = a)
  r2 = cvoer.rat.chisq(n=n, B=B, m=m, a = a)
  result = cbind(r1,r2)
  colnames(result) = c("normal","chi sqaure")
  rownames(result) = c("normal","normal.left","normal.right","precentile","precentile.left","precentile.right","basic","basic.left","basic.right")
  x_html=knitr::kable((result),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)
  }
compare(n=100,B=1000,m=5000,a=0.9)
compare(n=100,B=1000,m=5000,a=0.95)




## -----------------------------------------------------------------------------
set.seed(123)
library(boot)
x=as.numeric(iris[1:50, 2])
y=as.numeric(iris[1:50, 3])
R <- 999;z <- c(x,y);K <- 1:100;n<-length(x)
reps <- numeric(R);t0 <-cor.test(x,y,method="spearman")$estimate
for (i in 1:R) {
k <- sample(K, size = n, replace = FALSE)
x1 <- z[k];y1 <- z[-k] #complement of x1
reps[i] <- cor.test(x1,y1,method="spearman")$estimate
}
p.permutation <-mean(abs(c(t0, reps)) >= abs(t0))

pcortest=cor.test(x,y,method="spearman")$p.value
result=cbind(p.permutation,pcortest)
x_html=knitr::kable((result),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)

## -----------------------------------------------------------------------------


library(energy)

library(Ball)
library(RANN)
set.seed(1234)
m <- 1e2; k<-3; p<-2; mu <- 0; set.seed(12345)
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}


eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
x <- matrix(rnorm(n1*p,0,1.7),ncol=p);
y <- cbind(rnorm(n2),rnorm(n2,mean=mu,1.5));
z <- rbind(x,y)
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*1235)$p.value
}
alpha <- 0.1;
pow <- colMeans(p.values<alpha)
NNpower=pow[1]
energypower=pow[2]
bdpower=pow[3]
result=cbind(NNpower,energypower,bdpower)
x_html=knitr::kable((result),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)


## -----------------------------------------------------------------------------
library(energy)
library(Ball)
library(RANN)
set.seed(124)
m <- 1e2; k<-3; p<-2; mu <- .5; set.seed(12345)
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}


eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
x <- matrix(rnorm(n1*p,0,1.5),ncol=p);
y <- cbind(rnorm(n2,0.2,1.3),rnorm(n2,mean=mu,1.4));
z <- rbind(x,y)
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}
alpha <- 0.1;
pow <- colMeans(p.values<alpha)
NNpower=pow[1]
energypower=pow[2]
bdpower=pow[3]
result=cbind(NNpower,energypower,bdpower)
x_html=knitr::kable((result),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)

## -----------------------------------------------------------------------------
library(energy)
library(Ball)
library(RANN)
set.seed(1234)
m <- 1e2; k<-3; p<-2; mu <- .5; set.seed(12345)
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}
rbimodel=function(n,mean1,sd1,mean2,sd2,pro){
   p=sample(c(1,0),n,replace = TRUE,prob=c(pro,1-pro))
  x=p*rnorm(n,mean1,sd1)+(1-p)*rnorm(n,mean2,sd2)
  return(x)
}

eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
x <- matrix(rt(n1*p,df=1),ncol=p);
y <- cbind(rbimodel(n1,-1,1,1,2,0.5),rbimodel(n2,-1,1,1,2,0.5))
z <- rbind(x,y)
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}
alpha <- 0.1;
pow <- colMeans(p.values<alpha)
NNpower=pow[1]
energypower=pow[2]
bdpower=pow[3]
result=cbind(NNpower,energypower,bdpower)
x_html=knitr::kable((result),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)

## -----------------------------------------------------------------------------
library(energy)
library(Ball)
library(RANN)
set.seed(1234)
m <- 1e2; k<-3; p<-2; mu <- .5; set.seed(12345)
n1 <- 100;n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}
rbimodel=function(n,mean1,sd1,mean2,sd2,pro){
   p=sample(c(1,0),n,replace = TRUE,prob=c(pro,1-pro))
  x=p*rnorm(n,mean1,sd1)+(1-p)*rnorm(n,mean2,sd2)
  return(x)
}

eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)
for(i in 1:m){
x <- matrix(rnorm(n1*p,0,1.5),ncol=p);
y <- cbind(rnorm(n2,0.2,1.3),rnorm(n2,mean=mu,1.4));
z <- rbind(x,y)
p.values[i,1] <- eqdist.nn(z,N,k)$p.value
p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
}
alpha <- 0.1;
pow <- colMeans(p.values<alpha)
NNpower=pow[1]
energypower=pow[2]
bdpower=pow[3]
result=cbind(NNpower,energypower,bdpower)
x_html=knitr::kable((result),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)

## -----------------------------------------------------------------------------
set.seed(123)
 Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
        psi <- as.matrix(psi)
        n <- ncol(psi)
        k <- nrow(psi)

        psi.means <- rowMeans(psi)     #row means
        B <- n * var(psi.means)        #between variance est.
        psi.w <- apply(psi, 1, "var")  #within variances
        W <- mean(psi.w)               #within est.
        v.hat <- W*(n-1)/n + B/n+(B/(n*k))     #upper variance est.
        r.hat <- v.hat / W             #G-R statistic
        return(r.hat)
        }
rw.Metropolis <- function(n, sigma, x0, N) {
       # n: degree of freedom of t distribution
       # sigma:  standard variance of proposal distribution N(xt,sigma)
       # x0: initial value
       # N: size of random numbers required.
        x <- numeric(N)
        x[1] <- x0
        u <- runif(N)
        k <- 0
        for (i in 2:N) {
            y <- rnorm(1, x[i-1], sigma)
                if (u[i] <= (dt(y, n) / dt(x[i-1], n)))
                    x[i] <- y  
                else {
                    x[i] <- x[i-1]
                    k <- k + 1
                }
            }
        return(list(x=x, k=k))
}
N<- 8000     #length of chains
n=1
b <- 1000       #burn-in length
x0 <- c(-8, -5, 5, 8)
k=4
sigma=c(1,2.5,0.5,4)
X <- matrix(0, nrow=k, ncol=N)
 par(mfrow=c(1,1)) 
for(j in 1:4){
      for (i in 1:k){
        X[i, ] <- rw.Metropolis(n,sigma[j], x0[i],N)$x}
    #trace plots
     plot(1:N,X[1,],type="l",ylim=c(-30,30),xlab=bquote(sigma == .(round(sigma[j],3))),
             ylab="X")
     lines(1:N,X[2,],type="l",col=2)
     lines(1:N,X[3,],type="l",col=3)
     lines(1:N,X[4,],type="l",col=4)}

## -----------------------------------------------------------------------------
xt<- matrix(0, nrow=k, ncol=N)
for (i in 1:k){
        xt[i, ] <- rw.Metropolis(n,1, x0[i],N)$x}
psi <- t(apply(xt, 1, cumsum))
for (i in 1:nrow(psi)){psi[i,] <- psi[i,] / (1:ncol(psi))}

print(Gelman.Rubin(psi))

## -----------------------------------------------------------------------------
par(mfrow=c(1,1)) #reset to default
 a<- c(.05,seq(.1,.9,.1),.95)
    Q<- qt(a,n)
    
    mc<-cbind(xt[2,1001:N])
    Qrw<- apply(mc, 2,function(x) quantile(x,a))                 
    x_html=knitr::kable((round(cbind(Q, Qrw), 3)),"html")
   kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)

## -----------------------------------------------------------------------------
set.seed(123) 
N <- 12000 
index=1:N
burn <- 1000            #burn-in length
x0=rbind(c(5,0.6),c(6,0.7),c(7,0.8),c(8,0.5))
fxy=function(N,x0){
x<- matrix(0, N, 2)
x[1,]=x0
for (i in 2:N) {
       
         x[i,1]=rbinom(1,size=16,prob=x[i-1,2])
        x[i,2] = rbeta(1, shape1=x[i,1]+2, shape2=16-x[i,1]+4)
}
return(x)}
b <- burn + 1
xt=fxy(N,x0[3,])  
xt.b <- fxy(N,x0[3,])[b:N, ]
colMeans(xt.b)
    cov(xt.b)
    cor(xt.b)

    plot(xt, main="", cex=.5, xlab=bquote(X),
         ylab=bquote(Y))

## -----------------------------------------------------------------------------
xt1=rbind(fxy(N,x0[1,])[,1],fxy(N,x0[2,])[,1] ,fxy(N,x0[3,])[,1] ,fxy(N,x0[4,])[,1]  )
xt2=rbind(fxy(N,x0[1,])[,2],fxy(N,x0[2,])[,2] ,fxy(N,x0[3,])[,2] ,fxy(N,x0[4,])[,2]  )
psi1 <- t(apply(xt1, 1, cumsum))
psi2 <- t(apply(xt2, 1, cumsum))
for (i in 1:nrow(psi1)){psi1[i,] <- psi1[i,] / (1:ncol(psi1))
psi2[i,] <- psi2[i,] / (1:ncol(psi2))}

print(Gelman.Rubin(psi1))
print(Gelman.Rubin(psi2))

erg.mean<-function(x){ # compute ergodic mean 
        n<-length(x)
        result<-cumsum(x)/cumsum(rep(1,n))
      }
  
   par( mfrow=c(1,1))
   
   index<-1:N
   index2<-(b+1):N

   plot(index, xt[index,1], type='l', ylab='Values of x', xlab='Iterations', main='(a) Trace Plot of x')
   plot(index, xt[index,2], type='l', ylab='Values of y', xlab='Iterations', main='(b) Trace Plot of y')

   ergtheta0<-erg.mean( xt[index,1] )
   ergtheta02<-erg.mean( xt[index2,1] )
   ylims0<-range( c(ergtheta0,ergtheta02) )

   ergtheta1<-erg.mean( xt[index,2] )
   ergtheta12<-erg.mean( xt[index2,2] )
   ylims1<-range( c(ergtheta1,ergtheta12) )

   step<-10
   index3<-seq(1,N,step)
   index4<-seq(b+1,N,step)

   plot(index3 , ergtheta0[index3], type='l', ylab='Values of x', xlab='Iterations', main='(c) Ergodic Mean Plot of x')
   lines(index4, ergtheta02[index4-b], col=2, lty=2)

   plot(index3, ergtheta1[index3], type='l', ylab='Values of y', xlab='Iterations', main='(d) Ergodic Mean Plot of y' )
   lines(index4, ergtheta12[index4-b], col=2, lty=2)

   acf(xt[index2,1], main='Autocorrelations Plot for x')
   acf(xt[index2,2], main='Autocorrelations Plot for y') 






## -----------------------------------------------------------------------------
kth.term=function(a,k){
  if(k!=0){d=length(a)
  a1=as.numeric(t(a)%*%a)  ##calculate ||a||^2
  s1=(k+1)*log(a1)
  s2=lgamma((d+1)/2)
  s3=lgamma(k+3/2)
  s4=sum(log(1:k))
  s5=k*log(2)
  s6=log(2*k+1)
  s7=log(2*k+2)
  s8=lgamma(k+d/2+1)
  s=s1+s2+s3-s4-s5-s6-s7-s8
  return((-1)^(k%%2)*exp(s))}
  else{return(as.numeric(t(a)%*%a)/2*gamma((length(a)+1)/2)*gamma(3/2)/gamma(length(a)/2+1))}
}

## -----------------------------------------------------------------------------
paste("When k =5, a=c(1,2) the result is",kth.term(c(1,2),5))
paste("When k =30, a=c(1,2) the result is",kth.term(c(1,2),30))

## -----------------------------------------------------------------------------
sum.term=function(a){
  s=0
  k=0
  while(abs(kth.term(a,k))>10^(-20)){
    s=s+kth.term(a,k)
    k=k+1
  }
  return(s)
}
paste("When a=c(1,2) the result is",sum.term(c(1,2)))

## -----------------------------------------------------------------------------
f=function(k,a){
  s1.1=lgamma(k/2)+log(2)
  s1.2=log(sqrt(pi*(k-1)))
  s1.3=lgamma((k-1)/2)
  c1=sqrt(a^2*(k-1)/(k-a^2))
  f1=function(u){return((1+u^2/(k-1))^(-k/2))}
  i1=integrate(f1,lower=0,upper=c1)$value
  s2.1=lgamma((k+1)/2)+log(2)
  s2.2=log(sqrt(pi*k))
  s2.3=lgamma(k/2)
  c2=sqrt(a^2*k/(k+1-a^2))
  f2=function(u){return((1+u^2/k)^(-(k+1)/2))}
  i2=integrate(f2,lower=0,upper=c2)$value
  return(exp(s1.1-s1.2-s1.3)*i1-exp(s2.1-s2.2-s2.3)*i2)
}



## -----------------------------------------------------------------------------
K =4:25
root1=numeric(22)
root2=numeric(22)
root3=numeric(22)
for(i in 1:length(K)){
  root1[i]=uniroot(function(x) f(k=K[i],x) ,c(1,min(sqrt(K[i])-0.01,2)))$root
  root2[i]=-root1[i]
}
result=cbind(K,root1,root2,root3)
x_html=knitr::kable((result),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)


## -----------------------------------------------------------------------------
S = function(a, k){
1 - pt(sqrt((a^2 *k)/(k+1-a^2)), df = k)
}
f = function(a, k){
S(a, k-1)-S(a, k)
}
A = function(k){
uniroot(function(x) S(x,k-1)-S(x,k) , c(1,min(sqrt(k)-0.01,2)))$root
}
compare.11.4 = rep(0,length(K))
for(i in 1:length(K)){compare.11.4[i] =A(K[i])}
result2=cbind(compare.11.4,root1)
x_html=knitr::kable((result2),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "hover",full_width = TRUE,font_size=12)

## -----------------------------------------------------------------------------
trims <- c(0, 0.1, 0.2, 0.5)
x <- rcauchy(100)
lapply(trims, function(trim) mean(x, trim = trim))
lapply(trims, mean, x = x)

## -----------------------------------------------------------------------------
formulas = list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt 
)
rsq <- function(mod) summary(mod)$r.squared
result3 = lapply(formulas, function(f) lm(data = mtcars, formula = f))

lapply(result3,rsq)


## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})
result4=lapply(bootstraps,function(data)lm(data=data,formula=mpg ~ disp))
lapply(result4,rsq)

## -----------------------------------------------------------------------------
colsd=function(x){
  locus=1:length(x)
  vapply(locus,function(i)sd(x[[i]]),numeric(1))
}
colsd(mtcars)

## -----------------------------------------------------------------------------
colsdmix=function(x){
  locus1=1:length(x)
  locus2=vapply(locus1,function(i)is.numeric(x[[i]]),numeric(1))
  locus3=seq(along=locus2)[locus2==TRUE]
  vapply(locus3,function(i)sd(x[[i]]),numeric(1))
}
list2=list(list("a","b","c"),c(rnorm(10)),c(rbinom(10,5,0.5)))
colsdmix(list2)

## -----------------------------------------------------------------------------
library(parallel)
## repeat Page 204 Question 5
cl <- makeCluster(mc <- getOption("cl.cores", 4))
rsq <- function(mod) summary(mod)$r.squared
clusterExport(cl=cl, varlist=c("rsq"))
##redo Page 204 Question 3
parSapply(cl, bootstraps, function(x) rsq(lm(data = x, formula = mpg ~ disp)))
##Page 204 Question 4
parSapply(cl, formulas, function(f) rsq(lm(data =mtcars, formula = f)))

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)
set.seed(123)
cppFunction('NumericMatrix gibbsC(int N, int thin, int n, int a, int b) {
  NumericMatrix mat(N, 2);
  double x = 0, y = 0;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rbinom(1, n, y)[0];
      y = rbeta(1, (x+a), (n-x+b))[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}')
gibbsR <- function(N, thin,n,a,b) {
  mat <- matrix(nrow = N, ncol = 2)
  x <- y <- 0
  for (i in 1:N) {
    for (j in 1:thin) {
      x <- rbinom(1, n, y)
      y <- rbeta(1, (x+a), (n-x+b))
    }
    mat[i, ] <- c(x, y)
  }
  mat
}


gibR=gibbsR(1000,10,20,4,2)
gibC=gibbsC(1000,10,20,4,2)

qqplot(gibR[,1],gibC[,1])
qqplot(gibR[,2],gibC[,2])
ts <- microbenchmark(gibbR=gibbsR(1000,10,20,4,2),
gibbC=gibbsC(1000,10,20,4,2))
summary(ts)[,c(1,3,5,6)]

