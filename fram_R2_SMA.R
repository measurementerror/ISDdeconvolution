#########################################
##framingham heart study 
##framingham data from deconvolve package
#########################################

##11/6/23
library(deconvolve)
library(ggplot2)
library(decon)

kn<-function(x,h,sigma){
  (1/sqrt(2*pi))*exp((-1/2)*x^2)*(1-((sigma^2)/(2*h^2))*(x^2-1))
}

dens<-function(x,data,h,sigma){
  w<-0
  for (i in 1:length(data)){
    w<-w+(1/(length(data)*h))*kn((data[i]-x)/h,h,sigma)
  }
  return (w)
}

####naive ISD estimator with gaussian kernel
kernel<-function(u){
  1/(sqrt(2*pi))*exp(-0.5*u^2)
}

denOR<-function(x,data,h){
  w<-0
  for (i in 1:length(data)){
    w<-w+(1/(length(data)*h))*kernel((data[i]-x)/h)
  }
  return (w)
}

band<-function(data){
  length(data)^(-1/3)
}

###Stefanski Carroll bias-corrected estimator
firstde<-function(x){
  exp(-x)/(exp(-x)+1)^2
}
secde<-function(x){
  (-(exp(x)-1)*exp(x))/(exp(x)+1)^3
}

## data framingham from deconvolve package
View(framingham)
head(framingham)
attach(framingham)

## nested case control sample
Fram.exit <- rep(1,nrow(framingham))
fram1 <- data.frame(framingham,Fram.exit)

set.seed(243)
ncc <- ccwc(entry=0,exit=Fram.exit,fail=FIRSTCHD,controls=5,data=fram1,include=list(CHOLEST3,SBP21,SBP22,FIRSTCHD),
            match=list(AGE,SMOKE))


## error density check
err <- ncc$SBP21  - (ncc$SBP21 + ncc$SBP22)/2

qqplot.data <- function (vec) # argument: vector of numbers
{
  # following four lines from base R's qqline()
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  
  d <- data.frame(resids = vec)
  
  ggplot(d, aes(sample = resids)) + stat_qq() + geom_abline(slope = slope, intercept = int) + 
    labs(title="Normal Q-Q plot of the measurement error")
  
}

qqplot.data(err)
shapiro.test(err)


## parameter estimation 

con0 <- subset(ncc,FIRSTCHD == 0)
cas0 <- subset(ncc,FIRSTCHD == 1)


n0 <- nrow(con0)
n1 <- nrow(cas0)

## Based on measurement 1
z0 <- con0$SBP21
z1 <- cas0$SBP21

a <- max(quantile(z0,.1),quantile(z1,.1))
b <- min(quantile(z0,.9),quantile(z1,.9))

#a <- quantile(c(z0,z1),0.2)
#b <- quantile(c(z0,z1),0.8)
#a <- median(z0)
#b <- median(z1)



s <- sd(err)
##proposed estimator based on deconvolution
h0<-bw.dboot2(z0,sig=s,error="laplacian")
h1<-bw.dboot2(z1,sig=s,error="laplacian")

seq<-seq(a,b,by=0.01)
mid<-0.5*(seq[2:length(seq)]+seq[1:(length(seq)-1)])
width<-seq[2:length(seq)]-seq[1:(length(seq)-1)]

f<-rep(0,length(mid))
d0<-rep(0,length(mid))
d1<-rep(0,length(mid))

p0<-rep(0,length(mid))
p1<-rep(0,length(mid))

for (i in 1:length(mid)){
  p0[i]<-dens(mid[i],z0,h0,s)
  p1[i]<-dens(mid[i],z1,h1,s)
  
  d0[i]<-ifelse(p0[i]>0,p0[i],1/(n0))
  d1[i]<-ifelse(p1[i]>0,p1[i],1/(n1))
  
  f[i]<-(log(d1[i])-log(d0[i]))*(mid[i]-(a+b)/2)
  
}
B<-sum(f*width)
beta_bp <-B/((b-a)^3/12)

#### naive ISD estimator
t0<-sd(z0)*band(z0)
t1<-sd(z1)*band(z1)
d<-rep(0,length(mid))
for (i in 1:length(mid)){
  d[i]<-(log(denOR(mid[i],z1,t1))-log(denOR(mid[i],z0,t0)))*(mid[i]-(a+b)/2)
}
A<-sum(d*width)
beta_bp_naive <-A/((b-a)^3/12)


###MLE
y0<-rep(0,n0)
y1<-rep(1,n1)

zb<-c(z0,z1)
yb<-c(y0,y1)

beta_bp_mle <-as.numeric(glm(yb~zb,family = "binomial")$coefficients[2])

##Stefanski-Carroll bias-corrected estimator
s.n <-function(x){
  v<-rep(0,length(zb))
  for (i in 1:length(zb)){
    v[i]<-firstde(x*zb[i])*(zb[i]^2)
  }
  (1/(length(zb)))*sum(v)
}
j.n1<-(-1/(2*length(zb)))*sum(secde(beta_bp_mle*zb)*zb*beta_bp_mle)
j.n2<-(-1/length(zb))*sum(firstde(beta_bp_mle*zb))
B.n<-1/s.n(beta_bp_mle)*(j.n1+j.n2)
beta_bp_bc<-(1-(s^2)*B.n)*beta_bp_mle


c(beta_bp,beta_bp_naive,beta_bp_mle,beta_bp_bc)

## Based on measurement 2
z0 <- con0$SBP22
z1 <- cas0$SBP22

a <- max(quantile(z0,.1),quantile(z1,.1))
b <- min(quantile(z0,.9),quantile(z1,.9))

#a <- quantile(c(z0,z1),0.2)
#b <- quantile(c(z0,z1),0.8)
#a <- median(z0)
#b <- median(z1)



s <- sd(err)
##proposed estimator based on deconvolution
h0<-bw.dboot2(z0,sig=s,error="laplacian")
h1<-bw.dboot2(z1,sig=s,error="laplacian")

seq<-seq(a,b,by=0.01)
mid<-0.5*(seq[2:length(seq)]+seq[1:(length(seq)-1)])
width<-seq[2:length(seq)]-seq[1:(length(seq)-1)]

f<-rep(0,length(mid))
d0<-rep(0,length(mid))
d1<-rep(0,length(mid))

p0<-rep(0,length(mid))
p1<-rep(0,length(mid))

for (i in 1:length(mid)){
  p0[i]<-dens(mid[i],z0,h0,s)
  p1[i]<-dens(mid[i],z1,h1,s)
  
  d0[i]<-ifelse(p0[i]>0,p0[i],1/(n0))
  d1[i]<-ifelse(p1[i]>0,p1[i],1/(n1))
  
  f[i]<-(log(d1[i])-log(d0[i]))*(mid[i]-(a+b)/2)
  
}
B<-sum(f*width)
beta_bp <-B/((b-a)^3/12)

#### naive ISD estimator
t0<-sd(z0)*band(z0)
t1<-sd(z1)*band(z1)
d<-rep(0,length(mid))
for (i in 1:length(mid)){
  d[i]<-(log(denOR(mid[i],z1,t1))-log(denOR(mid[i],z0,t0)))*(mid[i]-(a+b)/2)
}
A<-sum(d*width)
beta_bp_naive <-A/((b-a)^3/12)


###MLE
y0<-rep(0,n0)
y1<-rep(1,n1)

zb<-c(z0,z1)
yb<-c(y0,y1)

beta_bp_mle <-as.numeric(glm(yb~zb,family = "binomial")$coefficients[2])

##Stefanski-Carroll bias-corrected estimator
s.n <-function(x){
  v<-rep(0,length(zb))
  for (i in 1:length(zb)){
    v[i]<-firstde(x*zb[i])*(zb[i]^2)
  }
  (1/(length(zb)))*sum(v)
}
j.n1<-(-1/(2*length(zb)))*sum(secde(beta_bp_mle*zb)*zb*beta_bp_mle)
j.n2<-(-1/length(zb))*sum(firstde(beta_bp_mle*zb))
B.n<-1/s.n(beta_bp_mle)*(j.n1+j.n2)
beta_bp_bc<-(1-(s^2)*B.n)*beta_bp_mle


c(beta_bp,beta_bp_naive,beta_bp_mle,beta_bp_bc)

