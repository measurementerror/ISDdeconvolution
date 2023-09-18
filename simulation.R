library(nimble) ##double exponential distribution package
library(decon)

###Fan & Truong (1993) and Qin & Feng (2002)
###deconvolution density estimator with gausian kernel + double exp error
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

###simulation

beta_0<-2
size<-500
beta<-rep(0,size)
beta_naive <-rep(0,size)
beta_mle<-rep(0,size)
beta_bc<-rep(0,size)


n0<-100
n1<-100
s<-0.5
rt<-sqrt(s^2/2)

# rate0<-3
# rate1<-rate0-beta_0

size <- 500
set.seed(234)

## Normal covariate
for (k in 1:size){
  
  u0<-rdexp(n0,0,rt) #error 
  z0<-rnorm(n0,0,1) #true unobserved
  x0<-u0+z0 #observed
  u1<-rdexp(n1,0,rt) 
  z1<-rnorm(n1,beta_0,1)
  x1<-u1+z1 
  
  
  ##proposed estimator based on deconvolution
  h0<-bw.dboot2(x0,sig=s,error="laplacian")
  h1<-bw.dboot2(x1,sig=s,error="laplacian")
  
  
  a<-mean(x0)
  b<-mean(x1)
  
  seq<-seq(a,b,by=0.01)
  mid<-0.5*(seq[2:length(seq)]+seq[1:(length(seq)-1)])
  width<-seq[2:length(seq)]-seq[1:(length(seq)-1)]
  
  f<-rep(0,length(mid))
  d0<-rep(0,length(mid))
  d1<-rep(0,length(mid))
  
  p0<-rep(0,length(mid))
  p1<-rep(0,length(mid))
  
  for (i in 1:length(mid)){
    p0[i]<-dens(mid[i],x0,h0,s)
    p1[i]<-dens(mid[i],x1,h1,s)
    
    d0[i]<-ifelse(p0[i]>0,p0[i],1/(n0))
    d1[i]<-ifelse(p1[i]>0,p1[i],1/(n1))
    
    f[i]<-(log(d1[i])-log(d0[i]))*(mid[i]-(a+b)/2)
    
  }
  B<-sum(f*width)
  beta[k]<-B/((b-a)^3/12)
  
  #### naive ISD estimator
  t0<-band(x0)
  t1<-band(x1)
  d<-rep(0,length(mid))
  for (i in 1:length(mid)){
    d[i]<-(log(denOR(mid[i],x1,t1))-log(denOR(mid[i],x0,t0)))*(mid[i]-(a+b)/2)
  }
  A<-sum(d*width)
  beta_naive[k]<-A/((b-a)^3/12)
  
  ###MLE
  y0<-rep(0,n0)
  y1<-rep(1,n1)
  
  xb<-c(x0,x1)
  yb<-c(y0,y1)
  
  beta_mle[k]<-as.numeric(glm(yb~xb,family = "binomial")$coefficients[2])
  
  ##Stefanski-Carroll bias-corrected estimator
  s.n <-function(x){
    v<-rep(0,length(xb))
    for (i in 1:length(xb)){
      v[i]<-firstde(x*xb[i])*(xb[i]^2)
    }
    (1/(length(xb)))*sum(v)
  }
  j.n1<-(-1/(2*length(xb)))*sum(secde(betahat[k]*xb)*xb*betahat[k])
  j.n2<-(-1/length(xb))*sum(firstde(betahat[k]*xb))
  B.n<-1/s.n(betahat[k])*(j.n1+j.n2)
  beta_bc[k]<-(1-(s^2)*B.n)*betahat[k]
  
  
}

result <- cbind(beta,beta_naive,beta_mle,beta_bc)