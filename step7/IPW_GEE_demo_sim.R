##############################
# demo of IPW-GEE approaches
##############################
set.seed(1234)
library(geepack)
library(CRTgeeDR)
library(lme4)

# calculate the true effect
true<-function(n,m.bar=50){
  # cluster size
  m<-round(runif(n,m.bar-10,m.bar+10))
  N<-sum(m)
  cluster<-rep(1:n,m)
  
  # outcome generation
  t<-rep(0,n)
  t[sample(1:n,n/2)]<-1
  t.rep<-rep(t,m)
  alpha<-rep(rnorm(n,0,0.18),m)
  x<-alpha+rnorm(N,0,0.5)
  delta<-rep(rnorm(n,0,0.2),m)
  p.out<-plogis(1.36*t.rep+x+delta)
  y<-rbinom(N,1,p.out)
  
  # analysis of complete data
  model<-geeglm(y~t.rep,family=binomial,id=cluster,corstr = "exchangeable")
  return(summary(model)$coef[2,1])
}
# true value
delta0<-true(n=1000,m.bar=50)

# simulation parameters
nsim<-200
n<-50
m.bar<-50

# matrices to save output
EST<-SE<-matrix(NA,nsim,4)
colnames(EST)<-colnames(SE)<-
  c("IPW1","IPW2","IPW1_CLU","IPW2_CLU")

# main function for simulations
for(i in 1:nsim){
  # cluster size
  m<-round(runif(n,m.bar-10,m.bar+10))
  N<-sum(m)
  cluster<-rep(1:n,m)
  
  # outcome generation
  t<-rep(0,n)
  t[sample(1:n,n/2)]<-1
  t.rep<-rep(t,m)
  alpha<-rep(rnorm(n,0,0.18),m)
  x<-alpha+rnorm(N,0,0.5)
  delta<-rep(rnorm(n,0,0.2),m)
  p.out<-plogis(1.36*t.rep+x+delta)
  y<-rbinom(N,1,p.out)
  
  # missingness indicator
  theta<-rep(rnorm(n,0,sqrt(4)),m)
  p.mis<-plogis(-1.8+t.rep+x+theta)
  p.obs<-1-p.mis
  R<-rbinom(N,1,p.obs)
  
  # summary of missingness
  # mean(p.mis)
  # tapply(R,cluster,sum)/m
  
  # estimate IP weights
  mism<-glm(R~t.rep+x,family=binomial(link="logit"))
  mism.cl<-glmer(R~t.rep+x+(1|cluster),family=binomial)
  w<-1/predict(mism,type="response")
  w.cl<-1/predict(mism.cl,type="response")
  
  # are the weights correctly estimated?
  # plot(1/w[R==1]~p.obs[R==1])
  # plot(1/w.cl[R==1]~p.obs[R==1])
  
  # observed data rows
  y.obs<-y[R==1]
  t.rep.obs<-t.rep[R==1]
  x.obs<-x[R==1]
  w.obs<-w[R==1]
  w.cl.obs<-w.cl[R==1]
  cluster.obs<-cluster[R==1]
  data.obs<-data.frame(y.obs,t.rep.obs,x.obs,w.obs,w.cl.obs,cluster.obs)
  
  # geepack analysis and geeCRTdr analysis
  m1<-geeglm(y.obs~t.rep.obs,family=binomial,id=cluster.obs,corstr = "exchangeable",weights=w.obs)
  m2<-geeDREstimation(y.obs~t.rep.obs,id="cluster.obs",data=data.obs,nameTRT = "t.rep.obs",
                        family=binomial,corstr = "exchangeable",weights=data.obs$w.obs)
  m.cl1<-geeglm(y.obs~t.rep.obs,family=binomial,id=cluster.obs,corstr = "exchangeable",weights=w.cl.obs)
  m.cl2<-geeDREstimation(y.obs~t.rep.obs,id="cluster.obs",data=data.obs,nameTRT = "t.rep.obs",
                        family=binomial,corstr = "exchangeable",weights=data.obs$w.cl.obs)
  ipw1<-summary(m1)
  ipw2<-summary(m2)
  ipw.cl1<-summary(m.cl1)
  ipw.cl2<-summary(m.cl2)
  
  # save results
  EST[i,1]<-ipw1$coef[2,1]
  EST[i,2]<-ipw2$beta[2]
  EST[i,3]<-ipw.cl1$coef[2,1]
  EST[i,4]<-ipw.cl2$beta[2]
  
  SE[i,1]<-ipw1$coef[2,2]
  SE[i,2]<-ipw2$se.robust[2]
  SE[i,3]<-ipw.cl1$coef[2,2]
  SE[i,4]<-ipw.cl2$se.robust[2]
  
  # loop control
  if(i%% 10==0)print(i)
}

# bias
colMeans(EST)-delta0
#        IPW1        IPW2    IPW1_CLU    IPW2_CLU 
# 0.005873991 0.003164442 0.024939635 0.023911156 

# comparison of sd
apply(EST,2,sd)
#      IPW1      IPW2  IPW1_CLU  IPW2_CLU 
# 0.1381772 0.1378128 0.1857456 0.1834562 

colMeans(SE)
#      IPW1      IPW2  IPW1_CLU  IPW2_CLU 
# 0.1293678 0.1290205 0.1697802 0.1678052 

# empirical coverage rates
colMeans(EST-1.96*SE<=delta0 & EST+1.96*SE>=delta0)
#     IPW1     IPW2 IPW1_CLU IPW2_CLU 
# 0.930    0.935    0.915    0.925 
