missing_icc=function(rho){
  pi=3.1415926
  sigma=(pi^2/3)/(1/rho-1)
  return(sigma)
}

missingicc=missing_icc(0.05)
# 0.004

library('lme4')
library(geepack)
library(CRTgeeDR)

myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

# expit function
expit=function(x){y=exp(x)/(1+exp(x));return(y)}

#

one_group3=function(mis,s,i,k,mm,seed=123,psi){
  set.seed(seed)
  if(s==1){b0=1;b1=1.36;b2=1}
  if(s==3){b0=1;b1=1.36;b2=0.588}
  ## parameters
  sigma_b=sqrt(0.2)
  mu_x=0
  sigma_alpha=sqrt(0.18)
  sigma_u=sqrt(3.37)
  
  phi=1
  x=matrix(0,k,mm+5*(mm/25))
  y=matrix(0,k,mm+5*(mm/25))
  pi=matrix(0,k,mm+5*(mm/25))
  cluster=matrix(0,k,mm+5*(mm/25))
  r=matrix(0,k,mm+5*(mm/25))
  mis_clu=matrix(0,k,mm+5*(mm/25))
  groups=matrix(i,k,mm+5*(mm/25))
  for(kk in 1:k){
    delta=rnorm(1,0,sigma_b)
    alpha=rnorm(1,mu_x,sigma_alpha)
    m=round(runif(1,mm-5*(mm/25),mm+5*(mm/25)))
    # print(m)
    cluster[kk,]=kk
    sigma_icc=missing_icc(0.05)
    clu=rnorm(1,0,sqrt(sigma_icc))
    # print(m)
    mis_clu[kk,]=clu
    for(j in 1:m){
      u=rnorm(1,0,sigma_u)
      x[kk,j]=u+alpha
      pi[kk,j]=expit(b0+b1*i+b2*x[kk,j]+delta)
      y[kk,j]=rbinom(1,1,pi[kk,j])
    }
  }
  #r=expit(psi+phi*x)
  # Different missingness generation models, 6 scenarios
  if(mis==1){r=expit(psi+phi*x)}
  if(mis==2){r=expit(psi+phi*x+groups)}
  if(mis==3){r=expit(psi+phi*x+groups+x*groups)}
  if(mis==4){r=expit(psi+phi*x+mis_clu)}
  if(mis==5){r=expit(psi+phi*x+groups+mis_clu)}
  if(mis==6){r=expit(psi+phi*x+groups+x*groups+mis_clu)}
  return(list(x=x,y=y,pi=pi,r=r,cluster=cluster))
}
#psi=-1.8 missing =30%
# combine two intervention arm
data_gene3=function(k,m,s,mis,seed=123,psi){
  set.seed(seed)
  # step1:
  a1=one_group3(i=1,k=k,m=m,seed=seed,mis=mis,s=s,psi=psi)
  a0=one_group3(i=0,k=k,m=m,seed=seed,mis=mis,s=s,psi=psi)
  # step2
  b1=data.frame(x=c(a1$x),y=c(a1$y),r=c(a1$r),cluster=c(a1$cluster))
  R=c()
  for(i in 1:dim(b1)[1]){
    R[i]=rbinom(1,1,b1$r[i])
  }
  b1$R=R
  b1$arm=1
  b1$cluster=b1$cluster+k
  b0=data.frame(x=c(a0$x),y=c(a0$y),r=c(a0$r),cluster=c(a0$cluster))
  R=c()
  for(i in 1:dim(b0)[1]){
    R[i]=rbinom(1,1,b0$r[i])
  }
  b0$R=R
  b0$arm=0
  # step3:
  data=rbind(b0,b1)
  data0=data[data$x!=0,]
  return(data0)
}



data_gen4=function(n,m.bar,icc,intercept,seed=123){
  set.seed(seed)
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
  missingicc=missing_icc(icc)
  theta<-rep(rnorm(n,0,sqrt(missingicc)),m)
  p.mis<-plogis(intercept+t.rep+x+theta)
  p.obs<-1-p.mis
  R<-rbinom(N,1,p.obs)
  temp=data.frame(y=y,x=x,arm=t.rep,delta=delta,R=R,cluster=cluster)
  return(temp)
}

mis_per=function(n,m.bar,icc,seed=123){
  data=c()
  for(i in seq(-3,3,0.1)){
    temp=data_gen4(n,m.bar,icc,intercept=i,seed)
    percent=sum(temp$R)/dim(temp)[1]
    data=rbind(data,c(percent,i))
  }
  data[,1]=abs(data[,1]-0.7)
  intercept=data[data[,1]==min(data[,1])][2]
  return(intercept)
}

mis_per_check=function(n,m,icc,intercept,seed=123){
  temp=data_gen4(n,m.bar,icc,intercept,seed)
  percent=sum(temp$R)/dim(temp)[1]
  print(percent)
}



for(icc in seq(0.1,0.9,0.1)){
  print('icc')
  print(icc)
  n=50
  m.bar=50
  
  res_est_ind=c();res_std_ind=c();res_warn_ind=c()
  res_est_ex=c();res_std_ex=c();res_warn_ex=c()
  res_clu_est_ind=c();res_clu_std_ind=c();res_clu_warn_ind=c()
  res_clu_est_ex=c();res_clu_std_ex=c();res_clu_warn_ex=c()
  ICCs=c()
  
  true_est_ex=c();true_std_ex=c();true_warn_ex=c()
  true_est_ind=c();true_std_ind=c();true_warn_ind=c()
  
  est=c();std=c();warn=c()
  
  for(time in 1:2){
    print(time)
    intercept=mis_per(n,m.bar,icc,seed=time)
    #mis_per_check(n,m.bar,icc,seed=time,intercept=intercept)
    
    d1=data_gen4(n,m.bar,icc,intercept,seed=time)
    d2=d1
    d2$y=ifelse(d2$R==0,NA,d2$y)
    d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=1-d2$R)
    
    logs=glm(missing ~ x+arm, data = d3,
             family = binomial(link='logit'))
    logs2=glmer(missing ~ x+arm+(1|cluster) , data = d3,
                family = binomial(link='logit'))
    
    weight=1/expit(predict(logs))
    weight2=1/expit(predict(logs2))
    
    d3$weight=weight
    d3$weight2=weight2
    
    ICCs=c(ICCs,icc)
    trues_ind=myTryCatch(geeDREstimation(formula=y~x+arm,
                                       id="cluster" , data = d1,nameY='y',
                                       nameTRT='arm',
                                       family =  binomial("logit"),
                                       corstr = "independence"))
    trues_ex=myTryCatch(geeDREstimation(formula=y~x+arm,
                                         id="cluster" , data = d1,nameY='y',
                                        nameTRT='arm',
                                         family =  binomial("logit"),
                                         corstr = "exchangeable"))
    
    if(is.null(trues_ind$value)==0){
      trues_ind_est=summary(trues_ind$value)$beta[3]
      trues_ind_std=summary(trues_ind$value)$se.robust[3]
      true_est_ind=c(true_est_ind,trues_ind_est)
      true_std_ind=c(true_std_ind,trues_ind_std)
      if(trues_ind$value$converged==0){true_warn_ind=c(true_warn_ind,time)}
      if(trues_ind$value$converged==1){true_warn_ind=c(true_warn_ind,0)}
    }
    if(is.null(trues_ex$value)==0){
      trues_ex_est=summary(trues_ex$value)$beta[3]
      trues_ex_std=summary(trues_ex$value)$se.robust[3]
      true_est_ex=c(true_est_ex,trues_ex_est)
      true_std_ex=c(true_std_ex,trues_ex_std)
      if(trues_ex$value$converged==0){true_warn_ex=c(true_warn_ex,time)}
      if(trues_ex$value$converged==1){true_warn_ex=c(true_warn_ex,0)}
    }
    
    if(is.null(trues_ind$value)==1){
      true_est_ind=c(true_est_ind,NA)
      true_std_ind=c(true_std_ind,NA)
      true_warn_ind=c(true_warn_ind,time)
    }
    if(is.null(trues_ex$value)==1){
      true_est_ex=c(true_est_ex,NA)
      true_std_ex=c(true_std_ex,NA)
      true_warn_ex=c(true_warn_ex,time)
    }
    
    ipw_ind=myTryCatch(geeDREstimation(formula=y~x+arm,
                               id="cluster" , data = d3,
                               nameMISS='missing',nameY='y',
                               nameTRT='arm',
                               weights = d3$weight,
                               family =  binomial("logit"),
                               corstr = "independence"))
    ipw_clu_ind=myTryCatch(geeDREstimation(formula=y~x+arm,
                               id="cluster" , data = d3,
                               nameMISS='missing',nameY='y',
                               nameTRT='arm',
                               weights = d3$weight2,
                               family =  binomial("logit"),
                               corstr = "independence"))
    ipw_ex=myTryCatch(geeDREstimation(formula=y~x+arm,
                               id="cluster" , data = d3,
                               nameMISS='missing',nameY='y',
                               nameTRT='arm',
                               weights = d3$weight,
                               family =  binomial("logit"),
                               corstr = "exchangeable"))
    ipw_clu_ex=myTryCatch(geeDREstimation(formula=y~x+arm,
                               id="cluster" , data = d3,
                               nameMISS='missing',nameY='y',
                               nameTRT='arm',
                               weights = d3$weight2,
                               family =  binomial("logit"),
                               corstr = "exchangeable"))
   
    if(is.null(ipw_ind$value)==0){
      ipw_ind_est=summary(ipw_ind$value)$beta[3]
      ipw_ind_std=summary(ipw_ind$value)$se.robust[3]
      res_est_ind=c(res_est_ind,ipw_ind_est)
      res_std_ind=c(res_std_ind,ipw_ind_std)
      if(ipw_ind$value$converged==0){res_warn_ind=c(res_warn_ind,time)}
      if(ipw_ind$value$converged==1){res_warn_ind=c(res_warn_ind,0)}
    }
    if(is.null(ipw_ex$value)==0){
      ipw_ex_est=summary(ipw_ex$value)$beta[3]
      ipw_ex_std=summary(ipw_ex$value)$se.robust[3]
      res_est_ex=c(res_est_ex,ipw_ex_est)
      res_std_ex=c(res_std_ex,ipw_ex_std)
      if(ipw_ex$value$converged==0){res_warn_ex=c(res_warn_ex,time)}
      if(ipw_ex$value$converged==1){res_warn_ex=c(res_warn_ex,0)}
    }
    if(is.null(ipw_clu_ind$value)==0){
      ipw_clu_ind_est=summary(ipw_clu_ind$value)$beta[3]
      ipw_clu_ind_std=summary(ipw_clu_ind$value)$se.robust[3]
      res_clu_est_ind=c(res_clu_est_ind,ipw_clu_ind_est)
      res_clu_std_ind=c(res_clu_std_ind,ipw_clu_ind_std)
      if(ipw_clu_ind$value$converged==0){res_clu_warn_ind=c(res_clu_warn_ind,time)}
      if(ipw_clu_ind$value$converged==1){res_clu_warn_ind=c(res_clu_warn_ind,0)}
    }
    if(is.null(ipw_clu_ex$value)==0){
      ipw_clu_ex_est=summary(ipw_clu_ex$value)$beta[3]
      ipw_clu_ex_std=summary(ipw_clu_ex$value)$se.robust[3]
      res_clu_est_ex=c(res_clu_est_ex,ipw_clu_ex_est)
      res_clu_std_ex=c(res_clu_std_ex,ipw_clu_ex_std)
      if(ipw_clu_ex$value$converged==0){res_clu_warn_ex=c(res_clu_warn_ex,time)}
      if(ipw_clu_ex$value$converged==1){res_clu_warn_ex=c(res_clu_warn_ex,0)}
    }

    if(is.null(ipw_ind$value)==1){
      res_est_ind=c(res_est_ind,NA)
      res_std_ind=c(res_std_ind,NA)
      res_warn_ind=c(res_warn_ind,time)
    }
    if(is.null(ipw_ex$value)==1){
      res_est_ex=c(res_est_ex,NA)
      res_std_ex=c(res_std_ex,NA)
      res_warn_ex=c(res_warn_ex,time)
    }
    if(is.null(ipw_clu_ind$value)==1){
      res_clu_est_ind=c(res_clu_est_ind,NA)
      res_clu_std_ind=c(res_clu_std_ind,NA)
      res_clu_warn_ind=c(res_clu_warn_ind,time)
    }
    if(is.null(ipw_clu_ex$value)==1){
      res_clu_est_ex=c(res_clu_est_ex,NA)
      res_clu_std_ex=c(res_clu_std_ex,NA)
      res_clu_warn_ex=c(res_clu_warn_ex,time)
    }
  }
  
  EST=data.frame(true_est_ind,true_est_ex,res_est_ind,res_est_ex,res_clu_est_ind,res_clu_est_ex,icc,time)
  STD=data.frame(true_std_ind,true_std_ex,res_std_ind,res_std_ex,res_clu_std_ind,res_clu_std_ex,icc,time)
  WARN=data.frame(true_warn_ind,true_warn_ex,res_warn_ind,res_warn_ex,res_clu_warn_ind,res_clu_warn_ex,icc,time)

  names(EST)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')
  names(STD)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')
  names(WARN)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')
  
  result=list(EST=EST,STD=STD,WARN=WARN)
  name_res=paste('ipw_result',icc,'.RData',sep='')
  save(result,file=name_res)
  
}


