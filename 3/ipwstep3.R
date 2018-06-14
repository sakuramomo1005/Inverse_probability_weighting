### renew:
## ICC effect 0.05,0.1,0.15,0.2,0.5
## Missing ICC: 0.1,0.15,0.2
## add standard MI
## using exchangable to analyze
## add random effect/cluster in missingness 
## add interactions 
## x, x+arm, x*arm
## truncating base on missing pendercentage,  Not weight itself
## do not round the weigth

library(jomo)
library(geepack)
library(lme4)
library(jomo)
library(CRTgeeDR)

#scenarios: 
# s=1: 
# mis=1 x
# mis=2 x+arm
# mis=3 x+arm+x:arm

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
one_group2=function(i,k,mm,seed=123,icc=0.05,s=1,mis=1){
  set.seed(seed)
  ## parameters
  if(s==1){b2=1;psi=-1.34}
  if(s==2){b2=0.588;psi=0.65}
  b0=1;b1=1.36
  sigma_b=sqrt(0.2)
  mu_x=0
  sigma_x=3.55
  sigma_alpha=sigma_x*icc
  sigma_u=sigma_x-sigma_alpha
  phi=1
  if(mm==25){
    x=matrix(0,k,mm+5)
    y=matrix(0,k,mm+5)
    pi=matrix(0,k,mm+5)
    cluster=matrix(0,k,mm+5)
    r=matrix(0,k,mm+5)
    mis_clu=matrix(0,k,mm+5)
    groups=matrix(i,k,mm+5)
  }
  if(mm==50){
    x=matrix(0,k,mm+10)
    y=matrix(0,k,mm+10)
    pi=matrix(0,k,mm+10)
    cluster=matrix(0,k,mm+10)
    r=matrix(0,k,mm+10)
    mis_clu=matrix(0,k,mm+10)
    groups=matrix(i,k,mm+10)
  }
  for(k in 1:k){
    delta=rnorm(1,0,sigma_b)
    alpha=rnorm(1,mu_x,sqrt(sigma_alpha))
    if(mm==25){m=round(runif(1,mm-5,mm+5))}
    if(mm==50){m=round(runif(1,mm-10,mm+10))}
    sigma_icc=0.4
    clu=rnorm(1,0,sqrt(sigma_icc))
    # print(m)
    cluster[k,]=k
    mis_clu[k,]=clu
    for(j in 1:m){
      u=rnorm(1,0,sqrt(sigma_u))
      x[k,j]=u+alpha
      pi[k,j]=expit(b0+b1*i+b2*x[k,j]+delta)
      y[k,j]=rbinom(1,1,pi[k,j])
    }
  }

  if(mis==1){r=expit(psi+phi*x)}
  if(mis==2){r=expit(psi+phi*x+groups)}
  if(mis==3){r=expit(psi+phi*x+groups+x*groups)}
  if(mis==4){r=expit(psi+phi*x+mis_clu)}
  if(mis==5){r=expit(psi+phi*x+groups+mis_clu)}
  if(mis==6){r=expit(psi+phi*x+groups+x*groups+mis_clu)}
  return(list(x=x,y=y,pi=pi,r=r,cluster=cluster))
}
data_gene2=function(k,m,seed=123,s=1,icc=0.05,mis=1){
  set.seed(seed)
  # step1:
  a1=one_group2(1,k,m,seed,icc,s,mis)
  a0=one_group2(0,k,m,seed,icc,s,mis)
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
mypool=function(mean0,sd0,num=5,print='no'){
  m=mean(mean0)
  v=mean(sd0)
  B=sd(mean0)
  v_hat=v+(1+1/num)*B
  l=m-1.96*v_hat
  u=m+1.96*v_hat
  if(print=='no'){
    return(list(mean=m,std=v_hat))
  }
  if(print=='yes'){
    print('mean (95% CI)')
    print(paste(round(m,2)," (",round(l,2),',',round(u,2),')',sep=''))
    return(list(mean=m,std=v_hat))
  }
}

miss_per=function(a){
  b=sum(a$R)/dim(a)[1]
  return(b)
}

missing_icc=function(icc){
  pi=3.142
  a=pi^2/(3*(1/icc-1))
  return(a)
}

ipw_gee_ana_clu=function(times,s=1,icc=0.05,mis=1){
  
## empty effects
for(emp in 1){
  true_est=c();true_warning=c()
  true_est2=c();true_warning2=c()
  
  uncra_est=c();uncra_std=c();uncra_warning=c()
  uncra_est2=c();uncra_std2=c();uncra_warning2=c()
  
  cra_est=c();cra_std=c();cra_warning=c()
  cra_est2=c();cra_std2=c();cra_warning2=c()
  
  ipwgeeno_est=c();ipwgeeno_std=c();ipwgeeno_warn=c()
  ipwgeeno_est2=c();ipwgeeno_std2=c();ipwgeeno_warn2=c()
  
  ipwcrtgeedr_est=c();ipwcrtgeedr_std=c();ipwcrtgee_warn=c()
  ipwcrtgeedr_est2=c();ipwcrtgeedr_std2=c();ipwcrtgee_warn2=c()
  
  ipwgeeyes_est=c();ipwgeeyes_std=c();ipwgeeyes_warn=c()
  ipwgeeyes_est2=c();ipwgeeyes_std2=c();ipwgeeyes_warn2=c()
  
  ipwcrtgeedr_yes_est=c();ipwcrtgeedr_yes_std=c();ipwcrtgee_yes_warn=c()
  ipwcrtgeedr_yes_est2=c();ipwcrtgeedr_yes_std2=c();ipwcrtgee_yes_warn2=c()
  
  data_tt=c()
  data_km=c()
  
  
  missing=c()
}

for(km in 1){
  if(km==1){k=25;m=25}
  if(km==2){k=50;m=25}
  if(km==3){k=25;m=50}
  if(km==4){k=50;m=50}
  print("***************")
  print("km")
  print(km)
  print('times')
  for(tt in 1:times){
    
    data_tt=c(data_tt,tt)
    data_km=c(data_km,km)
    missing=c(missing,miss_per(d1))
    
    print(tt)
    d1=data_gene2(k=k,m=m,seed=tt,s=s,icc=icc,mis=mis)
    d2=d1
    d2$y=ifelse(d2$R==1,NA,d2$y)
    d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=d2$R)
    
    ## weight without cluster effects 
    l1=glm(missing ~ x , data = d3, 
             family = binomial(link='logit'))
    warn1=myTryCatch(l1)$warning
    
    l2=glmer(missing ~ x+(1|cluster) , data = d3, 
             family = binomial(link='logit'))
    warn2=myTryCatch(l2)$warning
    
    if(is.null(warn1)*is.null(warn2)==1){
      true=myTryCatch(geese(formu,data=d1,id=cluster,
                            family = binomial(link='logit'),
                            corstr = 'independence'))
      true_est=c(true_est,summary(true$value)$mean['arm','estimate'])
      if(is.null(true$warning)==0){
        true_warning=c(true_warning,tt)}
      
      true2=myTryCatch(geese(formu,data=d1,id=cluster,
                             family = binomial(link='logit'),
                             corstr = 'exchangeable'))
      true_est2=c(true_est2,summary(true2$value)$mean['arm','estimate'])
      if(is.null(true2$warning)==0){
        true_warning2=c(true_warning2,tt)}
      
      ### unadjusted CRA
      d_nona=na.omit(d1)
      formu0=formula(y~arm)
      uncra=myTryCatch(geese(formu0,data=d_nona,id=cluster,
                             family = binomial(link='logit'),
                             corstr = 'independence'))
      uncra_est=c(uncra_est,summary(uncra$value)$mean['arm','estimate'])
      uncra_std=c(uncra_std,summary(uncra$value)$mean['arm','san.se'])
      if(is.null(uncra$warning)==0){
        uncra_warning=c(uncra_warning,tt)}
      
      uncra2=myTryCatch(geese(formu0,data=d_nona,id=cluster,
                              family = binomial(link='logit'),
                              corstr = 'exchangeable'))
      uncra_est2=c(uncra_est2,summary(uncra2$value)$mean['arm','estimate'])
      uncra_std2=c(uncra_std2,summary(uncra2$value)$mean['arm','san.se'])
      if(is.null(uncra2$warning)==0){
        uncra_warning2=c(uncra_warning2,tt)}
      
      ### adjusted CRA
      cra=myTryCatch(geese(formu,data=d_nona,id=cluster,
                           family = binomial(link='logit'),
                           corstr = 'independence'))
      cra_est=c(cra_est,summary(cra$value)$mean['arm','estimate'])
      cra_std=c(cra_std,summary(cra$value)$mean['arm','san.se'])
      if(is.null(cra$warning)==0){
        cra_warning=c(cra_warning,tt)}
      
      cra2=myTryCatch(geese(formu,data=d_nona,id=cluster,
                            family = binomial(link='logit'),
                            corstr = 'exchangeable'))
      cra_est2=c(cra_est2,summary(cra2$value)$mean['arm','estimate'])
      cra_std2=c(cra_std2,summary(cra2$value)$mean['arm','san.se'])
      if(is.null(cra2$warning)==0){
        cra_warning2=c(cra_warning2,tt)}
      
      weight=expit(predict(l1))
      weight=ifelse(weight<0.05,0.05,weight)
      d3$weight_no=1/weight
      d3$weight_no_round=round(1/weight)
      
      weight2=expit(predict(l2))
      weight2=ifelse(weight2<0.05,0.05,weight2)
      d3$weight_clu=1/weight2
      d3$weight_clu_round=round(1/weight2)
      d4=na.omit(d3)
      
      #ipw geese no cluster effects
      formu=formula(y~x+arm)
      ipwgeeno0=geese(formu,data=d4,id=cluster,
                      family = binomial(link='logit'),
                      weights = weight_no,
                      corstr = 'independence')
      tempgeeno_est=summary(ipwgeeno0)$mean['arm','estimate']
      tempgeeno_std=summary(ipwgeeno0)$mean['arm','san.se']
      ipwgeeno_est=c(ipwgeeno_est,tempgeeno_est)
      ipwgeeno_std=c(ipwgeeno_std,tempgeeno_std)
      
      ipwgeeno=myTryCatch(geese(formu,data=d4,id=cluster,
                                family = binomial(link='logit'),
                                weights = weight_no_round,
                                corstr = 'independence'))
      if(is.null(ipwgeeno$warning)==0){ipwgeeno_warn=c(ipwgeeno_warn,tt)}
      
      ipwgeeno2=geese(formu,data=d4,id=cluster,
                      family = binomial(link='logit'),
                      weights = weight_no,
                      corstr = 'exchangeable')
      tempgeeno_est2=summary(ipwgeeno2)$mean['arm','estimate']
      tempgeeno_std2=summary(ipwgeeno2)$mean['arm','san.se']
      ipwgeeno_est2=c(ipwgeeno_est2,tempgeeno_est2)
      ipwgeeno_std2=c(ipwgeeno_std2,tempgeeno_std2)
      
      ipwgeeno2=myTryCatch(geese(formu,data=d4,id=cluster,
                                 family = binomial(link='logit'),
                                 weights = weight_no_round,
                                 corstr = 'exchangeable'))
      if(is.null(ipwgeeno2$warning)==0){
        ipwgeeno_warn2=c(ipwgeeno_warn2,tt)}
      
      ### IPW-GEE CRTgeeDR, no cluster
      ipwcrtgeedr=myTryCatch(geeDREstimation(formula=y~x+arm,
                                             id="cluster" , data = d3,
                                             nameMISS='missing',nameY='y',
                                             nameTRT='arm',
                                             family =  binomial("logit"), 
                                             corstr = "independence",
                                             weights = d3$weight_no))
      tempcrtgee_est=summary(ipwcrtgeedr$value)$beta[3]
      tempcrtgee_std=summary(ipwcrtgeedr$value)$se.robust[3]
      ipwcrtgeedr_est=c(ipwcrtgeedr_est,tempcrtgee_est)
      ipwcrtgeedr_std=c(ipwcrtgeedr_std,tempcrtgee_std)
      if(ipwcrtgeedr$value$converged==0){ipwcrtgee_warn=c(ipwcrtgee_warn,tt)}
      
      ipwcrtgeedr2=myTryCatch(geeDREstimation(formula=y~x+arm,
                                              id="cluster" , data = d3,
                                              nameMISS='missing',nameY='y',
                                              nameTRT='arm',
                                              family =  binomial("logit"), 
                                              corstr = "exchangeable",
                                              weights = d3$weight_no))
      tempcrtgee_est2=summary(ipwcrtgeedr2$value)$beta[3]
      tempcrtgee_std2=summary(ipwcrtgeedr2$value)$se.robust[3]
      ipwcrtgeedr_est2=c(ipwcrtgeedr_est2,tempcrtgee_est2)
      ipwcrtgeedr_std2=c(ipwcrtgeedr_std2,tempcrtgee_std2)
      if(is.null(ipwcrtgeedr2$warning)==0){
        ipwcrtgee_warn2=c(ipwcrtgee_warn2,tt)}
     
     #ipw geese with cluster effects
      ipwgeeyes=geese(formu,data=d4,id=cluster,
                      family = binomial(link='logit'),
                      weights = weight_clu,
                      corstr = 'independence')
      
      tempgeeyes_est=summary(ipwgeeyes)$mean['arm','estimate']
      tempgeeyes_std=summary(ipwgeeyes)$mean['arm','san.se']
      ipwgeeyes_est=c(ipwgeeyes_est,tempgeeyes_est)
      ipwgeeyes_std=c(ipwgeeyes_std,tempgeeyes_std)
      
      ipwgeeyes=myTryCatch(geese(formu,data=d4,id=cluster,
                                 family = binomial(link='logit'),
                                 weights = weight_clu_round,
                                 corstr = 'independence'))
      if(is.null(ipwgeeyes$warning)==0){
        ipwgeeyes_warn=c(ipwgeeyes_warn,tt)}
      
      ipwgeeyes2=geese(formu,data=d4,id=cluster,
                       family = binomial(link='logit'),
                       weights = weight_clu,
                       corstr = 'exchangeable')
      tempgeeyes_est2=summary(ipwgeeyes2)$mean['arm','estimate']
      tempgeeyes_std2=summary(ipwgeeyes2)$mean['arm','san.se']
      ipwgeeyes_est2=c(ipwgeeyes_est2,tempgeeyes_est2)
      ipwgeeyes_std2=c(ipwgeeyes_std2,tempgeeyes_std2)
      
      ipwgeeyes2=myTryCatch(geese(formu,data=d4,id=cluster,
                                  family = binomial(link='logit'),
                                  weights = weight_clu_round,
                                  corstr = 'exchangeable'))
      if(is.null(ipwgeeyes2$warning)==0){
        ipwgeeyes_warn2=c(ipwgeeyes_warn2,tt)}
      
      ### IPW-GEE with cluster effect with package CRTgeeDR
      ipwcrtgeedr_yes=myTryCatch(geeDREstimation(formula=y~x+arm,
                                                 id="cluster" , data = d3,
                                                 weights = d3$weight_clu,
                                                 nameMISS='missing',nameY='y',
                                                 nameTRT='arm',
                                                 family =  binomial("logit"), 
                                                 corstr = "independence"))
      tempcrtgee_yes_est=summary(ipwcrtgeedr_yes$value)$beta[3]
      tempcrtgee_yes_std=summary(ipwcrtgeedr_yes$value)$se.robust[3]
      ipwcrtgeedr_yes_est=c(ipwcrtgeedr_yes_est,tempcrtgee_yes_est)
      ipwcrtgeedr_yes_std=c(ipwcrtgeedr_yes_std,tempcrtgee_yes_std)
      if(is.null(ipwcrtgeedr_yes$warning)==0){ipwcrtgee_yes_warn=c(ipwcrtgee_yes_warn,tt)}
      
      ipwcrtgeedr_yes2=myTryCatch(geeDREstimation(formula=y~x+arm,
                                                  id="cluster" , data = d3,
                                                  weights = d3$weight_clu,
                                                  nameMISS='missing',nameY='y',
                                                  nameTRT='arm',
                                                  family =  binomial("logit"), 
                                                  corstr = "exchangeable"))
       
      tempcrtgee_yes_est2=summary(ipwcrtgeedr_yes2$value)$beta[3]
      tempcrtgee_yes_std2=summary(ipwcrtgeedr_yes2$value)$se.robust[3]
      ipwcrtgeedr_yes_est2=c(ipwcrtgeedr_yes_est2,tempcrtgee_yes_est2)
      ipwcrtgeedr_yes_std2=c(ipwcrtgeedr_yes_std2,tempcrtgee_yes_std2)
      
      if(ipwcrtgeedr_yes2$value$converged==0){
        ipwcrtgee_yes_warn2=c(ipwcrtgee_yes_warn2,tt)}
    }
  }}

data_result=data.frame(true_ind=true_est,
                       true_exh=true_est2,
                       uncra_ind=uncra_est,
                       uncra_exh=uncra_est2,
                       cra_est_ind=cra_est,
                       cra_est_exh=cra_est2,
                       ipwgeeno_est_ind=ipwgeeno_est,
                       ipwgeeno_est_exh=ipwgeeno_est2,
                       ipwcrtgeedr_est_ind=ipwcrtgeedr_est,
                       ipwcrtgeedr_est_exh=ipwcrtgeedr_est2,
                       ipwgeeyes_est_ind=ipwgeeyes_est,
                       ipwgeeyes_est_exh=ipwgeeyes_est2,
                       ipwcrtgeedr_yes_est_ind=ipwcrtgeedr_yes_est,
                       ipwcrtgeedr_yes_est_exh=ipwcrtgeedr_yes_est2,
                       uncra_std_ind=uncra_std,
                       uncra_std_exh=uncra_std2,
                       cra_std_ind=cra_std,
                       cra_std_exh=cra_std2,
                       ipwgeeno_std_ind=ipwgeeno_std,
                       ipwgeeno_std_exh=ipwgeeno_std2,
                       ipwcrtgeedr_std_ind=ipwcrtgeedr_std,
                       ipwcrtgeedr_std_exh=ipwcrtgeedr_std2,
                       ipwgeeyes_std_ind=ipwgeeyes_std,
                       ipwgeeyes_std_exh=ipwgeeyes_std2,
                       ipwcrtgeedr_yes_std_ind=ipwcrtgeedr_yes_std,
                       ipwcrtgeedr_yes_std_exh=ipwcrtgeedr_yes_std2,
                       times=data_tt,
                       km=data_km,
                       missing=missing)
data_warn=list(true_warning=true_warning,
                     true_warning2=true_warning2,
                     uncra_warning=uncra_warning,
                     uncra_warning2=uncra_warning2,
                     cra_warning=cra_warning,
                     cra_warning2=cra_warning2,
                     ipwgeeno_warn=ipwgeeno_warn,
                     ipwgeeno_warn2=ipwgeeno_warn2,
                     ipwcrtgee_warn=ipwcrtgee_warn,
                     ipwcrtgee_warn2=ipwcrtgee_warn2,
               ipwgeeyes_warn=ipwgeeyes_warn,
               ipwgeeyes_warn2=ipwgeeyes_warn2,
               ipwcrtgee_yes_warn=ipwcrtgee_yes_warn,
               ipwcrtgee_yes_warn2=ipwcrtgee_yes_warn2)     
return(list(data_result=data_result,data_warn=data_warn))
}

d11=ipw_gee_ana_clu(200,s=1,icc=0.05,mis=4)
d22=ipw_gee_ana_clu(5,s=1,icc=0.05,mis=5)
d33=ipw_gee_ana_clu(200,s=1,icc=0.05,mis=6)

d44=ipw_gee_ana_clu(200,s=2,icc=0.05,mis=4)
d55=ipw_gee_ana_clu(200,s=2,icc=0.05,mis=5)
d66=ipw_gee_ana_clu(200,s=2,icc=0.05,mis=6)


one_group=function(i,k,mm,seed=123,icc=0.05,s=1,mis=1){
  set.seed(seed)
  ## parameters
  if(s==1){b2=1;psi=-1.34}
  if(s==2){b2=0.588;psi=0.65}
  b0=1;b1=1.36
  sigma_b=sqrt(0.2)
  mu_x=0
  sigma_x=3.55
  sigam_alpha=sigma_x*icc
  sigam_u=sigma_x-sigam_alpha
  phi=1
  if(mm==25){
    x=matrix(0,k,mm+5)
    y=matrix(0,k,mm+5)
    pi=matrix(0,k,mm+5)
    cluster=matrix(0,k,mm+5)
    r=matrix(0,k,mm+5)
  }
  if(mm==50){
    x=matrix(0,k,mm+10)
    y=matrix(0,k,mm+10)
    pi=matrix(0,k,mm+10)
    cluster=matrix(0,k,mm+10)
    r=matrix(0,k,mm+10)
  }
  for(k in 1:k){
    delta=rnorm(1,0,sigma_b)
    alpha=rnorm(1,mu_x,sigma_alpha)
    if(mm==25){m=round(runif(1,mm-5,mm+5))}
    if(mm==50){m=round(runif(1,mm-10,mm+10))}
    # print(m)
    cluster[k,]=k
    for(j in 1:m){
      u=rnorm(1,0,sigma_u)
      x[k,j]=u+alpha
      pi[k,j]=expit(b0+b1*i+b2*x[k,j]+delta)
      y[k,j]=rbinom(1,1,pi[k,j])
    }
  }
  if(mis==1){r=expit(psi+phi*x)}
  if(mis==2){r=expit(psi+phi*x+i)}
  if(mis==3){r=expit(psi+phi*x+i+x*i)}
  return(list(x=x,y=y,pi=pi,r=r,cluster=cluster))
}
data_gene=function(k,m,seed=123,s=1,icc=0.05,mis=1){
  set.seed(seed)
  # step1:
  a1=one_group(1,k,m,seed,icc,s,mis)
  a0=one_group(0,k,m,seed,icc,s,mis)
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
mypool=function(mean0,sd0,num=5,print='no'){
  m=mean(mean0)
  v=mean(sd0)
  B=sd(mean0)
  v_hat=v+(1+1/num)*B
  l=m-1.96*v_hat
  u=m+1.96*v_hat
  if(print=='no'){
    return(list(mean=m,std=v_hat))
  }
  if(print=='yes'){
    print('mean (95% CI)')
    print(paste(round(m,2)," (",round(l,2),',',round(u,2),')',sep=''))
    return(list(mean=m,std=v_hat))
  }
}

miss_per=function(a){
  b=sum(a$R)/dim(a)[1]
  return(b)
}

ipw_gee_ana=function(times,s=1,icc=0.05,mis=1){
  times=times
  ## empty effects
  for(emp in 1){
    true_est=c();true_warning=c()
    true_est2=c();true_warning2=c()
    
    uncra_est=c();uncra_std=c();uncra_warning=c()
    uncra_est2=c();uncra_std2=c();uncra_warning2=c()
    
    cra_est=c();cra_std=c();cra_warning=c()
    cra_est2=c();cra_std2=c();cra_warning2=c()
    
    ipwgeeno_est=c();ipwgeeno_std=c();ipwgeeno_warn=c()
    ipwgeeno_est2=c();ipwgeeno_std2=c();ipwgeeno_warn2=c()
    
    ipwcrtgeedr_est=c();ipwcrtgeedr_std=c();ipwcrtgee_warn=c()
    ipwcrtgeedr_est2=c();ipwcrtgeedr_std2=c();ipwcrtgee_warn2=c()
    
    ipwgeeyes_est=c();ipwgeeyes_std=c();ipwgeeyes_warn=c()
    ipwgeeyes_est2=c();ipwgeeyes_std2=c();ipwgeeyes_warn2=c()
    
    ipwcrtgeedr_yes_est=c();ipwcrtgeedr_yes_std=c();ipwcrtgee_yes_warn=c()
    ipwcrtgeedr_yes_est2=c();ipwcrtgeedr_yes_std2=c();ipwcrtgee_yes_warn2=c()
    
    data_tt=c()
    data_km=c()
    
    missing=c()
  }
  
  for(km in 1:2){
    if(km==1){k=25;m=25}
    if(km==2){k=50;m=25}
    if(km==3){k=25;m=50}
    if(km==4){k=50;m=50}
    print("***************")
    print("km")
    print(km)
    print('times')
    for(tt in 1:times){
      data_tt=c(data_tt,tt)
      data_km=c(data_km,km)
      print(tt)
      d1=data_gene(k=k,m=m,seed=tt,s=s,icc=icc,mis=mis)
      missing=c(missing,miss_per(d1))
      d2=d1
      d2$y=ifelse(d2$R==1,NA,d2$y)
      d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=d2$R)
      
      l=glm(missing ~ x , data = d3, 
            family = binomial(link='logit'))
      weight=expit(predict(l))
      weight=ifelse(weight<0.005,0.005,weight)
      d3$weight=1/weight
      d3$weight2=round(1/weight)
      d4=na.omit(d3)
      
      formu=formula(y~x+arm)
      true=myTryCatch(geese(formu,data=d1,id=cluster,
                            family = binomial(link='logit'),
                            corstr = 'independence'))
      true_est=c(true_est,summary(true$value)$mean['arm','estimate'])
      if(is.null(true$warning)==0){
        true_warning=c(true_warning,tt)}
      
      true2=myTryCatch(geese(formu,data=d1,id=cluster,
                             family = binomial(link='logit'),
                             corstr = 'exchangeable'))
      true_est2=c(true_est2,summary(true2$value)$mean['arm','estimate'])
      if(is.null(true2$warning)==0){
        true_warning2=c(true_warning2,tt)}
      
      ### unadjusted CRA
      d_nona=na.omit(d1)
      formu0=formula(y~arm)
      uncra=myTryCatch(geese(formu0,data=d_nona,id=cluster,
                             family = binomial(link='logit'),
                             corstr = 'independence'))
      uncra_est=c(uncra_est,summary(uncra$value)$mean['arm','estimate'])
      uncra_std=c(uncra_std,summary(uncra$value)$mean['arm','san.se'])
      if(is.null(uncra$warning)==0){
        uncra_warning=c(uncra_warning,tt)}
      
      uncra2=myTryCatch(geese(formu0,data=d_nona,id=cluster,
                              family = binomial(link='logit'),
                              corstr = 'exchangeable'))
      uncra_est2=c(uncra_est2,summary(uncra2$value)$mean['arm','estimate'])
      uncra_std2=c(uncra_std2,summary(uncra2$value)$mean['arm','san.se'])
      if(is.null(uncra2$warning)==0){
        uncra_warning2=c(uncra_warning2,tt)}
      
      ### adjusted CRA
      cra=myTryCatch(geese(formu,data=d_nona,id=cluster,
                           family = binomial(link='logit'),
                           corstr = 'independence'))
      cra_est=c(cra_est,summary(cra$value)$mean['arm','estimate'])
      cra_std=c(cra_std,summary(cra$value)$mean['arm','san.se'])
      if(is.null(cra$warning)==0){
        cra_warning=c(cra_warning,tt)}
      
      cra2=myTryCatch(geese(formu,data=d_nona,id=cluster,
                            family = binomial(link='logit'),
                            corstr = 'exchangeable'))
      cra_est2=c(cra_est2,summary(cra2$value)$mean['arm','estimate'])
      cra_std2=c(cra_std2,summary(cra2$value)$mean['arm','san.se'])
      if(is.null(cra2$warning)==0){
        cra_warning2=c(cra_warning2,tt)}
      
      ### IPW-no cluster effect, by using glm and geese
      ipwgeeno0=geese(formu,data=d4,id=cluster,
                      family = binomial(link='logit'),
                      weights = weight,
                      corstr = 'independence')
      tempgeeno_est=summary(ipwgeeno0)$mean['arm','estimate']
      tempgeeno_std=summary(ipwgeeno0)$mean['arm','san.se']
      
      ipwgeeno_est=c(ipwgeeno_est,tempgeeno_est)
      ipwgeeno_std=c(ipwgeeno_std,tempgeeno_std)
      
      ipwgeeno=myTryCatch(geese(formu,data=d4,id=cluster,
                                family = binomial(link='logit'),
                                weights = weight2,
                                corstr = 'independence'))
      if(is.null(ipwgeeno$warning)==0){ipwgeeno_warn=c(ipwgeeno_warn,tt)}
      
      ipwgeeno2=geese(formu,data=d4,id=cluster,
                      family = binomial(link='logit'),
                      weights = weight,
                      corstr = 'exchangeable')
      tempgeeno_est2=summary(ipwgeeno2)$mean['arm','estimate']
      tempgeeno_std2=summary(ipwgeeno2)$mean['arm','san.se']
      
      ipwgeeno_est2=c(ipwgeeno_est2,tempgeeno_est2)
      ipwgeeno_std2=c(ipwgeeno_std2,tempgeeno_std2)
      
      ipwgeeno2=myTryCatch(geese(formu,data=d4,id=cluster,
                                 family = binomial(link='logit'),
                                 weights = weight2,
                                 corstr = 'exchangeable'))
      if(is.null(ipwgeeno2$warning)==0){
        ipwgeeno_warn2=c(ipwgeeno_warn2,tt)}
      
      ### IPW-GEE with package CRTgeeDR
      ipwcrtgeedr=myTryCatch(geeDREstimation(formula=y~x+arm,
                                             id="cluster" , data = d3,
                                             nameMISS='missing',nameY='y',
                                             nameTRT='arm',
                                             family =  binomial("logit"), 
                                             corstr = "independence",
                                             model.weights=I(missing==1)~x))
      if(is.null(ipwcrtgeedr$value)==0){
        tempcrtgee_est=summary(ipwcrtgeedr$value)$beta[3]
        tempcrtgee_std=summary(ipwcrtgeedr$value)$se.robust[3]
        ipwcrtgeedr_est=c(ipwcrtgeedr_est,tempcrtgee_est)
        ipwcrtgeedr_std=c(ipwcrtgeedr_std,tempcrtgee_std)}
      if(ipwcrtgeedr$value$converged==0){
        ipwcrtgee_warn=c(ipwcrtgee_warn,tt)}
      if(is.null(ipwcrtgeedr$value)==1){
        ipwcrtgeedr_est=c(ipwcrtgeedr_est,NA)
        ipwcrtgeedr_std=c(ipwcrtgeedr_std,NA)}
      
      ipwcrtgeedr2=myTryCatch(geeDREstimation(formula=y~x+arm,
                                              id="cluster" , data = d3,
                                              nameMISS='missing',nameY='y',
                                              nameTRT='arm',
                                              family =  binomial("logit"), 
                                              corstr = "exchangeable",
                                              model.weights=I(missing==1)~x))
      
      if(is.null(ipwcrtgeedr2$value)==0){
        tempcrtgee_est2=summary(ipwcrtgeedr2$value)$beta[3]
        tempcrtgee_std2=summary(ipwcrtgeedr2$value)$se.robust[3]
        ipwcrtgeedr_est2=c(ipwcrtgeedr_est2,tempcrtgee_est2)
        ipwcrtgeedr_std2=c(ipwcrtgeedr_std2,tempcrtgee_std2)
        if(is.null(ipwcrtgeedr2$warning)==0){
          ipwcrtgee_warn2=c(ipwcrtgee_warn2,tt)}
      }
      if(is.null(ipwcrtgeedr2$value)==1){
        ipwcrtgeedr_est2=c(ipwcrtgeedr_est2,NA)
        ipwcrtgeedr_std2=c(ipwcrtgeedr_std2,NA)}
    }}
  
  data_result=data.frame(true_ind=true_est,
                         true_exh=true_est2,
                         uncra_ind=uncra_est,
                         uncra_exh=uncra_est2,
                         cra_est_ind=cra_est,
                         cra_est_exh=cra_est2,
                         ipwgeeno_est_ind=ipwgeeno_est,
                         ipwgeeno_est_exh=ipwgeeno_est2,
                         ipwcrtgeedr_est_ind=ipwcrtgeedr_est,
                         ipwcrtgeedr_est_exh=ipwcrtgeedr_est2,
                         uncra_std_ind=uncra_std,
                         uncra_std_exh=uncra_std2,
                         cra_std_ind=cra_std,
                         cra_std_exh=cra_std2,
                         ipwgeeno_std_ind=ipwgeeno_std,
                         ipwgeeno_std_exh=ipwgeeno_std2,
                         ipwcrtgeedr_std_ind=ipwcrtgeedr_std,
                         ipwcrtgeedr_std_exh=ipwcrtgeedr_std2,
                         times=data_tt,
                         km=data_km,
                         missing=missing)
  data_warn=list(true_warning=true_warning,
                 true_warning2=true_warning2,
                 uncra_warning=uncra_warning,
                 uncra_warning2=uncra_warning2,
                 cra_warning=cra_warning,
                 cra_warning2=cra_warning2,
                 ipwgeeno_warn=ipwgeeno_warn,
                 ipwgeeno_warn2=ipwgeeno_warn2,
                 ipwcrtgee_warn=ipwcrtgee_warn,
                 ipwcrtgee_warn2=ipwcrtgee_warn2,
                 ipwgeeyes_warn=ipwgeeyes_warn,
                 ipwgeeyes_warn2=ipwgeeyes_warn2,
                 ipwcrtgee_yes_warn=ipwcrtgee_yes_warn,
                 ipwcrtgee_yes_warn2=ipwcrtgee_yes_warn2)     
  return(list(data_result=data_result,data_warn=data_warn))
}

d1=ipw_gee_ana(100,s=1,icc=0.05,mis=1)
d2=ipw_gee_ana(100,s=1,icc=0.05,mis=2)
d3=ipw_gee_ana(100,s=1,icc=0.05,mis=3)

d4=ipw_gee_ana(100,s=2,icc=0.05,mis=1)
d5=ipw_gee_ana(100,s=2,icc=0.05,mis=2)
d6=ipw_gee_ana(100,s=2,icc=0.05,mis=3)

d7=ipw_gee_ana(100,s=1,icc=0.1,mis=1)
d8=ipw_gee_ana(100,s=1,icc=0.1,mis=2)
d9=ipw_gee_ana(100,s=1,icc=0.1,mis=3)

ipwcrtgee_warn2
ipwgeeyes_est
ipwgeeyes_std
ipwgeeyes_warn
ipwgeeyes_est2
ipwgeeyes_std2
ipwgeeyes_warn2
ipwcrtgeedr_yes_est
ipwcrtgeedr_yes_std
ipwcrtgee_yes_warn
ipwcrtgeedr_yes_est2
ipwcrtgeedr_yes_std2
ipwcrtgee_yes_warn2
data_tt
data_km
