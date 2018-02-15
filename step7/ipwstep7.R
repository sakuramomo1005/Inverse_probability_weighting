## from missingness ICC calculate variance
missing_icc=function(rho){
  pi=3.1415926
  sigma=(pi^2/3)/(1/rho-1)
  return(sigma)
}

# for example, calculate the variance of ICC 0.05
missingicc=missing_icc(0.05)


library(lme4)
library(geepack)
library(CRTgeeDR)

## function to catch errors and warns
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

# function for generate one intervention arm
# k cluster numbers; mm, cluster size; intercept: missingness generation intercpet; vars: uijl's variance
one_group3=function(mis=5,i,k,mm,icc,intercept,seed=123,vars){
  set.seed(seed)
  #if(s==1){b0=1;b1=1.36;b2=1;psi=-1.8}
  #if(s==3){b0=1;b1=1.36;b2=0.588;psi=-1.8}
  ## parameters
  b0=1;b1=1.36;b2=1;psi=intercept
  sigma_b=sqrt(0.2)
  mu_x=0
  sigma_alpha=sqrt(0.18)
  sigma_u=sqrt(vars) # change the variance 
  
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
    sigma_icc=missing_icc(icc) ## missingness ICC
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

# function that combines two intervention arm

data_gene3=function(k,m,s,mis=5,icc,intercept,seed=123,vars){
  set.seed(seed)
  # step1:
  a1=one_group3(i=1,k=k,mm=m,seed=seed,mis=mis,intercept=intercept,icc=icc,vars=vars)
  a0=one_group3(i=0,k=k,mm=m,seed=seed,mis=mis,intercept=intercept,icc=icc,vars=vars)
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

# function for calculate intercept in missingess generation model
mis_per=function(n,m.bar,icc,vars,seed=123){
  data=c()
  for(i in seq(-5,5,0.1)){
    temp=data_gene3(k=n,m=m.bar,icc=icc,intercept=i,seed=seed,vars=vars)
    percent=sum(temp$R)/dim(temp)[1]
    data=rbind(data,c(percent,i))
    #print(percent)
  }
  data[,1]=abs(data[,1]-0.3)
  intercept=data[data[,1]==min(data[,1])][2]
  return(intercept)
}

## function to change missing percentage
mis_per_check=function(vars,n,m,icc,intercept,seed=123){
  temp=data_gene3(k=n,m=m,icc=icc,intercept=intercept,seed,vars=vars)
  percent=sum(temp$R)/dim(temp)[1]
 # print(percent)
}


## calculate intercepts in different scenarions
temp1=c()
for(var in c(0.2,0.5,0.8,1.1)){
  n=50;m.bar=50;icc=0.5
  data=c()
  for(i in seq(-5,5,0.1)){
    temp=data_gene3(k=n,m=m.bar,icc=icc,intercept=i,seed=i,vars=var)
    percent=sum(temp$R)/dim(temp)[1]
    data=rbind(data,c(percent,i))
  }
  data[,1]=abs(data[,1]-0.3)
  intercept=data[data[,1]==min(data[,1])][2]
  print(data[data[,1]==min(data[,1])])
  print(intercept)
  temp1=rbind(temp1,c(var,intercept))
}

## calculate intercepts in different scenarions
temp2=c()
for(icc in seq(0.1,0.9,0.1)){
  var=0.2
  n=50;m.bar=50
  data=c()
  for(i in seq(-5,5,0.1)){
    temp=data_gene3(k=n,m=m.bar,icc=icc,intercept=i,seed=i,vars=var)
    percent=sum(temp$R)/dim(temp)[1]
    data=rbind(data,c(percent,i))
  }
  data[,1]=abs(data[,1]-0.3)
  intercept=data[data[,1]==min(data[,1])][2]
  print(data[data[,1]==min(data[,1])])
  print(intercept)
  temp2=rbind(temp2,c(icc,intercept))
}
temp1=data.frame(temp1)
names(temp1)=c('var','intercept')
temp2=data.frame(temp2)
names(temp2)=c('icc','intercept')


## simulations
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
  
  for(time in 1:1000){
    print('icc')
    print(icc)
    
    print(time)
    
    intercept=temp2[temp2$icc==icc,]$intercept
    mp=mis_per_check(n=50,m=50,icc=icc,seed=time,intercept=intercept,vars=0.2)
    
    d1=data_gene3(k=50,m=50,icc=icc,intercept=intercept,seed=time,vars=0.2)
    print(dim(d1))
    d2=d1
    d2$y=ifelse(d2$R==0,NA,d2$y)
    d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=d2$R)
    
    # weight calculation
    logs=glm(missing ~ x+arm, data = d3,
             family = binomial(link='logit'))
    logs2=glmer(missing ~ x+arm+(1|cluster) , data = d3,
                family = binomial(link='logit'))
    
    weight=predict(logs,type="response")
    weight2=predict(logs2,type="response")
    
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
    d3=na.omit(d3)
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
  
  mp=c(mp,intercept)
  names(EST)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')
  names(STD)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')
  names(WARN)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')
  
  result=list(EST=EST,STD=STD,WARN=WARN,mp=mp)
  name_res=paste('ipw_result554',icc,'.RData',sep='')
  save(result,file=name_res)
 
}

### readin the results and draw a table
setwd('')
icc=0.1
table_ind=c()
table_exch=c()
for(icc in seq(0.1,0.9,0.1)){
  names=paste('ipw_result54',icc,'.RData',sep='')
  load(names)
  print(result$mp)
  est=result$EST
 
  std=result$STD
  warn=result$WARN
  apply(warn[,1:6],2,sum)
  
  for(i in 1:6){
    est[warn[,i]>0,i]=NA
    std[warn[,i]>0,i]=NA
  }
  mean_value=apply(est,2,mean,na.rm=TRUE)
  mcsd_value=apply(est,2,sd,na.rm=TRUE)
  std_value=apply(std,2,mean,na.rm=TRUE)
  
  true_ind=mean_value[1]
  true_ex=mean_value[2]
  
  est_n=na.omit(est)
  std_n=na.omit(std)
  non_con=(400-dim(est_n)[1])/400
  cova=((est_n-1.96*std_n)<true_ind) & ((est_n+1.96*std_n)>true_ind)
  cov_ind=apply(cova,2,sum)/dim(est_n)[1]
  cova=((est_n-1.96*std_n)<true_ex) & ((est_n+1.96*std_n)>true_ex)
  cov_ex=apply(cova,2,sum)/dim(est_n)[1]
  
  ind=c(mean_value[c(3,5)]-true_ind,mcsd_value[c(3,5)],std_value[c(3,5)],cov_ind[c(3,5)],non_con)
  ind=round(ind,3)
  table_ind=rbind(table_ind,ind)
  exch=c(mean_value[c(4,6)]-true_ex,mcsd_value[c(4,6)],std_value[c(4,6)],cov_ind[c(4,6)],non_con)
  exch=round(exch,3)
  table_exch=rbind(table_exch,exch)
}

