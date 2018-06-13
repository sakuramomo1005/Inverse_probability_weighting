library(lme4)
library(geepack)
library(jomo)
library(CRTgeeDR)

errors=0
warn=0
errors_ipw1=0
warn_ipw1=0
errors_ipw2=0
warn_ipw2=0
errors_jomo=0
warn_jomo=0
warn_n=0
warn_time1=c()
warn_time2=c()
ipw_mm=c()
ipw_mm_sd=c()
warnss=c()


### data generation function
one_group=function(i,k,m,seed=123){
  set.seed(seed)
  x=matrix(NA,k,m)
  y=matrix(NA,k,m)
  pi=matrix(NA,k,m)
  
  for(k in 1:k){
    delta=rnorm(1,0,sigma_b)
    alpha=rnorm(1,mu_x,sigma_alpha)
    for(j in 1:m){
      u=rnorm(1,0,sigma_u)
      x[k,j]=u+alpha
      pi[k,j]=expit(b0+b1*i+b2*x[k,j]+delta)
      y[k,j]=rbinom(1,1,pi[k,j])
    }
  }
  return(list(x=x,y=y,pi=pi))
}

data_generation=function(k,m,seed=123,print=0){
  set.seed(seed)
  trt=one_group(1,k,m,seed)
  con=one_group(0,k,m,seed)
  X=rbind(trt$x,con$x)
  Y=rbind(trt$y,con$y)
  #  missing
  r=expit(psi+phi*X)
  Y_mis=Y
  R=matrix(NA,2*k,m)
  for(i in 1:(2*k)){
    for(j in 1:m){
      R[i,j]=rbinom(1,1,r[i,j])
      if(R[i,j]==1){Y_mis[i,j]=NA}
    }
  }
  # missing percentage:
  if(print==1){
    print(mean(r))
    print(mean(expit(psi+phi*trt$x)))
    print(mean(expit(psi+phi*con$x)))
    print(sum(is.na(Y_mis))/(k*m*2))
    mispre=sum(is.na(Y_mis))/(k*m*2)
  }
  X=trans(X,'X')
  Y=trans(Y,'Y')
  Y_mis=trans(Y_mis,'Y_mis')
  R=trans(R,'R')
  return(list(X=X,Y=Y,Y_mis=Y_mis,R=R,mispre=mispre))
}

### expit function 
expit=function(x){
  y=exp(x)/(1+exp(x))
  return(y)
}

trans=function(X,name){
  x1=as.vector(X)
  x2=matrix(x1,m,k*2,byrow = TRUE)
  X=matrix(x2,2*k*m,1,byrow = TRUE)
  X=data.frame(X)
  colnames(X)=name
  X$group=c(rep(1,k*m),rep(0,k*m))
  X$cluster=rep(1:(2*k),each=m)
  return(X)
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
analysis_results_mmi=function(mmi){
  m=c()
  std=c()
  #icc=c()
  for(i in 1:5){
    temp=mmi[mmi$Imputation==i,]
    #temp=temp[-1,]
    temp0=temp[,c('group','X','cluster','Y_mis')]
    formu=formula(Y_mis~X+group)
    temp0$Y_mis=as.numeric(temp0$Y_mis)
    temp0$Y_mis=temp0$Y_mis-1
    mmi2=geese(formu,data=temp0,id=cluster,
               family = binomial(link='logit'),
               corstr = 'independence')
    est=summary(mmi2)
    est_trt=est$mean['group','estimate']
    sd_trt=est$mean['group','san.se']
    m=c(m,est_trt)
    std=c(std,sd_trt)
  }
  return(mypool(m,std,print='yes'))
}

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

## parameters
b0=1;b1=1.36;b2=1
sigma_b=sqrt(0.2)
mu_x=0
sigma_alpha=sqrt(0.18)
sigma_u=sqrt(3.37)
psi=-1.34
phi=1
nburn=100
nbetween=1000
nimp=5

start_time=Sys.time()
for(km in 1:4){
  if(km==1){
    k=25
    m=25}
  if(km==2){
    k=50
    m=25}
  if(km==3){
    k=25
    m=50}
  if(km==4){
    k=50
    m=50}
### Data Generation
  
cca_mean=c()
cca_std=c()
ipw_mean=c()
ipw_std=c()
ipw_mean2=c()
ipw_std2=c()
jomo_mean=c()
jomo_std=c()
missingper=c()
error_cca=c()
warning_cca=c()
error_ipw1=c()
warning_ipw1=c()
error_ipw2_glmer=c()
warning_ipw2_glmer=c()
error_ipw2_gee=c()
warning_ipw2_gee=c()
error_mmi=c()
warning_mmi=c()
error_mmi_result=c()
warning_mmi_result=c()

for(times in 1:1000){
# generate the orginal dataset
print("*************")
print(times)
data_sim=data_generation(k,m,times+130,1)
missingper=c(missingper,data_sim$mispre)
## functions
# 1. pool function
## ipw-gee no cluster effects
#install.packages('geepack')
mis=data_sim$R
mis$X=data_sim$X$X
head(mis)
# estimate the model and store results in m
logs=glm(R ~ X , data = mis, 
         family = binomial(link='logit'))
logsum=summary(logs)
weight=expit(logsum$coefficients[,'Estimate'][1]+logsum$coefficients[,'Estimate'][2]*mis$X)
weight=matrix(weight,((2*k)*m),1)

org_data=data_sim$Y_mis
org_data$X=data_sim$X$X
org_data$weight=round(1/weight)
head(org_data)
cca=na.omit(org_data)
dim(cca)
formu=formula(Y_mis~X+group)
cca_adj=myTryCatch(
  geese(formu,data=cca,id=cluster,
              family = binomial(link='logit'),
              corstr = 'independence'))
error_cca=c(error_cca,cca_adj$error)
warning_cca=c(warning_cca,cca_adj$warning)
if(is.null(cca_adj$error)==0){errors=errors+1}
if(is.null(cca_adj$warning)==0){warn=warn+1}

cca_adj=cca_adj$value
est0=summary(cca_adj)
est_trt0=est0$mean['group','estimate']
sd_trt0=est0$mean['group','san.se']
cca_mean=c(cca_mean,est_trt0)
cca_std=c(cca_std,sd_trt0)

ipw=myTryCatch(
geese(formu,data=cca,id=cluster,
           family = binomial(link='logit'),
           weights = cca$weight,
           corstr = 'independence'))

error_ipw1=c(error_ipw1,ipw$error)
warning_ipw1=c(warning_ipw1,ipw$warning)
if(is.null(ipw$error)==0){errors_ipw1=errors_ipw1+1}
if(is.null(ipw$warning)==0){
  warn_ipw1=warn_ipw1+1
  warn_time1=c(warn_time1,times)}

ipw=ipw$value
est=summary(ipw)
est_trt=est$mean['group','estimate']
sd_trt=est$mean['group','san.se']
ipw_mean=c(ipw_mean,est_trt)
ipw_std=c(ipw_std,sd_trt)

## ipw-gee with cluster effects
mis2=data_sim$R
mis2$X=data_sim$X$X
head(mis2)
# estimate the model and store results in m
logs2 = myTryCatch(
  glmer(R ~ X+(1|cluster) , data = mis2,
              family = binomial(link='logit')))
error_ipw2_glmer=c(error_ipw2_glmer,logs2$error)
warning_ipw2_glmer=c(warning_ipw2_glmer,logs2$warning)

logs2=logs2$value
logsum2=summary(logs2)

ran=as.data.frame(VarCorr(logs2))['sdcor']
ran=as.numeric(ran)
rans=c()

for(i in 1:(2*k)){
  rans=c(rans,rnorm(1,0,sqrt(ran)))
}

rans2=rep(rans,each=m)
mis2$rans=rans2
head(mis2)
weight=expit(logsum$coefficients[,'Estimate'][1]+
               logsum$coefficients[,'Estimate'][2]*mis2$X+
               mis2$rans)

org_data2=data_sim$Y_mis
org_data2$X=data_sim$X$X
org_data2$weight=round(1/weight)
head(org_data2)
cca2=na.omit(org_data2)
head(cca2)

formu=formula(Y_mis~X+group)
ipw2=myTryCatch(
  geese(formu,data=cca2,id=cluster,
            family = binomial(link='logit'),
            #family=quasibinomial(),
            weights = cca2$weight,
            corstr = 'independence'))
error_ipw2_gee=c(error_ipw2_gee,ipw2$error)
warning_ipw2_gee=c(warning_ipw2_gee,ipw2$warning)
if(is.null(ipw2$error)==0){errors_ipw2=errors_ipw2+1}
if(is.null(ipw2$warning)==0){
  warn_ipw2=warn_ipw2+1
  warn_time2=c(warn_time2,times)}

ipw2=ipw2$value
est2=summary(ipw2)
est_trt2=est2$mean['group','estimate']
sd_trt2=est2$mean['group','san.se']

ipw_mean2=c(ipw_mean2,est_trt2)
ipw_std2=c(ipw_std2,sd_trt2)

## IPW-GEE using CRTgeeDR package
org_data=data_sim$Y_mis
org_data$X=data_sim$X$X
org_data$R=data_sim$R$R
names(org_data)=c('OUTCOME','TRT','CLUSTER','AGE','MISSING')
head(org_data)
ipwresults<-myTryCatch(geeDREstimation(formula=OUTCOME~TRT+AGE,
                                       id="CLUSTER" , data = org_data,
                                       family = "binomial", corstr = "independence",
                                       model.weights=I(MISSING==0)~AGE))
ipwresults=ipwresults$value
warnss=c(warnss,ipwresults$warning)
if(is.null(ipwresults$warning)==0){warn_n=warn_n+1}
res_m=summary(ipwresults)
ipw_mm=c(ipw_mm,res_m$beta[2])
ipw_mm_sd=c(ipw_mm_sd,res_m$se.model[2])
## MMI-GEE
mmidata=org_data[,1:4]
head(mmidata)
cluster=mmidata$cluster
mmidata$Y_mis=as.factor(mmidata$Y_mis)
mmi=myTryCatch(
  jomo(mmidata,nburn=nburn,clus=cluster,
         nbetween=nbetween,nimp=nimp,output = 0))
error_mmi=c(error_mmi,mmi$error)
warning_mmi=c(warning_mmi,mmi$warning)

mmi=mmi$value
m0=c()
std0=c()
#icc=c()
result_mmi=function(mmi){
  for(i in 1:5){
    temp=mmi[mmi$Imputation==i,]
    #temp=temp[-1,]
    temp0=temp[,c('group','X','cluster','Y_mis')]
    formu=formula(Y_mis~X+group)
    temp0$Y_mis=as.numeric(temp0$Y_mis)
    temp0$Y_mis=temp0$Y_mis-1
    mmi2=geese(formu,data=temp0,id=cluster,
               family = binomial(link='logit'),
               corstr = 'independence')
    est=summary(mmi2)
    est_trt=est$mean['group','estimate']
    sd_trt=est$mean['group','san.se']
    m0=c(m0,est_trt)
    std0=c(std0,sd_trt)
  }
  return(list(m0,std0))
}
rmmi=myTryCatch(
  result_mmi(mmi))
error_mmi_result=c(error_mmi_result,rmmi$error)
warning_mmi_result=c(warning_mmi_result,rmmi$warning)
if(is.null(rmmi$error)==0){errors_jomo=errors_jomo+1}
if(is.null(rmmi$warning)==0){warn_jomo=warn_jomo+1}

rmmi=rmmi$value
m0=rmmi[[1]]
std0=rmmi[[2]]
mmi_result=mypool(m0,std0,print='yes')
mmi_m=mmi_result$mean
mmi_s=mmi_result$std
jomo_mean=c(jomo_mean,mmi_m)
jomo_std=c(jomo_std,mmi_s)
}

end_time=Sys.time()
end_time-start_time

#save the results
final_data=list(
  cca_mean=cca_mean,
  cca_std=cca_std,
  ipw_mean=ipw_mean,
  ipw_std=ipw_std,
  ipw_mean2=ipw_mean2,
  ipw_std2=ipw_std2,
  jomo_mean= jomo_mean,
  jomo_std=jomo_std,
  missingper=missingper,
  error_cca=error_cca,
  warning_cca=warning_cca,
  error_ipw1=error_ipw1,
  warning_ipw1=warning_ipw1,
  error_ipw2_glmer=error_ipw2_glmer,
  warning_ipw2_glmer=warning_ipw2_glmer,
  error_ipw2_gee=error_ipw2_gee,
  warning_ipw2_gee=warning_ipw2_gee,
  error_mmi=error_mmi,
  warning_mmi=warning_mmi,
  error_mmi_result=error_mmi_result,
  warning_mmi_result=warning_mmi_result,
  warn_ipw1=warn_ipw1,
  warn_ipw2=warn_ipw2,
  warn_time1= warn_time1,
  warn_time2=warn_time2,
  ipw_mm=ipw_mm,
  ipw_mm_sd=ipw_mm_sd,
  warnss=warnss,
  warn_n=warn
)
if(km==1){save(final_data,file='ipw.RData')}
if(km==2){save(final_data,file='ipw2.RData')}
if(km==3){save(final_data,file='ipw3.RData')}
if(km==4){save(final_data,file='ipw4.RData')}
}






