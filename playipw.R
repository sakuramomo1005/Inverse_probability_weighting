## Y_{ijl}=\alpha_i+\tau_i \sigma_y X_{ijl}+\delta_{ij}+\epslion_{ijl}
library('lme4')
library('jomo')
library('mitml')
expit=function(t){
  out=exp(t)/(1+exp(t))
  return(out)
}

## 1. pool function
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
    temp=mitmlComplete(mmi,i)
   # temp$group=sim_data$group
    model=lmer(y~x+group+(1|cluster),temp)
    summary( model)
    est_trt=summary(model)$coefficients[,'Estimate']['group']
    est_trt=as.vector(est_trt)
    sd_trt=summary(model)$coefficients[,'Std. Error']['group']
    sd_trt=as.vector(sd_trt)
    m=c(m,est_trt)
    std=c(std,sd_trt)
  }
  return(mypool(m,std,print='yes'))
}

analysis_results_mmi2=function(mmi){
  m=c()
  std=c()
  #icc=c()
  for(i in 1:5){
    temp=mmi[mmi$Imputation==i,]
    temp=temp[-1,]
    temp$group=sim_data$group
    model=lmer(y~x+group+(1|cluster),temp)
    summary( model)
    est_trt=summary(model)$coefficients[,'Estimate']['group']
    est_trt=as.vector(est_trt)
    sd_trt=summary(model)$coefficients[,'Std. Error']['group']
    sd_trt=as.vector(sd_trt)
    m=c(m,est_trt)
    std=c(std,sd_trt)
  }
  return(mypool(m,std,print='yes'))
}

analysis_results_mmi3=function(mmi){
  m=c()
  std=c()
  #icc=c()
  for(i in 1:5){
    temp=mmi[mmi$Imputation==i,]
    #temp=temp[-1,]
    temp$group=sim_data$group
    model=lmer(y~x+group+(1|cluster),temp)
    summary( model)
    est_trt=summary(model)$coefficients[,'Estimate']['group']
    est_trt=as.vector(est_trt)
    sd_trt=summary(model)$coefficients[,'Std. Error']['group']
    sd_trt=as.vector(sd_trt)
    m=c(m,est_trt)
    std=c(std,sd_trt)
  }
  return(mypool(m,std,print='yes'))
}

start_time=Sys.time()
#sleep_for_a_minute()

k=50
m=25

ipw_mean=c()
ipw_std=c()

jomo_mean=c()
jomo_std=c()

jomo_no_m=c()
jomo_no_std=c()

mmi_mean_pan=c()
mmi_std_pan=c()

mmi_mean_pan2=c()
mmi_std_pan2=c()

mmi_mean_jomo=c()
mmi_std_jomo=c()

mmi_mean_jomo2=c()
mmi_std_jomo2=c()

datageneration=function(trt,m=25,k=50,b0=1,b1=1,a1=20,a3=5,rou=0.1,sigma_y=10,tau=0.5){
  ## generate 
  x=matrix(NA,k,m)
  y=matrix(NA,k,m)
  R=matrix(NA,k,m)
  u=matrix(NA,k,m)
  trt=trt
  for( i in 1:k){
    det=rnorm(1,0,sqrt(rou)*sigma_y)
    u[i,]=rnorm(1,0,1)
    for(j in 1:m){
      x[i,j]=rnorm(1,0,1)
      e=rnorm(1,0,sqrt(1-tau^2-rou)*sigma_y)
      y[i,j]=a1+tau*sigma_y*x[i,j]+a3*trt+u[i,j]+e
      R[i,j]=expit(b0+b1*x[i,j]+u[i,j])
    }
  }
  z=b0+b1*x+u
  pr = 1/(1+exp(-z))         # pass through an inv-logit function
  ind = rbinom(k*m,1,pr) 
  #mis=ifelse(R>0.6,1,0)
  data_con=data.frame(y= as.vector((y)),x=as.vector((x)),u= as.vector((u)),ind=as.vector((ind)),
                      R=as.vector((R)),cluster=rep(c(1:k),m),group=trt)
  return(data_con)
}

for(rs in 1:50){
trt=datageneration(trt=1)
con=datageneration(trt=0)
sim_data=rbind(trt,con)
sim_data[sim_data$group==0,]$cluster=sim_data[sim_data$group==0,]$cluster+50
head(sim_data)

### 
sum(sim_data$ind==1)/(k*m*2)
# estimate the model and store results in m
logs <- glmer(ind ~ x + (1|cluster), data = sim_data, family = binomial)
logsum=summary(logs)
ran=as.data.frame(VarCorr(logs))['sdcor']
ran=as.numeric(ran)
rans=c()
for(i in 1:m){
  rans=c(rans,rnorm(1,0,sqrt(ran)))
}
x=matrix(sim_data$x,(2*k),m)
rans=rep(rans,(2*k))
rans=matrix(rans,2*k,m)
weight=expit(logsum$coefficients[,'Estimate'][1]+logsum$coefficients[,'Estimate'][2]*x+rans)
weight=matrix(weight,((2*k)*m),1)

sim_data$weight=weight
cca=sim_data[sim_data$ind==1,]
dim(cca)

ipw=lmer(y~x+group+(1|cluster),data=cca,weights = weight)
est_u=as.data.frame(VarCorr(ipw))['sdcor'][1,1]
est_trt=summary(ipw)$coefficients[,'Estimate']['group']
est_trt=as.vector(est_trt)
sd_trt=summary(ipw)$coefficients[,'Std. Error']['group']
sd_trt=as.vector(sd_trt)
ipw_mean=c(ipw_mean,est_trt)
ipw_std=c(ipw_std,sd_trt)

sim_data_na=sim_data
sim_data_na$y[sim_data_na$ind==0]=NA
#head(sim_data_na,30)
#formula<-as.formula(y~x+group+(1|cluster))
#mmi<-jomo.lmer(formula,data=sim_data_na, level=type, nburn=10, nbetween=10)
#mmi=jomo(sim_data_na,clus=cluster,nburn=nburn,nbetween=nbetween,nimp=nimp,meth='random')
type=c(1,3,0,0,0,-2,0,0)
names(type)=colnames(sim_data_na)

type2=c(1,3,0,0,0,-2,2,0)
names(type2)=colnames(sim_data_na)

mmi_pan=panImpute(sim_data_na,type=type, n.burn=1000, n.iter=100, m=5)
mmi_jomo=jomoImpute(sim_data_na,type=type, n.burn=1000, n.iter=100, m=5)

mmi_pan2=panImpute(sim_data_na,type=type2, n.burn=1000, n.iter=100, m=5)
mmi_jomo2=jomoImpute(sim_data_na,type=type2, n.burn=1000, n.iter=100, m=5)

mmipan_result=analysis_results_mmi(mmi_pan)
mmijomo_result=analysis_results_mmi(mmi_jomo)

mmipan_result2=analysis_results_mmi(mmi_pan2)
mmijomo_result2=analysis_results_mmi(mmi_jomo2)

mmi_m=mmipan_result$mean
mmi_s=mmipan_result$std
mmi_mean_pan=c(mmi_mean_pan,mmi_m)
mmi_std_pan=c(mmi_std_pan,mmi_s)

mmi_m=mmipan_result2$mean
mmi_s=mmipan_result2$std
mmi_mean_pan2=c(mmi_mean_pan2,mmi_m)
mmi_std_pan2=c(mmi_std_pan2,mmi_s)

mmi_m=mmijomo_result$mean
mmi_s=mmijomo_result$std
mmi_mean_jomo=c(mmi_mean_jomo,mmi_m)
mmi_std_jomo=c(mmi_std_jomo,mmi_s)

mmi_m=mmijomo_result2$mean
mmi_s=mmijomo_result2$std
mmi_mean_jomo2=c(mmi_mean_jomo2,mmi_m)
mmi_std_jomo2=c(mmi_std_jomo2,mmi_s)


### jomo without cluster 
ttt0=sim_data_na[,c('y','x')]
mmi1=jomo(ttt0)
mmi_jomo0=analysis_results_mmi2(mmi1)
mmi_m=mmi_jomo0$mean
mmi_s=mmi_jomo0$std
jomo_no_m=c(jomo_no_m,mmi_m)
jomo_no_std=c(jomo_no_std,mmi_s)

### jomo with cluster
ttt=sim_data_na[,c('y','x','cluster')]
cluster=ttt$cluster
mmi2=jomo(ttt,clus=cluster)
mmi_result=analysis_results_mmi3(mmi2)

mmi_m=mmi_result$mean
mmi_s=mmi_result$std
jomo_mean=c(jomo_mean,mmi_m)
jomo_std=c(jomo_std,mmi_s)


print('repeat time')
print(rs)
}


end_time=Sys.time()
end_time-start_time

ipw_mean
ipw_std
ci(ipw_mean,ipw_std)
jomo_mean
jomo_std
ci(jomo_mean,jomo_std)
jomo_no_m
jomo_no_std
ci(jomo_no_m,jomo_no_std)
mmi_mean_pan
mmi_std_pan
ci(mmi_mean_pan,mmi_std_pan)
mmi_mean_pan2
mmi_std_pan2
ci(mmi_mean_pan2,mmi_std_pan2)
mmi_mean_jomo
mmi_std_jomo
ci(mmi_mean_jomo,mmi_std_jomo)
mmi_mean_jomo2
mmi_std_jomo2
ci(mmi_mean_jomo2,mmi_std_jomo2)
### converge
## ci
ci=function(m,std){
  print(paste('mean',round(mean(m),2)))
  print(paste('std',round(mean(std),2)))
  ci_l=m-1.96*std
  ci_u=m+1.96*std
  print(paste('95% ci:',"(",round(mean(ci_l),2),",",round(mean(ci_u),2),")",sep=''))
  num1=ifelse(ci_l<5,1,0)
  num2=ifelse(ci_u>5,1,0)
  n=num1+num2
  nn=sum(n==2)
  print('converage percent')
  print(nn/length(m))
}
