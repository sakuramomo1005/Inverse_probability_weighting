
### Functions:
###****************** Data Generation ***********************###

# generation one arm
# if cluster size mean = 25, then cluster size is generated  from UNI(20,30)
# if cluster size mean = 50, then cluster size is generated  from UNI(40,60) 

# mis: the missing scenarios
# s: the data generation scenraio, s=1 or 3
# i: the indicator of intervention arm
# k: number of clusters in each intervention arm
# mm: mean of number of individuals in each cluster
one_group3=function(mis,s,i,k,mm,seed=123){
  set.seed(seed)
  if(s==1){b0=1;b1=1.36;b2=1;psi=-1.34}
  if(s==3){b0=1;b1=1.36;b2=0.588;psi=-1.34}
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
    sigma_icc=0.4
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

# combine two intervention arm 
data_gene3=function(k,m,s,mis,seed=123){
  set.seed(seed)
  # step1:
  a1=one_group3(i=1,k=k,m=m,seed=seed,mis=mis,s=s)
  a0=one_group3(i=0,k=k,m=m,seed=seed,mis=mis,s=s)
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

# function to report warnings and errors
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

# function to calculate the converage
conv=function(means,sds,true){
  n=0
  temp=data.frame(means=means,sds=sds)
  temp=na.omit(temp)
  means=temp$means
  sds=temp$sds
  for(i in 1:length(means)){
    if(means[i]-1.96*sds[i]<true & means[i]+1.96*sds[i]> true){n=n+1}
  }
  return(n)
}

# function to calculate the missing percentage
miss_per=function(a){
  b=sum(a$R)/dim(a)[1]
  return(b)
}

## function for MI
# function to calculate results
# mmi is the results generated from MI
# org_data is the data set used to do MI
# num is the number of imputations
result_mmi=function(mmi,org_data,num=5){
  m0=c()
  std0=c()
  m02=c()
  std02=c()
  #icc=c()
  for(i in 1:num){
    #print(i)
    temp=mmi[mmi$Imputation==i,]
    #temp=temp[-1,]
    temp$cluster=org_data$cluster
    temp$arm=org_data$arm
    rownames(temp)=NULL
    temp0=temp[,c('arm','x','cluster','y')]
    formu=formula(y1~x+arm)
    temp0$y1=as.numeric(temp0$y)-1  # I have to change the factors to numerics, 
    # since geese cannot handle factors
    # calculate gee with independent working correlation matrix 
    mmi2=myTryCatch(geese(formu,data=temp0,id=cluster,
                          family = binomial(link='logit'),
                          corstr = 'independence'))
    # calculate gee with independent working exchangeable matrix 
    mmi22=myTryCatch(geese(formu,data=temp0,id=cluster,
                           family = binomial(link='logit'),
                           corstr = 'exchangeable'))
    if(is.null(mmi2$value)+is.null(mmi22$value)>0){
      m0=c(m0,NA)
      std0=c(std0,NA)
      m02=c(m02,NA)
      std02=c(std02,NA)
    }
    if(is.null(mmi2$value)+is.null(mmi22$value)==0){
      mmi2=mmi2$value;mmi22=mmi22$value
      est=summary(mmi2)
      est_trt=est$mean['arm','estimate']
      sd_trt=est$mean['arm','san.se']
      m0=c(m0,est_trt)
      std0=c(std0,sd_trt)
      
      est2=summary(mmi22)
      est_trt2=est2$mean['arm','estimate']
      sd_trt2=est2$mean['arm','san.se']
      m02=c(m02,est_trt2)
      std02=c(std02,sd_trt2)
    }
  }
  return(list(m0=m0,std0=std0,m02=m02,std02=std02))
  ## m0, std0 are the results of exchangeable and m02, std02 are the results of independent
}
mypool=function(mean0,sd0,num=5,print='no'){
  m=mean(mean0,na.rm=TRUE)
  v=mean(sd0,na.rm=TRUE)
  B=sd(mean0,na.rm=TRUE)
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


###****************** Data Analysis ***********************###

# generating 1000,000 samples to check the max weight
# x=alpha+u, alpha is cluster level and u is individual level 
try_u=rnorm(1000000,0,sqrt(3.37))
try_a=rnorm(1000,0,sqrt(0.18))
try_a2=rep(try_a,1000)
try_x=try_u+try_a2
## the max weight
1/expit(-1.34+min(try_x))

# based on our scenarios, generate k=25, m=25 to check the max weight
check_x1=one_group3(mis=1,s=1,i=1,k=25,mm=25,seed=123)
mean(check_x1$x)
sd(check_x1$x)
range(check_x1$x)
hist(check_x1$x)

## max weight
1/expit(-1.34+min(check_x1$x))



###************************ Data Analysis ***********************###
###************************      IPW      ***********************###

times=200

## gee with independent working correlation matrix
# weight method: truncation 1
for(s in c(1:3)){for(mis in 1:6){
  
  # variables to report mis, s, time, k and m, missing percentage
  ipwmis=c();ipws=c();ipwtime=c();ipwkm=c();missper=c()
  
  # variables to report true effect
  true_est=c();true_warning=c()
  
  # variables to reports unadjusted cca
  uncra_est=c();uncra_std=c();uncra_warning=c()
  
  # variables to reports adjusted cca
  cra_est=c();cra_std=c();cra_warning=c()
  
  # variables to reports ipw with packages lme4 and geepack, without cluster effects, truncation 1
  ipwgeeno_est=c();ipwgeeno_std=c();ipwgeeno_warn=c()
  
  # variables to reports ipw with packages CRTgeeDR, without cluster effects, truncation 1 
  ipwcrtgeedr_est=c();ipwcrtgeedr_std=c();ipwcrtgee_warn=c()
  
  
  # variables to reports ipw with packages lme4 and geepack, with cluster effects, truncation 1
  ipwgeeyes_est=c();ipwgeeyes_std=c();ipwgeeyes_warn=c()
  
  # variables to reports ipw with packages CRTgeeDR, with cluster effects, truncation 1
  ipwcrtgeedr_yes_est=c();ipwcrtgeedr_yes_std=c();ipwcrtgee_yes_warn=c()
  
  
  for(km in 1:4){
    
    if(km==1){k=25;m=25}
    if(km==2){k=50;m=25}
    if(km==3){k=25;m=50}
    if(km==4){k=50;m=50}
    
    for(tt in 1:times){
      print("***************")
      print('mis')
      print(mis)
      print("km")
      print(km)
      print('times')
      
      print(tt)
      d1=data_gene3(s=s,k=k,m=m,seed=tt,mis=mis)
      d2=d1
      d2$y=ifelse(d2$R==1,NA,d2$y)
      d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=d2$R)
      
      if(mis==1 | mis==4){
        logs=glm(missing ~ x , data = d3,
                 family = binomial(link='logit'))
        logs2=glmer(missing ~ x+(1|cluster) , data = d3,
                    family = binomial(link='logit'))
      }
      if(mis==2 | mis==5){
        logs=glm(missing ~ x+arm, data = d3,
                 family = binomial(link='logit'))
        logs2=glmer(missing ~ x+arm+(1|cluster) , data = d3,
                    family = binomial(link='logit'))
      }
      if(mis==3 | mis==6){
        logs=glm(missing ~ x+arm+x:arm, data = d3,
                 family = binomial(link='logit'))
        logs2=glmer(missing ~ x+arm+x:arm+(1|cluster) , data = d3,
                    family = binomial(link='logit'))
      }
      
      weight=expit(predict(logs))
      weight2=expit(predict(logs2))
      
      d3_weight=ifelse((1/weight)>1000,1000,(1/weight))
      d3$weight2=ifelse(1/weight2>1000,1000,1/weight2)
      
      # those two round weights are special for warning catching. 
      ## I put a round here, since without round, geese will generate warnings. We only want warnings that show non convergence.   
      d3$weight_round=ifelse(round(1/weight)>1000,1000,round(1/weight))
      d3$weight_round2=ifelse(round(1/weight2)>1000,1000,round(1/weight2))
      
      d4=na.omit(d3)
      
      if(mis==1|mis==4){
        ipwcrtgeedr=myTryCatch(geeDREstimation(formula=y~x+arm,
                                               id="cluster" , data = d3,
                                               nameMISS='missing',nameY='y',
                                               nameTRT='arm',
                                               family =  binomial("logit"),
                                               corstr = "independence",
                                               model.weights=I(missing==1)~x))
      }
      if(mis==2|mis==5){
        ipwcrtgeedr=myTryCatch(geeDREstimation(formula=y~x+arm,
                                               id="cluster" , data = d3,
                                               nameMISS='missing',nameY='y',
                                               nameTRT='arm',
                                               family =  binomial("logit"),
                                               corstr = "independence",
                                               model.weights=I(missing==1)~x+arm))
      }
      if(mis==3| mis==6){
        ipwcrtgeedr=myTryCatch(geeDREstimation(formula=y~x+arm,
                                               id="cluster" , data = d3,
                                               nameMISS='missing',nameY='y',
                                               nameTRT='arm',
                                               family =  binomial("logit"),
                                               corstr = "independence",
                                               model.weights=I(missing==1)~x*arm))
      }
      
      
      ipwcrtgeedr_yes=myTryCatch(geeDREstimation(formula=y~x+arm,
                                                 id="cluster" , data = d3,
                                                 nameMISS='missing',nameY='y',
                                                 nameTRT='arm',
                                                 family =  binomial("logit"),
                                                 corstr = "independence",
                                                 weights = d3$weight2))
      
      if(is.null(ipwcrtgeedr$value)+is.null(ipwcrtgeedr0$value)+is.null(ipwcrtgeedr_yes$value)==0){
        ipwtime=c(ipwtime,tt)
        ipwkm=c(ipwkm,km)
        ipwmis=c(ipwmis,mis)
        ipws=c(ipws,s)
        missper=c(missper,miss_per(d1)) 
        
        #d3$y=as.factor(d3$y)
        ### True values
        formu=formula(y~x+arm)
        true=myTryCatch(geese(formu,data=d1,id=cluster,
                              family = binomial(link='logit'),
                              corstr = 'independence'))
        true_est=c(true_est,summary(true$value)$mean['arm','estimate'])
        if(is.null(true$warning)==0){
          true_warning=c(true_warning,tt)}
        if(is.null(true$warning)!=0){
          true_warning=c(true_warning,0)}
        
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
        if(is.null(uncra$warning)!=0){
          uncra_warning=c(uncra_warning,0)}
        
        ### adjusted CRA
        cra=myTryCatch(geese(formu,data=d_nona,id=cluster,
                             family = binomial(link='logit'),
                             corstr = 'independence'))
        cra_est=c(cra_est,summary(cra$value)$mean['arm','estimate'])
        cra_std=c(cra_std,summary(cra$value)$mean['arm','san.se'])
        if(is.null(cra$warning)==0){
          cra_warning=c(cra_warning,tt)}
        if(is.null(cra$warning)!=0){
          cra_warning=c(cra_warning,0)}
        
        ### IPW-no cluster effect, by using glm and geese
        ipwgeeno=myTryCatch(geese(formu,data=d4,id=cluster,
                                  family = binomial(link='logit'),
                                  weights = d4$weight,
                                  corstr = 'independence'))
        tempgeeno_est=summary(ipwgeeno$value)$mean['arm','estimate']
        tempgeeno_std=summary(ipwgeeno$value)$mean['arm','san.se']
        ipwgeeno_est=c(ipwgeeno_est,tempgeeno_est)
        ipwgeeno_std=c(ipwgeeno_std,tempgeeno_std)
        ipwgeeno=myTryCatch(geese(formu,data=d4,id=cluster,
                                  family = binomial(link='logit'),
                                  weights = d4$weight_round,
                                  corstr = 'independence'))
        if(is.null(ipwgeeno$warning)==0){ipwgeeno_warn=c(ipwgeeno_warn,tt)}
        if(is.null(ipwgeeno$warning)!=0){ipwgeeno_warn=c(ipwgeeno_warn,0)}
        
        ### IPW-GEE with package CRTgeeDR
        tempcrtgee_est=summary(ipwcrtgeedr$value)$beta[3]
        tempcrtgee_std=summary(ipwcrtgeedr$value)$se.robust[3]
        ipwcrtgeedr_est=c(ipwcrtgeedr_est,tempcrtgee_est)
        ipwcrtgeedr_std=c(ipwcrtgeedr_std,tempcrtgee_std)
        if(ipwcrtgeedr$value$converged==0){ipwcrtgee_warn=c(ipwcrtgee_warn,tt)}
        if(ipwcrtgeedr$value$converged==1){ipwcrtgee_warn=c(ipwcrtgee_warn,0)}
        
        ### IPW-WITH cluster effect, by using glm and geese
        formu=formula(y~x+arm)
        ipwgeeyes=myTryCatch(geese(formu,data=d4,id=cluster,
                                   family = binomial(link='logit'),
                                   weights = d4$weight2,
                                   corstr = 'independence'))
        tempgeeyes_est=summary(ipwgeeyes$value)$mean['arm','estimate']
        tempgeeyes_std=summary(ipwgeeyes$value)$mean['arm','san.se']
        ipwgeeyes_est=c(ipwgeeyes_est,tempgeeyes_est)
        ipwgeeyes_std=c(ipwgeeyes_std,tempgeeyes_std)
        ipwgeeyes=myTryCatch(geese(formu,data=d4,id=cluster,
                                   family = binomial(link='logit'),
                                   weights = d4$weight_round2,
                                   corstr = 'independence'))
        if(is.null(ipwgeeyes$warning)==0){ipwgeeyes_warn=c(ipwgeeyes_warn,tt)}
        if(is.null(ipwgeeyes$warning)!=0){ipwgeeyes_warn=c(ipwgeeyes_warn,0)}
        
        ### IPW-GEE with cluster effect with package CRTgeeDR
        
        tempcrtgee_yes_est=summary(ipwcrtgeedr_yes$value)$beta[3]
        tempcrtgee_yes_std=summary(ipwcrtgeedr_yes$value)$se.robust[3]
        ipwcrtgeedr_yes_est=c(ipwcrtgeedr_yes_est,tempcrtgee_yes_est)
        ipwcrtgeedr_yes_std=c(ipwcrtgeedr_yes_std,tempcrtgee_yes_std)
        if(ipwcrtgeedr_yes$value$converged==0){ipwcrtgee_yes_warn=c(ipwcrtgee_yes_warn,tt)}
        if(ipwcrtgeedr_yes$value$converged==1){ipwcrtgee_yes_warn=c(ipwcrtgee_yes_warn,0)}
      }
    }}
  
  data_mean=data.frame(uncra_est=uncra_est,
                       cra_est=cra_est,
                       ipwgeeno_est=ipwgeeno_est,
                       ipwcrtgeedr_est=ipwcrtgeedr_est,
                       ipwgeeyes_est=ipwgeeyes_est,
                       ipwcrtgeedr_yes_est=ipwcrtgeedr_yes_est,
                       ipwtime=ipwtime,
                       ipwk=ipwkm,true_est=true_est,missper=missper)
  data_sd=data.frame(uncra_std=uncra_std,
                     cra_std=cra_std,
                     ipwgeeno_std=ipwgeeno_std,
                     ipwcrtgeedr_std=ipwcrtgeedr_std,
                     ipwgeeyes_std=ipwgeeyes_std,
                     ipwcrtgeedr_yes_std=ipwcrtgeedr_yes_std,
                     ipwtime=ipwtime,
                     ipwk=ipwkm)
  data_warn=data.frame(uncra_warn=uncra_warning,
                       cra_warn=cra_warning,
                       ipwgeeno_warn=ipwgeeno_warn,
                       ipwcrtgeedr_warn=ipwcrtgee_warn,
                       ipwgeeyes_warn=ipwgeeyes_warn,
                       ipwcrtgeedr_yes_warn=ipwcrtgee_yes_warn,
                       ipwtime=ipwtime,
                       ipwk=ipwkm)
  
  fff=list(data_mean=data_mean,data_sd=data_sd,data_warn=data_warn)
  file_name=paste('ipws3m1ind',s,mis,icc*100,'.RData',sep='')
  save(fff,file=file_name)
}}

## the same with previous but exchangeable working correlation matrix
# weight method: truncation 1
for(s in c(1:3)){
  for(mis in 1:6){
    
    ipwmis=c(); ipws=c() ;ipwtime=c() ; ipwkm=c() ;missper=c()
    
    true_est=c();true_warning=c()
    
    uncra_est=c();uncra_std=c();uncra_warning=c()
    
    cra_est=c();cra_std=c();cra_warning=c()
    
    ipwgeeno_est=c();ipwgeeno_std=c();ipwgeeno_warn=c()
    
    ipwcrtgeedr_est=c();ipwcrtgeedr_std=c();ipwcrtgee_warn=c()
    
    ipwgeeyes_est=c();ipwgeeyes_std=c();ipwgeeyes_warn=c()
    
    ipwcrtgeedr_yes_est=c();ipwcrtgeedr_yes_std=c();ipwcrtgee_yes_warn=c()
    
    
    for(km in 1:4){
      if(km==1){k=25;m=25}
      if(km==2){k=50;m=25}
      if(km==3){k=25;m=50}
      if(km==4){k=50;m=50}
      
      for(tt in 1:times){
        print("***************")
        print('mis')
        print(mis)
        print("km")
        print(km)
        print('times')
        print(tt)
        d1=data_gene3(s=s,k=k,m=m,seed=tt,mis=mis)
        d2=d1
        d2$y=ifelse(d2$R==1,NA,d2$y)
        d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=d2$R)
        #d3$y=as.factor(d3$y)
        ### True values
        if(mis==1 | mis==4){
          logs=glm(missing ~ x , data = d3,
                   family = binomial(link='logit'))
          logs2=glmer(missing ~ x+(1|cluster) , data = d3,
                      family = binomial(link='logit'))
        }
        if(mis==2 | mis==5){
          logs=glm(missing ~ x+arm, data = d3,
                   family = binomial(link='logit'))
          logs2=glmer(missing ~ x+arm+(1|cluster) , data = d3,
                      family = binomial(link='logit'))
        }
        if(mis==3 | mis==6){
          logs=glm(missing ~ x+arm+x:arm, data = d3,
                   family = binomial(link='logit'))
          logs2=glmer(missing ~ x+arm+x:arm+(1|cluster) , data = d3,
                      family = binomial(link='logit'))
        }
        
        weight=expit(predict(logs))
        weight2=expit(predict(logs2))
        
        d3$weight=ifelse((1/weight)>1000,1000,(1/weight))
        d3$weight2=ifelse(1/weight2>1000,1000,1/weight2)
        
        d3$weight_round=ifelse(round(1/weight)>1000,1000,round(1/weight))
        d3$weight_round2=ifelse(round(1/weight2)>1000,1000,round(1/weight2))
        
        d4=na.omit(d3)
        
        ipwcrtgeedr=myTryCatch(geeDREstimation(formula=y~x+arm,
                                               id="cluster" , data = d3,
                                               nameMISS='missing',nameY='y',
                                               nameTRT='arm',
                                               weights = d3$weight,
                                               family =  binomial("logit"),
                                               corstr = "exchangeable"))
        
        ipwcrtgeedr_yes=myTryCatch(geeDREstimation(formula=y~x+arm,
                                                   id="cluster" , data = d3,
                                                   nameMISS='missing',nameY='y',
                                                   nameTRT='arm',
                                                   family =  binomial("logit"),
                                                   corstr = "exchangeable",
                                                   weights = d3$weight2))
        
        if(is.null(ipwcrtgeedr$value)+is.null(ipwcrtgeedr_yes$value)==0){
          formu=formula(y~x+arm)
          ipwtime=c(ipwtime,tt)
          ipwkm=c(ipwkm,km)
          ipwmis=c(ipwmis,mis)
          ipws=c(ipws,s)
          missper=c(missper,miss_per(d1)) 
          true=myTryCatch(geese(formu,data=d1,id=cluster,
                                family = binomial(link='logit'),
                                corstr = 'exchangeable'))
          true_est=c(true_est,summary(true$value)$mean['arm','estimate'])
          if(is.null(true$warning)==0){
            true_warning=c(true_warning,tt)}
          if(is.null(true$warning)!=0){
            true_warning=c(true_warning,0)}
          ### unadjusted CRA
          d_nona=na.omit(d1)
          formu0=formula(y~arm)
          uncra=myTryCatch(geese(formu0,data=d_nona,id=cluster,
                                 family = binomial(link='logit'),
                                 corstr = 'exchangeable'))
          uncra_est=c(uncra_est,summary(uncra$value)$mean['arm','estimate'])
          uncra_std=c(uncra_std,summary(uncra$value)$mean['arm','san.se'])
          if(is.null(uncra$warning)==0){
            uncra_warning=c(uncra_warning,tt)}
          if(is.null(uncra$warning)!=0){
            uncra_warning=c(uncra_warning,0)}
          ### adjusted CRA
          cra=myTryCatch(geese(formu,data=d_nona,id=cluster,
                               family = binomial(link='logit'),
                               corstr = 'exchangeable'))
          cra_est=c(cra_est,summary(cra$value)$mean['arm','estimate'])
          cra_std=c(cra_std,summary(cra$value)$mean['arm','san.se'])
          if(is.null(cra$warning)==0){
            cra_warning=c(cra_warning,tt)}
          if(is.null(cra$warning)!=0){
            cra_warning=c(cra_warning,0)}
          ### IPW-no cluster effect, by using glm and geese
          
          ipwgeeno=myTryCatch(geese(formu,data=d4,id=cluster,
                                    family = binomial(link='logit'),
                                    weights = d4$weight,
                                    corstr = 'exchangeable'))
          tempgeeno_est=summary(ipwgeeno$value)$mean['arm','estimate']
          tempgeeno_std=summary(ipwgeeno$value)$mean['arm','san.se']
          ipwgeeno_est=c(ipwgeeno_est,tempgeeno_est)
          ipwgeeno_std=c(ipwgeeno_std,tempgeeno_std)
          ipwgeeno=myTryCatch(geese(formu,data=d4,id=cluster,
                                    family = binomial(link='logit'),
                                    weights = d4$weight_round,
                                    corstr = 'exchangeable'))
          if(is.null(ipwgeeno$warning)==0){ipwgeeno_warn=c(ipwgeeno_warn,tt)}
          if(is.null(ipwgeeno$warning)!=0){ipwgeeno_warn=c(ipwgeeno_warn,0)}
          
          ### IPW-GEE with package CRTgeeDR
          
          tempcrtgee_est=summary(ipwcrtgeedr$value)$beta[3]
          tempcrtgee_std=summary(ipwcrtgeedr$value)$se.robust[3]
          ipwcrtgeedr_est=c(ipwcrtgeedr_est,tempcrtgee_est)
          ipwcrtgeedr_std=c(ipwcrtgeedr_std,tempcrtgee_std)
          if(ipwcrtgeedr$value$converged==0){ipwcrtgee_warn=c(ipwcrtgee_warn,tt)}
          if(ipwcrtgeedr$value$converged==1){ipwcrtgee_warn=c(ipwcrtgee_warn,0)}
          
          ### IPW-WITH cluster effect, by using glm and geese
          formu=formula(y~x+arm)
          ipwgeeyes=myTryCatch(geese(formu,data=d4,id=cluster,
                                     family = binomial(link='logit'),
                                     weights = d4$weight2,
                                     corstr = 'exchangeable'))
          tempgeeyes_est=summary(ipwgeeyes$value)$mean['arm','estimate']
          tempgeeyes_std=summary(ipwgeeyes$value)$mean['arm','san.se']
          ipwgeeyes_est=c(ipwgeeyes_est,tempgeeyes_est)
          ipwgeeyes_std=c(ipwgeeyes_std,tempgeeyes_std)
          ipwgeeyes=myTryCatch(geese(formu,data=d4,id=cluster,
                                     family = binomial(link='logit'),
                                     weights = d4$weight_round2,
                                     corstr = 'exchangeable'))
          
          if(is.null(ipwgeeyes$warning)==0){ipwgeeyes_warn=c(ipwgeeyes_warn,tt)}
          if(is.null(ipwgeeyes$warning)!=0){ipwgeeyes_warn=c(ipwgeeyes_warn,0)}
          
          ### IPW-GEE with cluster effect with package CRTgeeDR
          
          tempcrtgee_yes_est=summary(ipwcrtgeedr_yes$value)$beta[3]
          tempcrtgee_yes_std=summary(ipwcrtgeedr_yes$value)$se.robust[3]
          ipwcrtgeedr_yes_est=c(ipwcrtgeedr_yes_est,tempcrtgee_yes_est)
          ipwcrtgeedr_yes_std=c(ipwcrtgeedr_yes_std,tempcrtgee_yes_std)
          if(ipwcrtgeedr_yes$value$converged==0){ipwcrtgee_yes_warn=c(ipwcrtgee_yes_warn,tt)}
          if(ipwcrtgeedr_yes$value$converged==1){ipwcrtgee_yes_warn=c(ipwcrtgee_yes_warn,0)}
        }}
    }
    data_mean=data.frame(uncra_est=uncra_est,
                         cra_est=cra_est,
                         ipwgeeno_est=ipwgeeno_est,
                         ipwcrtgeedr_est=ipwcrtgeedr_est,
                         ipwgeeyes_est=ipwgeeyes_est,
                         ipwcrtgeedr_yes_est=ipwcrtgeedr_yes_est,
                         ipwtime=ipwtime,
                         ipwk=ipwkm,true_est=true_est,missper=missper)
    data_sd=data.frame(uncra_std=uncra_std,
                       cra_std=cra_std,
                       ipwgeeno_std=ipwgeeno_std,
                       ipwcrtgeedr_std=ipwcrtgeedr_std,
                       ipwgeeyes_std=ipwgeeyes_std,
                       ipwcrtgeedr_yes_std=ipwcrtgeedr_yes_std,
                       ipwtime=ipwtime,
                       ipwk=ipwkm,true_est=true_est)
    data_warn=data.frame(uncra_warn=uncra_warning,
                         cra_warn=cra_warning,
                         ipwgeeno_warn=ipwgeeno_warn,
                         ipwcrtgeedr_warn=ipwcrtgee_warn,
                         ipwgeeyes_warn=ipwgeeyes_warn,
                         ipwcrtgeedr_yes_warn=ipwcrtgee_yes_warn,
                         ipwtime=ipwtime,
                         ipwk=ipwkm,true_warning=true_warning)
    fff=list(data_mean=data_mean,data_sd=data_sd,data_warn=data_warn)
    file_name=paste('ipws3m1_ex',s,mis,icc*100,'.RData',sep='')
    save(fff,file=file_name)
  }
}

## gee with independent working correlation matrix
# weight method: truncation 2 and stabilization
for(s in c(1,3)){
  for(mis in 1:6){
    
    # variables to report mis, s, time, k and m, missing percentage
    ipwmis=c();ipws=c();ipwtime=c();ipwkm=c()
    missper=c()
    
    # variables to report true effect
    true_est=c();true_warning=c()
    
    # variables to reports unadjusted cca
    uncra_est=c();uncra_std=c();uncra_warning=c()
    
    # variables to reports adjusted cca
    cra_est=c();cra_std=c();cra_warning=c()
    
    #### weight method: truncation, wihtin 5th and 95th percentiles.
    # variables to reports ipw with packages lme4 and geepack, without cluster effects, truncation
    ipwgeeno_est=c();ipwgeeno_std=c();ipwgeeno_warn=c()
    
    # variables to reports ipw with packages CRTgeeDR, without cluster effects, truncation 
    ipwcrtgeedr_est=c();ipwcrtgeedr_std=c();ipwcrtgee_warn=c()
    
    # variables to reports ipw with packages lme4 and geepack, with cluster effects, truncation 
    ipwgeeyes_est=c();ipwgeeyes_std=c();ipwgeeyes_warn=c()
    
    # variables to reports ipw with packages CRTgeeDR, with cluster effects, truncation 
    ipwcrtgeedr_yes_est=c();ipwcrtgeedr_yes_std=c();ipwcrtgee_yes_warn=c()
    
    #### weight method: stabilization
    # variables to reports ipw with packages lme4 and geepack, without cluster effects, stabilization
    ipwgeeno_ests=c();ipwgeeno_stds=c();ipwgeeno_warns=c()
    
    # variables to reports ipw with packages CRTgeeDR, without cluster effects, stabilization 
    ipwcrtgeedr_ests=c();ipwcrtgeedr_stds=c();ipwcrtgee_warns=c()
    
    # variables to reports ipw with packages lme4 and geepack, with cluster effects, stabilization 
    ipwgeeyes_ests=c();ipwgeeyes_stds=c();ipwgeeyes_warns=c()
    
    # variables to reports ipw with packages CRTgeeDR, with cluster effects, stabilization
    ipwcrtgeedr_yes_ests=c();ipwcrtgeedr_yes_stds=c();ipwcrtgee_yes_warns=c()
    
    for(km in 1:4){
      if(km==1){k=25;m=25}
      if(km==2){k=50;m=25}
      if(km==3){k=25;m=50}
      if(km==4){k=50;m=50}
      
      for(tt in 1:times){
        print("***************")
        print('s')
        print(s)
        print('mis')
        print(mis)
        print("km")
        print(km)
        print('times')
        print(tt)
        
        d1=data_gene3(s=s,k=k,m=m,seed=tt,mis=mis)
        d2=d1
        d2$y=ifelse(d2$R==1,NA,d2$y)
        d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=d2$R)
        #d3$y=as.factor(d3$y)
        ### True values
        if(mis==1 | mis==4){
          logs=glm(missing ~ x , data = d3,
                   family = binomial(link='logit'))
          logs2=glmer(missing ~ x+(1|cluster) , data = d3,
                      family = binomial(link='logit'))
        }
        if(mis==2 | mis==5){
          logs=glm(missing ~ x+arm, data = d3,
                   family = binomial(link='logit'))
          logs2=glmer(missing ~ x+arm+(1|cluster) , data = d3,
                      family = binomial(link='logit'))
        }
        if(mis==3 | mis==6){
          logs=glm(missing ~ x+arm+x:arm, data = d3,
                   family = binomial(link='logit'))
          logs2=glmer(missing ~ x+arm+x:arm+(1|cluster) , data = d3,
                      family = binomial(link='logit'))
        }
        
        ###### calculation of weights ######
        # truncation
        weight=1/expit(predict(logs))
        weight2=1/expit(predict(logs2))
        
        tsw = ifelse(weight < quantile(weight, probs=.05), quantile(weight, probs=.05), weight)
        tsw = ifelse(weight > quantile(weight, probs=.95), quantile(weight, probs=.95), tsw)
        
        tsw2 = ifelse(weight2 < quantile(weight2, probs=.05), quantile(weight2, probs=.05), weight2)
        tsw2 = ifelse(weight2 > quantile(weight2, probs=.95), quantile(weight2, probs=.95), tsw2)
        
        d3$weight=tsw
        d3$weight2=tsw2
        # those two round weights are special for warning catching. 
        ## I put a round here, since without round, geese will generate warnings. We only want warnings that show non convergence.   
        d3$weight_round=round(tsw)
        d3$weight_round2=round(tsw2)
        
        # stabilization
        # the numerator for stablilized weights
        num0 = predict(glm(missing ~ 1,data = d3,
                           family = binomial(link='logit')),type="response")
        # the propensity score
        ps=expit(predict(logs))
        ps2=expit(predict(logs2))
        
        sw=ifelse(d3$missing==1, num0/ps, (1-num0)/(1-ps))
        sw2=ifelse(d3$missing==1, num0/ps2, (1-num0)/(1-ps2))
        
        d3$sw=sw
        d3$sw2=sw2
        d3$sw_round=round(sw)
        d3$sw_round2=round(sw2)
        
        d4=na.omit(d3)
        # truncation
        ipwcrtgeedr=myTryCatch(geeDREstimation(formula=y~x+arm,
                                               id="cluster" , data = d3,
                                               nameMISS='missing',nameY='y',
                                               nameTRT='arm',
                                               weights = d3$weight,
                                               family =  binomial("logit"),
                                               corstr = "independence"))
        # truncation-for warning catching
        ipwcrtgeedr_yes=myTryCatch(geeDREstimation(formula=y~x+arm,
                                                   id="cluster" , data = d3,
                                                   nameMISS='missing',nameY='y',
                                                   nameTRT='arm',
                                                   family =  binomial("logit"),
                                                   corstr = "independence",
                                                   weights = d3$weight2))
        
        # stabilization
        ipwcrtgeedrs=myTryCatch(geeDREstimation(formula=y~x+arm,
                                                id="cluster" , data = d3,
                                                nameMISS='missing',nameY='y',
                                                nameTRT='arm',
                                                weights = d3$sw,
                                                family =  binomial("logit"),
                                                corstr = "independence"))
        
        # stabilization-for warning catching
        ipwcrtgeedr_yess=myTryCatch(geeDREstimation(formula=y~x+arm,
                                                    id="cluster" , data = d3,
                                                    nameMISS='missing',nameY='y',
                                                    nameTRT='arm',
                                                    family =  binomial("logit"),
                                                    corstr = "independence",
                                                    weights = d3$sw2))
        
        if(is.null(ipwcrtgeedr$value)+is.null(ipwcrtgeedr_yes$value)+
           is.null(ipwcrtgeedrs$value)+is.null(ipwcrtgeedr_yess$value)==0){
          ipwtime=c(ipwtime,tt)
          ipwkm=c(ipwkm,km)
          ipwmis=c(ipwmis,mis)
          ipws=c(ipws,s)
          missper=c(missper,miss_per(d1)) 
          
          formu=formula(y~x+arm)
          true=myTryCatch(geese(formu,data=d1,id=cluster,
                                family = binomial(link='logit'),
                                corstr = 'independence'))
          true_est=c(true_est,summary(true$value)$mean['arm','estimate'])
          if(is.null(true$warning)==0){
            true_warning=c(true_warning,tt)}
          if(is.null(true$warning)!=0){
            true_warning=c(true_warning,0)}
          
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
          if(is.null(uncra$warning)!=0){
            uncra_warning=c(uncra_warning,0)}
          
          ### adjusted CRA
          cra=myTryCatch(geese(formu,data=d_nona,id=cluster,
                               family = binomial(link='logit'),
                               corstr = 'independence'))
          cra_est=c(cra_est,summary(cra$value)$mean['arm','estimate'])
          cra_std=c(cra_std,summary(cra$value)$mean['arm','san.se'])
          if(is.null(cra$warning)==0){
            cra_warning=c(cra_warning,tt)}
          if(is.null(cra$warning)!=0){
            cra_warning=c(cra_warning,0)}
          
          ### IPW-no cluster effect, by using glm and geese
          ipwgeeno=myTryCatch(geese(formu,data=d4,id=cluster,
                                    family = binomial(link='logit'),
                                    weights = d4$weight,
                                    corstr = 'independence'))
          tempgeeno_est=summary(ipwgeeno$value)$mean['arm','estimate']
          tempgeeno_std=summary(ipwgeeno$value)$mean['arm','san.se']
          ipwgeeno_est=c(ipwgeeno_est,tempgeeno_est)
          ipwgeeno_std=c(ipwgeeno_std,tempgeeno_std)
          ipwgeeno=myTryCatch(geese(formu,data=d4,id=cluster,
                                    family = binomial(link='logit'),
                                    weights = d4$weight_round,
                                    corstr = 'independence'))
          if(is.null(ipwgeeno$warning)==0){ipwgeeno_warn=c(ipwgeeno_warn,tt)}
          if(is.null(ipwgeeno$warning)!=0){ipwgeeno_warn=c(ipwgeeno_warn,0)}
          
          ### IPW-GEE with package CRTgeeDR
          tempcrtgee_est=summary(ipwcrtgeedr$value)$beta[3]
          tempcrtgee_std=summary(ipwcrtgeedr$value)$se.robust[3]
          ipwcrtgeedr_est=c(ipwcrtgeedr_est,tempcrtgee_est)
          ipwcrtgeedr_std=c(ipwcrtgeedr_std,tempcrtgee_std)
          if(ipwcrtgeedr$value$converged==0){ipwcrtgee_warn=c(ipwcrtgee_warn,tt)}
          if(ipwcrtgeedr$value$converged==1){ipwcrtgee_warn=c(ipwcrtgee_warn,0)}
          
          ### IPW-WITH cluster effect, by using glm and geese
          formu=formula(y~x+arm)
          ipwgeeyes=myTryCatch(geese(formu,data=d4,id=cluster,
                                     family = binomial(link='logit'),
                                     weights = d4$weight2,
                                     corstr = 'independence'))
          tempgeeyes_est=summary(ipwgeeyes$value)$mean['arm','estimate']
          tempgeeyes_std=summary(ipwgeeyes$value)$mean['arm','san.se']
          ipwgeeyes_est=c(ipwgeeyes_est,tempgeeyes_est)
          ipwgeeyes_std=c(ipwgeeyes_std,tempgeeyes_std)
          ipwgeeyes=myTryCatch(geese(formu,data=d4,id=cluster,
                                     family = binomial(link='logit'),
                                     weights = d4$weight_round2,
                                     corstr = 'independence'))
          
          if(is.null(ipwgeeyes$warning)==0){ipwgeeyes_warn=c(ipwgeeyes_warn,tt)}
          if(is.null(ipwgeeyes$warning)!=0){ipwgeeyes_warn=c(ipwgeeyes_warn,0)}
          
          ### IPW-GEE with cluster effect with package CRTgeeDR
          tempcrtgee_yes_est=summary(ipwcrtgeedr_yes$value)$beta[3]
          tempcrtgee_yes_std=summary(ipwcrtgeedr_yes$value)$se.robust[3]
          ipwcrtgeedr_yes_est=c(ipwcrtgeedr_yes_est,tempcrtgee_yes_est)
          ipwcrtgeedr_yes_std=c(ipwcrtgeedr_yes_std,tempcrtgee_yes_std)
          if(ipwcrtgeedr_yes$value$converged==0){ipwcrtgee_yes_warn=c(ipwcrtgee_yes_warn,tt)}
          if(ipwcrtgeedr_yes$value$converged==1){ipwcrtgee_yes_warn=c(ipwcrtgee_yes_warn,0)}
          
          
          
          ### stabilization
          ipwgeenos=myTryCatch(geese(formu,data=d4,id=cluster,
                                     family = binomial(link='logit'),
                                     weights = d4$sw,
                                     corstr = 'independence'))
          tempgeeno_est=summary(ipwgeenos$value)$mean['arm','estimate']
          tempgeeno_std=summary(ipwgeenos$value)$mean['arm','san.se']
          ipwgeeno_ests=c(ipwgeeno_ests,tempgeeno_est)
          ipwgeeno_stds=c(ipwgeeno_stds,tempgeeno_std)
          ipwgeenos=myTryCatch(geese(formu,data=d4,id=cluster,
                                     family = binomial(link='logit'),
                                     weights = d4$sw_round,
                                     corstr = 'independence'))
          if(is.null(ipwgeenos$warning)==0){ipwgeeno_warns=c(ipwgeeno_warns,tt)}
          if(is.null(ipwgeenos$warning)!=0){ipwgeeno_warns=c(ipwgeeno_warns,0)}
          
          ### IPW-GEE with package CRTgeeDR
          
          tempcrtgee_est=summary(ipwcrtgeedrs$value)$beta[3]
          tempcrtgee_std=summary(ipwcrtgeedrs$value)$se.robust[3]
          ipwcrtgeedr_ests=c(ipwcrtgeedr_ests,tempcrtgee_est)
          ipwcrtgeedr_stds=c(ipwcrtgeedr_stds,tempcrtgee_std)
          if(ipwcrtgeedrs$value$converged==0){ipwcrtgee_warns=c(ipwcrtgee_warns,tt)}
          if(ipwcrtgeedrs$value$converged==1){ipwcrtgee_warns=c(ipwcrtgee_warns,0)}
          
          ### IPW-WITH cluster effect, by using glm and geese
          formu=formula(y~x+arm)
          ipwgeeyess=myTryCatch(geese(formu,data=d4,id=cluster,
                                      family = binomial(link='logit'),
                                      weights = d4$sw2,
                                      corstr = 'independence'))
          tempgeeyes_est=summary(ipwgeeyess$value)$mean['arm','estimate']
          tempgeeyes_std=summary(ipwgeeyess$value)$mean['arm','san.se']
          ipwgeeyes_ests=c(ipwgeeyes_ests,tempgeeyes_est)
          ipwgeeyes_stds=c(ipwgeeyes_stds,tempgeeyes_std)
          ipwgeeyess=myTryCatch(geese(formu,data=d4,id=cluster,
                                      family = binomial(link='logit'),
                                      weights = d4$sw_round2,
                                      corstr = 'independence'))
          
          if(is.null(ipwgeeyess$warning)==0){ipwgeeyes_warns=c(ipwgeeyes_warns,tt)}
          if(is.null(ipwgeeyess$warning)!=0){ipwgeeyes_warns=c(ipwgeeyes_warns,0)}
          
          ### IPW-GEE with cluster effect with package CRTgeeDR
          
          tempcrtgee_yes_est=summary(ipwcrtgeedr_yess$value)$beta[3]
          tempcrtgee_yes_std=summary(ipwcrtgeedr_yess$value)$se.robust[3]
          ipwcrtgeedr_yes_ests=c(ipwcrtgeedr_yes_ests,tempcrtgee_yes_est)
          ipwcrtgeedr_yes_stds=c(ipwcrtgeedr_yes_stds,tempcrtgee_yes_std)
          if(ipwcrtgeedr_yess$value$converged==0){ipwcrtgee_yes_warns=c(ipwcrtgee_yes_warns,tt)}
          if(ipwcrtgeedr_yess$value$converged==1){ipwcrtgee_yes_warns=c(ipwcrtgee_yes_warns,0)}
        }}
    }
    data_mean=data.frame(uncra_est=uncra_est,
                         cra_est=cra_est,
                         ipwgeeno_est=ipwgeeno_est,
                         ipwcrtgeedr_est=ipwcrtgeedr_est,
                         ipwgeeyes_est=ipwgeeyes_est,
                         ipwcrtgeedr_yes_est=ipwcrtgeedr_yes_est,
                         ipwgeeno_ests=ipwgeeno_ests,
                         ipwcrtgeedr_ests=ipwcrtgeedr_ests,
                         ipwgeeyes_ests=ipwgeeyes_ests,
                         ipwcrtgeedr_yes_ests=ipwcrtgeedr_yes_ests,
                         ipwtime=ipwtime,
                         ipwk=ipwkm,true_est=true_est,missper=missper)
    data_sd=data.frame(uncra_std=uncra_std,
                       cra_std=cra_std,
                       ipwgeeno_std=ipwgeeno_std,
                       ipwcrtgeedr_std=ipwcrtgeedr_std,
                       ipwgeeyes_std=ipwgeeyes_std,
                       ipwcrtgeedr_yes_std=ipwcrtgeedr_yes_std,
                       ipwgeeno_stds=ipwgeeno_stds,
                       ipwcrtgeedr_stds=ipwcrtgeedr_stds,
                       ipwgeeyes_stds=ipwgeeyes_stds,
                       ipwcrtgeedr_yes_stds=ipwcrtgeedr_yes_stds,
                       ipwtime=ipwtime,
                       ipwk=ipwkm,true_est=true_est)
    data_warn=data.frame(uncra_warn=uncra_warning,
                         cra_warn=cra_warning,
                         ipwgeeno_warn=ipwgeeno_warn,
                         ipwcrtgeedr_warn=ipwcrtgee_warn,
                         ipwgeeyes_warn=ipwgeeyes_warn,
                         ipwcrtgeedr_yes_warn=ipwcrtgee_yes_warn,
                         ipwgeeno_warns=ipwgeeno_warns,
                         ipwcrtgeedr_warns=ipwcrtgee_warns,
                         ipwgeeyes_warns=ipwgeeyes_warns,
                         ipwcrtgeedr_yes_warns=ipwcrtgee_yes_warns,
                         ipwtime=ipwtime,
                         ipwk=ipwkm,true_warning=true_warning)
    fff=list(data_mean=data_mean,data_sd=data_sd,data_warn=data_warn)
    file_name=paste('ipws3m3ind',s,mis,icc*100,'.RData',sep='')
    save(fff,file=file_name)
  }
}

## the same with previous but exchangeable working correlation matrix
# weight method: truncation 2 and stabilization
for(s in c(1,3)){
  for(mis in 1:6){
    
    ipwmis=c();ipws=c();ipwtime=c();ipwkm=c();missper=c()
    
    true_est=c();true_warning=c()
    
    uncra_est=c();uncra_std=c();uncra_warning=c()
    cra_est=c();cra_std=c();cra_warning=c()
    
    ipwgeeno_est=c();ipwgeeno_std=c();ipwgeeno_warn=c()
    ipwcrtgeedr_est=c();ipwcrtgeedr_std=c();ipwcrtgee_warn=c()
    ipwgeeyes_est=c();ipwgeeyes_std=c();ipwgeeyes_warn=c()
    ipwcrtgeedr_yes_est=c();ipwcrtgeedr_yes_std=c();ipwcrtgee_yes_warn=c()
    
    ipwgeeno_ests=c();ipwgeeno_stds=c();ipwgeeno_warns=c()
    ipwcrtgeedr_ests=c();ipwcrtgeedr_stds=c();ipwcrtgee_warns=c()
    ipwgeeyes_ests=c();ipwgeeyes_stds=c();ipwgeeyes_warns=c()
    ipwcrtgeedr_yes_ests=c();ipwcrtgeedr_yes_stds=c();ipwcrtgee_yes_warns=c()
    
    for(km in 1:4){
      if(km==1){k=25;m=25}
      if(km==2){k=50;m=25}
      if(km==3){k=25;m=50}
      if(km==4){k=50;m=50}
      
      for(tt in 1:times){
        print("***************")
        print('s')
        print(s)
        print('mis')
        print(mis)
        print("km")
        print(km)
        print('times')
        print(tt)
        
        d1=data_gene3(s=s,k=k,m=m,seed=tt,mis=mis)
        d2=d1
        d2$y=ifelse(d2$R==1,NA,d2$y)
        d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=d2$R)
        #d3$y=as.factor(d3$y)
        ### True values
        if(mis==1 | mis==4){
          logs=glm(missing ~ x , data = d3,
                   family = binomial(link='logit'))
          logs2=glmer(missing ~ x+(1|cluster) , data = d3,
                      family = binomial(link='logit'))
        }
        if(mis==2 | mis==5){
          logs=glm(missing ~ x+arm, data = d3,
                   family = binomial(link='logit'))
          logs2=glmer(missing ~ x+arm+(1|cluster) , data = d3,
                      family = binomial(link='logit'))
        }
        if(mis==3 | mis==6){
          logs=glm(missing ~ x+arm+x:arm, data = d3,
                   family = binomial(link='logit'))
          logs2=glmer(missing ~ x+arm+x:arm+(1|cluster) , data = d3,
                      family = binomial(link='logit'))
        }
        
        # truncation
        weight=1/expit(predict(logs))
        weight2=1/expit(predict(logs2))
        
        tsw = ifelse(weight < quantile(weight, probs=.05), quantile(weight, probs=.05), weight)
        tsw = ifelse(weight > quantile(weight, probs=.95), quantile(weight, probs=.95), tsw)
        
        tsw2 = ifelse(weight2 < quantile(weight2, probs=.05), quantile(weight2, probs=.05), weight2)
        tsw2 = ifelse(weight2 > quantile(weight2, probs=.95), quantile(weight2, probs=.95), tsw2)
        
        d3$weight=tsw
        d3$weight2=tsw2
        d3$weight_round=round(tsw)
        d3$weight_round2=round(tsw2)
        
        # stabilization
        num0 = predict(glm(missing ~ 1,data = d3,
                           family = binomial(link='logit')),type="response")
        ps=expit(predict(logs))
        ps2=expit(predict(logs2))
        
        sw=ifelse(d3$missing==1, num0/ps, (1-num0)/(1-ps))
        sw2=ifelse(d3$missing==1, num0/ps2, (1-num0)/(1-ps2))
        
        d3$sw=sw
        d3$sw2=sw2
        d3$sw_round=round(sw)
        d3$sw_round2=round(sw2)
        
        d4=na.omit(d3)
        # truncation
        ipwcrtgeedr=myTryCatch(geeDREstimation(formula=y~x+arm,
                                               id="cluster" , data = d3,
                                               nameMISS='missing',nameY='y',
                                               nameTRT='arm',
                                               weights = d3$weight,
                                               family =  binomial("logit"),
                                               corstr = "exchangeable"))
        # truncation
        ipwcrtgeedr_yes=myTryCatch(geeDREstimation(formula=y~x+arm,
                                                   id="cluster" , data = d3,
                                                   nameMISS='missing',nameY='y',
                                                   nameTRT='arm',
                                                   family =  binomial("logit"),
                                                   corstr = "exchangeable",
                                                   weights = d3$weight2))
        
        # stabilization
        ipwcrtgeedrs=myTryCatch(geeDREstimation(formula=y~x+arm,
                                                id="cluster" , data = d3,
                                                nameMISS='missing',nameY='y',
                                                nameTRT='arm',
                                                weights = d3$sw,
                                                family =  binomial("logit"),
                                                corstr = "exchangeable"))
        
        # stabilization
        ipwcrtgeedr_yess=myTryCatch(geeDREstimation(formula=y~x+arm,
                                                    id="cluster" , data = d3,
                                                    nameMISS='missing',nameY='y',
                                                    nameTRT='arm',
                                                    family =  binomial("logit"),
                                                    corstr = "exchangeable",
                                                    weights = d3$sw2))
        
        if(is.null(ipwcrtgeedr$value)+is.null(ipwcrtgeedr_yes$value)+
           is.null(ipwcrtgeedrs$value)+is.null(ipwcrtgeedr_yess$value)==0){
          ipwtime=c(ipwtime,tt)
          ipwkm=c(ipwkm,km)
          ipwmis=c(ipwmis,mis)
          ipws=c(ipws,s)
          missper=c(missper,miss_per(d1)) 
          
          formu=formula(y~x+arm)
          true=myTryCatch(geese(formu,data=d1,id=cluster,
                                family = binomial(link='logit'),
                                corstr = 'exchangeable'))
          true_est=c(true_est,summary(true$value)$mean['arm','estimate'])
          if(is.null(true$warning)==0){
            true_warning=c(true_warning,tt)}
          if(is.null(true$warning)!=0){
            true_warning=c(true_warning,0)}
          ### unadjusted CRA
          d_nona=na.omit(d1)
          formu0=formula(y~arm)
          uncra=myTryCatch(geese(formu0,data=d_nona,id=cluster,
                                 family = binomial(link='logit'),
                                 corstr = 'exchangeable'))
          uncra_est=c(uncra_est,summary(uncra$value)$mean['arm','estimate'])
          uncra_std=c(uncra_std,summary(uncra$value)$mean['arm','san.se'])
          if(is.null(uncra$warning)==0){
            uncra_warning=c(uncra_warning,tt)}
          if(is.null(uncra$warning)!=0){
            uncra_warning=c(uncra_warning,0)}
          ### adjusted CRA
          cra=myTryCatch(geese(formu,data=d_nona,id=cluster,
                               family = binomial(link='logit'),
                               corstr = 'exchangeable'))
          cra_est=c(cra_est,summary(cra$value)$mean['arm','estimate'])
          cra_std=c(cra_std,summary(cra$value)$mean['arm','san.se'])
          if(is.null(cra$warning)==0){
            cra_warning=c(cra_warning,tt)}
          if(is.null(cra$warning)!=0){
            cra_warning=c(cra_warning,0)}
          ### IPW-no cluster effect, by using glm and geese
          
          ipwgeeno=myTryCatch(geese(formu,data=d4,id=cluster,
                                    family = binomial(link='logit'),
                                    weights = d4$weight,
                                    corstr = 'exchangeable'))
          tempgeeno_est=summary(ipwgeeno$value)$mean['arm','estimate']
          tempgeeno_std=summary(ipwgeeno$value)$mean['arm','san.se']
          ipwgeeno_est=c(ipwgeeno_est,tempgeeno_est)
          ipwgeeno_std=c(ipwgeeno_std,tempgeeno_std)
          ipwgeeno=myTryCatch(geese(formu,data=d4,id=cluster,
                                    family = binomial(link='logit'),
                                    weights = d4$weight_round,
                                    corstr = 'exchangeable'))
          if(is.null(ipwgeeno$warning)==0){ipwgeeno_warn=c(ipwgeeno_warn,tt)}
          if(is.null(ipwgeeno$warning)!=0){ipwgeeno_warn=c(ipwgeeno_warn,0)}
          
          ### IPW-GEE with package CRTgeeDR
          
          tempcrtgee_est=summary(ipwcrtgeedr$value)$beta[3]
          tempcrtgee_std=summary(ipwcrtgeedr$value)$se.robust[3]
          ipwcrtgeedr_est=c(ipwcrtgeedr_est,tempcrtgee_est)
          ipwcrtgeedr_std=c(ipwcrtgeedr_std,tempcrtgee_std)
          if(ipwcrtgeedr$value$converged==0){ipwcrtgee_warn=c(ipwcrtgee_warn,tt)}
          if(ipwcrtgeedr$value$converged==1){ipwcrtgee_warn=c(ipwcrtgee_warn,0)}
          
          ### IPW-WITH cluster effect, by using glm and geese
          formu=formula(y~x+arm)
          ipwgeeyes=myTryCatch(geese(formu,data=d4,id=cluster,
                                     family = binomial(link='logit'),
                                     weights = d4$weight2,
                                     corstr = 'exchangeable'))
          tempgeeyes_est=summary(ipwgeeyes$value)$mean['arm','estimate']
          tempgeeyes_std=summary(ipwgeeyes$value)$mean['arm','san.se']
          ipwgeeyes_est=c(ipwgeeyes_est,tempgeeyes_est)
          ipwgeeyes_std=c(ipwgeeyes_std,tempgeeyes_std)
          ipwgeeyes=myTryCatch(geese(formu,data=d4,id=cluster,
                                     family = binomial(link='logit'),
                                     weights = d4$weight_round2,
                                     corstr = 'exchangeable'))
          
          if(is.null(ipwgeeyes$warning)==0){ipwgeeyes_warn=c(ipwgeeyes_warn,tt)}
          if(is.null(ipwgeeyes$warning)!=0){ipwgeeyes_warn=c(ipwgeeyes_warn,0)}
          
          ### IPW-GEE with cluster effect with package CRTgeeDR
          
          tempcrtgee_yes_est=summary(ipwcrtgeedr_yes$value)$beta[3]
          tempcrtgee_yes_std=summary(ipwcrtgeedr_yes$value)$se.robust[3]
          ipwcrtgeedr_yes_est=c(ipwcrtgeedr_yes_est,tempcrtgee_yes_est)
          ipwcrtgeedr_yes_std=c(ipwcrtgeedr_yes_std,tempcrtgee_yes_std)
          if(ipwcrtgeedr_yes$value$converged==0){ipwcrtgee_yes_warn=c(ipwcrtgee_yes_warn,tt)}
          if(ipwcrtgeedr_yes$value$converged==1){ipwcrtgee_yes_warn=c(ipwcrtgee_yes_warn,0)}
          
          
          
          ### stabilization
          ipwgeenos=myTryCatch(geese(formu,data=d4,id=cluster,
                                     family = binomial(link='logit'),
                                     weights = d4$sw,
                                     corstr = 'exchangeable'))
          tempgeeno_est=summary(ipwgeenos$value)$mean['arm','estimate']
          tempgeeno_std=summary(ipwgeenos$value)$mean['arm','san.se']
          ipwgeeno_ests=c(ipwgeeno_ests,tempgeeno_est)
          ipwgeeno_stds=c(ipwgeeno_stds,tempgeeno_std)
          ipwgeenos=myTryCatch(geese(formu,data=d4,id=cluster,
                                     family = binomial(link='logit'),
                                     weights = d4$sw_round,
                                     corstr = 'exchangeable'))
          if(is.null(ipwgeenos$warning)==0){ipwgeeno_warns=c(ipwgeeno_warns,tt)}
          if(is.null(ipwgeenos$warning)!=0){ipwgeeno_warns=c(ipwgeeno_warns,0)}
          
          ### IPW-GEE with package CRTgeeDR
          
          tempcrtgee_est=summary(ipwcrtgeedrs$value)$beta[3]
          tempcrtgee_std=summary(ipwcrtgeedrs$value)$se.robust[3]
          ipwcrtgeedr_ests=c(ipwcrtgeedr_ests,tempcrtgee_est)
          ipwcrtgeedr_stds=c(ipwcrtgeedr_stds,tempcrtgee_std)
          if(ipwcrtgeedrs$value$converged==0){ipwcrtgee_warns=c(ipwcrtgee_warns,tt)}
          if(ipwcrtgeedrs$value$converged==1){ipwcrtgee_warns=c(ipwcrtgee_warns,0)}
          
          ### IPW-WITH cluster effect, by using glm and geese
          formu=formula(y~x+arm)
          ipwgeeyess=myTryCatch(geese(formu,data=d4,id=cluster,
                                      family = binomial(link='logit'),
                                      weights = d4$sw2,
                                      corstr = 'exchangeable'))
          tempgeeyes_est=summary(ipwgeeyess$value)$mean['arm','estimate']
          tempgeeyes_std=summary(ipwgeeyess$value)$mean['arm','san.se']
          ipwgeeyes_ests=c(ipwgeeyes_ests,tempgeeyes_est)
          ipwgeeyes_stds=c(ipwgeeyes_stds,tempgeeyes_std)
          ipwgeeyess=myTryCatch(geese(formu,data=d4,id=cluster,
                                      family = binomial(link='logit'),
                                      weights = d4$sw_round2,
                                      corstr = 'exchangeable'))
          
          if(is.null(ipwgeeyess$warning)==0){ipwgeeyes_warns=c(ipwgeeyes_warns,tt)}
          if(is.null(ipwgeeyess$warning)!=0){ipwgeeyes_warns=c(ipwgeeyes_warns,0)}
          
          ### IPW-GEE with cluster effect with package CRTgeeDR
          
          tempcrtgee_yes_est=summary(ipwcrtgeedr_yess$value)$beta[3]
          tempcrtgee_yes_std=summary(ipwcrtgeedr_yess$value)$se.robust[3]
          ipwcrtgeedr_yes_ests=c(ipwcrtgeedr_yes_ests,tempcrtgee_yes_est)
          ipwcrtgeedr_yes_stds=c(ipwcrtgeedr_yes_stds,tempcrtgee_yes_std)
          if(ipwcrtgeedr_yess$value$converged==0){ipwcrtgee_yes_warns=c(ipwcrtgee_yes_warns,tt)}
          if(ipwcrtgeedr_yess$value$converged==1){ipwcrtgee_yes_warns=c(ipwcrtgee_yes_warns,0)}
        }}
    }
    data_mean=data.frame(uncra_est=uncra_est,
                         cra_est=cra_est,
                         ipwgeeno_est=ipwgeeno_est,
                         ipwcrtgeedr_est=ipwcrtgeedr_est,
                         ipwgeeyes_est=ipwgeeyes_est,
                         ipwcrtgeedr_yes_est=ipwcrtgeedr_yes_est,
                         ipwgeeno_ests=ipwgeeno_ests,
                         ipwcrtgeedr_ests=ipwcrtgeedr_ests,
                         ipwgeeyes_ests=ipwgeeyes_ests,
                         ipwcrtgeedr_yes_ests=ipwcrtgeedr_yes_ests,
                         ipwtime=ipwtime,
                         ipwk=ipwkm,true_est=true_est,missper=missper)
    data_sd=data.frame(uncra_std=uncra_std,
                       cra_std=cra_std,
                       ipwgeeno_std=ipwgeeno_std,
                       ipwcrtgeedr_std=ipwcrtgeedr_std,
                       ipwgeeyes_std=ipwgeeyes_std,
                       ipwcrtgeedr_yes_std=ipwcrtgeedr_yes_std,
                       ipwgeeno_stds=ipwgeeno_stds,
                       ipwcrtgeedr_stds=ipwcrtgeedr_stds,
                       ipwgeeyes_stds=ipwgeeyes_stds,
                       ipwcrtgeedr_yes_stds=ipwcrtgeedr_yes_stds,
                       ipwtime=ipwtime,
                       ipwk=ipwkm,true_est=true_est)
    data_warn=data.frame(uncra_warn=uncra_warning,
                         cra_warn=cra_warning,
                         ipwgeeno_warn=ipwgeeno_warn,
                         ipwcrtgeedr_warn=ipwcrtgee_warn,
                         ipwgeeyes_warn=ipwgeeyes_warn,
                         ipwcrtgeedr_yes_warn=ipwcrtgee_yes_warn,
                         ipwgeeno_warns=ipwgeeno_warns,
                         ipwcrtgeedr_warns=ipwcrtgee_warns,
                         ipwgeeyes_warns=ipwgeeyes_warns,
                         ipwcrtgeedr_yes_warns=ipwcrtgee_yes_warns,
                         ipwtime=ipwtime,
                         ipwk=ipwkm,true_warning=true_warning)
    fff=list(data_mean=data_mean,data_sd=data_sd,data_warn=data_warn)
    file_name=paste('ipws3m3_ex',s,mis,icc*100,'.RData',sep='')
    save(fff,file=file_name)
  }
}


############ Calculate the results
data_final4=c()
for(mis in 1:6){
  nnames=paste('C:/Users/apple/Documents/ipws3m1ind',s,mis,'5.RData',sep='')
  load(nnames)
  data_mean=fff$data_mean
  data_sd=fff$data_sd
  data_warn=fff$data_warn
  
  for(kk in 1:4){
    # in kk th scenario, (km=1,2,3,4)
    data_mean_k=data_mean[data_mean$ipwk==kk,]
    data_sd_k=data_sd[data_sd$ipwk==kk,]
    data_warn_k=data_warn[data_warn$ipwk==kk,]
    
    for(par in c(1:6,9)){
      #print(par)
      data_mean_k[,par]=ifelse(data_warn_k[,par]!=0,NA,data_mean_k[,par])
      #data_mean_k[,par]=ifelse(abs(data_mean_k[,par])>3,NA,data_mean_k[,par])
      data_sd_k[,par]=ifelse(data_warn_k[,par]!=0,NA,data_sd_k[,par])
      #data_sd_k[,par]=ifelse(abs(data_sd_k[,par])>3,NA,data_sd_k[,par])
    }
    
    ms=apply(data_mean_k,2,mean,na.rm=TRUE);ms=round(ms,3)
    mcsd=apply(data_mean_k,2,sd,na.rm=TRUE);mcsd=round(mcsd,3)
    sds=apply(data_sd_k,2,mean,na.rm=TRUE);sds=round(sds,3)
    
    convs=c()
    unc_time=c()
    
    trues=mean(data_mean_k$true_est,na.rm=TRUE)
    misp=mean(data_mean_k$missper,na.rm=TRUE)
    
    for(par in 1:6){
      ns=sum(is.na(data_mean_k[,par])==0)
      convs=c(convs,round(conv(data_mean_k[,par],data_sd_k[,par],trues)/ns,2))
      unc_time=c(unc_time,sum(is.na(data_mean_k[,par]),na.rm=TRUE))
    }
    data_final3=data.frame(mis=mis,k=kk,trues=trues,misp=misp,ms=ms[1:6],
                           mcsd=mcsd[1:6],sds=sds[1:6],convs=convs, unc_time=unc_time)
    data_final4=rbind(data_final4,data_final3)
  }
}
write.csv(data_final4,paste('datas3m1ind',s,icc,'.csv',sep=''))

data_final4=c()
for(mis in 1:6){
  nnames=paste('C:/Users/apple/Documents/ipws3m1_ex',s,mis,'5.RData',sep='')
  load(nnames)
  data_mean=fff$data_mean
  data_sd=fff$data_sd
  data_warn=fff$data_warn
  
  for(kk in 1:4){
    data_mean_k=data_mean[data_mean$ipwk==kk,]
    data_sd_k=data_sd[data_sd$ipwk==kk,]
    data_warn_k=data_warn[data_warn$ipwk==kk,]
    
    for(par in c(1:6,9)){
      #print(par)
      data_mean_k[,par]=ifelse(data_warn_k[,par]!=0,NA,data_mean_k[,par])
      data_mean_k[,par]=ifelse(abs(data_mean_k[,par])>3,NA,data_mean_k[,par])
      data_sd_k[,par]=ifelse(data_warn_k[,par]!=0,NA,data_sd_k[,par])
      data_sd_k[,par]=ifelse(abs(data_sd_k[,par])>3,NA,data_sd_k[,par])
    }
    ms=apply(data_mean_k,2,mean,na.rm=TRUE);ms=round(ms,3)
    mcsd=apply(data_mean_k,2,sd,na.rm=TRUE);mcsd=round(mcsd,3)
    sds=apply(data_sd_k,2,mean,na.rm=TRUE);sds=round(sds,3)
    convs=c()
    unc_time=c()
    trues=mean(data_mean_k$true_est,na.rm=TRUE)
    misp=mean(data_mean_k$missper,na.rm=TRUE)
    for(par in 1:6){
      ns=sum(is.na(data_mean_k[,par])==0)
      convs=c(convs,round(conv(data_mean_k[,par],data_sd_k[,par],trues)/ns,2))
      unc_time=c(unc_time,sum(is.na(data_mean_k[,par]),na.rm=TRUE))
    }
    data_final3=data.frame(mis=mis,k=kk,trues=trues,misp=misp,ms=ms[1:6],
                           mcsd=mcsd[1:6],sds=sds[1:6],convs=convs, unc_time=unc_time)
    data_final4=rbind(data_final4,data_final3)
  }
}
write.csv(data_final4,paste('datas3m1_ex',s,icc,'.csv',sep=''))

data_final4=c()
for(mis in 1:6){
  nnames=paste('C:/Users/apple/Documents/ipws3m3ind',s,mis,'5.RData',sep='')
  load(nnames)
  data_mean=fff$data_mean
  data_sd=fff$data_sd
  data_warn=fff$data_warn
  
  for(kk in 1:4){
    data_mean_k=data_mean[data_mean$ipwk==kk,]
    data_sd_k=data_sd[data_sd$ipwk==kk,]
    data_warn_k=data_warn[data_warn$ipwk==kk,]
    
    for(par in c(1:10,13)){
      #print(par)
      data_mean_k[,par]=ifelse(data_warn_k[,par]!=0,NA,data_mean_k[,par])
      data_mean_k[,par]=ifelse(abs(data_mean_k[,par])>3,NA,data_mean_k[,par])
      data_sd_k[,par]=ifelse(data_warn_k[,par]!=0,NA,data_sd_k[,par])
      data_sd_k[,par]=ifelse(abs(data_sd_k[,par])>3,NA,data_sd_k[,par])
    }
    ms=apply(data_mean_k,2,mean,na.rm=TRUE);ms=round(ms,3)
    mcsd=apply(data_mean_k,2,sd,na.rm=TRUE);mcsd=round(mcsd,3)
    sds=apply(data_sd_k,2,mean,na.rm=TRUE);sds=round(sds,3)
    convs=c()
    unc_time=c()
    trues=mean(data_mean_k$true_est,na.rm=TRUE)
    misp=mean(data_mean_k$missper,na.rm=TRUE)
    for(par in 1:10){
      ns=sum(is.na(data_mean_k[,par])==0)
      convs=c(convs,round(conv(data_mean_k[,par],data_sd_k[,par],trues)/ns,2))
      unc_time=c(unc_time,sum(is.na(data_mean_k[,par]),na.rm=TRUE))
    }
    data_final3=data.frame(mis=mis,k=kk,trues=trues,misp=misp,ms=ms[1:10],
                           mcsd=mcsd[1:10],sds=sds[1:10],convs=convs, unc_time=unc_time)
    data_final4=rbind(data_final4,data_final3)
  }
}
write.csv(data_final4,paste('datas3m3ind',s,icc,'.csv',sep=''))

data_final4=c()
for(mis in 1:6){
  nnames=paste('C:/Users/apple/Documents/ipws3m3_ex',s,mis,'5.RData',sep='')
  load(nnames)
  data_mean=fff$data_mean
  data_sd=fff$data_sd
  data_warn=fff$data_warn
  
  for(kk in 1:4){
    data_mean_k=data_mean[data_mean$ipwk==kk,]
    data_sd_k=data_sd[data_sd$ipwk==kk,]
    data_warn_k=data_warn[data_warn$ipwk==kk,]
    
    for(par in c(1:10,13)){
      #print(par)
      data_mean_k[,par]=ifelse(data_warn_k[,par]!=0,NA,data_mean_k[,par])
      data_mean_k[,par]=ifelse(abs(data_mean_k[,par])>3,NA,data_mean_k[,par])
      data_sd_k[,par]=ifelse(data_warn_k[,par]!=0,NA,data_sd_k[,par])
      data_sd_k[,par]=ifelse(abs(data_sd_k[,par])>3,NA,data_sd_k[,par])
    }
    ms=apply(data_mean_k,2,mean,na.rm=TRUE);ms=round(ms,3)
    mcsd=apply(data_mean_k,2,sd,na.rm=TRUE);mcsd=round(mcsd,3)
    sds=apply(data_sd_k,2,mean,na.rm=TRUE);sds=round(sds,3)
    convs=c()
    unc_time=c()
    trues=mean(data_mean_k$true_est,na.rm=TRUE)
    misp=mean(data_mean_k$missper,na.rm=TRUE)
    for(par in 1:10){
      ns=sum(is.na(data_mean_k[,par])==0)
      convs=c(convs,round(conv(data_mean_k[,par],data_sd_k[,par],trues)/ns,2))
      unc_time=c(unc_time,sum(is.na(data_mean_k[,par]),na.rm=TRUE))
    }
    data_final3=data.frame(mis=mis,k=kk,trues=trues,misp=misp,ms=ms[1:10],
                           mcsd=mcsd[1:10],sds=sds[1:10],convs=convs, unc_time=unc_time)
    data_final4=rbind(data_final4,data_final3)
  }
}
write.csv(data_final4,paste('datas3m3_ex',s,icc,'.csv',sep=''))


# convergence table

setwd('c:\\Users\\apple\\Documents')

ex1=read.csv('datas3m1_ex30.05.csv')
ex1=ex1[ex1$k==3,]
ex3=read.csv('datas3m3_ex30.05.csv')
ex3=ex3[ex3$k==3,]

rownames(ex1)=NULL
rownames(ex3)=NULL
ex1=ex1[,c('X','unc_time')]
ex3=ex3[,c('X','unc_time')]

m1=c()
m2=c()
m3=c()

for(i in 1:6){
  m1=rbind(m1,ex1[(3+(i-1)*6):(6+(i-1)*6),]) 
  m2=rbind(m2,ex3[(3+(i-1)*10):(6+(i-1)*10),])
  m3=rbind(m3,ex3[(7+(i-1)*10):(10+(i-1)*10),])
}

nam=c('IPW1','IPW2','IPW1_CLU','IPW2_CLU')
m1$X=rep(nam,6)
m2$X=rep(nam,6)
m3$X=rep(nam,6)

ind1=read.csv('datas3m1ind30.05.csv')
ind1=ind1[ind1$k==3,]
ind3=read.csv('datas3m3ind30.05.csv')
ind3=ind3[ind3$k==3,]

rownames(ind1)=NULL
rownames(ind3)=NULL
ind1=ind1[,c('X','unc_time')]
ind3=ind3[,c('X','unc_time')]

M1=c()
M2=c()
M3=c()

for(i in 1:6){
  M1=rbind(M1,ind1[c((3+(i-1)*7):(4+(i-1)*7),(6+(i-1)*7):(7+(i-1)*7)),]) 
  M2=rbind(M2,ind3[(3+(i-1)*10):(6+(i-1)*10),])
  M3=rbind(M3,ind3[(7+(i-1)*10):(10+(i-1)*10),])
}

nam=c('IPW1','IPW2','IPW1_CLU','IPW2_CLU')
M1$X=rep(nam,6)
M2$X=rep(nam,6)
M3$X=rep(nam,6)

rownames(M1)=NULL
rownames(M2)=NULL

unc=cbind(M1,M2,M3,m1,m2,m3)
unc=unc[,c(1,2,seq(4,12,2))]

names(unc)=c('Methods','M1','M2','M3','M1','M2','M3')


###************************ Data Analysis ***********************###
###************************      MI      ***********************###

# for multiple imputation, we considered standard MI and MMI, by using independent/ exchangeable working correlation matrix
for(mis in 1:6){
  
  ### Data Generation
  error_mmi_result=c();warning_mmi_result=c()
  jomo_mean=c();jomo_std=c()
  jomo_mean_ex=c();jomo_std_ex=c()
  errors_jomo=c();warn_jomo=c()
  error_mmi_result2=c();warning_mmi_result2=c()
  jomo_mean2=c();jomo_std2=c()
  jomo_mean2_ex=c();jomo_std2_ex=c()
  errors_jomo2=c();warn_jomo2=c()
  
  nimp=15
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
    for(tt in 1:times){
      print('****************')
      print('mis')
      print(mis)
      print(km)
      print(tt)
      d1=data_gene3(k=k,m=m,s=s,mis=mis,seed=tt)
      d2=d1
      d2$km=km
      d2$times=tt
      d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm)
      d3$y=ifelse(d2$R==1,NA,d3$y)
      d3$y=as.factor(d3$y)
      
      # standard MI
      d5=d3[,-3]
      d6=jomo(d5,nimp=nimp)
      rmmi2=myTryCatch(
        result_mmi(d6,org_data=d3,num=nimp))
      
      # MMI
      if(mis==1){
        tempd=d3[,c('y','x')]
        d4=jomo(tempd,nimp=nimp)
      }
      if(mis==2){
        tempd=d3[,c('y','x','arm')]
        d4=jomo(tempd,nimp=nimp)
      }
      if(mis==3){
        d3$armx=d3$x*d3$arm
        tempd=d3[,c('y','x','arm','armx')]
        d4=jomo(tempd,nimp=nimp)
      }
      if(mis==4){
        tempd=d3[,c('y','x','cluster')]
        d4=jomo(tempd,clus = d3$cluster,nimp=nimp)
      }
      if(mis==5){
        tempd=d3[,c('y','x','arm','cluster')]
        d4=jomo(tempd,clus = d3$cluster,nimp=nimp)
      }
      if(mis==6){
        d3$armx=d3$x*d3$arm
        tempd=d3[,c('y','x','arm','armx','cluster')]
        d4=jomo(tempd,clus = d3$cluster,nimp=nimp)
      }
      rmmi0=myTryCatch(
        result_mmi(d4,org_data=d3,num=nimp))
      
      error_mmi_result=c(error_mmi_result,rmmi0$error)
      warning_mmi_result=c(warning_mmi_result,rmmi0$warning)
      if(is.null(rmmi0$error)==0){errors_jomo=c(errors_jomo,tt)}
      if(is.null(rmmi0$warning)==0){warn_jomo=c(warn_jomo,tt)}
      
      rmmi=rmmi0$value
      m0=rmmi$m0
      std0=rmmi$std0
      mmi_result=mypool(m0,std0,num=nimp,print='yes')
      mmi_m=mmi_result$mean
      mmi_s=mmi_result$std
      jomo_mean=c(jomo_mean,mmi_m)
      jomo_std=c(jomo_std,mmi_s)
      
      m02=rmmi$m02
      std02=rmmi$std02
      mmi_result_ex=mypool(m02,std02,num=nimp,print='yes')
      mmi_m_ex=mmi_result_ex$mean
      mmi_s_ex=mmi_result_ex$std
      jomo_mean_ex=c(jomo_mean_ex,mmi_m_ex)
      jomo_std_ex=c(jomo_std_ex,mmi_s_ex)
      
      error_mmi_result2=c(error_mmi_result2,rmmi2$error)
      warning_mmi_result2=c(warning_mmi_result2,rmmi2$warning)
      if(is.null(rmmi2$error)==0){errors_jomo2=c(errors_jomo2,tt)}
      if(is.null(rmmi2$warning)==0){warn_jomo2=c(warn_jomo2,tt)}
      
      rmmi2=rmmi2$value
      m2=rmmi2$m0
      std2=rmmi2$std0
      mmi_result2=mypool(m02,std02,num=nimp,print='yes')
      mmi_m2=mmi_result2$mean
      mmi_s2=mmi_result2$std
      jomo_mean2=c(jomo_mean2,mmi_m2)
      jomo_std2=c(jomo_std2,mmi_s2)
      
      m2_ex=rmmi2$m02
      std2_ex=rmmi2$std02
      mmi_result2_ex=mypool(m2_ex,std2_ex,num=nimp,print='yes')
      mmi_m2_ex=mmi_result2_ex$mean
      mmi_s2_ex=mmi_result2_ex$std
      jomo_mean2_ex=c(jomo_mean2_ex,mmi_m2_ex)
      jomo_std2_ex=c(jomo_std2_ex,mmi_s2_ex)
      
    }}
  
  final_data=list(error_mmi_result=error_mmi_result,
                  warning_mmi_result=warning_mmi_result,
                  jomo_mean=jomo_mean,
                  jomo_std=jomo_std,
                  jomo_mean_ex=jomo_mean_ex,
                  jomo_std_ex=jomo_std_ex,
                  errors_jomo=errors_jomo,
                  warn_jomo= warn_jomo,
                  error_mmi_result2=error_mmi_result2,
                  warning_mmi_result2=warning_mmi_result2,
                  jomo_mean2=jomo_mean2,
                  jomo_std2=jomo_std2,
                  jomo_mean2_ex=jomo_mean2_ex,
                  jomo_std2_ex=jomo_std2_ex,
                  errors_jomo2=errors_jomo2,
                  warn_jomo2= warn_jomo2)
  
  #fff3=list(data_mean=data_mean,data_sd=data_sd,data_warn=data_warn)
  file_name=paste('ipw8mmismisicc',s,mis,icc*100,'.RData',sep='')
  save(final_data,file=file_name)
}

############ Calculate the results
data_mmi=c()
for(mis in 1:6){
  nnames=paste('C:/Users/apple/Documents/ipw8mmismisicc1',mis,'5.RData',sep='')
  load(nnames)
  
  # MMI
  mmi_m_km=matrix(0,4,100)
  mmi_sd_km=matrix(0,4,100)
  
  mmi_m_km_ex=matrix(0,4,100)
  mmi_sd_km_ex=matrix(0,4,100)
  
  smi_m_km=matrix(0,4,100)
  smi_sd_km=matrix(0,4,100)
  
  smi_m_km_ex=matrix(0,4,100)
  smi_sd_km_ex=matrix(0,4,100)
  
  convs=c()
  
  for(i in 1:4){
    mmi_m_km[i,]=final_data$jomo_mean[(1+100*(i-1)):(i*100)]
    mmi_sd_km[i,]=final_data$jomo_std[(1+100*(i-1)):(i*100)]
    mmi_m_km_ex[i,]=final_data$jomo_mean_ex[(1+100*(i-1)):(i*100)]
    mmi_sd_km_ex[i,]=final_data$jomo_std_ex[(1+100*(i-1)):(i*100)]
  }
  for(i in 1:4){
    smi_m_km[i,]=final_data$jomo_mean2[(1+100*(i-1)):(i*100)]
    smi_sd_km[i,]=final_data$jomo_std2[(1+100*(i-1)):(i*100)]
    smi_m_km_ex[i,]=final_data$jomo_mean2_ex[(1+100*(i-1)):(i*100)]
    smi_sd_km_ex[i,]=final_data$jomo_std2_ex[(1+100*(i-1)):(i*100)]
  }
  
  mmi_m_km=ifelse(abs(mmi_m_km)>3,NA,mmi_m_km)
  mmi_m_km_ex=ifelse(abs(mmi_m_km_ex)>3,NA,mmi_m_km_ex)
  smi_m_km=ifelse(abs(smi_m_km)>3,NA,smi_m_km)
  smi_m_km_ex=ifelse(abs(smi_m_km_ex)>3,NA,smi_m_km_ex)
  mmi_sd_km=ifelse(abs(mmi_sd_km)>3,NA,mmi_sd_km)
  mmi_sd_km_ex=ifelse(abs(mmi_sd_km_ex)>3,NA,mmi_sd_km_ex)
  smi_sd_km=ifelse(abs(smi_sd_km)>3,NA,smi_sd_km)
  smi_sd_km_ex=ifelse(abs(smi_sd_km_ex)>3,NA,smi_sd_km_ex)
  
  mmi_m=apply(mmi_m_km,1,mean,na.rm=TRUE)
  mmi_mcsd=apply(mmi_m_km,1,sd,na.rm=TRUE)
  mmi_sd=apply(mmi_sd_km,1,mean,na.rm=TRUE)
  mmi_m_ex=apply(mmi_m_km_ex,1,mean,na.rm=TRUE)
  mmi_mcsd_ex=apply(mmi_m_km,1,sd,na.rm=TRUE)
  mmi_sd_ex=apply(mmi_sd_km_ex,1,mean,na.rm=TRUE)
  
  smi_m=apply(smi_m_km,1,mean,na.rm=TRUE)
  smi_mcsd=apply(smi_m_km,1,sd,na.rm=TRUE)
  smi_sd=apply(smi_sd_km,1,mean,na.rm=TRUE)
  smi_m_ex=apply(smi_m_km_ex,1,mean,na.rm=TRUE)
  smi_mcsd_ex=apply(smi_m_km,1,sd,na.rm=TRUE)
  smi_sd_ex=apply(smi_sd_km_ex,1,mean,na.rm=TRUE)
  
  mi_m=c(mmi_m,mmi_m_ex,smi_m,smi_m_ex)
  mi_sd=c(mmi_sd,mmi_sd_ex,smi_sd,smi_sd_ex)
  mi_mcsd=c(mmi_mcsd,mmi_mcsd_ex,smi_mcsd,smi_mcsd_ex)
  
  trues=c(1.316,1.318,1.377,1.377)
  
  
  for(i in 1:4){convs=c(convs,(conv(mmi_m_km[i,],mmi_sd_km[i,],trues[i])))}
  for(i in 1:4){convs=c(convs,(conv(mmi_m_km_ex[i,],mmi_sd_km_ex[i,],trues[i])))}
  
  for(i in 1:4){convs=c(convs,(conv(smi_m_km[i,],smi_sd_km[i,],trues[i])))}
  for(i in 1:4){convs=c(convs,(conv(smi_m_km_ex[i,],smi_sd_km_ex[i,],trues[i])))}
  
  data_temp=data.frame(mi_m=mi_m,mi_mcsd=mi_mcsd,mi_sd=mi_sd,cov=convs,k=rep(1:4,4),mis=rep(mis,16))
  data_mmi=rbind(data_mmi,data_temp)
}
rnames=c('m','m','m','m','m_ex','m_ex','m_ex','m_ex',
         's','s','s','s','s_ex','s_ex','s_ex','s_ex')
rnames=rep(rnames,6)
data_mmi$name=rnames
write.csv(data_mmi,'data_mmi.csv')


mis=1
d1=data_gene3(k=25,m=25,s=1,mis=1,psi=-1.8)
d2=d1
d2$y=ifelse(d2$R==1,NA,d2$y)
d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=d2$R)

logs=glm(missing ~ x+arm, data = d3,
         family = binomial(link='logit'))
weight=1/expit(predict(logs))

histrv=hist(weight)

histrv$breaks
a=histrv$breaks
n=length(a)
a1=a[-n]
a2=a[-1]
breaks=paste(a1,"-",a2,sep='')
breaks=c(breaks,'total')
clusize=as.data.frame(breaks)
counts=histrv$counts
counts=c(counts,sum(counts))
clusize$counts=counts

weights=weight[weight<500]
histrv=hist(weights)

histrv$breaks
a=histrv$breaks
n=length(a)
a1=a[-n]
a2=a[-1]
breaks=paste(a1,"-",a2,sep='')
breaks=c(breaks,'total')
clusize=as.data.frame(breaks)
counts=histrv$counts
counts=c(counts,sum(counts))
clusize$counts=counts

