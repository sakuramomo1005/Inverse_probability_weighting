library('CRTgeeDR')

## CRTgeeDR package example:
data(data.sim)
#### IPW GEE
head(data.sim)
sum(data.sim$MISSING)/dim(data.sim)[1]

# in this section, MISSING = 1 means missing and MISSING = 0 means observed.
ipwresults<-geeDREstimation(formula=OUTCOME~TRT,
                            id="CLUSTER" , data = data.sim,
                            family = "binomial", corstr = "independence",
                            model.weights=I(MISSING==0)~TRT*AGE) # missing = 0, observed
summary(ipwresults)
ipwresults2<-geeDREstimation(formula=OUTCOME~TRT,
                            id="CLUSTER" , data = data.sim,
                            family = "binomial", corstr = "independence",
                            sandwich.nuisance=TRUE,
                            model.weights=I(MISSING==0)~TRT*AGE) # missing = 0, observed
summary(ipwresults2)



### in our scenarios:
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

## from missingness ICC calculate variance
missing_icc=function(rho){
  pi=3.1415926
  sigma=(pi^2/3)/(1/rho-1)
  return(sigma)
}

# for example, calculate the variance of ICC 0.05
missingicc=missing_icc(0.05)

# function for calculate intercept in missingess generation model
mis_per=function(n,m,icc,vars,seed=123){
  data=c()
  for(i in seq(-5,5,0.1)){
    temp=data_gene5(mis=5,k=n,m=m,icc=icc,intercept=i,seed=seed,vars=vars)
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
  temp=data_gene5(mis=5,k=n,m=m,icc=icc,intercept=intercept,seed,vars=vars)
  percent=sum(temp$R)/dim(temp)[1]
  print(percent)
}
################################ 1. data generation:
################################
# function for generate one intervention arm
# k cluster numbers; mm, cluster size; intercept: missingness generation intercpet; vars: uijl's variance
one_group5=function(mis=5,i,k,m,icc,intercept,vars,seed=123){
  set.seed(seed)
  b0=1;b1=1.36;b2=1;psi=intercept
  sigma_b=sqrt(0.2)
  mu_x=0
  sigma_alpha=sqrt(0.18)
  sigma_u=sqrt(vars) # change the variance 
  sigma_icc=missing_icc(icc) ## missingness ICC
  
  phi=1
  x=matrix(0,k,m+5*(m/25))  ## k clusters in each arm; max 60 in each cluster 
  y=matrix(0,k,m+5*(m/25))  ## k clusters in each arm; max 60 in each cluster 
  pi=matrix(0,k,m+5*(m/25))
  cluster=matrix(0,k,m+5*(m/25))
  r=matrix(0,k,m+5*(m/25))
  delta0=matrix(0,k,m+5*(m/25))
  mis_clu=matrix(0,k,m+5*(m/25))
  groups=matrix(i,k,m+5*(m/25))
  
  for(kk in 1:k){
    delta=rnorm(1,0,sigma_b)        #cluster level effects
    alpha=rnorm(1,mu_x,sigma_alpha) #cluster level effects
    
    m_size=round(runif(1,m-5*(m/25),m+5*(m/25))) # cluster size
    cluster[kk,]=kk # cluster marker 
    delta0[kk,]=delta
    
    # cluster effects
    clu=rnorm(1,0,sqrt(sigma_icc))
    mis_clu[kk,]=clu
    
    for(j in 1:m_size){
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
  return(list(x=c(x),y=c(y),pi=c(pi),r=c(r),cluster=c(cluster),groups=c(groups),delta=c(delta0)))
}

# function that combines two intervention arm
data_gene5=function(k,m,mis=5,icc,intercept,vars,seed=123){
  set.seed(seed)
  # step1:
  a1=one_group5(i=1,k=k,m=m,seed=seed,mis=mis,intercept=intercept,icc=icc,vars=vars)
  a0=one_group5(i=0,k=k,m=m,seed=seed,mis=mis,intercept=intercept,icc=icc,vars=vars)
  # step2
  b1=data.frame(x=c(a1$x),y=c(a1$y),r=c(a1$r),cluster=c(a1$cluster),arm=(a1$group),delta=(a1$delta))
  R=c()
  for(i in 1:dim(b1)[1]){
    R[i]=rbinom(1,1,b1$r[i])
  }
  b1$R=R
  b1$cluster=b1$cluster+k
  b0=data.frame(x=c(a0$x),y=c(a0$y),r=c(a0$r),cluster=c(a0$cluster),arm=(a0$group),delta=(a0$delta))
  R=c()
  for(i in 1:dim(b0)[1]){
    R[i]=rbinom(1,1,b0$r[i])
  }
  b0$R=R
  # step3:
  data=rbind(b0,b1)
  data0=data[data$x!=0,]
  return(data0)
}

### determin how to write the formula in model.weights
ipwtry=geeDREstimation(formula=y~arm+x,
                        nameTRT = "arm",nameMISS = "R", nameY = "y",
                             id="cluster" , data = temp,
                             family = "binomial", corstr = "independence",
                             sandwich.nuisance=TRUE,
                             model.weights=I(R==0)~arm+x) # R = 0, observed
summary(ipwtry)

## model.weights cannot add cluster effects 
ipwtry2=geeDREstimation(formula=y~arm+x,
                        nameTRT = "arm",nameMISS = "R", nameY = "y",
                        id="cluster" , data = temp,
                        family = "binomial", corstr = "independence",
                        sandwich.nuisance=TRUE,
                        model.weights=I(R==0)~arm+x+(1|cluster)) # R = 0, observed
summary(ipwtry2)


vars=0.2
## simulations
for(icc in seq(0.1,0.9,0.1)){
  print('icc')
  print(icc)
  n=50
  m=50
  
  res_est_ind1=c();res_std_ind1=c();res_warn_ind1=c()
  res_est_ex1=c();res_std_ex1=c();res_warn_ex1=c()
  res_clu_est_ind1=c();res_clu_std_ind1=c();res_clu_warn_ind1=c()
  res_clu_est_ex1=c();res_clu_std_ex1=c();res_clu_warn_ex1=c()
  
  res_est_ind2=c();res_std_ind2=c();res_warn_ind2=c()
  res_est_ex2=c();res_std_ex2=c();res_warn_ex2=c()
  res_clu_est_ind2=c();res_clu_std_ind2=c();res_clu_warn_ind2=c()
  res_clu_est_ex2=c();res_clu_std_ex2=c();res_clu_warn_ex2=c()
  
  res_est_ind3=c();res_std_ind3=c();res_warn_ind3=c()
  res_est_ex3=c();res_std_ex3=c();res_warn_ex3=c()
  res_clu_est_ind3=c();res_clu_std_ind3=c();res_clu_warn_ind3=c()
  res_clu_est_ex3=c();res_clu_std_ex3=c();res_clu_warn_ex3=c()
  
  ICCs=c()
  TIMEs=c()
  
  
  true_est_ex=c();true_std_ex=c();true_warn_ex=c()
  true_est_ind=c();true_std_ind=c();true_warn_ind=c()
  
  est=c();std=c();warn=c()
  
  #intercept=temp2[temp2$icc==icc,]$intercept
  #mp=mis_per_check(n=50,m=50,icc=icc,seed=time,intercept=intercept,vars=0.2)
  
  intercept=mis_per(n=50,m=50,icc=icc,vars=vars)
  print(paste('intercept',intercept))
  mp=mis_per_check(vars=vars,n=50,m=50,icc=icc,intercept=intercept)
  print(mp)
  
  for(times in 1:500){
    print(paste('intercept',intercept,'mp',mp))
    print(paste('icc',icc,'times',times))
                     
    TIMEs=c(TIMEs,times)
    
    d1=data_gene5(k=50,m=50,icc=icc,intercept=intercept,seed=times,vars=vars)
    d2=d1;d2$y=ifelse(d2$R==1,NA,d2$y)
    d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=d2$R)
    
    # weight calculation
    logs=glm(missing ~ x+arm, data = d3,
             family = binomial(link='logit'))
    logs2=glmer(missing ~ x+arm+(1|cluster) , data = d3,
                family = binomial(link='logit'))
    
    weight1=1/predict(logs,type="response")
    weight2=1/predict(logs2,type="response")
    
   # weight11=1/expit(predict(logs))
   # weight22=1/expit(predict(logs2)) different methods have the same results. 
                     
    d3$weight=weight1
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
    
    # ipw_ind1,ipw with independent working correlation matrix, no cluster effects
    # ipw_ind1,by using statement weigths= ..., 
    ipw_ind1=myTryCatch(geeDREstimation(formula=y~x+arm,
                                       id="cluster" , data = d3,
                                       nameMISS='missing',nameY='y',
                                       nameTRT='arm',
                                       weights = d3$weight,
                                       family =  binomial("logit"),
                                       corstr = "independence"))
    # ipw_ind1,ipw with independent working correlation matrix, no cluster effects
    # ipw_ind1,by using statement weigths= ...,  added sandwich.nuisance=TRUE,
    ipw_ind2=myTryCatch(geeDREstimation(formula=y~x+arm,
                                        id="cluster" , data = d3,
                                        nameMISS='missing',nameY='y',
                                        nameTRT='arm',
                                        weights = d3$weight,
                                        sandwich.nuisance=TRUE,
                                        family =  binomial("logit"),
                                        corstr = "independence"))
    # ipw_ind1,ipw with independent working correlation matrix, no cluster effects
    # ipw_ind1,by using statement model.weigths 
    ipw_ind3=myTryCatch(geeDREstimation(formula=y~x+arm,
                           nameTRT = "arm",nameMISS = "missing", nameY = "y",
                           id="cluster" , data = d3,
                           family =  binomial("logit"), corstr = "independence",
                           sandwich.nuisance=TRUE,
                           model.weights=I(missing==0)~arm+x)) # missing = 0, observed
    
    # ipw_clu_ind,ipw with independent working correlation matrix, cluster effects
    # by using statement weigths=
    ipw_clu_ind1=myTryCatch(geeDREstimation(formula=y~x+arm,
                                           id="cluster" , data = d3,
                                           nameMISS='missing',nameY='y',
                                           nameTRT='arm',
                                           weights = d3$weight2,
                                           family =  binomial("logit"),
                                           corstr = "independence"))
    # add  sandwich.nuisance = TRUE
    ipw_clu_ind2=myTryCatch(geeDREstimation(formula=y~x+arm,
                                           id="cluster" , data = d3,
                                           nameMISS='missing',nameY='y',
                                           nameTRT='arm',
                                           weights = d3$weight2,
                                           sandwich.nuisance = TRUE,
                                           family =  binomial("logit"),
                                           corstr = "independence"))
    # ipw_clu_ind,ipw with independent working correlation matrix, cluster effects
    # by using statement model.weights
    ipw_clu_ind3=myTryCatch(geeDREstimation(formula=y~x+arm,
                                           id="cluster" , data = d3,
                                           nameMISS='missing',nameY='y',
                                           nameTRT='arm',
                                           sandwich.nuisance=TRUE,
                                           model.weights=I(missing==0)~arm+x,
                                           family =  binomial("logit"),
                                           corstr = "independence")) 

# exchangeable

    # ipw_ex1,ipw with independent working correlation matrix, no cluster effects
    # ipw_ex1,by using statement weigths= ..., 
    ipw_ex1=myTryCatch(geeDREstimation(formula=y~x+arm,
                                    id="cluster" , data = d3,
                                    nameMISS='missing',nameY='y',
                                    nameTRT='arm',
                                    weights = d3$weight,
                                    family =  binomial("logit"),
                                    corstr = "exchangeable"))
    # ipw_ex2,ipw with independent working correlation matrix, no cluster effects
    # ipw_ex21,by using statement weigths= ...,  added sandwich.nuisance=TRUE,
    ipw_ex2=myTryCatch(geeDREstimation(formula=y~x+arm,
                                    id="cluster" , data = d3,
                                    nameMISS='missing',nameY='y',
                                    nameTRT='arm',
                                    weights = d3$weight,
                                    sandwich.nuisance=TRUE,
                                    family =  binomial("logit"),
                                    corstr = "exchangeable"))
    # ipw_ind1,ipw with independent working correlation matrix, no cluster effects
    # ipw_ind1,by using statement model.weigths 
    ipw_ex3=myTryCatch(geeDREstimation(formula=y~x+arm,
                                    nameTRT = "arm",nameMISS = "missing", nameY = "y",
                                    id="cluster" , data = d3,
                                    family =  binomial("logit"), 
                                    corstr = "exchangeable",
                                    sandwich.nuisance=TRUE,
                                    model.weights=I(missing==0)~arm+x)) # missing = 0, observed

    # ipw_clu_ex1,ipw with independent working correlation matrix, cluster effects
    # by using statement weigths=
    ipw_clu_ex1=myTryCatch(geeDREstimation(formula=y~x+arm,
                                        id="cluster" , data = d3,
                                        nameMISS='missing',nameY='y',
                                        nameTRT='arm',
                                        weights = d3$weight2,
                                        family =  binomial("logit"),
                                        corstr = "exchangeable"))
    # add  sandwich.nuisance = TRUE
    ipw_clu_ex2=myTryCatch(geeDREstimation(formula=y~x+arm,
                                        id="cluster" , data = d3,
                                        nameMISS='missing',nameY='y',
                                        nameTRT='arm',
                                        weights = d3$weight2,
                                        sandwich.nuisance = TRUE,
                                        family =  binomial("logit"),
                                        corstr = "exchangeable"))
    # ipw_clu_ex3,ipw with independent working correlation matrix, cluster effects
    # by using statement model.weights
    ipw_clu_ex3=myTryCatch(geeDREstimation(formula=y~x+arm,
                                        id="cluster" , data = d3,
                                        nameMISS='missing',nameY='y',
                                        nameTRT='arm',
                                        sandwich.nuisance=TRUE,
                                        model.weights=I(missing==0)~arm+x,
                                        family =  binomial("logit"),
                                         corstr = "exchangeable"))


    
    if(is.null(ipw_ind1$value)==0){
      ipw_ind_est1=summary(ipw_ind1$value)$beta[3]
      ipw_ind_std1=summary(ipw_ind1$value)$se.robust[3]
      res_est_ind1=c(res_est_ind1,ipw_ind_est1)
      res_std_ind1=c(res_std_ind1,ipw_ind_std1)
      if(ipw_ind1$value$converged==0){res_warn_ind1=c(res_warn_ind1,time)}
      if(ipw_ind1$value$converged==1){res_warn_ind1=c(res_warn_ind1,0)}
    }
    if(is.null(ipw_ex1$value)==0){
      ipw_ex_est1=summary(ipw_ex1$value)$beta[3]
      ipw_ex_std1=summary(ipw_ex1$value)$se.robust[3]
      res_est_ex1=c(res_est_ex1,ipw_ex_est1)
      res_std_ex1=c(res_std_ex1,ipw_ex_std1)
      if(ipw_ex1$value$converged==0){res_warn_ex1=c(res_warn_ex1,time)}
      if(ipw_ex1$value$converged==1){res_warn_ex1=c(res_warn_ex1,0)}
    }
    if(is.null(ipw_clu_ind1$value)==0){
      ipw_clu_ind_est1=summary(ipw_clu_ind1$value)$beta[3]
      ipw_clu_ind_std1=summary(ipw_clu_ind1$value)$se.robust[3]
      res_clu_est_ind1=c(res_clu_est_ind1,ipw_clu_ind_est1)
      res_clu_std_ind1=c(res_clu_std_ind1,ipw_clu_ind_std1)
      if(ipw_clu_ind1$value$converged==0){res_clu_warn_ind1=c(res_clu_warn_ind1,time)}
      if(ipw_clu_ind1$value$converged==1){res_clu_warn_ind1=c(res_clu_warn_ind1,0)}
    }
    if(is.null(ipw_clu_ex1$value)==0){
      ipw_clu_ex_est1=summary(ipw_clu_ex1$value)$beta[3]
      ipw_clu_ex_std1=summary(ipw_clu_ex1$value)$se.robust[3]
      res_clu_est_ex1=c(res_clu_est_ex1,ipw_clu_ex_est1)
      res_clu_std_ex1=c(res_clu_std_ex1,ipw_clu_ex_std1)
      if(ipw_clu_ex1$value$converged==0){res_clu_warn_ex1=c(res_clu_warn_ex1,time)}
      if(ipw_clu_ex1$value$converged==1){res_clu_warn_ex1=c(res_clu_warn_ex1,0)}
    }
    
    if(is.null(ipw_ind1$value)==1){
      res_est_ind1=c(res_est_ind1,NA)
      res_std_ind1=c(res_std_ind1,NA)
      res_warn_ind1=c(res_warn_ind1,time)
    }
    if(is.null(ipw_ex1$value)==1){
      res_est_ex1=c(res_est_ex1,NA)
      res_std_ex1=c(res_std_ex1,NA)
      res_warn_ex1=c(res_warn_ex1,time)
    }
    if(is.null(ipw_clu_ind1$value)==1){
      res_clu_est_ind1=c(res_clu_est_ind1,NA)
      res_clu_std_ind1=c(res_clu_std_ind1,NA)
      res_clu_warn_ind1=c(res_clu_warn_ind1,time)
    }
    if(is.null(ipw_clu_ex1$value)==1){
      res_clu_est_ex1=c(res_clu_est_ex1,NA)
      res_clu_std_ex1=c(res_clu_std_ex1,NA)
      res_clu_warn_ex1=c(res_clu_warn_ex1,time)
    }
    
    
    if(is.null(ipw_ind2$value)==0){
      ipw_ind_est2=summary(ipw_ind2$value)$beta[3]
      ipw_ind_std2=summary(ipw_ind2$value)$se.robust[3]
      res_est_ind2=c(res_est_ind2,ipw_ind_est2)
      res_std_ind2=c(res_std_ind2,ipw_ind_std2)
      if(ipw_ind2$value$converged==0){res_warn_ind2=c(res_warn_ind2,time)}
      if(ipw_ind2$value$converged==1){res_warn_ind2=c(res_warn_ind2,0)}
    }
    if(is.null(ipw_ex2$value)==0){
      ipw_ex_est2=summary(ipw_ex2$value)$beta[3]
      ipw_ex_std2=summary(ipw_ex2$value)$se.robust[3]
      res_est_ex2=c(res_est_ex2,ipw_ex_est2)
      res_std_ex2=c(res_std_ex2,ipw_ex_std2)
      if(ipw_ex2$value$converged==0){res_warn_ex2=c(res_warn_ex2,time)}
      if(ipw_ex2$value$converged==1){res_warn_ex2=c(res_warn_ex2,0)}
    }
    if(is.null(ipw_clu_ind2$value)==0){
      ipw_clu_ind_est2=summary(ipw_clu_ind2$value)$beta[3]
      ipw_clu_ind_std2=summary(ipw_clu_ind2$value)$se.robust[3]
      res_clu_est_ind2=c(res_clu_est_ind2,ipw_clu_ind_est2)
      res_clu_std_ind2=c(res_clu_std_ind2,ipw_clu_ind_std2)
      if(ipw_clu_ind2$value$converged==0){res_clu_warn_ind2=c(res_clu_warn_ind2,time)}
      if(ipw_clu_ind2$value$converged==1){res_clu_warn_ind2=c(res_clu_warn_ind2,0)}
    }
    if(is.null(ipw_clu_ex2$value)==0){
      ipw_clu_ex_est2=summary(ipw_clu_ex2$value)$beta[3]
      ipw_clu_ex_std2=summary(ipw_clu_ex2$value)$se.robust[3]
      res_clu_est_ex2=c(res_clu_est_ex2,ipw_clu_ex_est2)
      res_clu_std_ex2=c(res_clu_std_ex2,ipw_clu_ex_std2)
      if(ipw_clu_ex2$value$converged==0){res_clu_warn_ex2=c(res_clu_warn_ex2,time)}
      if(ipw_clu_ex2$value$converged==1){res_clu_warn_ex2=c(res_clu_warn_ex2,0)}
    }
    
    if(is.null(ipw_ind2$value)==1){
      res_est_ind2=c(res_est_ind2,NA)
      res_std_ind2=c(res_std_ind2,NA)
      res_warn_ind2=c(res_warn_ind2,time)
    }
    if(is.null(ipw_ex2$value)==1){
      res_est_ex2=c(res_est_ex2,NA)
      res_std_ex2=c(res_std_ex2,NA)
      res_warn_ex2=c(res_warn_ex2,time)
    }
    if(is.null(ipw_clu_ind2$value)==1){
      res_clu_est_ind2=c(res_clu_est_ind2,NA)
      res_clu_std_ind2=c(res_clu_std_ind2,NA)
      res_clu_warn_ind2=c(res_clu_warn_ind2,time)
    }
    if(is.null(ipw_clu_ex2$value)==1){
      res_clu_est_ex2=c(res_clu_est_ex2,NA)
      res_clu_std_ex2=c(res_clu_std_ex2,NA)
      res_clu_warn_ex2=c(res_clu_warn_ex2,time)
    }
    
    
    if(is.null(ipw_ind3$value)==0){
      ipw_ind_est3=summary(ipw_ind3$value)$beta[3]
      ipw_ind_std3=summary(ipw_ind3$value)$se.robust[3]
      res_est_ind3=c(res_est_ind3,ipw_ind_est3)
      res_std_ind3=c(res_std_ind3,ipw_ind_std3)
      if(ipw_ind3$value$converged==0){res_warn_ind3=c(res_warn_ind3,time)}
      if(ipw_ind3$value$converged==1){res_warn_ind3=c(res_warn_ind3,0)}
    }
    if(is.null(ipw_ex3$value)==0){
      ipw_ex_est3=summary(ipw_ex3$value)$beta[3]
      ipw_ex_std3=summary(ipw_ex3$value)$se.robust[3]
      res_est_ex3=c(res_est_ex3,ipw_ex_est3)
      res_std_ex3=c(res_std_ex3,ipw_ex_std3)
      if(ipw_ex3$value$converged==0){res_warn_ex3=c(res_warn_ex3,time)}
      if(ipw_ex3$value$converged==1){res_warn_ex3=c(res_warn_ex3,0)}
    }
    if(is.null(ipw_clu_ind3$value)==0){
      ipw_clu_ind_est3=summary(ipw_clu_ind3$value)$beta[3]
      ipw_clu_ind_std3=summary(ipw_clu_ind3$value)$se.robust[3]
      res_clu_est_ind3=c(res_clu_est_ind3,ipw_clu_ind_est3)
      res_clu_std_ind3=c(res_clu_std_ind3,ipw_clu_ind_std3)
      if(ipw_clu_ind3$value$converged==0){res_clu_warn_ind3=c(res_clu_warn_ind3,time)}
      if(ipw_clu_ind3$value$converged==1){res_clu_warn_ind3=c(res_clu_warn_ind3,0)}
    }
    if(is.null(ipw_clu_ex3$value)==0){
      ipw_clu_ex_est3=summary(ipw_clu_ex3$value)$beta[3]
      ipw_clu_ex_std3=summary(ipw_clu_ex3$value)$se.robust[3]
      res_clu_est_ex3=c(res_clu_est_ex3,ipw_clu_ex_est3)
      res_clu_std_ex3=c(res_clu_std_ex3,ipw_clu_ex_std3)
      if(ipw_clu_ex3$value$converged==0){res_clu_warn_ex3=c(res_clu_warn_ex3,time)}
      if(ipw_clu_ex3$value$converged==1){res_clu_warn_ex3=c(res_clu_warn_ex3,0)}
    }
    
    if(is.null(ipw_ind3$value)==1){
      res_est_ind3=c(res_est_ind3,NA)
      res_std_ind3=c(res_std_ind3,NA)
      res_warn_ind3=c(res_warn_ind3,time)
    }
    if(is.null(ipw_ex3$value)==1){
      res_est_ex3=c(res_est_ex3,NA)
      res_std_ex3=c(res_std_ex3,NA)
      res_warn_ex3=c(res_warn_ex3,time)
    }
    if(is.null(ipw_clu_ind3$value)==1){
      res_clu_est_ind3=c(res_clu_est_ind3,NA)
      res_clu_std_ind3=c(res_clu_std_ind3,NA)
      res_clu_warn_ind3=c(res_clu_warn_ind3,time)
    }
    if(is.null(ipw_clu_ex3$value)==1){
      res_clu_est_ex3=c(res_clu_est_ex3,NA)
      res_clu_std_ex3=c(res_clu_std_ex3,NA)
      res_clu_warn_ex3=c(res_clu_warn_ex3,time)
    }
    
  }
  
  EST1=data.frame(true_est_ind,true_est_ex,res_est_ind1,res_est_ex1,res_clu_est_ind1,res_clu_est_ex1,icc,TIMEs)
  STD1=data.frame(true_std_ind,true_std_ex,res_std_ind1,res_std_ex1,res_clu_std_ind1,res_clu_std_ex1,icc,TIMEs)
  WARN1=data.frame(true_warn_ind,true_warn_ex,res_warn_ind1,res_warn_ex1,res_clu_warn_ind1,res_clu_warn_ex1,icc,TIMEs)

  EST2=data.frame(true_est_ind,true_est_ex,res_est_ind2,res_est_ex2,res_clu_est_ind2,res_clu_est_ex2,icc,TIMEs)
  STD2=data.frame(true_std_ind,true_std_ex,res_std_ind2,res_std_ex2,res_clu_std_ind2,res_clu_std_ex2,icc,TIMEs)
  WARN2=data.frame(true_warn_ind,true_warn_ex,res_warn_ind2,res_warn_ex2,res_clu_warn_ind2,res_clu_warn_ex2,icc,TIMEs)
  
  EST3=data.frame(true_est_ind,true_est_ex,res_est_ind3,res_est_ex3,res_clu_est_ind3,res_clu_est_ex3,icc,TIMEs)
  STD3=data.frame(true_std_ind,true_std_ex,res_std_ind3,res_std_ex3,res_clu_std_ind3,res_clu_std_ex3,icc,TIMEs)
  WARN3=data.frame(true_warn_ind,true_warn_ex,res_warn_ind3,res_warn_ex3,res_clu_warn_ind3,res_clu_warn_ex3,icc,TIMEs)
  
  mp=c(mp,intercept)
  names(EST1)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')
  names(STD1)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')
  names(WARN1)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')

  names(EST2)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')
  names(STD2)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')
  names(WARN2)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')

  names(EST3)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')
  names(STD3)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')
  names(WARN3)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')
  
  result=list(EST1=EST1,STD1=STD1,WARN1=WARN1,
              EST2=EST2,STD2=STD2,WARN2=WARN2,
              EST3=EST3,STD3=STD3,WARN3=WARN3,mp=mp)
  name_res=paste('ipw_result_stand1',icc,vars,'.RData',sep='')
  save(result,file=name_res)
}


### add expression 4
vars=0.2
################################ simulations
################################

for(icc in seq(0.1,0.9,0.1)){
  print('icc')
  print(icc)
  n=50
  m=50
  
  res_est_ind3=c();res_std_ind3=c();res_warn_ind3=c()
  res_est_ex3=c();res_std_ex3=c();res_warn_ex3=c()
  res_clu_est_ind3=c();res_clu_std_ind3=c();res_clu_warn_ind3=c()
  res_clu_est_ex3=c();res_clu_std_ex3=c();res_clu_warn_ex3=c()
  
  ICCs=c()
  TIMEs=c()
  
  missings=c()
  
  true_est_ex=c();true_std_ex=c();true_warn_ex=c()
  true_est_ind=c();true_std_ind=c();true_warn_ind=c()
  
  est=c();std=c();warn=c()
  
  #intercept=temp2[temp2$icc==icc,]$intercept
  #mp=mis_per_check(n=50,m=50,icc=icc,seed=time,intercept=intercept,vars=0.2)
  
  intercept=mis_per(n=50,m=50,icc=icc,vars=vars)
  print(paste('intercept',intercept))
  mp=mis_per_check(vars=vars,n=50,m=50,icc=icc,intercept=intercept)
  print(mp)
  
  for(times in 1:500){
    print(paste('intercept',intercept,'mp',mp))
    print(paste('icc',icc,'times',times))
    
    TIMEs=c(TIMEs,times)
    
    d1=data_gene5(k=50,m=50,icc=icc,intercept=intercept,seed=times,vars=vars)
    d2=d1;d2$y=ifelse(d2$R==1,NA,d2$y)
    d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=d2$R)
    print(sum(is.na(d3$y))/dim(d3)[1])
    missings=c(missings,sum(is.na(d3$y))/dim(d3)[1])
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
    
    ipw_ind4=myTryCatch(geeDREstimation(formula=y~x+arm,
                                        nameTRT = "arm",nameMISS = "missing", nameY = "y",
                                        id="cluster" , data = d3,
                                        family =  binomial("logit"), corstr = "independence",
                                        sandwich.nuisance=TRUE,
                                        model.weights=I(missing==1)~arm+x)) 
    
    ipw_clu_ind4=myTryCatch(geeDREstimation(formula=y~x+arm,
                                            id="cluster" , data = d3,
                                            nameMISS='missing',nameY='y',
                                            nameTRT='arm',
                                            sandwich.nuisance=TRUE,
                                            model.weights=I(missing==1)~arm+x,
                                            family =  binomial("logit"),
                                            corstr = "independence")) 
    
    # exchangeable
    ipw_ex4=myTryCatch(geeDREstimation(formula=y~x+arm,
                                       nameTRT = "arm",nameMISS = "missing", nameY = "y",
                                       id="cluster" , data = d3,
                                       family =  binomial("logit"), 
                                       corstr = "exchangeable",
                                       sandwich.nuisance=TRUE,
                                       model.weights=I(missing==1)~arm+x)) # missing = 0, observed
    
    ipw_clu_ex4=myTryCatch(geeDREstimation(formula=y~x+arm,
                                           id="cluster" , data = d3,
                                           nameMISS='missing',nameY='y',
                                           nameTRT='arm',
                                           sandwich.nuisance=TRUE,
                                           model.weights=I(missing==1)~arm+x,
                                           family =  binomial("logit"),
                                           corstr = "exchangeable"))
    
    
    if(is.null(ipw_ind4$value)==0){
      ipw_ind_est3=summary(ipw_ind4$value)$beta[3]
      ipw_ind_std3=summary(ipw_ind4$value)$se.robust[3]
      res_est_ind3=c(res_est_ind3,ipw_ind_est3)
      res_std_ind3=c(res_std_ind3,ipw_ind_std3)
      if(ipw_ind4$value$converged==0){res_warn_ind3=c(res_warn_ind3,times)}
      if(ipw_ind4$value$converged==1){res_warn_ind3=c(res_warn_ind3,0)}
    }
    if(is.null(ipw_ex4$value)==0){
      ipw_ex_est3=summary(ipw_ex4$value)$beta[3]
      ipw_ex_std3=summary(ipw_ex4$value)$se.robust[3]
      res_est_ex3=c(res_est_ex3,ipw_ex_est3)
      res_std_ex3=c(res_std_ex3,ipw_ex_std3)
      if(ipw_ex4$value$converged==0){res_warn_ex3=c(res_warn_ex3,times)}
      if(ipw_ex4$value$converged==1){res_warn_ex3=c(res_warn_ex3,0)}
    }
    if(is.null(ipw_clu_ind4$value)==0){
      ipw_clu_ind_est3=summary(ipw_clu_ind4$value)$beta[3]
      ipw_clu_ind_std3=summary(ipw_clu_ind4$value)$se.robust[3]
      res_clu_est_ind3=c(res_clu_est_ind3,ipw_clu_ind_est3)
      res_clu_std_ind3=c(res_clu_std_ind3,ipw_clu_ind_std3)
      if(ipw_clu_ind4$value$converged==0){res_clu_warn_ind3=c(res_clu_warn_ind3,times)}
      if(ipw_clu_ind4$value$converged==1){res_clu_warn_ind3=c(res_clu_warn_ind3,0)}
    }
    if(is.null(ipw_clu_ex4$value)==0){
      ipw_clu_ex_est3=summary(ipw_clu_ex4$value)$beta[3]
      ipw_clu_ex_std3=summary(ipw_clu_ex4$value)$se.robust[3]
      res_clu_est_ex3=c(res_clu_est_ex3,ipw_clu_ex_est3)
      res_clu_std_ex3=c(res_clu_std_ex3,ipw_clu_ex_std3)
      if(ipw_clu_ex4$value$converged==0){res_clu_warn_ex3=c(res_clu_warn_ex3,times)}
      if(ipw_clu_ex4$value$converged==1){res_clu_warn_ex3=c(res_clu_warn_ex3,0)}
    }
    
    if(is.null(ipw_ind4$value)==1){
      res_est_ind3=c(res_est_ind3,NA)
      res_std_ind3=c(res_std_ind3,NA)
      res_warn_ind3=c(res_warn_ind3,time)
    }
    if(is.null(ipw_ex4$value)==1){
      res_est_ex3=c(res_est_ex3,NA)
      res_std_ex3=c(res_std_ex3,NA)
      res_warn_ex3=c(res_warn_ex3,time)
    }
    if(is.null(ipw_clu_ind4$value)==1){
      res_clu_est_ind3=c(res_clu_est_ind3,NA)
      res_clu_std_ind3=c(res_clu_std_ind3,NA)
      res_clu_warn_ind3=c(res_clu_warn_ind3,time)
    }
    if(is.null(ipw_clu_ex4$value)==1){
      res_clu_est_ex3=c(res_clu_est_ex3,NA)
      res_clu_std_ex3=c(res_clu_std_ex3,NA)
      res_clu_warn_ex3=c(res_clu_warn_ex3,time)
    }
    
  }
  
  EST3=data.frame(true_est_ind,true_est_ex,res_est_ind3,res_est_ex3,res_clu_est_ind3,res_clu_est_ex3,icc,TIMEs)
  STD3=data.frame(true_std_ind,true_std_ex,res_std_ind3,res_std_ex3,res_clu_std_ind3,res_clu_std_ex3,icc,TIMEs)
  WARN3=data.frame(true_warn_ind,true_warn_ex,res_warn_ind3,res_warn_ex3,res_clu_warn_ind3,res_clu_warn_ex3,icc,TIMEs)
  
  mp=c(mp,intercept)
  names(EST3)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')
  names(STD3)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')
  names(WARN3)=c('true_ind','true_ex','ipw_ind','ipw_ex','ipw_clu_ind','ipw_ex_ind','icc','time')
  
  result=list(EST3=EST3,STD3=STD3,WARN3=WARN3,mp=mp)
  name_res=paste('ipw_result_stand1_3',icc,vars,'.RData',sep='')
  save(result,file=name_res)
}



setwd('c://R//ipw//step8')

### read in the results and draw tables
setwd('')
icc=0.1
vars=0.2

setwd("c:/R/ipw/step8/r")
table_ind1=c()
table_exch1=c()
for(icc in seq(0.1,0.9,0.1)){
  names=paste('ipw_result_stand11',icc,vars,'.RData',sep='')
  load(names)
  print(result$mp)
  mean(result$missings)
  names(est)[6]='ipw_clu_ex'
  est=result$EST1
  std=result$STD1
  warn=result$WARN1
  apply(warn[,1:6],2,sum)
  print(mean(result$missings))
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
  non_con=(200-dim(est_n)[1])/200
  cova=((est_n-1.96*std_n)<true_ind) & ((est_n+1.96*std_n)>true_ind)
  cov_ind=apply(cova,2,sum)/dim(est_n)[1]
  cova=((est_n-1.96*std_n)<true_ex) & ((est_n+1.96*std_n)>true_ex)
  cov_ex=apply(cova,2,sum)/dim(est_n)[1]
  
  ind=c(mean_value[c(3,5)]-true_ind,mcsd_value[c(3,5)],std_value[c(3,5)],cov_ind[c(3,5)],non_con)
  ind=round(ind,3)
  table_ind1=rbind(table_ind1,ind)
  exch=c(mean_value[c(4,6)]-true_ex,mcsd_value[c(4,6)],std_value[c(4,6)],cov_ind[c(4,6)],non_con)
  exch=round(exch,3)
  table_exch1=rbind(table_exch1,exch)
}

table_ind2=c()
table_exch2=c()
for(icc in seq(0.1,0.6,0.1)){
  names=paste('ipw_result_stand11',icc,vars,'.RData',sep='')
  load(names)
  print(result$mp)
  est=result$EST2
  std=result$STD2
  warn=result$WARN2
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
  non_con=(500-dim(est_n)[1])/500
  cova=((est_n-1.96*std_n)<true_ind) & ((est_n+1.96*std_n)>true_ind)
  cov_ind=apply(cova,2,sum)/dim(est_n)[1]
  cova=((est_n-1.96*std_n)<true_ex) & ((est_n+1.96*std_n)>true_ex)
  cov_ex=apply(cova,2,sum)/dim(est_n)[1]
  
  ind=c(mean_value[c(3,5)]-true_ind,mcsd_value[c(3,5)],std_value[c(3,5)],cov_ind[c(3,5)],non_con)
  ind=round(ind,3)
  table_ind2=rbind(table_ind2,ind)
  exch=c(mean_value[c(4,6)]-true_ex,mcsd_value[c(4,6)],std_value[c(4,6)],cov_ind[c(4,6)],non_con)
  exch=round(exch,3)
  table_exch2=rbind(table_exch2,exch)
}

table_ind3=c()
table_exch3=c()
for(icc in seq(0.1,0.9,0.1)){
  names=paste('ipw_result_stand11',icc,vars,'.RData',sep='')
  load(names)
  print(result$mp)
  est=result$EST3
  std=result$STD3
  warn=result$WARN3
  apply(warn[,1:6],2,sum)
  print(mean(result$missings))
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
  non_con=(500-dim(est_n)[1])/500
  cova=((est_n-1.96*std_n)<true_ind) & ((est_n+1.96*std_n)>true_ind)
  cov_ind=apply(cova,2,sum)/dim(est_n)[1]
  cova=((est_n-1.96*std_n)<true_ex) & ((est_n+1.96*std_n)>true_ex)
  cov_ex=apply(cova,2,sum)/dim(est_n)[1]
  
  ind=c(mean_value[c(3,5)]-true_ind,mcsd_value[c(3,5)],std_value[c(3,5)],cov_ind[c(3,5)],non_con)
  ind=round(ind,3)
  table_ind3=rbind(table_ind3,ind)
  exch=c(mean_value[c(4,6)]-true_ex,mcsd_value[c(4,6)],std_value[c(4,6)],cov_ind[c(4,6)],non_con)
  exch=round(exch,3)
  table_exch3=rbind(table_exch3,exch)
}

table_ind4=c()
table_exch4=c()
for(icc in seq(0.1,0.9,0.1)){
  names=paste('ipw_result_stand1_3',icc,vars,'.RData',sep='')
  load(names)
  print(result$mp)
  est=result$EST3
  std=result$STD3
  warn=result$WARN3
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
  non_con=(500-dim(est_n)[1])/500
  cova=((est_n-1.96*std_n)<true_ind) & ((est_n+1.96*std_n)>true_ind)
  cov_ind=apply(cova,2,sum)/dim(est_n)[1]
  cova=((est_n-1.96*std_n)<true_ex) & ((est_n+1.96*std_n)>true_ex)
  cov_ex=apply(cova,2,sum)/dim(est_n)[1]
  
  ind=c(mean_value[c(3,5)]-true_ind,mcsd_value[c(3,5)],std_value[c(3,5)],cov_ind[c(3,5)],non_con)
  ind=round(ind,3)
  table_ind4=rbind(table_ind4,ind)
  exch=c(mean_value[c(4,6)]-true_ex,mcsd_value[c(4,6)],std_value[c(4,6)],cov_ind[c(4,6)],non_con)
  exch=round(exch,3)
  table_exch4=rbind(table_exch4,exch)
}




##errors:
icc=0.9;intercept=-4;times=500;vars=0.2
d1=data_gene5(k=50,m=50,icc=icc,intercept=intercept,seed=times,vars=vars)
d2=d1;d2$y=ifelse(d2$R==1,NA,d2$y)
d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=d2$R)
head(d3)
myTryCatch(geeDREstimation(formula=y~arm+x,
                nameTRT = "arm",nameMISS = "missing", nameY = "y",
                id="cluster" , data = d3,
                family =  binomial("logit"), corstr = "independence",
                sandwich.nuisance=TRUE,
                model.weights=I(missing==0)~arm+x))

icc=0.1;intercept=-1.4;times=1;vars=0.2
d1=data_gene5(k=50,m=50,icc=icc,intercept=intercept,seed=times,vars=vars)
d2=d1;d2$y=ifelse(d2$R==1,NA,d2$y)
d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=d2$R)
head(d3)
geeDREstimation(formula=y~arm+x,
                nameTRT = "arm",nameMISS = "missing", nameY = "y",
                id="cluster" , data = d3,
                family =  binomial("logit"), corstr = "independence",
                sandwich.nuisance=TRUE,
                model.weights=I(missing==0)~arm+x)


icc=0.1;intercept=-1.4;vars=0.2;times=1
d1=data_gene5(k=50,m=50,icc=icc,intercept=intercept,seed=times,vars=vars)
d2=d1;d2$y=ifelse(d2$R==1,NA,d2$y)
d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=d2$R)

logs=glm(missing ~ x+arm, data = d3,
         family = binomial(link='logit'))

weight1=1/predict(logs,type="response")
d3$weight=weight1

geeDREstimation(formula=y~x+arm,
                id="cluster" , data = d3,
                nameMISS='missing',nameY='y',
                nameTRT='arm',
                sandwich.nuisance=TRUE,
                model.weights=I(missing==0)~arm+x,
                family =  binomial("logit"),
                corstr = "independence")

geeDREstimation(formula=y~x+arm,
                id="cluster" , data = d3,
                nameMISS='missing',nameY='y',
                nameTRT='arm',
                sandwich.nuisance=TRUE,
                model.weights=I(missing==1)~arm+x,
                family =  binomial("logit"),
                corstr = "independence")
