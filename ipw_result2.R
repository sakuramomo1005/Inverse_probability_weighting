####### compare geepack and CRTgeeDR
k=25;m=25;tt=1
compare_data=data_gene(k,m,tt)
head(compare_data,3)
## using CRTgeeDR package 
result1=geeDREstimation(formula=y~x+arm,
                        id="cluster" , data = compare_data,
                        nameMISS='missing',nameY='y',
                        nameTRT='arm',
                        family=binomial("logit"), 
                        corstr = "independence",
                        model.weights=I(R==1)~x)
summary(result1)
## method 2
print("calculate the weight")
w1=glm(R ~ x , data = compare_data, 
       family = binomial(link='logit'))
weight=expit(predict(w1))
compare_data$weight=round(1/weight)
compare_data=na.omit(compare_data)
print("fit the gee model")
result2=geese(formula=y~x+arm,data=compare_data,id=cluster,
              family = binomial(link='logit'),
              weights = weight,
              corstr = 'independence')
summary(result2)


library(jomo)
library(geepack)
library(lme4)
library(jomo)

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
one_group=function(i,k,mm,seed=123){
  set.seed(seed)
  ## parameters
  b0=1;b1=1.36;b2=1
  sigma_b=sqrt(0.2)
  mu_x=0
  sigma_alpha=sqrt(0.18)
  sigma_u=sqrt(3.37)
  psi=-1.34
  phi=1
  x=matrix(0,k,mm+5)
  y=matrix(0,k,mm+5)
  pi=matrix(0,k,mm+5)
  cluster=matrix(0,k,mm+5)
  r=matrix(0,k,mm+5)
  for(k in 1:k){
    delta=rnorm(1,0,sigma_b)
    alpha=rnorm(1,mu_x,sigma_alpha)
    m=round(runif(1,mm-5,mm+5))
   # print(m)
    cluster[k,]=k
    for(j in 1:m){
      u=rnorm(1,0,sigma_u)
      x[k,j]=u+alpha
      pi[k,j]=expit(b0+b1*i+b2*x[k,j]+delta)
      y[k,j]=rbinom(1,1,pi[k,j])
    }
  }
  r=expit(psi+phi*x)
  return(list(x=x,y=y,pi=pi,r=r,cluster=cluster))
}
data_gene=function(k,m,seed=123){
  set.seed(seed)
  # step1:
  a1=one_group(1,k,m,seed)
  a0=one_group(0,k,m,seed)
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

true_est=c();true_warning=c()
uncra_est=c();uncra_std=c();uncra_warning=c()
cra_est=c();cra_std=c();cra_warning=c()
ipwgeeno_est=c();ipwgeeno_std=c();ipwgeeno_warn=c()
ipwcrtgeedr_est=c();ipwcrtgeedr_std=c();ipwcrtgee_warn=c()
ipwgeeyes_est=c();ipwgeeyes_std=c();ipwgeeyes_warn=c()
ipwcrtgeedr_yes_est=c();ipwcrtgeedr_yes_std=c();ipwcrtgee_yes_warn=c()



times=1000
for(km in 1:4){
  if(km==1){k=25;m=25}
  if(km==2){k=50;m=25}
  if(km==3){k=25;m=50}
  if(km==4){k=50;m=50}
  print("***************")
  print("km")
  print(km)
  print('times')
  for(tt in 1:times){
    print(tt)
    d1=data_gene(k,m,tt)
    d2=d1
    d2$y=ifelse(d2$R==1,NA,d2$y)
    d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=d2$R)
    #d3$y=as.factor(d3$y)
    ### True values
    formu=formula(y~x+arm)
    true=myTryCatch(geese(formu,data=d1,id=cluster,
                          family = binomial(link='logit'),
                          corstr = 'independence'))
    true_est=c(true_est,summary(true$value)$mean['arm','estimate'])
    if(is.null(true$warning)==0){
      true_warning=c(true_warning,tt)}
    
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
    
    ### adjusted CRA
    cra=myTryCatch(geese(formu,data=d_nona,id=cluster,
                           family = binomial(link='logit'),
                           corstr = 'independence'))
    cra_est=c(cra_est,summary(cra$value)$mean['arm','estimate'])
    cra_std=c(cra_std,summary(cra$value)$mean['arm','san.se'])
    if(is.null(cra$warning)==0){
      cra_warning=c(cra_warning,tt)}
    
    ### IPW-no cluster effect, by using glm and geese
    logs=glm(missing ~ x , data = d3, 
             family = binomial(link='logit'))
    logsum=summary(logs)
    weight=expit(predict(logs))
    d33=d3
    d33$weight=round(1/weight)
    d4=na.omit(d33)
    ipwgeeno=myTryCatch(geese(formu,data=d4,id=cluster,
             family = binomial(link='logit'),
             weights = weight,
             corstr = 'independence'))
    tempgeeno_est=summary(ipwgeeno$value)$mean['arm','estimate']
    tempgeeno_std=summary(ipwgeeno$value)$mean['arm','san.se']
    ipwgeeno_est=c(ipwgeeno_est,tempgeeno_est)
    ipwgeeno_std=c(ipwgeeno_std,tempgeeno_std)
    if(is.null(ipwgeeno$warning)==0){ipwgeeno_warn=c(ipwgeeno_warn,tt)}
    
    ### IPW-GEE with package CRTgeeDR
    ipwcrtgeedr=myTryCatch(geeDREstimation(formula=y~x+arm,
                                id="cluster" , data = d3,
                                nameMISS='missing',nameY='y',
                                nameTRT='arm',
                                family =  binomial("logit"), 
                                corstr = "independence",
                                model.weights=I(missing==1)~x))
    tempcrtgee_est=summary(ipwcrtgeedr$value)$beta[3]
    tempcrtgee_std=summary(ipwcrtgeedr$value)$se.robust[3]
    ipwcrtgeedr_est=c(ipwcrtgeedr_est,tempcrtgee_est)
    ipwcrtgeedr_std=c(ipwcrtgeedr_std,tempcrtgee_std)
    if(is.null(ipwcrtgeedr$warning)==0){ipwcrtgee_warn=c(ipwcrtgee_warn,tt)}
  
    ### IPW-WITH cluster effect, by using glm and geese
    logs=glmer(missing ~ x+(1|cluster) , data = d3, 
             family = binomial(link='logit'))
    logsum=summary(logs)
    weight=expit(predict(logs))
    d3$weight=round(1/weight)
    d4=na.omit(d3)
    formu=formula(y~x+arm)
    ipwgeeyes=myTryCatch(geese(formu,data=d4,id=cluster,
                              family = binomial(link='logit'),
                              weights = weight,
                              corstr = 'independence'))
    tempgeeyes_est=summary(ipwgeeyes$value)$mean['arm','estimate']
    tempgeeyes_std=summary(ipwgeeyes$value)$mean['arm','san.se']
    ipwgeeyes_est=c(ipwgeeyes_est,tempgeeyes_est)
    ipwgeeyes_std=c(ipwgeeyes_std,tempgeeyes_std)
    if(is.null(ipwgeeyes$warning)==0){ipwgeeyes_warn=c(ipwgeeyes_warn,tt)}
    
    ### IPW-GEE with cluster effect with package CRTgeeDR
    ipwcrtgeedr_yes=myTryCatch(geeDREstimation(formula=y~x+arm,
                                           id="cluster" , data = d3,
                                           nameMISS='missing',nameY='y',
                                           nameTRT='arm',
                                           family =  binomial("logit"), 
                                           corstr = "independence",
                                           model.weights=I(missing==1)~x+(1|cluster)))
    tempcrtgee_yes_est=summary(ipwcrtgeedr_yes$value)$beta[3]
    tempcrtgee_yes_std=summary(ipwcrtgeedr_yes$value)$se.robust[3]
    ipwcrtgeedr_yes_est=c(ipwcrtgeedr_yes_est,tempcrtgee_yes_est)
    ipwcrtgeedr_yes_std=c(ipwcrtgeedr_yes_std,tempcrtgee_yes_std)
    if(is.null(ipwcrtgeedr_yes$warning)==0){ipwcrtgee_yes_warn=c(ipwcrtgee_yes_warn,tt)}
    }}


##  deal with the results
#(ipwgeeno_warn)
#[1]  12  13  14  51  58  62  63  71 108 115 131 135 139 194 225 233 283
#[18] 288 309 314 339 346 378 383 390 413 415 430 475 502 504 517 519 524
#[35] 541 572 581 601 633 646 691 720 729 731 750 751 771 774 802 803 829
#[52] 838 865 885 927 933 945  13  27  37  58  71 108 135 208 233 250 288
#[69] 295 309 332 356 390 415 430 485 489 504 517 541 601 624 633 711 720
#[86] 731 774 786 803 829 927 933 987   6  32  87 153 221 254 289 390 407
#[103] 415 427 436 442 453 457 488 496 517 546 662 731 761 767 780 821 829
#[120] 835 899 906 913 920 921 929 944 962  42 121 177 185 221 241 390 401
#[137] 403 407 415 427 434 453 484 485 549 662 731 761 808 920 923 949 961

## UNCOVERGENCE
# round 1:
uncon1=c(12,13,14,51,58,62,63,71,108,115,131,135,139,194,225,233,283,
        288,309,314,339,346,378,383,390,413,415,430,475,502,504,517,519,524,
        541,572,581,601,633,646,691,720,729,731,750,751,771,774,802,803,829,
        838,865,885,927,933,945)
uncon2=c(13,27,37,58,71,108,135,208,233,250,288,
         295,309,332,356,390,415,430,485,489,504,517,541,601,624,633,711,720,
         731,774,786,803,829,927,933,987)
uncon3=c(6,32,87,153,221,254,289,390,407,
         415,427,436,442,453,457,488,496,517,546,662,731,761,767,780,821,829,
         835,899,906,913,920,921,929,944,962)
uncon4=c(42,121,177,185,221,241,390,401,
         403,407,415,427,434,453,484,485,549,662,731,761,808,920,923,949,961)


uncon2=uncon2+1000
uncon3=uncon3+2000
uncon4=uncon4+3000

unconyes1=c(12,13,14,51,58,62,63,71,87,108,115,131,135,139,194,225,233,
            283,288,309,314,339,346,378,383,390,413,415,430,475,502,504,517,519,
            524,541,546,547,572,581,601,633,646,691,715,720,729,731,750,751,755,
            771,774,802,803,829,838,865,885,927,933,937,945)
unconyes2=c(13,27,37,58,71,
            108,135,208,233,250,288,295,309,332,356,390,415,430,485,489,504,517,
            541,601,624,633,711,720,731,774,786,803,829,927,933,987)
unconyes3=c(6,32,87,
            153,221,254,289,390,407,415,427,436,442,453,457,488,496,517,546,662,
            731,761,767,780,821,829,835,899,906,913,920,921,929,944,962)
unconyes4=c(42,121,
            177,185,221,241,390,401,403,407,415,427,434,453,484,485,549,662,731,
            761,808,920,923,949,961)

unconyes2=unconyes2+1000
unconyes3=unconyes3+2000
unconyes4=unconyes4+3000

ipwgeeno_est[c(uncon1,uncon2,uncon3,uncon4)]
ipwgeeyes_est[c(unconyes1,unconyes2,unconyes3,unconyes4)]

## delete the unconverge
mean(ipwgeeno_est[-c(uncon1,uncon2,uncon3,uncon4)])

## 25 25
length(uncon1)
length(unconyes1)
true_est1=true_est[1:1000]
uncra_est1=uncra_est[1:1000]
cra_est1=cra_est[1:1000]
ipwgeeno_est1=ipwgeeno_est[1:1000]
ipwgeeno_est1=ipwgeeno_est1[-c(uncon1)]
ipwcrtgeedr_est1=ipwcrtgeedr_est[1:1000]
ipwgeeyes_est1=ipwgeeyes_est[1:1000]
ipwgeeyes_est1=ipwgeeyes_est1[-c(unconyes1)]
ipwcrtgeedr_yes_est1=ipwcrtgeedr_yes_est[1:1000]
mean(true_est1)
mean(uncra_est1)
mean(cra_est1)
mean(ipwgeeno_est1)
mean(ipwcrtgeedr_est1)
mean(ipwgeeyes_est1)
mean(ipwcrtgeedr_yes_est1)

uncra_std1=uncra_std[1:1000]
cra_std1=cra_std[1:1000]
ipwgeeno_std1=ipwgeeno_std[1:1000]
ipwgeeno_std1=ipwgeeno_std1[-c(uncon1)]
ipwcrtgeedr_std1=ipwcrtgeedr_std[1:1000]
ipwgeeyes_std1=ipwgeeyes_std[1:1000]
ipwgeeyes_std1=ipwgeeyes_std1[-c(unconyes1)]
ipwcrtgeedr_yes_std1=ipwcrtgeedr_yes_std[1:1000]
for(iii in 1){
  print(mean(uncra_std1))
  print(mean(cra_std1))
  print(mean(ipwgeeno_std1))
  print(mean(ipwcrtgeedr_std1))
  print(mean(ipwgeeyes_std1))
  print(mean(ipwcrtgeedr_yes_std1))
}

## 50 25
length(uncon2)
length(unconyes2)
true_est2=true_est[1001:2000]
uncra_est2=uncra_est[1001:2000]
cra_est2=cra_est[1001:2000]
ipwgeeno_est2=ipwgeeno_est[1001:2000]
ipwgeeno_est2=ipwgeeno_est2[-c(uncon2)]
ipwcrtgeedr_est2=ipwcrtgeedr_est[1001:2000]
ipwgeeyes_est2=ipwgeeyes_est[1001:2000]
ipwgeeyes_est2=ipwgeeyes_est2[-c(unconyes2)]
ipwcrtgeedr_yes_est2=ipwcrtgeedr_yes_est[1001:2000]
mean(true_est2)
mean(uncra_est2)
mean(cra_est2)
mean(ipwgeeno_est2)
mean(ipwcrtgeedr_est2)
mean(ipwgeeyes_est2)
mean(ipwcrtgeedr_yes_est2)

uncra_std2=uncra_std[1001:2000]
cra_std2=cra_std[1001:2000]
ipwgeeno_std2=ipwgeeno_std[1001:2000]
ipwgeeno_std2=ipwgeeno_std2[-c(uncon2)]
ipwcrtgeedr_std2=ipwcrtgeedr_std[1001:2000]
ipwgeeyes_std2=ipwgeeyes_std[1001:2000]
ipwgeeyes_std2=ipwgeeyes_std2[-c(unconyes2)]
ipwcrtgeedr_yes_std2=ipwcrtgeedr_yes_std[1001:2000]
for(iii in 1){
  print(mean(uncra_std2))
  print(mean(cra_std2))
  print(mean(ipwgeeno_std2))
  print(mean(ipwcrtgeedr_std2))
  print(mean(ipwgeeyes_std2))
  print(mean(ipwcrtgeedr_yes_std2))
}

## 25 50
length(uncon3)
length(unconyes3)
true_est3=true_est[2001:3000]
uncra_est3=uncra_est[2001:3000]
cra_est3=cra_est[2001:3000]
ipwgeeno_est3=ipwgeeno_est[2001:3000]
ipwgeeno_est3=ipwgeeno_est3[-c(uncon3)]
ipwcrtgeedr_est3=ipwcrtgeedr_est[2001:3000]
ipwgeeyes_est3=ipwgeeyes_est[2001:3000]
ipwgeeyes_est3=ipwgeeyes_est3[-c(unconyes3)]
ipwcrtgeedr_yes_est3=ipwcrtgeedr_yes_est[2001:3000]
mean(true_est3)
mean(uncra_est3)
mean(cra_est3)
mean(ipwgeeno_est3)
mean(ipwcrtgeedr_est3)
mean(ipwgeeyes_est3)
mean(ipwcrtgeedr_yes_est3)

uncra_std3=uncra_std[2001:3000]
cra_std3=cra_std[2001:3000]
ipwgeeno_std3=ipwgeeno_std[2001:3000]
ipwgeeno_std3=ipwgeeno_std3[-c(uncon3)]
ipwcrtgeedr_std3=ipwcrtgeedr_std[2001:3000]
ipwgeeyes_std3=ipwgeeyes_std[2001:3000]
ipwgeeyes_std3=ipwgeeyes_std3[-c(unconyes3)]
ipwcrtgeedr_yes_std3=ipwcrtgeedr_yes_std[2001:3000]
for(iii in 1){
print(mean(uncra_std3))
print(mean(cra_std3))
print(mean(ipwgeeno_std3))
print(mean(ipwcrtgeedr_std3))
print(mean(ipwgeeyes_std3))
print(mean(ipwcrtgeedr_yes_std3))}

## 50 50
length(uncon4)
length(unconyes4)
true_est4=true_est[3001:4000]
uncra_est4=uncra_est[3001:4000]
cra_est4=cra_est[3001:4000]
ipwgeeno_est4=ipwgeeno_est[3001:4000]
ipwgeeno_est4=ipwgeeno_est4[-c(uncon4)]
ipwcrtgeedr_est4=ipwcrtgeedr_est[3001:4000]
ipwgeeyes_est4=ipwgeeyes_est[3001:4000]
ipwgeeyes_est4=ipwgeeyes_est4[-c(unconyes4)]
ipwcrtgeedr_yes_est4=ipwcrtgeedr_yes_est[3001:4000]
mean(true_est4)
mean(uncra_est4)
mean(cra_est4)
mean(ipwgeeno_est4)
mean(ipwcrtgeedr_est4)
mean(ipwgeeyes_est4)
mean(ipwcrtgeedr_yes_est4)

uncra_std4=uncra_std[3001:4000]
cra_std4=cra_std[3001:4000]
ipwgeeno_std4=ipwgeeno_std[3001:4000]
ipwgeeno_std4=ipwgeeno_std4[-c(uncon4)]
ipwcrtgeedr_std4=ipwcrtgeedr_std[3001:4000]
ipwgeeyes_std4=ipwgeeyes_std[3001:4000]
ipwgeeyes_std4=ipwgeeyes_std4[-c(unconyes4)]
ipwcrtgeedr_yes_std4=ipwcrtgeedr_yes_std[3001:4000]
for(iii in 1){
print(mean(uncra_std4))
print(mean(cra_std4))
print(mean(ipwgeeno_std4))
print(mean(ipwcrtgeedr_std4))
print(mean(ipwgeeyes_std4))
print(mean(ipwcrtgeedr_yes_std4))}

###### converge percent
converge=function(m,std,truth){
  n=length(m)
  covn=0
  nn=sum(((m-std*1.96)<truth) & ((m+1.96*std) >truth))
  return(nn/n)
}

for(iii in 1){
  print(converge(uncra_est1,uncra_std1,mean(true_est1)))
  print(converge(cra_est1,cra_std1,mean(true_est1)))
  print(converge(ipwgeeno_est1, ipwgeeno_std1,mean(true_est1)))
  print(converge(ipwcrtgeedr_est1, ipwcrtgeedr_std1,mean(true_est1)))
  print(converge(ipwgeeyes_est1, ipwgeeyes_std1,mean(true_est1)))
  print(converge(ipwcrtgeedr_est1, ipwcrtgeedr_std1,mean(true_est1)))
}
for(iii in 2){
  print(converge(uncra_est2,uncra_std2,mean(true_est2)))
  print(converge(cra_est2,cra_std2,mean(true_est2)))
  print(converge(ipwgeeno_est2, ipwgeeno_std2,mean(true_est2)))
  print(converge(ipwcrtgeedr_est2, ipwcrtgeedr_std2,mean(true_est2)))
  print(converge(ipwgeeyes_est2, ipwgeeyes_std2,mean(true_est2)))
  print(converge(ipwcrtgeedr_est2, ipwcrtgeedr_std2,mean(true_est2)))
}
for(iii in 3){
  print(converge(uncra_est3,uncra_std3,mean(true_est3)))
  print(converge(cra_est3,cra_std3,mean(true_est3)))
  print(converge(ipwgeeno_est3, ipwgeeno_std3,mean(true_est3)))
  print(converge(ipwcrtgeedr_est3, ipwcrtgeedr_std3,mean(true_est3)))
  print(converge(ipwgeeyes_est3, ipwgeeyes_std3,mean(true_est3)))
  print(converge(ipwcrtgeedr_est3, ipwcrtgeedr_std3,mean(true_est3)))
}
for(iii in 4){
  print(converge(uncra_est4,uncra_std4,mean(true_est4)))
  print(converge(cra_est4,cra_std4,mean(true_est4)))
  print(converge(ipwgeeno_est4, ipwgeeno_std4,mean(true_est4)))
  print(converge(ipwcrtgeedr_est4, ipwcrtgeedr_std4,mean(true_est4)))
  print(converge(ipwgeeyes_est4, ipwgeeyes_std4,mean(true_est4)))
  print(converge(ipwcrtgeedr_est4, ipwcrtgeedr_std4,mean(true_est4)))
}

mean(true_est)
mean(uncra_est)
mean(cra_est)
mean(ipwgeeno_est)
mean(ipwcrtgeedr_est)
mean(ipwgeeyes_est)
mean(ipwcrtgeedr_yes_est)

mean(uncra_std)
mean(cra_std)
mean(ipwgeeno_std)
mean(ipwcrtgeedr_std)
mean(ipwgeeyes_std)
mean(ipwcrtgeedr_yes_std)


### cluster effect
miss_per=c()
means=c()
stds=c()
mins=c()
maxs=c()
times=1000
for(km in 1:4){
  if(km==1){k=25;m=25}
  if(km==2){k=50;m=25}
  if(km==3){k=25;m=50}
  if(km==4){k=50;m=50}
  print("***************")
  print("km")
  print(km)
  print('times')
  for(tt in 1:times){
    print(tt)
    d1=data_gene(k,m,tt)
    miss_per=c(miss_per,sum(d1$R)/(k*m*2))
    mins=c(mins,min(table(d1$cluster)))
    maxs=c(maxs,max(table(d1$cluster)))
    means=c(means,mean(table(d1$cluster)))
    stds=c(stds,sd(table(d1$cluster)))
  }
}

miss_per1=miss_per[1:1000]
mins1=mins[1:1000]
maxs1=maxs[1:1000]
means1=means[1:1000]
stds1=stds[1:1000]

mean(miss_per1)
min(mins1)
max(maxs1)
mean(means1)
mean(stds1)

miss_per2=miss_per[1001:2000]
mins2=mins[1001:2000]
maxs2=maxs[1001:2000]
means2=means[1001:2000]
stds2=stds[1001:2000]

mean(miss_per2)
min(mins2)
mean(maxs2)
mean(means2)
mean(stds2)

miss_per3=miss_per[2001:3000]
mins3=mins[2001:3000]
maxs3=maxs[2001:3000]
means3=means[2001:3000]
stds3=stds[2001:3000]

mean(miss_per3)
min(mins3)
max(maxs3)
mean(means3)
mean(stds3)

miss_per4=miss_per[3001:4000]
mins4=mins[3001:4000]
maxs4=maxs[3001:4000]
means4=means[3001:4000]
stds4=stds[3001:4000]

mean(miss_per4)
min(mins4)
max(maxs4)
mean(means4)
mean(stds4)


#### rerun the code, with the changes of WEIGHT

ipwgeeno_est=c();ipwgeeno_std=c();ipwgeeno_warn=c()
ipwgeeyes_est=c();ipwgeeyes_std=c();ipwgeeyes_warn=c()

times=1000
for(km in 1:4){
  if(km==1){k=25;m=25}
  if(km==2){k=50;m=25}
  if(km==3){k=25;m=50}
  if(km==4){k=50;m=50}
  print("***************")
  print("km")
  print(km)
  print('times')
  for(tt in 1:times){
    print(tt)
    d1=data_gene(k,m,tt)
    d2=d1
    d2$y=ifelse(d2$R==1,NA,d2$y)
    d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=d2$R)
    #d3$y=as.factor(d3$y)
    ### True values
    formu=formula(y~x+arm)
    ### IPW-no cluster effect, by using glm and geese
    logs=glm(missing ~ x , data = d3, 
             family = binomial(link='logit'))
    logsum=summary(logs)
    weight=expit(predict(logs))
    d33=d3
    d33$weight=round(1/weight)
    d33$weight=ifelse(d33$weight>50,50,d33$weight)
    d4=na.omit(d33)
    ipwgeeno=myTryCatch(geese(formu,data=d4,id=cluster,
                              family = binomial(link='logit'),
                              weights = weight,
                              corstr = 'independence'))
    tempgeeno_est=summary(ipwgeeno$value)$mean['arm','estimate']
    tempgeeno_std=summary(ipwgeeno$value)$mean['arm','san.se']
    ipwgeeno_est=c(ipwgeeno_est,tempgeeno_est)
    ipwgeeno_std=c(ipwgeeno_std,tempgeeno_std)
    if(is.null(ipwgeeno$warning)==0){ipwgeeno_warn=c(ipwgeeno_warn,tt)}
    
    ### IPW-WITH cluster effect, by using glm and geese
    logs=glmer(missing ~ x+(1|cluster) , data = d3, 
               family = binomial(link='logit'))
    logsum=summary(logs)
    weight=expit(predict(logs))
    d3$weight=round(1/weight)
    d3$weight=ifelse(d3$weight>50,50,d3$weight)
    d4=na.omit(d3)
    formu=formula(y~x+arm)
    ipwgeeyes=myTryCatch(geese(formu,data=d4,id=cluster,
                               family = binomial(link='logit'),
                               weights = weight,
                               corstr = 'independence'))
    tempgeeyes_est=summary(ipwgeeyes$value)$mean['arm','estimate']
    tempgeeyes_std=summary(ipwgeeyes$value)$mean['arm','san.se']
    ipwgeeyes_est=c(ipwgeeyes_est,tempgeeyes_est)
    ipwgeeyes_std=c(ipwgeeyes_std,tempgeeyes_std)
    if(is.null(ipwgeeyes$warning)==0){ipwgeeyes_warn=c(ipwgeeyes_warn,tt)}
  }}


ipwgeeno_est1=ipwgeeno_est[1:1000]
ipwgeeyes_est1=ipwgeeyes_est[1:1000]
ipwgeeno_std1=ipwgeeno_std[1:1000]
ipwgeeyes_std1=ipwgeeyes_std[1:1000]
for(iii in 1){
  print(mean(ipwgeeno_est1))
  print(mean(ipwgeeyes_est1))
  print(mean(ipwgeeno_std1))
  print(mean(ipwcrtgeedr_std1))
  print(converge(m=ipwgeeno_est1, std=ipwgeeno_std1,mean(true_est1)))
  print(converge(m=ipwgeeyes_est1, std=ipwgeeyes_std1,mean(true_est1)))
}

ipwgeeno_est2=ipwgeeno_est[1001:2000]
ipwgeeyes_est2=ipwgeeyes_est[1001:2000]
ipwgeeno_std2=ipwgeeno_std[1001:2000]
ipwgeeyes_std2=ipwgeeyes_std[1001:2000]
for(iii in 2){
  print(mean(ipwgeeno_est2))
  print(mean(ipwgeeyes_est2))
  print(mean(ipwgeeno_std2))
  print(mean(ipwcrtgeedr_std2))
  print(converge(m=ipwgeeno_est2, std=ipwgeeno_std2,mean(true_est2)))
  print(converge(m=ipwgeeyes_est2, std=ipwgeeyes_std2,mean(true_est2)))
}

ipwgeeno_est3=ipwgeeno_est[2001:3000]
ipwgeeyes_est3=ipwgeeyes_est[2001:3000]
ipwgeeno_std3=ipwgeeno_std[2001:3000]
ipwgeeyes_std3=ipwgeeyes_std[2001:3000]
for(iii in 3){
  print(mean(ipwgeeno_est3))
  print(mean(ipwgeeyes_est3))
  print(mean(ipwgeeno_std3))
  print(mean(ipwcrtgeedr_std3))
  print(converge(m=ipwgeeno_est3, std=ipwgeeno_std3,mean(true_est3)))
  print(converge(m=ipwgeeyes_est3, std=ipwgeeyes_std3,mean(true_est3)))
}

ipwgeeno_est4=ipwgeeno_est[3001:4000]
ipwgeeyes_est4=ipwgeeyes_est[3001:4000]
ipwgeeno_std4=ipwgeeno_std[3001:4000]
ipwgeeyes_std4=ipwgeeyes_std[3001:4000]
for(iii in 4){
  print(mean(ipwgeeno_est4))
  print(mean(ipwgeeyes_est4))
  print(mean(ipwgeeno_std4))
  print(mean(ipwcrtgeedr_std4))
  print(converge(m=ipwgeeno_est4, std=ipwgeeno_std4,mean(true_est4)))
  print(converge(m=ipwgeeyes_est4, std=ipwgeeyes_std4,mean(true_est4)))
}

### the results for MMI 
### I seperated this part since it ran really slowly

result_mmi=function(mmi){
  m0=c()
  std0=c()
  #icc=c()
  for(i in 1:5){
    temp=mmi[mmi$Imputation==i,]
    #temp=temp[-1,]
    temp0=temp[,c('arm','x','cluster','y')]
    formu=formula(y~x+arm)
    temp0$y1=as.numeric(temp0$y)-1
    # temp0$=temp0$Y_mis-1
    mmi2=geese(formu,data=temp0,id=cluster,
               family = binomial(link='logit'),
               corstr = 'independence')
    est=summary(mmi2)
    est_trt=est$mean['arm','estimate']
    sd_trt=est$mean['arm','san.se']
    m0=c(m0,est_trt)
    std0=c(std0,sd_trt)
  }
  return(list(m0,std0))
}

error_mmi_result=c()
warning_mmi_result=c()

jomo_mean=c()
jomo_std=c()

errors_jomo=c()
warn_jomo=c()

times=1000

for(km in 1:1){
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
    d1=data_gene(k,m,times)
    d2=d1
    d2$y=ifelse(d2$R==1,NA,d2$y)
    d3=data.frame(y=d2$y,x=d2$x,cluster=d2$cluster,arm=d2$arm,missing=d2$R)
    #d3$y=as.factor(d3$y)
    d4=jomo(d3,clus = d3$cluster,nimp=15)
    rmmi=myTryCatch(
      result_mmi(mmi))
    error_mmi_result=c(error_mmi_result,rmmi$error)
    warning_mmi_result=c(warning_mmi_result,rmmi$warning)
    if(is.null(rmmi$error)==0){errors_jomo=c(errors_jomo,tt)}
    if(is.null(rmmi$warning)==0){warn_jomo=c(warn_jomo,tt)}
    rmmi=rmmi$value
    m0=rmmi[[1]]
    std0=rmmi[[2]]
    mmi_result=mypool(m0,std0,print='yes')
    mmi_m=mmi_result$mean
    mmi_s=mmi_result$std
    jomo_mean=c(jomo_mean,mmi_m)
    jomo_std=c(jomo_std,mmi_s)
  }}

mean(jomo_mean[1:1000])
mean(jomo_mean[1001:2000])
mean(jomo_mean[2001:3000])
mean(jomo_mean[3001:4000])

mean(jomo_std[1:1000])
mean(jomo_std[1001:2000])
mean(jomo_std[2001:3000])
mean(jomo_std[3001:4000])

###### converge percent
converge=function(m,std,truth){
  n=length(m)
  covn=0
  nn=sum(((m-std*1.96)<truth) & ((m+1.96*std) >truth))
  return(nn/n)
}

converge(jomo_mean[1:1000],jomo_std[1:1000],1.325)
converge(jomo_mean[1001:2000],jomo_std[1001:2000],1.319)
converge(jomo_mean[2001:3000],jomo_std[2001:3000],1.323)
converge(jomo_mean[3001:4000],jomo_std[3001:4000],1.317)