cca_r=c()
ccastd_r=c()
mmi_r=c()
mmistd_r=c()
ipw_r=c()
ipwstd_r=c()
ipw_r2=c()
ipw_r2_std=c()
ipw_r3=c()
ipw_r3_std=c()
cov_cca=c()
cov_ipw1=c()
cov_ipw2=c()
cov_ipw3=c()
cov_mmi=c()

ipw_mean=c()
ipw_std=c()
ipw_mean2=c()
ipw_std2=c()

warning_ipw1=c()
warning_ipw2_glmer=c()
warning_ipw2_gee=c()
warning_mmi=c()
warning_mmi_result=c()

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

t=1000
eff0=c(1.326051, 1.318658, 1.319951, 1.318619)

for(km in 1:2){
  if(km==1){
    k=25
    m=25
    load("ipw.RData")}
  if(km==2){
    k=50
    m=25
    load("ipw2.RData")}
  if(km==3){
    k=25
    m=50
    load("ipw3.RData")}
  if(km==4){
    k=50
    m=50
    load("ipw4.RData")}
  
  eff=eff0[km]
  
  cca_mean=final_data[[1]]
  cca_std=final_data[[2]]
  ipw_mean=final_data[[3]]
  ipw_std=final_data[[4]]
  ipw_mean2=final_data[[5]]
  ipw_std2=final_data[[6]]
  jomo_mean=final_data[[7]]
  jomo_std=final_data[[8]]
  missingper=final_data[[9]]
  error_cca=final_data[[10]]
  warning_cca=final_data[[11]]
  error_ipw1=final_data[[12]]
  warning_ipw1=final_data[[13]]
  error_ipw2_glmer=final_data[[14]]
  warning_ipw2_glmer=final_data[[15]]
  error_ipw2_gee=final_data[[16]]
  warning_ipw2_gee=final_data[[17]]
  error_mmi=final_data[[18]]
  warning_mmi=final_data[[19]]
  error_mmi_result=final_data[[20]]
  warning_mmi_result=final_data[[21]]
  
  
  ### check errors:
  length(error_cca)
  length(error_ipw1)
  length(error_ipw2_glmer)
  length(error_ipw2_gee)
  length(error_mmi)
  length(error_mmi_result)
  
  ### check warnings:
  length(warning_cca)
  length(warning_ipw1) ##112
  length(warning_ipw2_glmer)
  length(warning_ipw2_gee) ##118
  length(warning_mmi)
  length(warning_mmi_result)
  
  ### average mean:
  mean(missingper)
  
  cca=mean(cca_mean)
  ccastd=mean(cca_std)
  cca_r=c(cca_r,cca)
  ccastd_r=c(ccastd_r,ccastd)
  ### converage:
  converge=function(mean,std){
    a=sum((mean-1.96*std)<eff & (mean+1.96*std) >eff)
    return(a/1000)
  }
  
  cca_co2525=converge(cca_mean,cca_std)
  cov_cca=c(cov_cca,cca_co2525)
  
  mmi=mean(jomo_mean)
  mmistd=mean(jomo_std)
  mmi_r=c(mmi_r,mmi)
  mmistd_r=c(mmistd_r,mmistd)
  mmi_co2525=converge(jomo_mean,jomo_std)
  cov_mmi=c(cov_mmi,mmi_co2525)
  
  ### ipw
  warn_time1=final_data['warn_time1']
  warn_ipw1=final_data['warn_ipw1']
  warn_time2=final_data['warn_ipw2']
  warn_ipw2=final_data['warn_ipw2']
  
  warns=c(warn_time1,warn_time2)
  warns=unique(warns)
  length(warns)
  ipw_mean_new=ipw_mean[setdiff(c(1:t),warns)]
  ipw=mean(ipw_mean_new)
  ipw_std_new=ipw_std[setdiff(c(1:t),warns)]
  ipwstd=mean(ipw_std_new)
  
  ipw_r=c(ipw_r,ipw)
  ipwstd_r=c(ipwstd_r,ipwstd)
  
  ipw1_co2525=converge(ipw_mean_new,ipw_std_new)
  cov_ipw1=c(cov_ipw1,ipw1_co2525)
  ## IPW with cluster 
  ipw_mean2_new=ipw_mean2[setdiff(c(1:t),warns)]
  ipw_cluster=mean(ipw_mean2_new)
  ipw_std2_new=ipw_std2[setdiff(c(1:t),warns)]
  ipw_std_cluster=mean(ipw_std2_new)
  
  ipw_r2=c(ipw_r2,ipw_cluster)
  ipw_r2_std=c(ipw_r2_std,ipw_std_cluster)
  
  ipw2_co2525=converge(ipw_mean2_new,ipw_std2_new)
  cov_ipw2=c(cov_ipw2,ipw2_co2525)
  
  ### IPW-GEE with package CRTgeeDR
  ipw_mm=final_data['ipw_mm']
  ipw_mm_sd=final_data['ipw_mm_sd']
  ipw_r3=c(ipw_r3,mean(ipw_mm))
  ipw_r3_std=c(ipw_r3_std,mean(ipw_mm_sd))
  cov_ipw3=c(cov_ipw3,converge(ipw_mm,ipw_mm_sd))
}

