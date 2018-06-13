##### generate the true effects of GEE

eff=c() ## true effects
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
for(times in 1:1000){
  ## full data effect:
  print("************")
  print(times)
  data_sim=data_generation(k,m,times+130,1)
  R=data_sim$R
  X=data_sim$X
  full=cbind(data_sim$X,data_sim$Y$Y)
  names(full)[4]='Y'
  formu=formula(Y~X+group)
  ful_res=geese(formu,data=full,id=cluster,
                family = binomial(link='logit'),
                corstr = 'independence')
  logsum_f=summary(ful_res)
  eff_gee=logsum_f$mean['group','estimate']
  eff=c(eff,eff_gee)
}
}