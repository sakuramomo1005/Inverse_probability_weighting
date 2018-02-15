### generating the dataset

N=200 # cluster number 
M=10  # cluster size

lambda=5
eij=rnorm(1,0,1)

## linkage for model 3, the inverse logit function
h=function(x){
  y=exp(x)/(1+exp(x))
  return(y)
}

## scenario 1
## MCAR
b0=1;b1=0;delta=0
## MAR
b0=0;b1=0.5;delta=0
## Cluster-specific nonignorable 1
b0=0;b1=0.5;delta=5
## cluster-specific nonignorable 2
b0=0;b1=0.5;delta=5

## Model 5:
tau=1
eij=rnorm(1,0,1)
ai=rnorm(1,0,1)
ui=rnorm(1,0,tau)
vi=ai+delta*ui

x1ij=rnorm(1,2,1)
x1ij=ifelse(x1ij>4,4,x1ij)
x1ij=ifelse(x1ij<0,0,x1ij)
xij=x1ij

### original data genetation function
ori_gen=function(s){
  ## MCAR
  if(s==1){b0=1;b1=0;delta=0}
  ## MAR
  if(s==2){b0=0;b1=0.5;delta=0}
  ## Cluster-specific nonignorable 1
  if(s==3){b0=0;b1=0.5;delta=5}
  ## cluster-specific nonignorable 2
  if(s==4){b0=0;b1=0.5;delta=5}
  
  yij=xij*lambda+vi+eij
  y=matrix(NA,N,M)
  x=matrix(NA,N,M)
  PR=matrix(NA,N,M)
  for(i in 1:N){
    ai=rnorm(1,0,1)
    ui=rnorm(1,0,tau)
    vi=ai+delta*ui
    for(j in 1:M){
      eij=rnorm(1,0,1)
      x1ij=rnorm(1,2,1)
      x1ij=ifelse(x1ij>4,4,x1ij)
      x1ij=ifelse(x1ij<0,0,x1ij)
      x[i,j]=x1ij
      PR[i,j]=h(b0+x[i,j]*b1+ui)
      y[i,j]=x[i,j]*lambda+vi+eij
    }
  }
  return(list(y=y,x=x,PR=PR))
}

sim_data1=ori_gen(1)
y=as.data.frame(sim_data1$y)
#y$num=c(1:200)

x=as.data.frame(sim_data1$x)
#x$num=c(1:200)

pr=as.data.frame(sim_data1$PR)
#pr$num=c(1:200)


### sampling 
# a
num1=sample(c(1:200),50,replace = TRUE)
new_y=c()
for(i in num1){
  new_y=rbind(new_y,y[i,])
}
rownames(new_y)=c(1:50)
# b
new_y1=c()
for(i in 1:50){
  n1=sample(c(1:200),1,replace = TRUE)
  for(j in 1:10){
    m1=sample(c(1:10),5,replace=TRUE)
  }
  temp=y[n1,m1]
  colnames(temp)=c(1:5)
  #temp$num=n1
  new_y1=rbind(new_y1,temp)
}
new_y2=c()
for(i in 1:50){
  n1=sample(c(1:200),1,replace = TRUE)
  for(j in 1:10){
    m1=sample(c(1:10),5,replace=TRUE)
  }
  temp=y[n1,m1]
  colnames(temp)=c(1:5)
  #temp$num=n1
  new_y2=rbind(new_y2,temp)
}
new_yy=cbind(new_y1,new_y2)
colnames(new_yy)=colnames(new_y)
newy=rbind(new_y,new_yy)
rownames(newy)=c(1:100)
head(newy)
