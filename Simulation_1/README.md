## We demonstrate one setting of simulation 1 here

### Setting:
- S: 20
- n: 50,000
- $\rho$: 0.9
- $h^2$: 1.25%

### Create output directory and load packages

```applescript
dir.create('~/output/',recursive=TRUE)
setwd("~/output/")
rm(list = ls())
library(susieR)
library(BGData)
library(BGLR)
```

### Genotype generation functions

```applescript
perm=function(x,prop=1,n=length(x),...){
  q=floor(n*prop)
  y=x
  if(q<n){
    tmp=sample(1:n,size=q,...)
    y[tmp]=sample(x[tmp],size=q)
  }
  return(y)
}

getBlock=function(n,p,freq=0.2,shape1=2,shape2=.5,replace=T){
  Z=matrix(nrow=n,ncol=p,NA)
  Z[,1]=rbinom(size=2,n=n,prob=freq)
  
  for(i in 2:p){
    Z[,i]=perm(Z[,i-1],prop=rbeta(n=1,shape1=shape1,shape2=shape2),replace=replace)
  }
  return(Z)
  
}

```

### Setting the parameters to generate the data

```applescript
setwd('~/output/')
jobID = 1
set.seed(19092264+jobID)
shape1_par=0.5 
h2_par= 0.0125
SSize=50000
nSNP_par=525
S=20
Oracle=S+10
Dis_z=25
if(S==20){
  if(Dis_z == 25){
  QTL=seq(25,500,by=25)
  }else{
    QTL=seq(from=20,by=Dis_z,length.out=S)
  }
}else if(S==5){
  if(Dis_z == 100){
    QTL=seq(100,500,by=100)
  print(QTL)
  }else{
    QTL=seq(from=20,by=Dis_z,length.out=S)
  }
}

Data_prep = function(shape1,h2,SSize,nSNP,QTL,shape2=3){
  n = SSize
  nQTL = length(QTL)  
  blsize = nSNP
  X = array(NA,dim = c(n,nSNP))
  X[,1:blsize] = getBlock(n,blsize,shape1=shape1,shape2=shape2,replace=T)
  beta = rnorm(nQTL) + 1
  signal = X[,QTL]%*%beta
  vG=var(signal)
  vE=vG*(1-h2)/h2
  error=rnorm(n=n,sd=sqrt(vE))
  y=signal+error
  rm(list = c("beta","signal","vG","vE","error"))
  X = matrix(as.double(X),nrow=SSize)  
  return(list(X,y))
}

master_input = Data_prep(shape1_par,h2_par,SSize,nSNP_par,QTL)
X = master_input[[1]]
y = master_input[[2]]
```
## Run SuSiE

```applescript
message('Running Susie ...')
fit.susie = susie(X,y,L=Oracle)
B = t(susie_get_posterior_samples(fit.susie,15000)$b)
SuSiE = function(X,B,threshold=a){
    out_list = susie_get_cs(fit.susie,X,coverage=(1-threshold))$cs
    output = data.frame(cluster_id = seq_len(length(out_list)),clusters = rep(NA,length(out_list))
              ,cPIP = rep(NA,length(out_list)), threshold = rep(threshold,length(out_list)))
    for(i in seq_len(length(out_list))){
      output$clusters[i] = paste(paste0('',out_list[[i]]),collapse=',')
      output$cPIP[i] = mean(apply(as.matrix(B[,out_list[[i]]])!=0,1,any))
    }
    return(output)
}
threshold = c(0,0.02,0.05,0.1,0.2)
output = SuSiE(X=X,B=B,threshold=threshold[1])
for(a in threshold[2:length(threshold)]){tmp=SuSiE(X=X,B=B,threshold=a);output=rbind(output,tmp)}
output$method = rep('Susie',nrow(output))
#print(dim(B))
samples = list(susie=B)
```
## Run BGLR
```applescript
message('Running BGLR ...')
Fit.BGLR=function(X,y){
    thin = 5
    prob_inc = Oracle/nSNP_par
    burnIn = 2500
    nIter = 15000
    fm = BLRXy(y=y,ETA=list(list(X=X,model='BayesC',probIn=prob_inc,counts=110,saveEffects=TRUE)),nIter=nIter,burnIn=burnIn,verbose=FALSE)
    samples = readBinMat('ETA_1_b.bin')
    return(samples)
}
B = Fit.BGLR(X=X,y=y)
for(i in 2:4){
  B=rbind(B,Fit.BGLR(X=X,y=y))
}
samples$BVS = B
```
## Perform hierarchical clustering on the predictors
```applescript
X=preprocess(X,scale=TRUE,center=TRUE)
XX = round(abs(crossprod(X)/SSize),2)
p = ncol(XX)
DIS = matrix(0,nrow=p,ncol=p)
for(a in 1:(p-1)){
  for(b in (a+1):p){
    DIS[a,b] = sqrt(XX[a,a]+XX[b,b]-2*XX[a,b])
    DIS[b,a] = DIS[a,b]	
  }
}
r = hclust(as.dist(DIS))
merge_mt = r$merge
```
## Run BHHT with Bayesian FDR error control
```applescript
Bayes = function(threshold,nSNP,merge_mt,B,prune=0){
    tmp = Bayes_HHT(alpha=threshold,nSNP=nSNP,merge_mt=merge_mt,B=B,type='Bayes',output_type='list')
    output = data.frame(cluster_id = seq_len(length(tmp)),clusters = rep(NA,length(tmp)),cPIP = rep(NA,length(tmp))
              , threshold = rep(threshold,length(tmp)), method = rep('Bayes',length(tmp)))
    for(i in seq_len(length(tmp))){
      tmp1 = tmp[[i]]
      id_prune = which(apply(as.matrix(B[,tmp1])!=0,2,mean)<=prune)
      tmp2 = if(length(id_prune)==0){tmp1}else{tmp1[-id_prune]}
      output$clusters[i] = paste(paste0('',tmp2),collapse=',')
      output$cPIP[i] = mean(apply(as.matrix(B[,tmp2])!=0,1,any))
    }
    return(output)
}
message('Running BHHT Bayes ...')
output = rbind(output,Bayes(threshold=threshold[1],nSNP=p,merge_mt=merge_mt,B=B,prune=0))
for(a in threshold[2:length(threshold)]){tmp=Subfam(threshold=a,nSNP=p,merge_mt=merge_mt,B=B,prune=0);output=rbind(output,tmp)}  
```
## Run SNP level testing
```applescript
Individual = function(threshold,nSNP,B){
    tmp = Bayes_HHT(alpha=threshold,nSNP=nSNP,merge_mt=merge_mt,B=B,type='Ind_lvl',output_type='list')
    output = data.frame(cluster_id = seq_len(length(tmp)),clusters = rep(NA,length(tmp)),cPIP = rep(NA,length(tmp))
              , threshold = rep(threshold,length(tmp)), method = rep('Ind_lvl',length(tmp)))
    for(i in seq_len(length(tmp))){
      output$clusters[i] = paste(paste0('',tmp[[i]]),collapse=',')
      output$cPIP[i] = mean(apply(as.matrix(B[,tmp[[i]]])!=0,1,any))
    }
    return(output)
}
message('Running Ind lvl ...')
output = rbind(output,Individual(threshold=threshold[1],nSNP=p,B=B))
for(a in threshold[2:length(threshold)]){tmp=Individual(threshold=a,nSNP=p,B=B);output=rbind(output,tmp)}
```
## Saving the output
```applescript
s = getwd()
args = paste0("n=",SSize,"_p=",nSNP_par,"_S=",S,"_shape1_par=",shape1_par,"_h2=",h2_par,"_Dis_z=",Dis_z)
dir.create(paste0(s,'/Setting_',args,"/Rep_",jobID),recursive=TRUE)
setwd(paste0(s,'/Setting_',args,"/Rep_",jobID))
save(samples,file = "samples.RData")
write.table(output,file = paste0("output.txt"),row.names=FALSE)
```
## Preparing plot data
```applescript
output_dir = '~/output/'
rep = 1
if(shape1_par==0.2){cor=0.99}else if(shape1_par==0.5){cor=0.9}
res_disc_FDR_susie = array(NA,dim=c(5,100,rep,1))
res_disc_FDR_MRHT = array(NA,dim=c(5,100,rep,1))
res_disc_FDR_ind = array(NA,dim=c(5,100,rep,1))
res_disc_BFDR_susie = array(NA,dim=c(5,100,rep,1))
res_disc_BFDR_MRHT = array(NA,dim=c(5,100,rep,1))
res_disc_BFDR_ind = array(NA,dim=c(5,100,rep,1))
res_disc_power_susie = array(NA,dim=c(5,100,rep,1))
res_disc_power_MRHT = array(NA,dim=c(5,100,rep,1))
res_disc_power_ind = array(NA,dim=c(5,100,rep,1))
id = c(1:rep)

for(i in id){
  message(i)
  setwd(paste0(output_dir,'Setting_',args,"/Rep_",i))
  output = read.table('output.txt',header=TRUE)
  ### Calculating power and FDR ar fixed resolution for SuSiE
  out_list = lapply(threshold,function(a){tmp=as.character(output$clusters[(output$method=='Susie'&output$threshold==a)]);
                  lapply(seq_along(tmp),function(i) as.numeric(strsplit(tmp[i],',')[[1]]))})
  c_loc_fdr = lapply(threshold,function(a){tmp=output$clusters[(output$method=='Susie'&output$threshold==a)];
                  tmp1=output$cPIP[(output$method=='Susie'&output$threshold==a)];
                  lapply(seq_along(tmp),function(i) if(tmp[i]==''){0}else{(1-as.numeric(tmp1[i]))})})
  clsize_out_list = lapply(out_list,function(a) lapply(a, length))

  tmp = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);sum(unlist(sapply(out_list[[k]][a],function(x) !any(QTL%in%x))))}))
  tmp1 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);sum(unlist(sapply(out_list[[k]][a],function(x) (QTL%in%x))))}))
  tmp3 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0.001}else{max(sum(unlist(clsize_out_list[[k]])<=j),0.001)}))
  tmp4 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);mean(unlist(c_loc_fdr[[k]][a]))}))
  res_disc_BFDR_susie[,,i,1] = tmp4
  res_disc_FDR_susie[,,i,1] = tmp/tmp3
  res_disc_power_susie[,,i,1] = tmp1/length(QTL)
  print("#### susie done")

  ### Calculating power and FDR at fixed resolution for BHHT: Bayes FDR
  out_list = lapply(threshold,function(a){tmp=as.character(output$clusters[(output$method=='Bayes'&output$threshold==a)]);
                  lapply(seq_along(tmp),function(i) as.numeric(strsplit(tmp[i],',')[[1]]))})
  c_loc_fdr = lapply(threshold,function(a){tmp=output$clusters[(output$method=='Bayes'&output$threshold==a)];
                  tmp1=output$cPIP[(output$method=='Bayes'&output$threshold==a)];
                  lapply(seq_along(tmp),function(i) if(tmp[i]==''){0}else{(1-as.numeric(tmp1[i]))})})
  clsize_out_list = lapply(out_list,function(a) lapply(a, length))

  tmp = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);sum(unlist(sapply(out_list[[k]][a],function(x) !any(QTL%in%x))))}))
  tmp1 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);sum(unlist(sapply(out_list[[k]][a],function(x) (QTL%in%x))))}))
  tmp3 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0.001}else{max(sum(unlist(clsize_out_list[[k]])<=j),0.001)}))
  tmp4 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);mean(unlist(c_loc_fdr[[k]][a]))}))
  res_disc_BFDR_MRHT[,,i,1] = tmp4
  res_disc_FDR_MRHT[,,i,1] = tmp/tmp3
  res_disc_power_MRHT[,,i,1] = tmp1/length(QTL)
  print("#### MRHT done")

  ### Calculating power and FDR ar fixed resolution for Ind lvl testing

  out_list = lapply(threshold,function(a){tmp=as.character(output$clusters[(output$method=='Ind_lvl'&output$threshold==a)]);
                  lapply(seq_along(tmp),function(i) as.numeric(strsplit(tmp[i],',')[[1]]))})
  c_loc_fdr = lapply(threshold,function(a){tmp=output$clusters[(output$method=='Ind_lvl'&output$threshold==a)];
                  tmp1=output$cPIP[(output$method=='Ind_lvl'&output$threshold==a)];
                  lapply(seq_along(tmp),function(i) if(tmp[i]==''){0}else{(1-as.numeric(tmp1[i]))})})
  clsize_out_list = lapply(out_list,function(a) lapply(a, length))

  tmp = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);sum(unlist(sapply(out_list[[k]][a],function(x) !any(QTL%in%x))))}))
  tmp1 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);sum(unlist(sapply(out_list[[k]][a],function(x) (QTL%in%x))))}))
  tmp3 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0.001}else{max(sum(unlist(clsize_out_list[[k]])<=j),0.001)}))
  tmp4 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);mean(unlist(c_loc_fdr[[k]][a]))}))
  res_disc_BFDR_ind[,,i,1] = tmp4
  res_disc_FDR_ind[,,i,1] = tmp/tmp3
  res_disc_power_ind[,,i,1] = tmp1/length(QTL)
  print("#### ind done")
}
f1 = function(x) mean(x,na.rm=T)
f2 = function(x) sd(x,na.rm=T)
Data = data.frame(FDR = c(c(apply(res_disc_FDR_susie,c(1,2,4),f1)),c(apply(res_disc_FDR_MRHT,c(1,2,4),f1))
                    ,c(apply(res_disc_FDR_ind,c(1,2,4),f1))),
                    Power = c(c(apply(res_disc_power_susie,c(1,2,4),f1)),c(apply(res_disc_power_MRHT,c(1,2,4),f1))
                    ,c(apply(res_disc_power_ind,c(1,2,4),f1))),
                    BFDR = c(c(apply(res_disc_BFDR_susie,c(1,2,4),f1)),c(apply(res_disc_BFDR_MRHT,c(1,2,4),f1))
                    ,c(apply(res_disc_BFDR_ind,c(1,2,4),f1))),
                    FDR_sd = c(c(apply(res_disc_FDR_susie,c(1,2,4),f2))/sqrt(rep),
                    c(apply(res_disc_FDR_MRHT,c(1,2,4),f2))/sqrt(rep),
                    c(apply(res_disc_FDR_ind,c(1,2,4),f2))/sqrt(rep)),
                    Power_sd = c(c(apply(res_disc_power_susie,c(1,2,4),f2))/sqrt(rep),
                    c(apply(res_disc_power_MRHT,c(1,2,4),f2))/sqrt(rep),
                    c(apply(res_disc_power_ind,c(1,2,4),f2))/sqrt(rep)),
                    BFDR_sd = c(c(apply(res_disc_BFDR_susie,c(1,2,4),f2))/sqrt(rep),
                    c(apply(res_disc_BFDR_MRHT,c(1,2,4),f2))/sqrt(rep),
                    c(apply(res_disc_BFDR_ind,c(1,2,4),f2))/sqrt(rep)))
Data = data.frame(Data,FDR_LB = Data$FDR - 1.96*Data$FDR_sd,FDR_UB = Data$FDR + 1.96*Data$FDR_sd,Power_LB = Data$Power - 1.96*Data$Power_sd,Power_UB = Data$Power + 1.96*Data$Power_sd,BFDR_LB = Data$BFDR - 1.96*Data$BFDR_sd,BFDR_UB = Data$BFDR + 1.96*Data$BFDR_sd)
label_model = c(rep("SuSiE",(500*1)),rep("DS-BFDR",(500*1)),rep("SNP-PIP",(500*1)))
label_res = rep(gl(100,5,labels=sapply(1:100,function(i) paste0("Resolution: ",i))),(3*1)) 
label_n = rep(gl(1,500,labels = SSize),3)
label_S = rep(gl(1,500,labels = paste0("S:",length(QTL))),3)
label_r = rep(gl(1,(500*3),labels = paste0("r:",cor)),1)
Data = data.frame(Data,Method = label_model,Res = label_res,n = label_n, S = label_S, r = label_r)
```
