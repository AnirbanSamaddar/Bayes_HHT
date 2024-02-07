# dir.create('~/output/',recursive=TRUE)
# setwd("~/output/")
setwd("/mnt/home/samadda1/Bayes_Zoom/output/")
rm(list = ls())

####Package dependency
library(susieR)
library(BGData)

###cmd read
args <- commandArgs(trailingOnly = TRUE)
print(args)
for(s in args){eval(parse(text=s))}

#jobID = jobID + 1 ##modification req


###Data generation and model fitting

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


####Sequential run for 10 replications 
seq_run = function(id){
  #setwd('~/output/')
  setwd("/mnt/home/samadda1/Bayes_Zoom/output/")
  jobID = 20*(jobID-1)+id
  ###For replication of results
  set.seed(19092264+jobID)
  shape1_par=shape1_par  # r ~ 0.9
  h2_par= h2_par
  SSize=n
  nSNP_par=p
  Oracle = S+10	
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

  
  Data_prep = function(shape1_par,h2_par,SSize,nSNP_par,QTL){
    
    Data_gen = function(shape1=shape1_par,shape2=3,h2=h2_par,n=SSize,b=1,nSNP=nSNP_par,nQTL=c(50),repl=1){
 
      blsize = floor(nSNP/b)
      X = array(NA,dim = c(n,nSNP))
      for(j in 1:b){
        X[,(blsize*(j-1)+1):(blsize*j)] = getBlock(n,blsize,shape1=shape1,shape2=shape2,replace=T)
      }
      
      QTL_1=QTL
      beta = rnorm(sum(nQTL)) + 1
      signal = X[,QTL_1]%*%beta
      vG=var(signal)
      vE=vG*(1-h2)/h2
      error=rnorm(n=n,sd=sqrt(vE))
      y=signal+error
      rm(list = c("beta","signal","vG","vE","error"))  
      X = matrix(as.double(X),nrow=SSize)
      
      
      return(list(X,y))
    } 
    
    L = Data_gen(nQTL=length(QTL))
    
    return(L)
  }
  
  master_input = Data_prep(shape1_par,h2_par,SSize,nSNP_par,QTL)
  X = master_input[[1]]
  y = master_input[[2]]
  
  
  ###Standard susie
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
  
  
  #### Fit BGLR
  message('Running BGLR ...')
  Fit.BGLR=function(X,y){
    dyn.load("/mnt/research/quantgen/projects/BGLR/crossproducts/source/sampler.so")
    source("/mnt/research/quantgen/projects/BGLR/crossproducts/source/BGLR_cross.R")
    thin = 5
    prob_inc = Oracle/nSNP_par
    burnIn = 2500
    nIter = 15000
    fm = BGLR_cross(y=y,X=X,probIn=prob_inc,counts=110,nIter=nIter,burnIn=burnIn,thin=thin,verbose = F)
    index =seq(from=1,to=(nIter-burnIn),by=thin)
    samples = t(fm$samples$BETA)[-c(1:burnIn),][index,]
    return(samples)
  }
  B = lapply(1:4,function(i) Fit.BGLR(X=X,y=y))
  B = rbind(B[[1]],B[[2]],B[[3]],B[[4]])
  samples$BVS = B
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
  #print(DIS)
  r = hclust(as.dist(DIS))
  merge_mt = r$merge

  #### Run group-wise inference
  source("https://raw.githubusercontent.com/AnirbanSamaddar/Bayes_HHT/main/src/Multi-resolution_Test_function.R")
  Node_wise = function(threshold,nSNP,merge_mt,B,prune=0){
    tmp = Method_Node_wise(alpha=threshold,nSNP=nSNP,merge_mt=merge_mt,B=B)[[1]]
    output = data.frame(cluster_id = seq_len(length(tmp)),clusters = rep(NA,length(tmp)),cPIP = rep(NA,length(tmp))
              , threshold = rep(threshold,length(tmp)), method = rep('Node_wise',length(tmp)))
    for(i in seq_len(length(tmp))){
      tmp1 = tmp[[i]]
      id_prune = which(apply(as.matrix(B[,tmp1])!=0,2,mean)<=prune)
      tmp2 = if(length(id_prune)==0){tmp1}else{tmp1[-id_prune]}
      output$clusters[i] = paste(paste0('',tmp2),collapse=',')
      output$cPIP[i] = mean(apply(as.matrix(B[,tmp2])!=0,1,any))
    }
    return(output)
  }  
   Bayes = function(threshold,nSNP,merge_mt,B,prune=0){
    tmp = Method_Bayes(alpha=threshold,nSNP=nSNP,merge_mt=merge_mt,B=B)[[1]]
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
  Subfam = function(threshold,nSNP,merge_mt,B,prune=0){
    tmp = Method_Subfam_test(alpha=threshold,nSNP=nSNP,merge_mt=merge_mt,B=B)[[1]]
    output = data.frame(cluster_id = seq_len(length(tmp)),clusters = rep(NA,length(tmp)),cPIP = rep(NA,length(tmp))
              , threshold = rep(threshold,length(tmp)), method = rep('Subfam',length(tmp)))
    for(i in seq_len(length(tmp))){
      tmp1 = tmp[[i]]
      id_prune = which(apply(as.matrix(B[,tmp1])!=0,2,mean)<=prune)
      tmp2 = if(length(id_prune)==0){tmp1}else{tmp1[-id_prune]}
      output$clusters[i] = paste(paste0('',tmp2),collapse=',')
      output$cPIP[i] = mean(apply(as.matrix(B[,tmp2])!=0,1,any))
    }
    return(output)
  }
  Individual = function(threshold,nSNP,B){
    tmp = Method_ind_lvl(alpha=threshold,nSNP=nSNP,B=B)[[1]]
    output = data.frame(cluster_id = seq_len(length(tmp)),clusters = rep(NA,length(tmp)),cPIP = rep(NA,length(tmp))
              , threshold = rep(threshold,length(tmp)), method = rep('Ind_lvl',length(tmp)))
    for(i in seq_len(length(tmp))){
      output$clusters[i] = paste(paste0('',tmp[[i]]),collapse=',')
      output$cPIP[i] = mean(apply(as.matrix(B[,tmp[[i]]])!=0,1,any))
    }
    return(output)
  }


  threshold = c(0,0.02,0.05,0.1,0.2)
  message('Running MRI Node_wise ...')
  output = rbind(output,Node_wise(threshold=threshold[1],nSNP=p,merge_mt=merge_mt,B=B,prune=0))
  for(a in threshold[2:length(threshold)]){tmp=Node_wise(threshold=a,nSNP=p,merge_mt=merge_mt,B=B,prune=0);output=rbind(output,tmp)}
  message('Running MRI Bayes ...')
  output = rbind(output,Bayes(threshold=threshold[1],nSNP=p,merge_mt=merge_mt,B=B,prune=0))
  for(a in threshold[2:length(threshold)]){tmp=Bayes(threshold=a,nSNP=p,merge_mt=merge_mt,B=B,prune=0);output=rbind(output,tmp)}
  message('Running MRI Subfam ...')
  output = rbind(output,Subfam(threshold=threshold[1],nSNP=p,merge_mt=merge_mt,B=B,prune=0))
  for(a in threshold[2:length(threshold)]){tmp=Subfam(threshold=a,nSNP=p,merge_mt=merge_mt,B=B,prune=0);output=rbind(output,tmp)}
  message('Running Ind lvl ...')
  output = rbind(output,Individual(threshold=threshold[1],nSNP=p,B=B))
  for(a in threshold[2:length(threshold)]){tmp=Individual(threshold=a,nSNP=p,B=B);output=rbind(output,tmp)}


####saving the output
  s = getwd()
  dir.create(paste0(s,'/Setting_',paste(args[2:length(args)],collapse='_'),"/Rep_",jobID),recursive=TRUE)
  setwd(paste0(s,'/Setting_',paste(args[2:length(args)],collapse='_'),"/Rep_",jobID))
  save(samples,file = "samples.RData")
  write.table(output,file = paste0("output.txt"),row.names=FALSE)
}

system.time(for(id in 18:20){
  message(paste0('Running Rep ',(20*(jobID-1)+id)))
  seq_run(id)
})

quit(save="no")
