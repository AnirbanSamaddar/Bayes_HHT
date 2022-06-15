rm(list = ls())
library(BGLR)
library(data.table)
library(ggplot2)
library(susieR)
library(ggpubr)
library(grid)

args <- commandArgs(trailingOnly = TRUE)
print(args)
if (length(args) < 1) {
  stop("Usage: Simulation_study.R core")
}
core = eval(parse(text=args[1])); print(core)


############# Genotype generating function ############# 
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


source('Multi-resolution_Test_function.R')
output_dir = '~/output'
dir.create(output_dir)
############# Function for MC replications. The variable 'core' simulates different settings. ############# 
seq_run = function(id,core = 1){
  setwd(output_dir)
  jobID = id
  ### For replication of results
  set.seed(19092264+jobID)
  if(core == 1){   ### cor = 0.9 h2 = 1.25% S = 20 n=10K
    shape1_par=0.5            # r ~ 0.9. Decreasing this means more correlation.
    h2_par=5/400              # h2 ~ 1.25%. Increasing this means higher signal-to-noise ratio.
    SSize=10000               # Sample size.
    nSNP_par=525              # Number of non-zero effects.
    Oracle = 30	              # Hyperparameter for SuSiE and Bayesian Spike-and-Slab model.
    QTL=seq(25,500,by=25)     # Position of the non-zero effects.
  }else if(core == 2){   ### cor = 0.9 h2 = 1.25% S = 20 n=50K
    shape1_par=0.5  # r ~ 0.9
    h2_par=5/400
    SSize=50000
    nSNP_par=525
    Oracle = 30	
    QTL=seq(25,500,by=25)
  }
  
  ### Preparing data
  Data_prep = function(shape1,h2,SSize,nSNP_par,QTL){
    X = getBlock(SSize,nSNP_par,shape1=shape1,shape2=3,replace=T)
    nQTL = length(QTL)
    beta = rnorm(nQTL) + 1
    signal = X[,QTL]%*%beta
    vG=var(signal)
    vE=vG*(1-h2)/h2
    error=rnorm(n=SSize,sd=sqrt(vE))
    y=signal+error
    rm(list = c("beta","signal","vG","vE","error"))  
    X = matrix(as.double(X),nrow=SSize)     
    return(list(X,y))  
  }
  master_input = Data_prep(shape1_par,h2_par,SSize,nSNP_par,QTL)
  
  
  ### Running SuSiE
  message("SuSiE starting")
  X = master_input[[1]]
  y = master_input[[2]]
  fit.susie = susie(X,y,L=Oracle)
  threshold = c(0,0.02,0.05,0.1,0.2)
  out_list = lapply(threshold,function(x) susie_get_cs(fit.susie,X,coverage=(1-x))$cs)
  B = t(susie_get_posterior_samples(fit.susie,15000)$b)
  samples = list(susie=B)

  ### Saving the output
  s=getwd()
  dir.create(paste0(s,"/Rep_",jobID,"_setting_",core))
  setwd(paste0(s,"/Rep_",jobID,"_setting_",core))
  save(out_list,file="out_list_susie.RData")
  

  #### Running BGLR model
  message("BGLR starting")
  nSNP = dim(master_input[[1]])[2]
  Fit.BGLR=function(X,y){

	thin = 5
	prob_inc = Oracle/nSNP
	burnIn = 2500
	nIter = 15000

  priors=list(list(model="BayesC",probIn=prob_inc,counts=110,saveEffects=TRUE))
  XX = crossprod(X)
  idPriors = rep(1,dim(XX)[1])
	Xy=as.vector(crossprod(X,y))
  fm = BLRCross(y=y,XX=XX,Xy=Xy,priors=priors,idPriors = idPriors,nIter=nIter,burnIn=burnIn,thin=thin,verbose = F)
  samples=readBinMat('ETA_1_b.bin')

  return(samples)
  }
  s = getwd()
  L = master_input
  B = lapply(1:4,function(i) Fit.BGLR(L[[1]],L[[2]]))
  B = rbind(B[[1]],B[[2]],B[[3]],B[[4]])
  samples$BVS = B
  save(samples,file = "samples.RData")

  ### Performing Hierarchical clustering (w euclidean distance) over the inputs   
  XX = crossprod(L[[1]])
  p = ncol(XX)
  DIS = matrix(0,nrow=p,ncol=p)
  for(a in 1:(p-1)){
  	for(b in (a+1):p){
		DIS[a,b] = sqrt(XX[a,a]+XX[b,b]-2*XX[a,b])
		DIS[b,a] = DIS[a,b]	
	}
  }
  r = hclust(as.dist(DIS))
  
  ### function for running the three hypothesis tests -- MRI: Bayes FDR, MRI: Subfam FDR, Ind lvl 
  Run = function(s,alpha_in=threshold){

    ### Run Algorithm 1: Discovery set FDR control
    message("Algorithm 1 starting")
    merge_mt = r$merge
    nSNP = ncol(B)
    out_list = lapply(1:length(alpha_in),function(i) Method_Bayes(alpha = alpha_in[i],nSNP,merge_mt,B)[[1]])

    ### Saving the output
    save(out_list,file=paste0(s,"/out_list_MRHT.RData"))
     
    
    ### Run Algorithm 2: Subfamily-wise FDR control
    message("Algorithm 2 starting")
    out_list = lapply(1:length(alpha_in),function(i) Method_Subfam_test(alpha = alpha_in[i],nSNP,merge_mt,B)[[1]])

    ### Saving the output
    save(out_list,file=paste0(s,"/out_list_subfam_test.RData"))


    ### Run supplementary Algorithm 1: Indivial level testing
    message("supplementary Algorithm 1 starting")
    out_list = lapply(1:length(alpha_in),function(i) Method_ind_lvl(alpha = alpha_in[i],nSNP,B)[[1]])

    ### Saving the output
    save(out_list,file=paste0(s,"/out_list_ind_effect.RData"))
  }
  
  ### Running the three hypothesis tests
  s = getwd()
  tmp1 = Run(s,threshold)
}

############# Sequentially running 10 MC reps #############
system.time(for(id in 1:10){
  seq_run(id,core=1)
})


quit(save='no')