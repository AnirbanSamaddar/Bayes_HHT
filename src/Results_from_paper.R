rm(list = ls())
set.seed(195021)
library(BGLR)
library(data.table)
library(ggplot2)
library(susieR)

data(mice)
dir='http://mtweb.cs.ucl.ac.uk/mus/www/GSCAN/HS_GENOTYPES/'
files=paste0(dir,'chr',c(1:19,'X'),'.map')
MAP=fread(files[1],data.table=FALSE)
for(i in 2:length(files)){
                MAP=rbind(MAP,fread(files[i],data.table=FALSE))
}
SNP=sapply(FUN=function(x)x[[1]],strsplit(colnames(mice.X),split='_'))
MAP=MAP[MAP$V2%in%SNP,]
X=mice.X[,SNP%in%MAP$V2]
SNP=SNP[SNP%in%MAP$V2]

X=scale(X[,1:500])
 
h2=0.1
QTL=seq(from=50,to=450,by=100)
nQTL=length(QTL); n=nrow(X)
b=rep(1,nQTL)*sqrt(h2/nQTL)
signal=X[,QTL]%*%b
error=rnorm(n,sd=sqrt(var(signal)/h2))
y=signal+error

LP=list(  list(X=X,model='BayesC',probIn=1/10,counts=12,saveEffects=TRUE))
fm=BGLR(y=y,ETA=LP,nIter=5000,burnIn=1000,thin=5,verbose=F)

B=readBinMat('ETA_1_b.bin')

Dataset = data.frame(SNP_id = c(1:500), bp_id = MAP$V3[1:500],prob_inc = fm$ETA[[1]]$d)
id_prob_inc = sapply(QTL,function(x){tmp = c(0*c(1:(x-11)),Dataset$prob_inc[(x-10):(x+10)]
                                                            ,0*c((x+11):500)); return(which.max(tmp))})

t = 0.95
f = function(x){
  cuttoff_id = FALSE
  for(i in 20:0){
    tmp = mean(apply(as.matrix(B[,(x-i):(x+i)])!=0,1,any))
    if(tmp<t){
      cuttoff_id = TRUE
      break
    }
  }
  return(list(tmp,ifelse(cuttoff_id,min((i+1),20),i)))
}
grp_prob_inc = sapply(id_prob_inc,function(x) f(x)[[1]])
grp_size = sapply(id_prob_inc,function(x) f(x)[[2]])
grp_size_iter = matrix(0,nrow=length(grp_size),ncol = (20-min(grp_size)+1))
grp_prob_inc_iter = grp_size_iter

for(i in 20:min(grp_size)){
  tmp = sapply(1:length(grp_size),function(j) max(i,grp_size[j]))
  grp_size_iter[,(20-i+1)] = tmp
  tmp1 = sapply(1:length(id_prob_inc), function(j) mean(apply(as.matrix(B[,(id_prob_inc[j]-tmp[j]):(id_prob_inc[j]+tmp[j])])!=0,1,any)))
  grp_prob_inc_iter[,(20-i+1)] = tmp1
}

grp_dataset = data.frame(start = c(id_prob_inc-grp_size_iter), end = c(id_prob_inc+grp_size_iter), jt_prob_inc = c(grp_prob_inc_iter), iteration = c(sapply(1:ncol(grp_size_iter),function(i) rep(i,length(grp_size)))))
grp_dataset$start_bp_id = Dataset$bp_id[grp_dataset$start]
grp_dataset$end_bp_id = Dataset$bp_id[grp_dataset$end]

png('motivation.png')
p2 = ggplot(Dataset,aes(x=bp_id,y=prob_inc)) + geom_rect(data=grp_dataset[91:95,] ,aes(x=NULL,y=NULL,xmin=start_bp_id,xmax=end_bp_id,ymin=0,ymax=jt_prob_inc),fill="lightblue") + geom_point(color = "gray40") + labs(x="SNP positions (in Mbp)",y="PIP") + ylim(0,1) + geom_hline(yintercept=0.95, linetype = "dashed", color="blue") + geom_vline(xintercept = Dataset$bp_id[QTL],linetype = "dashed",color = "blue") +theme_bw()
plot(p2)
dev.off()


cor_dataset = abs(cor(X))
tmp = c(1:500)
dataset2 = expand.grid(tmp,tmp)
dataset2$cor_v = sapply(1:nrow(dataset2),function(i) cor_dataset[dataset2[i,1],dataset2[i,2]])
colnames(dataset2) = c("SNP1","SNP2","Cor")
png('correlation_heatmap.png')
p = ggplot(dataset2,aes(SNP1,SNP2))+ geom_tile(aes(fill = Cor)) +
  scale_fill_gradient(low="yellow", high="red") + labs(x = "SNP id", y = "SNP id") + geom_vline(xintercept = dataset2$SNP1[1:500][QTL],linetype = "dashed",color = "blue")
plot(p)
dev.off()


dataset3 = data.frame(Iterations = rep(c(1:nrow(B)),2), Samples = abs(c(B[,249:250])), label = c(rep("SNP 249",nrow(B)),rep("SNP 250",nrow(B))))
png('trace_plot.png')
p = ggplot(dataset3,aes(x=Iterations,y=Samples,color = label)) + geom_point() + theme_bw() + labs(x = "Iteration", y = "Sample")
plot(p)
dev.off()



# Simulation study

###Data generating function

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

source('Mult-resolution_Test_function.R')
output_dir = '../output'
dir.create(output_dir)
setwd(output_dir)
seq_run = function(id,core = 1){
  jobID = id
  ###For replication of results
  set.seed(19092264+jobID)
  if(core == 1){   ### cor = 0.9 h2 = 5% S = 20 n=10K
    shape1_par=0.5  # r ~ 0.9
    h2_par=20/400
    SSize=10000
    nSNP_par=525
    Oracle = 30	
    QTL=seq(25,500,by=25)
  }else if(core == 2){   ### cor = 0.9 h2 = 5% S = 20 n=10K
    shape1_par=0.5  # r ~ 0.9
    h2_par=20/400
    SSize=50000
    nSNP_par=525
    Oracle = 30	
    QTL=seq(25,500,by=25)
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
  
  QTL_1_id = QTL
  
  
  ###Running standard susie
  X = master_input[[1]]
  y = master_input[[2]]
  fit.susie = susie(X,y,L=Oracle)
  threshold = c(0,0.02,0.05,0.1,0.2)
  out_list = lapply(threshold,function(x) susie_get_cs(fit.susie,X,coverage=(1-x))$cs)
  B = t(susie_get_posterior_samples(fit.susie,15000)$b)
  samples = list(susie=B)
  

  FDR_1 = rep(NA,5)
  Power_1 = rep(NA,5)
  BFDR_1 = rep(NA,5)
  cl_length = rep(NA,5)
  cl_loc_FDR = lapply(1:5, function(i) NA)
  
  for(i in 1:5){
    
    if(length(out_list[[i]])!=0){
      T = sum(sapply(out_list[[i]],function(s) ifelse(sum(QTL_1_id%in%s)>0,1,0)))
      
      FDR_1[i] = 1 - T/length(out_list[[i]])
      
      Power_1[i] = sum(unlist(out_list[[i]])%in%QTL_1_id)/length(QTL_1_id)
      
      coverage = sapply(out_list[[i]],function(x) mean(apply(as.matrix(B[,x])!=0,1,any)))
      
      BFDR_1[i] = 1 - mean(coverage)
      
      cl_length[i] = mean(sapply(out_list[[i]],length))
      
      cl_loc_FDR[[i]] = 1 - coverage
      
    }else{
      FDR_1[i] = 0
      Power_1[i] = 0
      BFDR_1[i] = 0
      cl_length[i] = 0
      cl_loc_FDR[[i]] = 0
    }
    
  }
  
  ans_susie = c(FDR_1,Power_1,BFDR_1,cl_length)
  s=getwd()
  dir.create(paste0(s,"/Rep_",jobID,"_setting_",core))
  setwd(paste0(s,"/Rep_",jobID,"_setting_",core))
  save(out_list,file="out_list_susie.RData")
  write(ans_susie,ncolumns = length(ans_susie),file = paste0("output_susie.txt"))
  save(cl_loc_FDR,file="lfdr_susie.RData")
  
  #### Fit BGLR
  nSNP = dim(master_input[[1]])[2]
  
  
  Fit.BGLR=function(X,y){

	thin = 5
	prob_inc = Oracle/nSNP
	burnIn = 2500
	nIter = 15000

    priors=list(list(model="BayesC",probIn=prob_inc,counts=110,saveEffects=TRUE))
    XX = crossprod(X)
    idPriors = rep(1,p)
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
  

  Run = function(s,alpha_in=threshold){

    #Run Algorithm 1: Discovery set FDR control
    merge_mt = r$merge
    nSNP = ncol(B)
    out_list = lapply(1:length(alpha_in),function(i) NA)
    t_min = rep(NA, length(alpha_in))
    FDR_1 = rep(NA,length(alpha_in))
    Power_1 = rep(NA,length(alpha_in))
    BFDR_1 = rep(NA,length(alpha_in))
    cl_length = rep(NA,length(alpha_in))
    cl_loc_FDR = lapply(1:5, function(i) NA)


    for(i in 1:length(alpha_in)){
      tmp <- Method_MRHT(alpha = alpha_in[i],nSNP,merge_mt,B)
      t_min[i] = tmp[[2]]
      out_list[[i]] = tmp[[1]]
      
      #### Function to calculate the Posterior prob. of H1 given a cluster
      pos_prob_of_alt = function(a) mean(apply(as.matrix(B[,a])!=0,1,any))
      prob_alt_out_list = unlist(lapply(out_list[[i]],pos_prob_of_alt))
      
      #### Calculating actual FDR, Power and predicted BFDR
      if(length(out_list[[i]][[1]])!=0){
        
        T = sum(sapply(out_list[[i]],function(s) ifelse(sum(QTL_1_id%in%s)>0,1,0)))
        
        FDR_1[i] = 1 - T/length(out_list[[i]])
        
        Power_1[i] = sum(unlist(out_list[[i]])%in%QTL_1_id)/length(QTL_1_id)
        
        BFDR_1[i] = 1 - sum(prob_alt_out_list)/length(out_list[[i]])
        
        cl_length[i] = mean(sapply(out_list[[i]],length))
        
        cl_loc_FDR[[i]] = 1 - prob_alt_out_list
        
      }else{
        FDR_1[i] = 0
        Power_1[i] = 0
        BFDR_1[i] = 0
        cl_length[i] = 0
        cl_loc_FDR[[i]] = 0
      }
    }
    ans_MRHT = c(FDR_1,Power_1,t_min,BFDR_1,cl_length)
    save(out_list,file=paste0(s,"/out_list_MRHT.RData"))
    save(cl_loc_FDR,file=paste0(s,"/lfdr_MRHT.RData"))
     
    
    #Run Algorithm 2: Subfamily-wise FDR control
    out_list = lapply(1:length(alpha_in),function(i) NA)
    FDR_1 = rep(NA,length(alpha_in))
    Power_1 = rep(NA,length(alpha_in))
    BFDR_1 = rep(NA,length(alpha_in))
    cl_length = rep(NA,length(alpha_in))
    cl_loc_FDR = lapply(1:5, function(i) NA)

    
    for(i in 1:length(alpha_in)){
      tmp <- Method_Subfam_test(alpha = alpha_in[i],nSNP,merge_mt,B)
      BFDR_1[i] = tmp[[2]]
      out_list[[i]] = tmp[[1]]
      
      #### Function to calculate the Posterior prob. of H1 given a cluster
      prob_alt_out_list = unlist(lapply(out_list[[i]],pos_prob_of_alt))


      #### Calculating actual FDR, Power and predicted BFDR
      if(length(out_list[[i]][[1]])!=0){
        
        T = sum(sapply(out_list[[i]],function(s) ifelse(sum(QTL_1_id%in%s)>0,1,0)))
        
        FDR_1[i] = 1 - T/length(out_list[[i]])
        
        Power_1[i] = sum(unlist(out_list[[i]])%in%QTL_1_id)/length(QTL_1_id)
        
        
        cl_length[i] = mean(sapply(out_list[[i]],length))
        
        cl_loc_FDR[[i]] = 1 - prob_alt_out_list

      }else{
        FDR_1[i] = 0
        Power_1[i] = 0
        BFDR_1[i] = 0
        cl_length[i] = 0
        cl_loc_FDR[[i]] = 0
      }
      
    }
    
    ans_subfam_test = c(FDR_1,Power_1,BFDR_1,cl_length)
    save(out_list,file=paste0(s,"/out_list_subfam_test.RData"))
    save(cl_loc_FDR,file=paste0(s,"/lfdr_subfam_test.RData"))


    #Run supplementary Algorithm 1: Indivial level testing
    out_list = lapply(1:length(alpha_in),function(i) NA)
    t_min = rep(NA, length(alpha_in))
    FDR_1 = rep(NA,length(alpha_in))
    Power_1 = rep(NA,length(alpha_in))
    BFDR_1 = rep(NA,length(alpha_in))
    cl_loc_FDR = lapply(1:5, function(i) NA)
    
    
    for(i in 1:length(alpha_in)){
      tmp <- Method_ind_lvl(alpha = alpha_in[i],nSNP,B)
      t_min[i] = tmp[[2]]
      out_list[[i]] = tmp[[1]]
      
      #### Function to calculate the Posterior prob. of H1 given a cluster
      prob_alt_out_list = unlist(lapply(out_list[[i]],pos_prob_of_alt))
      
      #### Calculating actual FDR, Power and predicted BFDR
      if(length(out_list[[i]])!=0){
        
        T = sum(sapply(out_list[[i]],function(s) ifelse(sum(QTL_1_id%in%s)>0,1,0)))
        
        FDR_1[i] = 1 - T/length(out_list[[i]])
        
        Power_1[i] = sum(unlist(out_list[[i]])%in%QTL_1_id)/length(QTL_1_id)
        
        BFDR_1[i] = 1 - sum(prob_alt_out_list)/length(out_list[[i]])
        
        cl_loc_FDR[[i]] = 1 - prob_alt_out_list

      }else{
        FDR_1[i] = 0
        Power_1[i] = 0
        BFDR_1[i] = 0
        cl_loc_FDR[[i]] = 0
      }
    }
    
    
    ans_ind_effect = c(FDR_1,Power_1,t_min,BFDR_1)
    save(out_list,file=paste0(s,"/out_list_ind_effect.RData"))
    save(cl_loc_FDR,file=paste0(s,"/lfdr_ind_effect.RData"))

    ###Print output
    write(ans_MRHT,ncolumns = length(ans_MRHT),file = paste0("output_MRHT.txt"))
    write(ans_ind_effect,ncolumns = length(ans_ind_effect),file = paste0("output_ind_effect.txt"))
    write(ans_subfam_test,ncolumns = length(ans_subfam_test),file = paste0("output_subfam_test.txt"))
    return(unlist(out_list))
  }
  
  s = getwd()
  tmp1 = Run(s,threshold)
  unlink(grep(pattern = "Rep_",list.files(s),value = T))
  
}


system.time(for(id in 1:10){
  seq_run(id,core=1)
})
