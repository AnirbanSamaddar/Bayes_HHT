## We demonstrate one setting of simulation 2

### Setting:
This is for demonstration. Please change the parameters in the below code to simulate another scenario.

- n: 10000
- S: 5
- $h^2$: 1.25%

```applescript
rm(list = ls())
n=10000
S=5
h2_par= 0.0125
```
### Create output directory and load packages

```applescript
dir.create('~/output/Simulation_2',recursive=TRUE)
setwd("~/output/Simulation_2")
library(susieR)
library(data.table)
library(BGData)
library(BGLR)
library(ggplot2)
```
### Generate the data

```applescript
setwd('~/output/Simulation_2/')
mc = 1
set.seed(19092264+mc)
message("Data loading ...")
chr = sample(1:22,size=1)
bed_file = paste0('path to UKB bed files','/chrom',chr,'.bed')
X = BEDMatrix(bed_file)
map_file = paste0('path to UKB map files','/chrom',chr,'.bim')
MAP = fread(map_file)

message("Data frame with variant sets ...")
TMP = data.frame(chr = chr, snp_id = MAP$V2, bp = MAP$V4, array = TRUE)

message('Create 100kbp block filter')
set.seed(19092264+(mc))
start = sample(which(TMP$bp>=(head(TMP$bp,1)+1e6) & TMP$bp<= (tail(TMP$bp,1)-1e6)),size=1)
TMP$block = FALSE
block = which(TMP$bp >= TMP$bp[start] & TMP$bp <= (TMP$bp[start] + 100e3))
TMP$block[block] = TRUE

message('Create QTL id')
set.seed(19092264+(mc))
cand_QTL = which(TMP$block & TMP$array)
split_vector = split(cand_QTL, cut(seq_along(cand_QTL), breaks = S))
QTLs = sapply(1:S,function(i){mid = floor(length(split_vector[[i]])/2);split_vector[[i]][mid]})
TMP$QTL_flag = FALSE
TMP$QTL_flag[QTLs] = TRUE

message('Create 500kbp flank filter')
start = head(which(TMP$block & TMP$array),1)
end = tail(which(TMP$block & TMP$array),1)
tmp1 = which(TMP$bp >= (TMP$bp[start] - 0.5e6) & TMP$bp <= TMP$bp[start])
tmp2 = which(TMP$bp >= (TMP$bp[end]) & TMP$bp <= (TMP$bp[end] + 0.5e6))
TMP$Flanks = FALSE
TMP$Flanks[c(tmp1,tmp2)] =TRUE

message("Create chunk (block + flank) ...")
set.seed(19092264+(mc))
sample_set = sample(1:nrow(X),size=n,replace=FALSE)
chunk = sort(which((TMP$array & (TMP$block + TMP$Flanks))))
core = sort(which((TMP$array & TMP$block)))
chunk_snps = X[sample_set,chunk]
chunk_snps = as.matrix(preprocess(chunk_snps, center = TRUE, impute = TRUE))
core_id = which(chunk%in%core)
qtls_id = which(core%in%which(TMP$QTL_flag))

message('Generate trait')
set.seed(19092264+(mc))
h2 = h2_par
Z = preprocess(X[sample_set,TMP$QTL_flag],center=TRUE,impute=TRUE)
effect_size = matrix(1,ncol=S)
signal = Z %*% t(effect_size)
vE = var(signal)*(1-h2)/h2
trait = signal + rnorm(nrow(signal),mean=0,sd=sqrt(vE))
print(var(signal)/var(trait))
```
## Run SuSiE

```applescript
set.seed(19092264+mc)
Oracle = S+10	
message('Running Susie ...')
XX = crossprod(chunk_snps)
Xy = t(chunk_snps) %*% trait
yy = sum(trait^2)
system.time(fit.susie <- susie_suff_stat(XtX=XX,Xty=Xy,yty=yy,n=n,L=Oracle))
B = t(susie_get_posterior_samples(fit.susie,15000)$b)
B = B[,core_id]

alphas = fit.susie$alpha[,core_id]
pairwise_cor = cor(chunk_snps[,core_id])
SuSiE = function(B,threshold=a){
    out_list = list()
    for(i in seq(Oracle)){
      # sort the alphas[i,] in decreasing order
      ordered_alphas = sort(alphas[i,], decreasing = TRUE)
      # get the indices of the ordered alphas
      ordered_alphas_index = order(alphas[i,], decreasing = TRUE)
      ordered_alphas_cumsum = cumsum(ordered_alphas)
      # get the minimum index of the ordered alphas that is greater than the (1-threshold)
      if(length(which(ordered_alphas_cumsum >= (1-threshold)))==0){
        ordered_alphas_index = NULL
      }else if(length(which(ordered_alphas_cumsum >= (1-threshold)))==length(core_id)){
        ordered_alphas_index = ordered_alphas_index[1]
      }else{
        tmp = which(ordered_alphas_cumsum < (1-threshold))
        ordered_alphas_index = ordered_alphas_index[seq((max(tmp)+1))]
      }
      if(length(ordered_alphas_index)>1){
        tmp = pairwise_cor[ordered_alphas_index,ordered_alphas_index]
        if(min(abs(tmp))<0.5){next}else{out_list[[i]] = ordered_alphas_index}
      }else{
        out_list[[i]] = ordered_alphas_index
      }
    }
    # remove the null elements from the list
    out_list = out_list[!sapply(out_list, is.null)]
    output = data.frame(cluster_id = seq_len(length(out_list)),clusters = rep(NA,length(out_list))
              ,cPIP = rep(NA,length(out_list))
              , threshold = rep(threshold,length(out_list)))
    for(i in seq_len(length(out_list))){
      output$clusters[i] = paste(paste0('',out_list[[i]]),collapse=',')
      output$cPIP[i] = mean(apply(as.matrix(B[,out_list[[i]]])!=0,1,any))
    }
    return(output)
}
threshold = c(0,0.02,0.05,0.1)
output = SuSiE(B=B,threshold=threshold[1])
for(a in threshold[2:length(threshold)]){tmp=SuSiE(B=B,threshold=a);output=rbind(output,tmp)}
output$method = rep('Susie',nrow(output))
samples = list(susie=B)
```
## Run BGLR
```applescript
message('Running BGLR ...')
set.seed(19092264+mc)
Fit.BGLR=function(X,y){
    thin = 5
    prob_inc = Oracle/ncol(X)
    burnIn = 2500
    nIter = 15000
    fm = BLRXy(y=y,ETA=list(list(X=X,model='BayesC',probIn=prob_inc,counts=110,saveEffects=TRUE)),nIter=nIter,burnIn=burnIn,verbose=FALSE)
    samples = readBinMat('ETA_1_b.bin')
    return(samples)
}
B = Fit.BGLR(X=chunk_snps,y=trait)[,core_id]
for(i in 2:4){
  B=rbind(B,Fit.BGLR(X=chunk_snps,y=trait)[,core_id])
}
samples$BVS = B
```
## Perform hierarchical clustering on the predictors
```applescript
XX = abs(XX[core_id,core_id]/n)
p = ncol(XX)
DIS = matrix(0,nrow=p,ncol=p)
for(a in 1:(p-1)){
    tryCatch({
      for(b in (a+1):p){
        DIS[a,b] = sqrt(XX[a,a]+XX[b,b]-2*XX[a,b])
        DIS[b,a] = DIS[a,b]	
      }
    },warning=function(w){print(paste0(w," at ",a," and ",b))})
}
r = hclust(as.dist(DIS))
merge_mt = r$merge
```
## Run BHHT with Bayesian FDR error control
```applescript
source("https://raw.githubusercontent.com/AnirbanSamaddar/Bayes_HHT/main/src/Multi-resolution_Test_function.R")
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
for(a in threshold[2:length(threshold)]){tmp=Bayes(threshold=a,nSNP=p,merge_mt=merge_mt,B=B,prune=0);output=rbind(output,tmp)}  
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
### Calculating the distances from QTLs
```applescript
fun = function(x){
  tmp=as.character(x)
  DS = as.numeric(strsplit(tmp,',')[[1]])
  D= as.matrix(dist(c(TMP$bp[TMP$QTL_flag],TMP$bp[core[DS]])))
  if(length(DS)==0){
    D = matrix(rep(1e50,length(qtls_id)),nrow=length(qtls_id))
  }else{
    D = D[1:length(qtls_id),(length(qtls_id)+1):ncol(D)]
  }
  #print(D)
  return(list(apply(as.matrix(D),1,min),apply(as.matrix(D),1,max)))
}
for(i in seq_len(nrow(output))){
  tmp = fun(output$clusters[i])
  output$min[i] = paste0(tmp[[1]],collapse=',')
  output$max[i] = paste0(tmp[[2]],collapse=',')
}
```
## Saving the output
```applescript
s = getwd()
args = paste0("n=",n,"_p=",p,"_S=",S,"_h2=",h2_par)
dir.create(paste0(s,'/Setting_',args,"/Rep_",mc),recursive=TRUE)
setwd(paste0(s,'/Setting_',args,"/Rep_",mc))
save(samples,file = "samples.RData")
write.table(output,file = paste0("output.txt"),row.names=FALSE)
```
## Preparing the plot data
```applescript
output_dir = '~/output/Simulation_2/'
rep = 1
resolution = c(0,10,100)
res_disc_FDR_susie = array(NA,dim=c(4,length(resolution),rep,1))
res_disc_FDR_MRHT = array(NA,dim=c(4,length(resolution),rep,1))
res_disc_FDR_ind = array(NA,dim=c(4,length(resolution),rep,1))
res_disc_BFDR_susie = array(NA,dim=c(4,length(resolution),rep,1))
res_disc_BFDR_MRHT = array(NA,dim=c(4,length(resolution),rep,1))
res_disc_BFDR_ind = array(NA,dim=c(4,length(resolution),rep,1))
res_disc_power_susie = array(NA,dim=c(4,length(resolution),rep,1))
res_disc_power_MRHT = array(NA,dim=c(4,length(resolution),rep,1))
res_disc_power_ind = array(NA,dim=c(4,length(resolution),rep,1))
id = c(1:rep)
QTL = qtls_id

for(i in id){
  message(i)
  setwd(paste0(output_dir,'Setting_',args,"/Rep_",i))
  output = read.table('output.txt',header=TRUE)
  ### Calculating power and FDR ar fixed resolution for SuSiE
  out_list = lapply(threshold,function(a){tmp=as.character(output$clusters[(output$method=='Susie'&output$threshold==a)]);
                  tmp = tmp[!duplicated(tmp)];lapply(seq_along(tmp),function(i) as.numeric(strsplit(tmp[i],',')[[1]]))})
  c_loc_fdr = lapply(threshold,function(a){tmp=output$clusters[(output$method=='Susie'&output$threshold==a)];tmp = tmp[!duplicated(tmp)];
                  tmp1=output$cPIP[(output$method=='Susie'&output$threshold==a)];
                  lapply(seq_along(tmp),function(i) if(tmp[i]==''){0}else{(1-as.numeric(tmp1[i]))})})
  dist_from_qtl = lapply(threshold,function(a){tmp1 = output$clusters[(output$method=='Susie'&output$threshold==a)];
                  tmp=as.character(output$min[(output$method=='Susie'&output$threshold==a)])[!duplicated(tmp1)];
                  lapply(seq_along(tmp),function(i) as.numeric(strsplit(tmp[i],',')[[1]]))})
  clsize_out_list = lapply(dist_from_qtl,function(a) lapply(a, min)) #lapply(out_list,function(a) lapply(a, length))

  tmp = sapply(1:length(resolution),function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])>resolution[j]*1e3);length(a)}))
  tmp1 = sapply(1:length(resolution),function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])<=resolution[j]*1e3);length(a)}))
  tmp3 = sapply(1:length(resolution),function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0.001}else{length(unlist(out_list[[k]]))}))
  tmp4 <- tryCatch({sapply(1:length(resolution),function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])>resolution[j]*1e3);if(length(a)==0){0}else{mean(unlist(c_loc_fdr[[k]][a]))}}))},
            warning = function(w) {
              message("Warning: ", w, "at i=",i)
              return(NA)
            })
  res_disc_BFDR_susie[,,i,1] = tmp4
  res_disc_FDR_susie[,,i,1] = tmp/tmp3
  res_disc_power_susie[,,i,1] = tmp1/S
  print("#### susie done")

  ### Calculating power and FDR at fixed resolution for BHHT: Bayes FDR
  out_list = lapply(threshold,function(a){tmp=as.character(output$clusters[(output$method=='Bayes'&output$threshold==a)]);
                  lapply(seq_along(tmp),function(i) as.numeric(strsplit(tmp[i],',')[[1]]))})
  c_loc_fdr = lapply(threshold,function(a){tmp=output$clusters[(output$method=='Bayes'&output$threshold==a)];
                  tmp1=output$cPIP[(output$method=='Bayes'&output$threshold==a)];
                  lapply(seq_along(tmp),function(i) if(tmp[i]==''){0}else{(1-as.numeric(tmp1[i]))})})
  dist_from_qtl = lapply(threshold,function(a){tmp=as.character(output$min[(output$method=='Bayes'&output$threshold==a)]);
                  lapply(seq_along(tmp),function(i) as.numeric(strsplit(tmp[i],',')[[1]]))})
  clsize_out_list = lapply(dist_from_qtl,function(a) lapply(a, min)) #lapply(out_list,function(a) lapply(a, length))

  tmp = sapply(1:length(resolution),function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])>resolution[j]*1e3);length(a)}))
  tmp1 = sapply(1:length(resolution),function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])<=resolution[j]*1e3);length(a)}))
  tmp3 = sapply(1:length(resolution),function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0.001}else{length(unlist(out_list[[k]]))}))
  tmp4 <- tryCatch({sapply(1:length(resolution),function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])>resolution[j]*1e3);if(length(a)==0){0}else{mean(unlist(c_loc_fdr[[k]][a]))}}))},
            warning = function(w) {
              message("Warning: ", w, "at i=",i)
              return(NA)
            })
  res_disc_BFDR_MRHT[,,i,1] = tmp4
  res_disc_FDR_MRHT[,,i,1] = tmp/tmp3
  res_disc_power_MRHT[,,i,1] = tmp1/S
  print("#### BHHT done")

  ### Calculating power and FDR ar fixed resolution for Ind lvl testing

  out_list = lapply(threshold,function(a){tmp=as.character(output$clusters[(output$method=='Ind_lvl'&output$threshold==a)]);
                  lapply(seq_along(tmp),function(i) as.numeric(strsplit(tmp[i],',')[[1]]))})
  c_loc_fdr = lapply(threshold,function(a){tmp=output$clusters[(output$method=='Ind_lvl'&output$threshold==a)];
                  tmp1=output$cPIP[(output$method=='Ind_lvl'&output$threshold==a)];
                  lapply(seq_along(tmp),function(i) if(tmp[i]==''){0}else{(1-as.numeric(tmp1[i]))})})
  dist_from_qtl = lapply(threshold,function(a){tmp=as.character(output$min[(output$method=='Ind_lvl'&output$threshold==a)]);
                  lapply(seq_along(tmp),function(i) as.numeric(strsplit(tmp[i],',')[[1]]))})
  clsize_out_list = lapply(dist_from_qtl,function(a) lapply(a, min)) #lapply(out_list,function(a) lapply(a, length))

  tmp = sapply(1:length(resolution),function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])>resolution[j]*1e3);length(a)}))
  tmp1 = sapply(1:length(resolution),function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])<=resolution[j]*1e3);length(a)}))
  tmp3 = sapply(1:length(resolution),function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0.001}else{length(unlist(out_list[[k]]))}))
  tmp4 <- tryCatch({sapply(1:length(resolution),function(j) sapply(1:length(out_list),function(k) if(length(unlist(out_list[[k]]))==0){0}else{a = which(unlist(clsize_out_list[[k]])>resolution[j]*1e3);if(length(a)==0){0}else{mean(unlist(c_loc_fdr[[k]][a]))}}))},
            warning = function(w) {
              message("Warning: ", w, "at i=",i)
              return(NA)
            })
  res_disc_BFDR_ind[,,i,1] = tmp4
  res_disc_FDR_ind[,,i,1] = tmp/tmp3
  res_disc_power_ind[,,i,1] = tmp1/S
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
label_model = c(rep("SuSiE",(4*length(resolution))),rep("DS-BFDR",(4*length(resolution))),rep("SNP-PIP",(4*length(resolution))))
label_res = rep(gl(length(resolution),4,labels=sapply(1:length(resolution),function(i) paste0("Resolution: ",resolution[i]," kbp"))),(3*1)) 
label_n = rep(gl(1,4*length(resolution)*3,labels = n),1)
label_S = rep(gl(1,4*length(resolution)*3,labels = paste0("S:",length(QTL))),1)
Data = data.frame(Data,Method = label_model,Res = label_res,n = label_n, S = label_S)
```
## Plotting Power vs FDR at resolution 0,10,100 kbp
```applescript
output_dir = '~/output/Simulation_2/'
figures_dir = 'figures/'
dir.create(paste0(output_dir,figures_dir,'/Setting_',args),recursive=TRUE)
setwd(paste0(output_dir,figures_dir,'/Setting_',args))
DATA = Data
maxClustSize=c(0,10,100)
alpha = c(0,0.02,0.05,0.1)
tmp = rep(gl(length(alpha),1,labels=alpha),nrow(DATA)/length(alpha))
DATA$alpha = tmp
for(i in 1:3){
    res=paste0('Resolution: ',maxClustSize[i],' kbp')
    DATA1=DATA[DATA$Res==res,]
    p=ggplot( DATA1,aes(x=FDR,y=Power))+
        geom_line(aes(size=Method,linetype=Method,color=Method),show.legend=TRUE)+ #scale_colour_manual(name="Method",values=c("red", "blue", "green"))+ 
        scale_size_manual(name="Method",values=c(1,0.5,0.5))+
        geom_point(aes(shape=alpha,color=Method),size=2) + #guides(color=FALSE) +
        labs(shape="Level") + scale_shape_manual(labels=c(0,0.02,0.05,0.1),values=c(21,22,23,24))+
        #scale_size_manual(name="Method",values=c(0.5,5,1),breaks=unique(DATA1$Method),labels=unique(DATA1$Method))+
        geom_errorbar(aes(ymin=Power-Power_sd, ymax=Power+Power_sd), width=.2)+ 
        geom_errorbar(aes(xmin=FDR-FDR_sd, xmax=FDR+FDR_sd))+ 
        xlab('Empirical FDR')+
        xlim(c(0,.1))+
        geom_vline(aes(xintercept=.05),linetype='dashed',col='grey29')+
        ggtitle(paste0(' 5 causal variants (Resolution ',maxClustSize[i],' kbp)'))+
        theme(legend.position = c(0.43, 0.20))+
        ylim(c(0,1))
    ggsave(file=paste0("power_fdr_plot_res",maxClustSize[i],"kbp.png"),plot = p,height=8,width = 10)
}
```
[Menue](https://github.com/AnirbanSamaddar/Bayes_HHT?tab=readme-ov-file#examples)
