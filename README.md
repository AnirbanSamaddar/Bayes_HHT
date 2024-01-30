# Bayesian Hierarchical Hypothesis Testing (BHHT)

## Example (Mice)
### Simulation
  - 300 SNPs
  - 5 QTLs more or less each in the center of separate LD-block
  - 10% heritability

```applescript
library(BGLR)
library(ggplot2)

set.seed(195021)
data(mice)
QTL=round(c(65.5,62.5*2,62.5*3.5,65.2*4.1))
p=300
X=scale(mice.X[,1:p],scale=FALSE,center=TRUE)
MAP=mice.map[1:p,]
signal=rowSums(X[,QTL])
error=rnorm(nrow(X))
error=error*sqrt(var(signal)*90/10)
y=signal+error
```
### Collecting posterior samples
  - Five chains, each 12,000 iterations, burn-in 2,000, and thinning interval 5.
  - BGLR software, model BayesC

```applescript
suppressMessages(
  fm<-BLRXy(y=y,ETA=list(list(X=X,model='BayesC',probIn=1/100,counts=1000,saveEffects=TRUE)),nIter=12000,burnIn=2000,verbose=FALSE)
)
B=readBinMat('ETA_1_b.bin')
for(i in 2:5){
  print(i)
	tmp=sample(1:p,size=p)
	suppressMessages(
 	  fm<-BLRXy(y=y,ETA=list(list(X=X[,tmp],model='BayesC',probIn=1/100,counts=1000,saveEffects=TRUE)),nIter=12000,burnIn=2000,verbose=FALSE)
	)
	B=rbind(B,readBinMat('ETA_1_b.bin')[,order(tmp)])
}
```
### Bayesian hierarchical hypothesis testing

```applescript
PIP=colMeans(B!=0)
source("~/Bayes_HHT/src/Multi-resolution_Test_function.R")
XX = abs(crossprod(X))
DIS = matrix(0,nrow=p,ncol=p)
for(a in 1:(p-1)){
	for(b in (a+1):p){
		DIS[a,b] = sqrt(XX[a,a]+XX[b,b]-2*XX[a,b])
		DIS[b,a] = DIS[a,b]	
  }
} 
r = hclust(as.dist(DIS))
merge_mt = r$merge
alpha_in = 0.05
Bayes = function(alpha,nSNP,merge_mt,B){
  tmp = Method_Bayes(alpha=alpha,nSNP=nSNP,merge_mt=merge_mt,B=B)[[1]]
  output = data.frame(cluster_id = seq_len(length(tmp)),clusters = rep(NA,length(tmp)),cPIP = rep(NA,length(tmp)))
  for(i in seq_len(length(tmp))){
    output$clusters[i] = paste(paste0('',tmp[[i]]),collapse=',')
    output$cPIP[i] = mean(apply(B[,tmp[[i]]]!=0,1,any))
    output$first[i]=min(tmp[[i]])
    output$last[i]=max(tmp[[i]])
    output$size[i]=length(tmp[[i]])
  }
  return(output)
}
DS_BHHT = Bayes(alpha=alpha_in,nSNP=p,merge_mt=merge_mt,B=B)
ds_bht=eval(parse(text=paste0('c(',sapply(X=DS_BHHT['clusters'],FUN=paste,collapse=','),')')))
PIP = apply(B!=0,2,mean)

DF=data.frame("Mbp"=MAP$mbp,PIP=PIP,"Variant Type"='Marker',check.names=F)
DF$'Variant Type'[QTL]='Causal'
DF$BHHT='Non Rej.'
DF$BHHT[ds_bht]='Discovery'

dir.create('~/Toy_example/',recursive=TRUE)
setwd('~/Toy_example/')
p2=ggplot(DF,aes(x=Mbp,y=PIP))+
  geom_point(aes(color=BHHT,size=BHHT,shape=`Variant Type`))+
  scale_shape_manual(values=c(16,1))+
  scale_size_manual(values=c(3,1.5))+
  ylim(c(0,1))+
  geom_hline(yintercept=c(.9),linetype='dashed')+
  geom_vline(xintercept=mice.map$mbp[QTL],linetype='dashed',color='grey50')

  for(i in 1:nrow(DS_BHHT)){
		p2=p2+annotate("rect", xmin = mice.map$mbp[DS_BHHT$first[i]]-.2, 
		             xmax = mice.map$mbp[DS_BHHT$last[i]]+.2, ymin =0, 
		             ymax = DS_BHHT$cPIP[i],alpha = .1,fill = "blue")
   }
dir.create('~/Toy_example/',recursive=TRUE)
setwd('~/Toy_example/')
ggsave('PIP_w_discovery_seg.png',width=10)
```

# Reference

# Work order

~~- Add the demonstrating example~~
~~- Change the test function~~
~~- Make necessary changes to the example~~
- Add Simulation 1 and Simulation 2 (using mice)
   
# Bayes_HHT
Implementation of Bayesian Hierarchical Hypothesis Testing 

In the simualation study, we have run each setting for 1000 MC reps. The following code is for demonstration where we run 10 MC reps for the setting -- ```n=10K S=20 r=0.9 h2=1.25%```.
```
$cd src
$core=1
$Rscript Simulation_study.R ${core}
```
The variable ```core``` defines the setting in ```Simulation_study.R``` file. Changing ```core=2``` will run another predefined setting where ```n=50K```. Please change the parameters in the lines 49-54 in ```Simulation_study.R``` file to run different settings. The output from the simulation study setting will be stored in the folder ```~\output```. 

To generate plots similar to Figure 3, please run the following code.
```
$cd src
$resolution='c(1,5)'
$core=1
$Rscript Make_plot.R ${resolution} ${core}
```
Note that to generate the power-FDR plot like Figure 3 one needs to fix the resolution. In the above code, change the array variable ```resolution``` to generate Figure 3 for different resolutions. The output will be stored in the folder ```$HOME'/figures/Setting_'${core}```.
