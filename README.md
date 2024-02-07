# Bayesian Hierarchical Hypothesis Testing (BHHT)


The function `Bayes_HHT()` takes posterior samples from a Bayesian variable selection model (`B`, effects in columns, samples in rows), a hierarchical clustering of the predictors (i.e., the columns of `X` in `y=Xb+e`), and a Bayesian FDR threshold (`alpha`) and returns credible sets, constiting on predictors that are jointly associated with the outcome. Further details about the methodology can be found in (Samaddar, Maiti, and de los Campos, 2023)[]. The following example demonstrate the use of the function.



## Demonstrating Example (Mice)

To illustrate the use of the function we simulate a phenotype using 300 SNPs, 5 with non-zero effects, and a 10% heritability (i.e., model R-squared).

### 1) Simulation
  - 300 SNPs
  - 5 QTLs more or less each in the center of separate LD-block
  - 10% heritability

```r
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
### 2) Collecting posterior samples
  - Five chains, each 12,000 iterations, burn-in 2,000, and thinning interval 5.
  - BGLR software, model BayesC

After running the code below, the matrix `B` contains the posterior samples (effects in columns, samples in rows).

```r
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
DS_BHHT = Bayes_HHT(alpha=alpha_in,nSNP=p,merge_mt=merge_mt,B=B,type='Bayes',output_type='table')ds_bht=eval(parse(text=paste0('c(',sapply(X=DS_BHHT['clusters'],FUN=paste,collapse=','),')')))
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
- Fill the texts in the landing page
- Modify the output of Bayes_HHT
   
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
