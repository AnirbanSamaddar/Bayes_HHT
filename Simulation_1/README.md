## We demonstrate one setting of simulation 1 here

### Setting:
- S: 20
- n: 50,000
- $\rho$: 0.9

### Create output directory and load packages

```applescript
dir.create('~/output/',recursive=TRUE)
setwd("~/output/")
rm(list = ls())
library(susieR)
library(BGData)
```

### Data generation

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
