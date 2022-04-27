################################ Test ####################################
rm(list=ls())
library(BGLR)
source("Multi-resolution_Test_function.R")

# Data generation
X = matrix(rnorm(1000000,5,1),10000,100)
X[,11] = X[,12]
beta = matrix(rep(1,5),ncol=1)
y = X[,seq(10,50,by=10)]%*%beta + rnorm(100,0,20)
p = ncol(X)

# Hyper-parameters    
thin = 5
prob_inc = 1/p
burnIn = 2500
nIter = 15000

# Fitting the model
priors=list(list(model="BayesC",probIn=prob_inc,counts=110,saveEffects=TRUE))
XX = crossprod(X)
idPriors = rep(1,p)
Xy=as.vector(crossprod(X,y))
fm = BLRCross(y=y,XX=XX,Xy=Xy,priors=priors,idPriors = idPriors,nIter=nIter,burnIn=burnIn,thin=thin,verbose = F)
B=readBinMat('ETA_1_b.bin')

# PIP plot
png('PIP.png')
plot(apply(B!=0,2,mean),cex=0.5)
dev.off()

# Calculating the euclidean distance matrix
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

# Bayes-FDR control
output_Bayes_FDR = Method_Bayes(alpha=alpha_in,nSNP=p,merge_mt=merge_mt,B=B)

# Subfamily-wise FDR control
output_subfam_FDR = Method_Subfam_test(alpha=alpha_in,nSNP=p,merge_mt=merge_mt,B=B)

print(output_Bayes_FDR)
print(output_subfam_FDR)


