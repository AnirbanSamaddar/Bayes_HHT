rm(list = ls())
set.seed(195021)
library(BGLR)
library(data.table)
library(ggplot2)
library(susieR)
library(ggpubr)
library(grid)

figures_dir = '~/figures/Motivation'
dir.create(figures_dir)

############# Simulating phenotype from mice genotypes ############# 
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


############# Fit the BGLR model #############
LP=list(list(X=X,model='BayesC',probIn=1/10,counts=12,saveEffects=TRUE))
fm=BGLR(y=y,ETA=LP,nIter=5000,burnIn=1000,thin=5,verbose=F)
B=readBinMat('ETA_1_b.bin')


############# Prepare data for Figure 1(a) #############
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


############# Figure 1(a) #############
png(paste0(figures_dir,'/Figure_1_a.png'))
p2 = ggplot(Dataset,aes(x=bp_id,y=prob_inc)) + geom_rect(data=grp_dataset[91:95,] ,aes(x=NULL,y=NULL,xmin=start_bp_id,xmax=end_bp_id,ymin=0,ymax=jt_prob_inc),fill="lightblue") + geom_point(color = "gray40") + labs(x="SNP positions (in Mbp)",y="PIP") + ylim(0,1) + geom_hline(yintercept=0.95, linetype = "dashed", color="blue") + geom_vline(xintercept = Dataset$bp_id[QTL],linetype = "dashed",color = "blue") +theme_bw()
plot(p2)
dev.off()

############# Figure 1(b) #############
cor_dataset = abs(cor(X))
tmp = c(1:500)
dataset2 = expand.grid(tmp,tmp)
dataset2$cor_v = sapply(1:nrow(dataset2),function(i) cor_dataset[dataset2[i,1],dataset2[i,2]])
colnames(dataset2) = c("SNP1","SNP2","Cor")
png(paste0(figures_dir,'/Figure_1_b.png'))
p = ggplot(dataset2,aes(SNP1,SNP2))+ geom_tile(aes(fill = Cor)) +
  scale_fill_gradient(low="yellow", high="red") + labs(x = "SNP id", y = "SNP id") + geom_vline(xintercept = dataset2$SNP1[1:500][QTL],linetype = "dashed",color = "blue")
plot(p)
dev.off()


############# Supplementary Figure 1 #############
dataset3 = data.frame(Iterations = rep(c(1:nrow(B)),2), Samples = abs(c(B[,249:250])), label = c(rep("SNP 249",nrow(B)),rep("SNP 250",nrow(B))))
png(paste0(figures_dir,'/Supp_Figure_1.png'))
p = ggplot(dataset3,aes(x=Iterations,y=Samples,color = label)) + geom_point() + theme_bw() + labs(x = "Iteration", y = "Sample")
plot(p)
dev.off()

quit(save='no')