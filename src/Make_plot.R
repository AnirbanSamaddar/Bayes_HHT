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
  stop("Usage: Make_plot.R resolution")
}
resolution = eval(parse(text = args[1])); print(resolution)


output_dir = '~/output'
figures_dir = '~/figures'
dir.create(figures_dir)

############# Preparing Figure 3 #############
Make_plot = function(core=1,rep=10,resolution = c(1,5)){

  if(core == 1){   ### cor = 0.9 h2 = 5% S = 20 n=10K
    ss = '10K'
    cor = 'r:0.9'
    QTL=seq(25,500,by=25)
  }else if(core == 2){   ### cor = 0.9 h2 = 5% S = 20 n=50K
    ss = '50K'
    cor = 'r:0.9'
    QTL=seq(25,500,by=25)
  }

  res_disc_FDR_susie = array(NA,dim=c(5,100,rep,1))
  res_disc_FDR_MRHT = array(NA,dim=c(5,100,rep,1))
  res_disc_FDR_subfam = array(NA,dim=c(5,100,rep,1))
  res_disc_FDR_ind = array(NA,dim=c(5,100,rep,1))
  res_disc_power_susie = array(NA,dim=c(5,100,rep,1))
  res_disc_power_MRHT = array(NA,dim=c(5,100,rep,1))
  res_disc_power_subfam = array(NA,dim=c(5,100,rep,1))
  res_disc_power_ind = array(NA,dim=c(5,100,rep,1))
  id = c(1:rep)


  for(i in id){

    setwd(paste0(output_dir,"/Rep_",i,"_setting_",core))

    load("out_list_susie.RData")
    out_list[sapply(out_list,is.null)] = list(integer(0))
    clsize_out_list = lapply(out_list,function(a) lapply(a, length))
    
    ### Calculating power and FDR ar fixed resolution for SuSiE

    tmp = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(out_list[[k]])==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);sum(unlist(sapply(out_list[[k]][a],function(x) !any(QTL%in%x))))}))
    tmp1 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(out_list[[k]])==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);sum(unlist(sapply(out_list[[k]][a],function(x) (QTL%in%x))))}))
    tmp3 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(out_list[[k]])==0){0.001}else{max(sum(unlist(clsize_out_list[[k]])<=j),0.001)}))
    res_disc_FDR_susie[,,i,1] = tmp/tmp3
    res_disc_power_susie[,,i,1] = tmp1/length(QTL)
    print("#### susie done")

    load("out_list_MRHT.RData")
    out_list[sapply(out_list,is.null)] = list(integer(0))
    clsize_out_list = lapply(out_list,function(a) lapply(a, length))

    ### Calculating power and FDR ar fixed resolution for MRI: Bayes FDR

    tmp = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(out_list[[k]][[1]])==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);sum(unlist(sapply(out_list[[k]][a],function(x) !any(QTL%in%x))))}))
    tmp1 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(out_list[[k]][[1]])==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);sum(unlist(sapply(out_list[[k]][a],function(x) (QTL%in%x))))}))
    tmp3 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(out_list[[k]][[1]])==0){0.001}else{max(sum(unlist(clsize_out_list[[k]])<=j),0.001)}))
    res_disc_FDR_MRHT[,,i,1] = tmp/tmp3
    res_disc_power_MRHT[,,i,1] = tmp1/length(QTL)
    print("#### MRHT done")

    load("out_list_subfam_test.RData")
    out_list[sapply(out_list,is.null)] = list(integer(0))
    clsize_out_list = lapply(out_list,function(a) lapply(a, length))

    ### Calculating power and FDR ar fixed resolution for MRI: Subfam FDR

    tmp = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(out_list[[k]][[1]])==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);sum(unlist(sapply(out_list[[k]][a],function(x) !any(QTL%in%x))))}))
    tmp1 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(out_list[[k]][[1]])==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);sum(unlist(sapply(out_list[[k]][a],function(x) (QTL%in%x))))}))
    tmp3 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(out_list[[k]][[1]])==0){0.001}else{max(sum(unlist(clsize_out_list[[k]])<=j),0.001)}))
    res_disc_FDR_subfam[,,i,1] = tmp/tmp3
    res_disc_power_subfam[,,i,1] = tmp1/length(QTL)
    print("#### subfam done")
    
    load("out_list_ind_effect.RData")
    out_list[sapply(out_list,is.null)] = list(integer(0))
    clsize_out_list = lapply(out_list,function(a) lapply(a, length))
    
    ### Calculating power and FDR ar fixed resolution for Ind lvl testing

    tmp = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(out_list[[k]])==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);sum(unlist(sapply(out_list[[k]][a],function(x) !any(QTL%in%x))))}))
    tmp1 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(out_list[[k]])==0){0}else{a = which(unlist(clsize_out_list[[k]])<=j);sum(unlist(sapply(out_list[[k]][a],function(x) (QTL%in%x))))}))
    tmp3 = sapply(1:100,function(j) sapply(1:length(out_list),function(k) if(length(out_list[[k]])==0){0.001}else{max(sum(unlist(clsize_out_list[[k]])<=j),0.001)}))
    res_disc_FDR_ind[,,i,1] = tmp/tmp3
    res_disc_power_ind[,,i,1] = tmp1/length(QTL)
    print("#### ind done")
  }
  print("#########Done#########")


  setwd(figures_dir)
  new_dir = paste0('Setting_',core)
  dir.create(new_dir)
  setwd(new_dir)

  ### Preparing data for Figure 3
  
  f1 = function(x) mean(x,na.rm=T)
  f2 = function(x) sd(x,na.rm=T)
  Data = data.frame(FDR = c(c(apply(res_disc_FDR_susie,c(1,2,4),f1)),c(apply(res_disc_FDR_MRHT,c(1,2,4),f1)),c(apply(res_disc_FDR_subfam,c(1,2,4),f1)),c(apply(res_disc_FDR_ind,c(1,2,4),f1))),Power = c(c(apply(res_disc_power_susie,c(1,2,4),f1)),c(apply(res_disc_power_MRHT,c(1,2,4),f1)),c(apply(res_disc_power_subfam,c(1,2,4),f1)),c(apply(res_disc_power_ind,c(1,2,4),f1))),FDR_sd = c(c(apply(res_disc_FDR_susie,c(1,2,4),f2))/sqrt(rep),c(apply(res_disc_FDR_MRHT,c(1,2,4),f2))/sqrt(rep),c(apply(res_disc_FDR_subfam,c(1,2,4),f2))/sqrt(rep),c(apply(res_disc_FDR_ind,c(1,2,4),f2))/sqrt(rep)),Power_sd = c(c(apply(res_disc_power_susie,c(1,2,4),f2))/sqrt(rep),c(apply(res_disc_power_MRHT,c(1,2,4),f2))/sqrt(rep),c(apply(res_disc_power_subfam,c(1,2,4),f2))/sqrt(rep),c(apply(res_disc_power_ind,c(1,2,4),f2))/sqrt(rep)))
  Data = data.frame(Data,FDR_LB = Data$FDR - 1.96*Data$FDR_sd,FDR_UB = Data$FDR + 1.96*Data$FDR_sd,Power_LB = Data$Power - 1.96*Data$Power_sd,Power_UB = Data$Power + 1.96*Data$Power_sd)
  label_model = c(rep("SuSiE",(500*1)),rep("MRI: Bayes FDR",(500*1)),rep("MRI: Subfam FDR",(500*1)),rep("Ind. effects test",(500*1)))
  label_res = rep(gl(100,5,labels=sapply(1:100,function(i) paste0("Resolution: ",i))),(4*1)) 
  label_n = rep(gl(1,500,labels = ss),4)
  label_S = rep(gl(1,500,labels = paste0("S:",length(QTL))),(1*4))
  label_r = rep(gl(1,(500*4),labels = cor),1)
  Data = data.frame(Data,Method = label_model,Res = label_res,n = label_n, S = label_S, r = label_r)

  ### Preparing Figure 3

  Data1 = Data[(Data$Res%in%paste0(rep("Resolution: ",length(resolution)),resolution) ),]
  p1 = ggplot(Data1,aes(x=FDR,y=Power,color = Method)) + facet_wrap(~S+r+Res) + geom_point() + geom_line() + geom_errorbar(aes(ymin=Power_LB,ymax=Power_UB),width=0.002) +
  labs(x="FDR",y="Power") + theme_bw()+ theme(legend.position='none',axis.title.x=element_blank(),axis.title.y=element_blank()) + 
  xlim(0,0.25) + geom_vline(xintercept=0.05,linetype="dashed",color="black") + 
  scale_x_continuous(breaks=c(0,0.05,0.15,0.25))
  return(p1)

}

############# Plotting Figure 3 #############
p1 = Make_plot(rep=10,resolution = resolution)
setwd(paste0(figures_dir,'/Setting_1'))
png("Power_FDR_cor_V_S.png")
figure = ggarrange(p1,nrow=1,ncol=1,common.legend=TRUE)
annotate_figure(figure, left = textGrob("Power", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                    bottom = textGrob("Empirical FDR", gp = gpar(cex = 1)))
dev.off()
