###Group-wise inference functions needed
###Multi-resolution test approach. Input: merge matrix, nSNP, alpha
Method_Bayes = function(alpha=0.05,nSNP=100,merge_mt,B){
  
  stages = merge_mt
  
  
  #### Function to recreate clusters from Merge matrix  
  f1 = function(m){
    
    resolutions_1 = list(NULL)
    children_id = rep(NA,(nSNP+nrow(m)))
    for(i in 1:nrow(m)){
      if(max(m[i,])<0){
        resolutions_1[[i]] = -c(m[i,1],m[i,2])
        children_id[nSNP+i] = paste(-c(m[i,1],m[i,2]),collapse=",")
      }else if(min(m[i,])<0 & 0<max(m[i,])){
        resolutions_1[[i]] = c(resolutions_1[[m[i,(m[i,]>0)]]],-m[i,(m[i,]<0)])
        tmp1 = which(m[i,]>0)
        tmp2 = which(m[i,]<0)
        children_id[nSNP+i] = paste(c((nSNP+m[i,tmp1]),-m[i,tmp2]),collapse=",")
      }else{
        resolutions_1[[i]] = c(resolutions_1[[m[i,1]]],resolutions_1[[m[i,2]]])
        children_id[nSNP+i] = paste(nSNP+m[i,],collapse=",")
      }
    }   
    
    a = lapply(1:nSNP,function(i){return(i)})
    b = lapply(1:(nSNP+nrow(m)),function(i){ifelse(i<=nSNP,return(a[[i]]),return(resolutions_1[[(i-nSNP)]]))})
    return(list(b,children_id))
  }
  
  tmp = f1(stages)
  resolutions = tmp[[1]]
  children_id = tmp[[2]]
  
  
  #### Function for finding the discovery set level BFDR while controlling local FDR of all clusters at level t
  
  
  Inc_prob_calc = function(l){
    
    inc_prob = unlist(lapply(l,function(a) mean(apply(as.matrix(B[,a])!=0,1,any))))
    return(inc_prob)
  }
  
  Inc_prob = Inc_prob_calc(resolutions)
  FDR_ctrl_res = resolutions
  
  f = function(t){
    
    #### Selecting the outer node clusters where FDR is controlled at level t
    S = which(Inc_prob>=t)
    Rejection = Inc_prob>=t
    FDR_ctrl_3 = list(NULL)
    tmp=1
    id = rep(NA,length(S))
    for(i in S){
      if(is.na(children_id[i])){
        FDR_ctrl_3[[tmp]] = FDR_ctrl_res[[i]]
        id[tmp] = i
        tmp = tmp+1
      }else{
        tmp1 = as.numeric(strsplit(children_id[i],",")[[1]])
        if(all(!Rejection[tmp1])){
          FDR_ctrl_3[[tmp]] = FDR_ctrl_res[[i]]
          id[tmp] = i
          tmp = tmp+1			
        }
      }	
    }
    
    #### Calculating the BFDR of the discovery set
    a = 1-Inc_prob[id[which(!is.na(id))]]
    FDR_bar = sum(a)/max(length(FDR_ctrl_3),1)
    return(FDR_bar)
  }
  
  #### Finding the optimal t value for selecting Discovery set level BFDR
  Ordered_inc_prob = sort(Inc_prob,decreasing=TRUE)
  Ordered_FDR_bar = rep(1.1,length(Ordered_inc_prob)) ##Since BFDR <= 1 
  for(i in seq_along(Ordered_inc_prob)){
    temp_1 = f(Ordered_inc_prob[i])
    if(temp_1<=alpha){
      Ordered_FDR_bar[i] = temp_1
    }else{
      Ordered_FDR_bar[i] = temp_1
    }
  }
  
  I = ifelse(sum(Ordered_FDR_bar>alpha)>0,(which.max(Ordered_FDR_bar>alpha)-1),length(Ordered_FDR_bar))
  t_min = ifelse(I>0,Ordered_inc_prob[I],1) ###If I=0 then no selection
  
  #### Selecting the outer node clusters where FDR is controlled at level t_min 
  S = which(Inc_prob>=t_min)
  Rejection = Inc_prob>=t_min
  FDR_ctrl_3 = list(NULL)
  tmp=1
  id = rep(NA,length(S))
  
  
  for(i in S){
    if(is.na(children_id[i])){
      FDR_ctrl_3[[tmp]] = FDR_ctrl_res[[i]]
      id[tmp] = i
      tmp = tmp+1
    }else{
      tmp1 = as.numeric(strsplit(children_id[i],",")[[1]])
      if(all(!Rejection[tmp1])){
        FDR_ctrl_3[[tmp]] = FDR_ctrl_res[[i]]
        id[tmp] = i
        tmp = tmp+1			
      }
    }	
  }
  
  #output = data.frame(cluster_id = seq_len(length(FDR_ctrl_3)),clusters = rep(NA,length(FDR_ctrl_3)))
  #for(i in seq_len(length(FDR_ctrl_3))){
  #  output$clusters[i] = paste(paste0('',FDR_ctrl_3[[i]]),collapse=',')
  #}
  #return(output)
  return(list(FDR_ctrl_3,t_min))
}

###Test for individual effects. Input: nSNP, alpha
Method_ind_lvl = function(alpha=0.05,nSNP=100,B){
  
  
  Inc_prob = unlist(lapply(1:nSNP,function(a) mean(B[,a]!=0)))
  
  f = function(t){
    
    S = which(Inc_prob>=t)
    Local_fdr_ctrl_res = c(1:nSNP)[S]
    
    #### Calculating the BFDR of the discovery set
    a = unlist(lapply(Local_fdr_ctrl_res,function(b) (1-Inc_prob[b])))
    FDR_bar = sum(a)/max(length(Local_fdr_ctrl_res),1)
    return(FDR_bar)
  }
  
  #### Finding the optimal t value for selecting Discovery set level BFDR
  Ordered_inc_prob = sort(Inc_prob,decreasing=TRUE)
  Ordered_FDR_bar = rep(1.1,length(Ordered_inc_prob)) ##Since BFDR <= 1 
  for(i in seq_along(Ordered_inc_prob)){
    temp_1 = f(Ordered_inc_prob[i])
    if(temp_1<=alpha){
      Ordered_FDR_bar[i] = temp_1
    }else{
      break
    }
  }
  
  I = ifelse(sum(Ordered_FDR_bar>alpha)>0,(which.max(Ordered_FDR_bar>alpha)-1),length(Ordered_FDR_bar))
  t_min = ifelse(I>0,Ordered_inc_prob[I],1) ###If I=0 then no selection
  
  
  #### Selecting the clusters where FDR is controlled at level t_min
  S = which(Inc_prob>=t_min)
  Local_fdr_ctrl_res = c(1:nSNP)[S]
  
  
  return(list(Local_fdr_ctrl_res,t_min,Ordered_FDR_bar))
}

Method_Subfam_test = function(alpha=0.05,nSNP=100,merge_mt,B){
  stages = merge_mt
  
  
  #### Function to recreate clusters from Merge matrix  
  f1 = function(m){
    
    resolutions_1 = list(NULL)
    children_id = rep(NA,(nSNP+nrow(m)))
    for(i in 1:nrow(m)){
      if(max(m[i,])<0){
        resolutions_1[[i]] = -c(m[i,1],m[i,2])
        children_id[nSNP+i] = paste(-c(m[i,1],m[i,2]),collapse=",")
      }else if(min(m[i,])<0 & 0<max(m[i,])){
        resolutions_1[[i]] = c(resolutions_1[[m[i,(m[i,]>0)]]],-m[i,(m[i,]<0)])
        tmp1 = which(m[i,]>0)
        tmp2 = which(m[i,]<0)
        children_id[nSNP+i] = paste(c((nSNP+m[i,tmp1]),-m[i,tmp2]),collapse=",")
      }else{
        resolutions_1[[i]] = c(resolutions_1[[m[i,1]]],resolutions_1[[m[i,2]]])
        children_id[nSNP+i] = paste(nSNP+m[i,],collapse=",")
      }
    }   
    
    a = lapply(1:nSNP,function(i){return(i)})
    b = lapply(1:(nSNP+nrow(m)),function(i){ifelse(i<=nSNP,return(a[[i]]),return(resolutions_1[[(i-nSNP)]]))})
    return(list(b,children_id))
  }
  
  tmp = f1(stages)
  resolutions = tmp[[1]]
  children_id = tmp[[2]]
  
  
  #### Function for finding the discovery set level BFDR while controlling local FDR of all clusters at level t
  
  
  Inc_prob_calc = function(l){
    
    inc_prob = unlist(lapply(l,function(a) mean(apply(as.matrix(B[,a])!=0,1,any))))
    return(inc_prob)
  }
  
  Inc_prob = Inc_prob_calc(resolutions)
  rejection_id = rep(0,length(resolutions))
  tested_id = rep(0,length(resolutions))
  
  
  #### Subfamilywise testing
  if(Inc_prob[length(Inc_prob)]>=(1-alpha)){rejection_id[length(rejection_id)]=1
  ; tested_id[length(tested_id)]=1}
  
  for(i in length(resolutions):1){
    if(rejection_id[i]!=0){
      tmp1 = as.numeric(strsplit(children_id[i],",")[[1]])
      a = Inc_prob[tmp1]
      ordered_a = sort(a,decreasing = TRUE)
      f = function(t){
        sub_fam_rejection_set = which(a>=t)
        sub_fam_BFDR = sum((1-a[sub_fam_rejection_set]))/max(length(sub_fam_rejection_set),1)
        return(sub_fam_BFDR)
      }
      Ordered_BFDR = sapply(ordered_a,f)
      I = ifelse(sum(Ordered_BFDR>alpha)>0,(which.max(Ordered_BFDR>alpha)-1),length(Ordered_BFDR))
      t_min = ifelse(I>0,ordered_a[I],1) ## if I = 0 no selection
      sub_fam_rejection = which(a>=t_min)
      rejection_id[tmp1[sub_fam_rejection]] = 1
      tested_id[tmp1[sub_fam_rejection]] = 1
      sub_fam_non_rejection = which(a<t_min)
      tested_id[tmp1[sub_fam_non_rejection]] = 1
    }
    if(all(rejection_id[(i-1):1]==0)){
      break
    }
  }
  
  ####Selection of outer nodes
  FDR_ctrl_res = resolutions
  S = which(rejection_id!=0)
  FDR_ctrl_3 = list(NULL)
  tmp=1
  id = rep(NA,length(S))
  
  
  for(i in S){
    if(is.na(children_id[i])){
      FDR_ctrl_3[[tmp]] = FDR_ctrl_res[[i]]
      id[tmp] = i
      tmp = tmp+1
    }else{
      tmp1 = as.numeric(strsplit(children_id[i],",")[[1]])
      if(all(rejection_id[tmp1]==0)){
        FDR_ctrl_3[[tmp]] = FDR_ctrl_res[[i]]
        id[tmp] = i
        tmp = tmp+1			
      }
    }	
  }
  
  #### Calculating the BFDR of the discovery set
  a = 1-Inc_prob[id[which(!is.na(id))]]
  FDR_bar = sum(a)/max(length(FDR_ctrl_3),1)

  #output = data.frame(cluster_id = seq_len(length(FDR_ctrl_3)),clusters = rep(NA,length(FDR_ctrl_3)))
  #for(i in seq_len(length(FDR_ctrl_3))){
  #  output$clusters[i] = paste(paste0('',FDR_ctrl_3[[i]]),collapse=',')
  #}
  #return(output)
  return(list(FDR_ctrl_3,FDR_bar))  
}
