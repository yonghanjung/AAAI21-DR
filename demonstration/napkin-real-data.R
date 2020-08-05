library(bnlearn)
library(MASS)
library(nnet)
library(mise)

returnUnique = function(OBS){
  W = OBS[,1] # High dim surrogate
  R = OBS[,2] # Cofounder 0-numCate
  X = OBS[,3]
  Y = OBS[,4]
  Wunique = unique(W)[order(unique(W))]
  Runique = unique(R)[order(unique(R))]
  Xunique = unique(X)[order(unique(X))]
  Yunique = unique(Y)[order(unique(Y))]
  return(list(Wunique,Runique,Xunique,Yunique))
}

goodSample = function(OBS,mysize){
  OBSuniqueList = returnUnique(OBS)
  trialnum = 0
  stopSwitch = TRUE
  while(1){
    stopSwitch = TRUE 
    trialnum = trialnum + 1 
    sampled_idx = sample(c(1:nrow(OBS)),size=mysize)
    OBS.sampled = OBS[sampled_idx,]
    OBSSampleuniqueList = returnUnique(OBS.sampled)
    for (idx in 1:4){
      if (identical(OBSSampleuniqueList,OBSuniqueList) == FALSE){
        stopSwitch = FALSE 
      }
    }
    if (stopSwitch == TRUE){
      break
    }
  }
  rownames(OBS.sampled) = c(1:nrow(OBS.sampled))
  return(OBS.sampled)
}


dataGen = function(seednum,N){
  set.seed(seednum)
  data(alarm)
  DATA = data.matrix(alarm) - 1 
  DATA = data.frame(DATA[,c('STKV','CCHL','HR','CO','BP','TPR','ANES')])
  
  # Prob.Snode.A_K = matrix(runif(length(unique(DATA$ANES)) * length(unique(DATA$STKV)),min=0,max=1)
  #                           ,nrow = length(unique(DATA$ANES)),ncol = length(unique(DATA$STKV)) )
  Prob.Snode.A_K = matrix(c(0.15,0.75,0.8,0.2,0.3,0.7),nrow = length(unique(DATA$ANES)),ncol = length(unique(DATA$STKV)))
  
  mytheta = mapply(function(anesval,stkvval){
    rowidx_Pi = anesval+1 
    colidx_Pi = stkvval+1
    result_val = Prob.Snode.A_K[rowidx_Pi,colidx_Pi]
    return(result_val)
  },DATA$ANES,DATA$STKV)
  
  DATA[,'mytheta'] = mytheta
  taking_idx = sapply(mytheta,function(p) rbinom(1,1,p))
  DATA$taking_idx = taking_idx
  sampled_df = subset(DATA,taking_idx==1)
  rownames(sampled_df) = c(1:nrow(sampled_df))
  summary(sampled_df)
  
  # Making Napkin
  colnames(sampled_df) = c("K","W","R","X",'Y',"TPR","A","taking_idx")
  OBS.Large = sampled_df[,c("W","R","X","Y")]
  if (N >= nrow(OBS.Large)){
    N = nrow(OBS.Large)
  }
  OBS = goodSample(OBS.Large,N)
  
  return(list(OBS.Large,OBS))
}

