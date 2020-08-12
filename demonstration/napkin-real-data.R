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

nonRandomSample = function(OBS,mysize){
  OBSuniqueList = returnUnique(OBS)
  # Enumerate all possible values of column
  tmp = c()
  for (idx in 1:length(OBSuniqueList)){
    tmp = append(tmp,list(OBSuniqueList[[idx]]))
  }
  allpossible = expand.grid(tmp)
  colnames(allpossible) = colnames(OBS)
  mycollect = c()
  for (rowidx in 1:nrow(allpossible)){
    filtered_OBS = subset(OBS, W == allpossible[rowidx,'W'] & R == allpossible[rowidx,'R'] & X == allpossible[rowidx,'X'] & Y == allpossible[rowidx,'Y'] )
    if (nrow(filtered_OBS) > 0){
      mycollect = rbind(mycollect,filtered_OBS[1,]) 
    }
  }
  if (nrow(mycollect) < mysize){
    for (rowidx in 1:(mysize - nrow(mycollect))){
      mycollect = rbind(mycollect,OBS[rowidx,])
    }
  }
  row.names(mycollect) = c(1:nrow(mycollect))
  return(mycollect)
}

goodSample = function(OBS,mysize){
  OBSuniqueList = returnUnique(OBS)
  trialnum = 0
  trialMax = 10
  stopSwitch = TRUE
  while(1){
    trialnum = trialnum + 1 
    stopSwitch = TRUE 
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
    if (trialnum > trialMax){
      OBS.sampled = nonRandomSample(OBS,mysize)
    }
  }
  rownames(OBS.sampled) = c(1:nrow(OBS.sampled))
  return(OBS.sampled)
}


dataGen = function(seednum,N,Nmax){
  set.seed(seednum)
  data(alarm)
  DATA = data.matrix(alarm) - 1 
  DATA = data.frame(DATA[,c('STKV','CCHL','HR','CO','BP','TPR','ANES')])
  
  # Prob.Snode.A_K = matrix(runif(length(unique(DATA$ANES)) * length(unique(DATA$STKV)),min=0,max=1)
  #                           ,nrow = length(unique(DATA$ANES)),ncol = length(unique(DATA$STKV)) )
  Prob.Snode.A_K = matrix(c(0.85,0.75,0.8,0.65,0.65,0.7),nrow = length(unique(DATA$ANES)),ncol = length(unique(DATA$STKV)))
  
  mytheta = mapply(function(anesval,stkvval){
    rowidx_Pi = anesval+1 
    colidx_Pi = stkvval+1
    result_val = Prob.Snode.A_K[rowidx_Pi,colidx_Pi]
    return(result_val)
  },DATA$ANES,DATA$STKV)
  
  iterMax = 10
  iteridx = 0
  while(1){
    iteridx = iteridx + 1 
    DATA[,'mytheta'] = mytheta
    taking_idx = sapply(mytheta,function(p) rbinom(1,1,p))
    DATA$taking_idx = taking_idx
    sampled_df = subset(DATA,taking_idx==1)
    if (nrow(sampled_df) > Nmax){
      break 
    }
    if (iteridx > iterMax){
      sampled_df =DATA
      break 
    }
  }
  rownames(sampled_df) = c(1:nrow(sampled_df))
  summary(sampled_df)
  
  # Making Napkin
  colnames(sampled_df) = c("K","W","R","X",'Y',"TPR","A","taking_idx")
  # OBS.Large = sampled_df[,c("W","R","X","Y")]
  OBS = sampled_df[,c("W","R","X","Y")]
  # if (N >= nrow(OBS.Large)){
  #   N = nrow(OBS.Large)
  # }
  return(OBS)
  # OBS = goodSample(OBS.Large,N)
  # return(list(OBS.Large,OBS))
}

