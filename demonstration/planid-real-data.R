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


dataSampling_S = function(Val1, Val2, Prob.S, DATA){
  copy_DATA = DATA 
  # Prob.S = matrix(runif(length(unique(Val1)) * length(unique(Val2)),min=0.5,max=1)
  #                           ,nrow = length(unique(Val1)),ncol = length(unique(Val2)) )
  
  mytheta = mapply(function(val1,val2){
    rowidx_Pi = val1+1
    colidx_Pi = val2+1
    result_val = Prob.S[rowidx_Pi,colidx_Pi]
    return(result_val)
  },Val1,Val2)
  copy_DATA[,'Snode'] = mytheta
  taking_idx = sapply(mytheta,function(p) rbinom(1,1,p))
  copy_DATA$taking_idx = taking_idx
  sampled_df = subset(copy_DATA,taking_idx==1)
  selected_DATA = DATA[rownames(sampled_df),]
  rownames(selected_DATA) = c(1:nrow(selected_DATA))
  return(selected_DATA)
}

dataGen = function(seednum,N){
  data(alarm)
  set.seed(seednum)
  DATA = data.matrix(alarm) - 1
  DATA = data.frame(DATA[,c('PMB','VMCH','SHNT','INT','DISC','VTUB','KINK','VLNG','SAO2','VALV','PVS','ACO2','ANES','CCHL','TPR')])
  
  # Conditioning on INT and PVS 
  INTunique = unique(DATA$INT)[order(unique(DATA$INT))]
  PVSunique = unique(DATA$PVS)[order(unique(DATA$PVS))]
  cval = c(0,0)
  numCval = 0
  for (c1val in INTunique){
    for (c2val in PVSunique){
      filteredDATA = subset(DATA, INT==c1val & PVS==c2val)
      if (nrow(filteredDATA) > numCval){
        numCval = nrow(filteredDATA)
        cval = c(c1val,c2val)
      }
    }
  }
  DATA = subset(DATA, INT==cval[1] & PVS==cval[2]) 
  
  # Generate S nodes 
  # Prob.S = matrix(runif(length(unique(Val1)) * length(unique(Val2)),min=0.5,max=1)
  #                           ,nrow = length(unique(Val1)),ncol = length(unique(Val2)) )
  
  # Prob.S.PMB.VMCH = matrix(runif(length(unique(DATA$PMB)) * length(unique(DATA$VMCH)),min=0.5,max=1),nrow = length(unique(DATA$PMB)),ncol = length(unique(DATA$VMCH)) )
  Prob.S.PMB.VMCH = matrix(c(0.65,0.55,0.8,0.82,0.76,0.53,0.9,0.95),nrow = length(unique(DATA$PMB)),ncol = length(unique(DATA$VMCH)) )
  
  # Prob.S.VMCH.SHNT = matrix(runif(length(unique(DATA$SHNT)) * length(unique(DATA$VMCH)),min=0.5,max=1),nrow = length(unique(DATA$SHNT)),ncol = length(unique(DATA$VMCH)) )
  Prob.S.VMCH.SHNT = matrix(c(0.85,0.94,0.8,0.65,0.76,0.53,0.6,0.55),nrow = length(unique(DATA$PMB)),ncol = length(unique(DATA$VMCH)) )
  
  # Prob.S.SHNT.KINK = matrix(runif(length(unique(DATA$SHNT)) * length(unique(DATA$KINK)),min=0.5,max=1),nrow = length(unique(DATA$SHNT)),ncol = length(unique(DATA$KINK)) )
  Prob.S.SHNT.KINK = matrix(c(0.95,0.65,0.75,0.85),nrow = length(unique(DATA$SHNT)),ncol = length(unique(DATA$KINK)) )
  
  # Prob.S.PMB.TPR = matrix(runif(length(unique(DATA$PMB)) * length(unique(DATA$TPR)),min=0.5,max=1),nrow = length(unique(DATA$PMB)),ncol = length(unique(DATA$TPR)) )
  Prob.S.PMB.TPR = matrix(c(0.65,0.75,0.85,0.75,0.7,0.6),nrow = length(unique(DATA$PMB)),ncol = length(unique(DATA$TPR)) )
  
  # Prob.S.DISC.ANES = matrix(runif(length(unique(DATA$DISC)) * length(unique(DATA$ANES)),min=0.5,max=1),nrow = length(unique(DATA$DISC)),ncol = length(unique(DATA$ANES)) )
  Prob.S.DISC.ANES = matrix(c(0.65,0.85,0.75,0.5),nrow = length(unique(DATA$SHNT)),ncol = length(unique(DATA$KINK)) )
  
  DATA = dataSampling_S(Val1 = DATA$PMB, Val2 = DATA$VMCH,Prob.S.PMB.VMCH,DATA) # S1
  DATA = dataSampling_S(Val1 = DATA$SHNT, Val2 = DATA$VMCH,Prob.S.VMCH.SHNT,DATA) # S2
  DATA = dataSampling_S(Val1 = DATA$SHNT, Val2 = DATA$KINK,Prob.S.SHNT.KINK,DATA) # S3
  DATA = dataSampling_S(Val1 = DATA$PMB, Val2 = DATA$TPR,Prob.S.PMB.TPR,DATA) # S4
  DATA = dataSampling_S(Val1 = DATA$DISC, Val2 = DATA$ANES,Prob.S.DISC.ANES,DATA) # S5
  
  # Hiding variables 
  OBS.Large = DATA[,c("SHNT","VTUB","SAO2","VLNG","CCHL")]
  colnames(OBS.Large) = c("X1","Z","R","X2","Y")
  if (N >= nrow(OBS.Large)){
    N = nrow(OBS.Large)
  }
  OBS = goodSample(OBS.Large,N)
  
  return(list(OBS.Large,OBS))
}






