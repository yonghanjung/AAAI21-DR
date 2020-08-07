library(bnlearn)
library(MASS)
library(nnet)
library(mise)

myreturnUnique = function(DATA,covariate_to_keep){
  returnList = vector("list",length(covariate_to_keep))
  idx = 1 
  for (covariate_name in covariate_to_keep){
    myCovariate = DATA[,covariate_name]
    returnList[[idx]] = unique(myCovariate)[order(unique(myCovariate))]
    idx = idx + 1 
  }
  return(returnList)
}

checkGoodSample = function(OBS,subOBS,covariate_to_keep){
  OBSunique = myreturnUnique(OBS,covariate_to_keep)
  SUBunique = myreturnUnique(subOBS,covariate_to_keep)
  stopSwitch = TRUE 
  for (idx in 1:length(covariate_to_keep)){
    if (identical(OBSunique[[idx]],SUBunique[[idx]])==FALSE){
      stopSwitch = FALSE
    }
  }
  return(stopSwitch)
}

goodSample = function(OBS,mysize){
  covariate_to_keep = colnames(OBS)
  while(1){
    stopSwitch = TRUE 
    sampled_idx = sample(c(1:nrow(OBS)),size=mysize)
    OBS.sampled = OBS[sampled_idx,]
    stopSwitch = checkGoodSample(OBS,OBS.sampled,covariate_to_keep)
    if (stopSwitch == TRUE){
      break
    }
  }
  rownames(OBS.sampled) = c(1:nrow(OBS.sampled))
  return(OBS.sampled)
}


dataSampling_S = function(Val1, Val2, Prob.S, myDATA, covariate_to_keep){
  copy_DATA = myDATA 
  
  mytheta = mapply(function(val1,val2){
    rowidx_Pi = val1+1
    colidx_Pi = val2+1
    result_val = Prob.S[rowidx_Pi,colidx_Pi]
    return(result_val)
  },Val1,Val2)
  
  while(1){
    copy_DATA = myDATA 
    copy_DATA[,'Snode'] = mytheta
    taking_idx = sapply(mytheta,function(p) rbinom(1,1,p))
    copy_DATA$taking_idx = taking_idx
    sampled_df = subset(copy_DATA,taking_idx==1)
    sampledUnique = myreturnUnique(sampled_df,covariate_to_keep)
    if(checkGoodSample(copy_DATA,sampled_df,covariate_to_keep)){
      break 
    }
  }
  selected_DATA = subset(myDATA,taking_idx==1)
  rownames(selected_DATA) = c(1:nrow(selected_DATA))
  return(selected_DATA)
}

dataGen = function(seednum,N,Nmax){
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
  covariate_to_keep = c("SHNT","VTUB","SAO2","VLNG","CCHL")
  DATA.unique = myreturnUnique(DATA,covariate_to_keep)

  # Generate S nodes 
  # Prob.S = matrix(runif(length(unique(Val1)) * length(unique(Val2)),min=0.5,max=1)
  #                           ,nrow = length(unique(Val1)),ncol = length(unique(Val2)) )
  
  # Prob.S.PMB.VMCH = matrix(runif(length(unique(DATA$PMB)) * length(unique(DATA$VMCH)),min=0.5,max=1),nrow = length(unique(DATA$PMB)),ncol = length(unique(DATA$VMCH)) )
  Prob.S.PMB.VMCH = matrix(c(0.75,0.85,0.6,0.82,0.76,0.53,0.7,0.95),nrow = length(unique(DATA$PMB)),ncol = length(unique(DATA$VMCH)) )
  
  # Prob.S.VMCH.SHNT = matrix(runif(length(unique(DATA$SHNT)) * length(unique(DATA$VMCH)),min=0.5,max=1),nrow = length(unique(DATA$SHNT)),ncol = length(unique(DATA$VMCH)) )
  Prob.S.VMCH.SHNT = matrix(c(0.85,0.54,0.8,0.65,0.76,0.53,0.82,0.85),nrow = length(unique(DATA$VMCH)),ncol = length(unique(DATA$SHNT)) )
  
  # Prob.S.SHNT.KINK = matrix(runif(length(unique(DATA$SHNT)) * length(unique(DATA$KINK)),min=0.5,max=1),nrow = length(unique(DATA$SHNT)),ncol = length(unique(DATA$KINK)) )
  Prob.S.SHNT.KINK = matrix(c(0.75,0.85,0.75,0.85),nrow = length(unique(DATA$SHNT)), ncol = length(unique(DATA$KINK)) )
  
  # Prob.S.PMB.TPR = matrix(runif(length(unique(DATA$PMB)) * length(unique(DATA$TPR)),min=0.5,max=1),nrow = length(unique(DATA$PMB)),ncol = length(unique(DATA$TPR)) )
  Prob.S.PMB.TPR = matrix(c(0.87,0.75,0.85,0.75,0.7,0.6),nrow = length(unique(DATA$PMB)),ncol = length(unique(DATA$TPR)) )
  
  # Prob.S.DISC.ANES = matrix(runif(length(unique(DATA$DISC)) * length(unique(DATA$ANES)),min=0.5,max=1),nrow = length(unique(DATA$DISC)),ncol = length(unique(DATA$ANES)) )
  Prob.S.DISC.ANES = matrix(c(0.89,0.85,0.75,0.5),nrow = length(unique(DATA$DISC)),ncol = length(unique(DATA$ANES)) )
  
  while(1){
    DATA_sampled = dataSampling_S(Val1 = DATA$PMB, Val2 = DATA$VMCH, Prob.S.PMB.VMCH,DATA, covariate_to_keep) # S1
    DATA_sampled = dataSampling_S(Val1 = DATA_sampled$VMCH, Val2 = DATA_sampled$SHNT, Prob.S.VMCH.SHNT, DATA_sampled, covariate_to_keep) # S2
    DATA_sampled = dataSampling_S(Val1 = DATA_sampled$SHNT, Val2 = DATA_sampled$KINK, Prob.S.SHNT.KINK, DATA_sampled, covariate_to_keep) # S3
    DATA_sampled = dataSampling_S(Val1 = DATA_sampled$PMB, Val2 = DATA_sampled$TPR, Prob.S.PMB.TPR, DATA_sampled, covariate_to_keep) # S4
    DATA_sampled = dataSampling_S(Val1 = DATA_sampled$DISC, Val2 = DATA_sampled$ANES, Prob.S.DISC.ANES, DATA_sampled, covariate_to_keep) # S5
    if (nrow(DATA_sampled) > Nmax){
      break 
    }
  }
  
  # Hiding variables 
  OBS.Large = DATA[,covariate_to_keep]
  colnames(OBS.Large) = c("X1","Z","R","X2","Y")
  if (N >= nrow(OBS.Large)){
    N = nrow(OBS.Large)
  }
  OBS = goodSample(OBS.Large,N)
  
  return(list(OBS.Large,OBS))
}






