library(bnlearn)
library(MASS)
library(nnet)
library(mise)

dataSampling_S = function(Val1, Val2, DATA){
  copy_DATA = DATA 
  Prob.S = matrix(runif(length(unique(Val1)) * length(unique(Val2)),min=0.5,max=1)
                            ,nrow = length(unique(Val1)),ncol = length(unique(Val2)) )
  
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

dataGen = function(seednum,N,Ninv){
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
  DATA = dataSampling_S(Val1 = DATA$PMB, Val2 = DATA$VMCH, DATA) # S1
  DATA = dataSampling_S(Val1 = DATA$SHNT, Val2 = DATA$VMCH, DATA) # S2
  DATA = dataSampling_S(Val1 = DATA$SHNT, Val2 = DATA$KINK, DATA) # S3
  DATA = dataSampling_S(Val1 = DATA$PMB, Val2 = DATA$TPR, DATA) # S4
  DATA = dataSampling_S(Val1 = DATA$DISC, Val2 = DATA$ANES, DATA) # S5
  
  # Hiding variables 
  OBS = DATA[,c("SHNT","VTUB","SAO2","VLNG","CCHL")]
  colnames(OBS) = c("X1","Z","R","X2","Y")
  return(list(OBS,OBS))
}






