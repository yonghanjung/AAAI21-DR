library(bnlearn)
library(MASS)
library(nnet)
library(mise)

dataGen = function(seednum,N,Ninv){
  set.seed(seednum)
  data(alarm)
  DATA = data.matrix(alarm) - 1 
  DATA = data.frame(DATA[,c('STKV','CCHL','HR','CO','BP','TPR','ANES')])
  
  Prob.Snode.A_K = matrix(runif(length(unique(DATA$ANES)) * length(unique(DATA$STKV)),min=0,max=1)
                            ,nrow = length(unique(DATA$ANES)),ncol = length(unique(DATA$STKV)) )
  
  # Prob.Snode.CO_HR = matrix(runif(length(unique(DATA$CO)) * length(unique(DATA$ANES)),min=0,max=1)
  #                           ,nrow = length(unique(DATA$CO)),ncol = length(unique(DATA$ANES)) )
  
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
  OBS = sampled_df[,c("W","R","X","Y")]
  return(list(OBS,OBS))
}

