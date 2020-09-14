# setwd("Example_1-napkin")
library(xgboost)
library(boot)
library(mise)
library(tictoc)
mise()

computePerformance = function(OBS,answer,prediction){
  Xunique = unique(OBS$X)[order(unique(OBS$X))]
  idx = 1 
  proportion_X = rep(0,length(Xunique))
  for (xval in Xunique){
    proportion_X[idx] = nrow(subset(OBS,X==xval))/nrow(OBS)
    idx = idx + 1 
  }
  return(sum(abs(answer-prediction)*proportion_X))
}

computePerformanceAbsolute = function(OBS,answer,prediction){
  return(mean(abs(answer-prediction),na.rm=T))
}

ResultingPerformanceTable_Absolute = function(PIanswer,DRanswer,answer,OBS){
  finalperformance = c()
  for (mismode in c(0,1,2)){
    performance_PI = computePerformanceAbsolute(OBS,answer,get(paste("PIanswer.",mismode,sep="")))
    performance_DR = computePerformanceAbsolute(OBS,answer,get(paste("DRanswer.",mismode,sep="")))
    # performance_WERM = computePerformanceAbsolute(OBS,answer,get(paste("WERManswer.",mismode,sep="")))
    myresult = c(performance_PI,performance_DR)
    finalperformance = rbind(finalperformance,myresult)
  }
  tmp_mat = matrix(round(finalperformance,3),ncol=2)
  colnames(tmp_mat) = c('Plug-in','DR')
  rownames(tmp_mat) = c('(abs) mis0','(abs) mis1','(abs) mis2')
  
  winner_array = c()
  for (mismode in c(0,1,2)){
    winner_name = colnames(tmp_mat)[which.min(tmp_mat[(mismode+1),])]
    winner_array = c(winner_array,winner_name)
  }
  tmp_mat = cbind(tmp_mat,winner_array)
  colnames(tmp_mat)[ncol(tmp_mat)] = "Winner"
  
  return(tmp_mat)
}

ResultingPerformanceTable_Weight = function(PIanswer,DRanswer,answer,OBS){
  finalperformance = c()
  for (mismode in c(0,1,2)){
    performance_PI = computePerformance(OBS,answer,get(paste("PIanswer.",mismode,sep="")))
    performance_DR = computePerformance(OBS,answer,get(paste("DRanswer.",mismode,sep="")))
    # performance_WERM = computePerformance(OBS,answer,get(paste("WERManswer.",mismode,sep="")))
    myresult = c(performance_PI,performance_DR)
    finalperformance = rbind(finalperformance,myresult)
  }
  tmp_mat = matrix(round(finalperformance,3),ncol=2)
  colnames(tmp_mat) = c('Plug-in','DR')
  rownames(tmp_mat) = c('(weight) mis0','(weight) mis1','(weight) mis2')
  
  winner_array = c()
  for (mismode in c(0,1,2)){
    winner_name = colnames(tmp_mat)[which.min(tmp_mat[(mismode+1),])]
    winner_array = c(winner_array,winner_name)
  }
  tmp_mat = cbind(tmp_mat,winner_array)
  colnames(tmp_mat)[ncol(tmp_mat)] = "Winner"
  
  return(tmp_mat)
  # print(paste("Mismode:",mismode))
  # print(tmp_mat)
  
  # print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))
}

################################
# Dataset Generation 
################################
source('planid-syn-data.R') # Napkin dataset generation code 
# source('napkin-WERM.R')
source('napkin-syn-plugin.R')
source('napkin-syn-DR.R')

N = 1000
Nintv = 1000000
D = 10 # Cardinality of W 
numCate = 2
C = numCate - 1

seednum = sample(1:1000000,1)
mytmp = dataGen(seednum,N,Nintv,D,C)
OBS = mytmp[[1]] # Observational dataset 
INTV = mytmp[[2]]

answer = rep(0,4)
idx = 1
yval = 1 
for (x1val in c(0,1)){
  for (x2val in c(0,1)){
    answer[idx] = nrow(subset(INTV, X1.intv==x1val & X2.intv==x2val & Y.intv==yval))/nrow(subset(INTV, X1.intv==x1val & X2.intv==x2val))
    idx = idx + 1   
  }
}

proportion_X = rep(0,length(Xunique)); idx = 1
for (xval in Xunique){
  proportion_X[idx] = nrow(subset(OBS,X==xval))/nrow(OBS)
  idx = idx + 1 
}
proportion_X = round(proportion_X,4)

for (mismode in c(0,1,2)){
  tic(); PIanswer = PlugInEstimator(OBS,D,mismode,seednum); toc(); print('Done: PIanswer')
  tic(); DRanswer = DREstimator(OBS,D,mismode,seednum); toc(); print('Done: DRanswer')
  
  assign(paste("PIanswer.",mismode,sep=""),PIanswer)
  assign(paste("DRanswer.",mismode,sep=""),DRanswer)
}
PIanswer = c(PIanswer.0,PIanswer.1,PIanswer.2)
DRanswer = c(DRanswer.0,DRanswer.1,DRanswer.2)

print("answer: Weighted")
tmp_mat = ResultingPerformanceTable_Weight(PIanswer,DRanswer,answer,OBS)
print(tmp_mat)

print("answer: Absolute")
tmp_mat = ResultingPerformanceTable_Absolute(PIanswer,DRanswer,answer,OBS)
print(tmp_mat)

