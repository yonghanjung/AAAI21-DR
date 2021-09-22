library(mise)
library(tictoc)
# mise()

# taskset -c 22-31 Rscript napkin-real-main.R 0 10000

suppressMessages(source('napkin-real-data.R'))
suppressMessages(source('napkin-real-plugin.R'))
suppressMessages(source('napkin-real-DR.R'))
suppressMessages(source('napkin-real-vDML.R'))
suppressMessages(source('napkin-real-WERM.R'))
suppressMessages(source('napkin-real-asBDNaive-groundtruth.R'))

args = commandArgs(trailingOnly = TRUE)

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

ResultingPerformanceTable_Absolute = function(PIanswer,DRanswer,WERManswer,answer,OBS){
  finalperformance = c()
  for (mismode in c(0,1,2)){
    performance_PI = computePerformanceAbsolute(OBS,answer,get(paste("PIanswer.",mismode,sep="")))
    performance_DR = computePerformanceAbsolute(OBS,answer,get(paste("DRanswer.",mismode,sep="")))
    performance_WERM = computePerformanceAbsolute(OBS,answer,get(paste("WERManswer.",mismode,sep="")))
    myresult = c(performance_PI,performance_DR,performance_WERM)
    finalperformance = rbind(finalperformance,myresult)
  }
  tmp_mat = matrix(round(finalperformance,3),ncol=3)
  colnames(tmp_mat) = c('Plug-in','DR','WERM')
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

ResultingPerformanceTable_Weight = function(PIanswer,DRanswer,WERManswer,vDMLanswer,answer,OBS){
  finalperformance = c()
  for (mismode in c(0,1,2)){
    performance_PI = computePerformance(OBS,answer,get(paste("PIanswer.",mismode,sep="")))
    performance_DR = computePerformance(OBS,answer,get(paste("DRanswer.",mismode,sep="")))
    performance_WERM = computePerformance(OBS,answer,get(paste("WERManswer.",mismode,sep="")))
    performance_vDML = computePerformance(OBS,answer,get(paste("vDMLanswer.",mismode,sep="")))
    myresult = c(performance_PI,performance_DR,performance_WERM,performance_vDML)
    finalperformance = rbind(finalperformance,myresult)
  }
  tmp_mat = matrix(round(finalperformance,3),ncol=4)
  colnames(tmp_mat) = c('Plug-in','DR','WERM','vDML')
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

# mismode = as.numeric(args[1])

N = as.numeric(args[1])
seednum = as.numeric(args[2])
# seednum = sample(1:10000000,1)
print(c(N,seednum))
# seednum = 4536437


# seednum = 123
# mismode = 0
# N = 10000

Nmax = 10000
tmp = dataGen(seednum,N,Nmax)
DATA = tmp[[1]]
OBS.Large = tmp[[2]]
OBS = tmp[[3]]

# tic(); answer.plugin = BDEstimator(DATA); toc(); print('Done: Answer 0')
tic(); answer = BDNaiveEstimator(DATA); toc(); print('Done: Answer')
Xunique = unique(OBS.Large$X)[order(unique(OBS.Large$X))]
proportion_X = rep(0,length(Xunique)); idx = 1
for (xval in Xunique){
    proportion_X[idx] = nrow(subset(OBS.Large,X==xval))/nrow(OBS.Large)
    idx = idx + 1 
}
proportion_X = round(proportion_X,4)

for (mismode in c(0,1,2)){
  tic(); PIanswer = PlugInEstimator(OBS,mismode,seednum); toc(); print('Done: PIanswer')
  tic(); DRanswer = DREstimator(OBS,mismode,seednum); toc(); print('Done: DRanswer')
  tic(); WERManswer = WERMEstimator(OBS,mismode,seednum); toc(); print('Done: WERManswer')  
  tic(); vDMLanswer = vDMLEstimator(OBS,mismode,seednum); toc(); print('Done: vDMLanswer')  
  
  assign(paste("PIanswer.",mismode,sep=""),PIanswer)
  assign(paste("DRanswer.",mismode,sep=""),DRanswer)
  assign(paste("vDMLanswer.",mismode,sep=""),vDMLanswer)
  assign(paste("WERManswer.",mismode,sep=""),WERManswer)
}
PIanswer = c(PIanswer.0,PIanswer.1,PIanswer.2)
DRanswer = c(DRanswer.0,DRanswer.1,DRanswer.2)
WERManswer = c(WERManswer.0,WERManswer.1,WERManswer.2)
vDMLanswer = c(vDMLanswer.0,vDMLanswer.1,vDMLanswer.2)

# asBDanswer = asBDEstimator(OBS,mismode,seednum)

print("answer: Weighted")
tmp_mat = ResultingPerformanceTable_Weight(PIanswer,DRanswer,WERManswer,vDMLanswer,answer,OBS)
print(tmp_mat)

# print("answer: Absolute")
# tmp_mat = ResultingPerformanceTable_Absolute(PIanswer,DRanswer,WERManswer,answer,OBS)
# print(tmp_mat)
