library(mise)
library(tictoc)
# mise()

# taskset -c 22-31 Rscript napkin-real-main.R 0 10000

suppressMessages(source('napkin-real-data.R'))
suppressMessages(source('napkin-real-plugin.R'))
suppressMessages(source('napkin-real-DR.R'))
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
  performance_PI = computePerformanceAbsolute(OBS,answer,PIanswer)
  performance_DR = computePerformanceAbsolute(OBS,answer,DRanswer)
  performance_WERM = computePerformanceAbsolute(OBS,answer,WERManswer)
  
  tmp_mat = matrix(round(c(performance_PI,performance_DR,performance_WERM),3),ncol=3)
  colnames(tmp_mat) = c('Plug-in','DR','WERM')
  rownames(tmp_mat) = 'Error'
  return(tmp_mat)
}

ResultingPerformanceTable = function(PIanswer,DRanswer,WERManswer,answer,OBS){
  performance_PI = computePerformance(OBS,answer,PIanswer)
  performance_DR = computePerformance(OBS,answer,DRanswer)
  performance_WERM = computePerformance(OBS,answer,WERManswer)
  
  tmp_mat = matrix(round(c(performance_PI,performance_DR,performance_WERM),3),ncol=3)
  colnames(tmp_mat) = c('Plug-in','DR','WERM')
  rownames(tmp_mat) = 'Error'
  return(tmp_mat)
  # print(paste("Mismode:",mismode))
  # print(tmp_mat)
  
  # print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))
}

mismode = as.numeric(args[1])

seednum = sample(1:10000000,1)
# seednum = 4536437
N = as.numeric(args[2])

seednum = 9281412
mismode = 2 
N = 11500

Nmax = 1000
tmp = dataGen(seednum,N,Nmax)
DATA = tmp[[1]]
OBS.Large = tmp[[2]]
OBS = tmp[[3]]

# tic(); answer.plugin = BDEstimator(DATA); toc(); print('Done: Answer 0')
tic(); answer = BDNaiveEstimator(DATA); toc(); print('Done: Answer')

tic(); PIanswer = PlugInEstimator(OBS,mismode,seednum); toc(); print('Done: PIanswer')
tic(); DRanswer = DREstimator(OBS,mismode,seednum); toc(); print('Done: DRanswer')
tic(); WERManswer = WERMEstimator(OBS,mismode,seednum); toc(); print('Done: WERManswer')
# asBDanswer = asBDEstimator(OBS,mismode,seednum)

print("answer: Weighted")
tmp_mat = ResultingPerformanceTable(PIanswer,DRanswer,WERManswer,answer,OBS)
print(tmp_mat)
print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))

print("answer: Absolute")
tmp_mat = ResultingPerformanceTable_Absolute(PIanswer,DRanswer,WERManswer,answer,OBS)
print(tmp_mat)
print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))