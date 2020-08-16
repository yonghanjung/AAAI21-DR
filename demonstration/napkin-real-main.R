library(mise)
library(tictoc)
# mise()

suppressMessages(source('napkin-real-data.R'))
suppressMessages(source('napkin-real-naive.R'))
suppressMessages(source('napkin-real-plugin.R'))
suppressMessages(source('napkin-real-DR.R'))
suppressMessages(source('napkin-real-WERM.R'))
suppressMessages(source('napkin-real-DR-naive.R'))
suppressMessages(source('napkin-real-asBD-groundtruth.R'))
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
N = as.numeric(args[2]); Nmax = 1000
tmp = dataGen(seednum,N,Nmax)
DATA = tmp[[1]]
OBS.Large = tmp[[2]]
OBS = tmp[[3]]

tic(); answer.plugin = BDEstimator(DATA); toc(); print('Done: Answer 0')
tic(); answer.naive = BDNaiveEstimator(DATA); toc(); print('Done: Answer 1')


tic(); PIanswer = PlugInEstimator(OBS,mismode,seednum); toc(); print('Done: PIanswer')
tic(); DRanswer = DREstimator(OBS,mismode,seednum); toc(); print('Done: DRanswer')
tic(); WERManswer = WERMEstimator(OBS,mismode,seednum); toc(); print('Done: WERManswer')
# asBDanswer = asBDEstimator(OBS,mismode,seednum)

print("answer: plugin")
tmp_mat = ResultingPerformanceTable(PIanswer,DRanswer,WERManswer,answer.plugin,OBS)
print(tmp_mat)
print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))

print("answer: naive")
tmp_mat = ResultingPerformanceTable(PIanswer,DRanswer,WERManswer,answer.naive,OBS)
print(tmp_mat)
print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))
