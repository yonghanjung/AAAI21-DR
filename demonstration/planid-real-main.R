library(mise)
library(tictoc)
# mise()

suppressMessages(source('planid-real-data.R'))
suppressMessages(source('planid-real-naive.R'))
suppressMessages(source('planid-real-plugin.R'))
suppressMessages(source('planid-real-WERM.R'))
suppressMessages(source('planid-real-DR.R'))
suppressMessages(source('planid-real-DR-naive.R'))
suppressMessages(source('planid-real-asBD-groundtruth.R'))
suppressMessages(source('planid-real-asBDNaive-groundtruth.R'))

args = commandArgs(trailingOnly = TRUE)

computePerformance = function(OBS,answer,prediction){
  X1unique = unique(OBS$X1)[order(unique(OBS$X1))]
  X2unique = unique(OBS$X2)[order(unique(OBS$X2))]
  idx = 1 
  proportion_X = rep(0,length(X1unique)*length(X2unique))
  for (x1val in X1unique){
    for (x2val in X2unique){
      proportion_X[idx] = nrow(subset(OBS,X1==x1val & X2==x2val))/nrow(OBS)
      idx = idx + 1 
    }
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
# OBS = OBS.Large
answer.plugin = BDEstimator(DATA)
answer.naive = BDNaiveEstimator(DATA)
# answer1 = NaiveEstimator(OBS.Large)
# answer2 = DRNaiveEstimator(OBS,0)
# answer = (answer1+answer2)/2

tic(); PIanswer = PlugInEstimator(OBS,mismode,seednum); toc()
# PIanswer.Large = PlugInEstimator(OBS.Large,mismode,seednum)
tic(); DRanswer = DREstimator(OBS,mismode,seednum); toc()
# DRanswer.Large = DREstimator(OBS.Large,mismode,seednum)
tic(); WERManswer = WERMEstimator(OBS,mismode,seednum); toc()
# WERManswer = rep(0.5,length(PIanswer))

print("answer: plugin")
tmp_mat = ResultingPerformanceTable(PIanswer,DRanswer,WERManswer,answer.plugin,OBS)
print(tmp_mat)
print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))

print("answer: naive")
tmp_mat = ResultingPerformanceTable(PIanswer,DRanswer,WERManswer,answer.naive,OBS)
print(tmp_mat)
print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))

# performance_PI = computePerformance(OBS,answer,PIanswer)
# performance_DR = computePerformance(OBS,answer,DRanswer)
# performance_WERM = computePerformance(OBS,answer,WERManswer)
# # performance_OBS = mean(abs(answer-obsans))

# tmp_mat = matrix(round(c(performance_PI,performance_DR,performance_WERM),3),ncol=3)
# colnames(tmp_mat) = c('Plug-in','DR','WERM')
# rownames(tmp_mat) = 'Error'
# print(paste("Mismode:",mismode))

