library(mise)
library(tictoc)
mise()

suppressMessages(source('planid-real-data.R'))
# suppressMessages(source('planid-real-naive.R'))
suppressMessages(source('planid-real-plugin.R'))
suppressMessages(source('planid-real-WERM.R'))
suppressMessages(source('planid-real-DR.R'))
# suppressMessages(source('planid-real-DR-naive.R'))
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

computePerformanceAbsolute = function(OBS,answer,prediction){
  return(mean(abs(answer-prediction),na.rm=T))
  # X1unique = unique(OBS$X1)[order(unique(OBS$X1))]
  # X2unique = unique(OBS$X2)[order(unique(OBS$X2))]
  # idx = 1 
  # proportion_X = rep(0,length(X1unique)*length(X2unique))
  # for (x1val in X1unique){
  #   for (x2val in X2unique){
  #     proportion_X[idx] = nrow(subset(OBS,X1==x1val & X2==x2val))/nrow(OBS)
  #     idx = idx + 1 
  #   }
  # }
  # return(sum(abs(answer-prediction)*proportion_X))
}

ResultingPerformanceTable = function(PIanswer,DRanswer,WERManswer,answer,OBS){
  performance_PI = computePerformance(OBS,answer,PIanswer)
  performance_DR = computePerformance(OBS,answer,DRanswer)
  performance_WERM = computePerformance(OBS,answer,WERManswer)

  tmp_mat = matrix(round(c(performance_PI,performance_DR,performance_WERM),3),ncol=3)
  colnames(tmp_mat) = c('Plug-in','DR','WERM')
  rownames(tmp_mat) = 'Error (weight)'
  return(tmp_mat)
  # print(paste("Mismode:",mismode))
  # print(tmp_mat)

  # print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))
}

ResultingPerformanceTable_Absolute = function(PIanswer,DRanswer,WERManswer,answer,OBS){
  performance_PI = computePerformanceAbsolute(OBS,answer,PIanswer)
  performance_DR = computePerformanceAbsolute(OBS,answer,DRanswer)
  performance_WERM = computePerformanceAbsolute(OBS,answer,WERManswer)

  tmp_mat = matrix(round(c(performance_PI,performance_DR,performance_WERM),3),ncol=3)
  colnames(tmp_mat) = c('Plug-in','DR','WERM')
  rownames(tmp_mat) = 'Error (abs)'
  return(tmp_mat)
  # print(paste("Mismode:",mismode))
  # print(tmp_mat)

  # print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))
}

mismode = as.numeric(args[1])

seednum = sample(1:10000000,1)
# seednum = 4536437
N = as.numeric(args[2])

# seednum = 12356
print(seednum)
mismode = 0
N = 10000

Nmax = 1000
tmp = dataGen(seednum,N,Nmax)
DATA = tmp[[1]]
OBS.Large = tmp[[2]]
OBS = tmp[[3]]

answer = BDNaiveEstimator(DATA)

X1unique = unique(OBS.Large$X1)[order(unique(OBS.Large$X1))]
X2unique = unique(OBS.Large$X2)[order(unique(OBS.Large$X2))]
proportion_X = rep(0,length(X1unique)*length(X2unique)); idx = 1
for (x1val in X1unique){
  for (x2val in X2unique){
    proportion_X[idx] = nrow(subset(OBS.Large,X1==x1val & X2==x2val))/nrow(OBS.Large)
    idx = idx + 1 
  }
}
proportion_X = round(proportion_X,4)

# answer = BDEstimator(DATA)


tic(); PIanswer = PlugInEstimator(OBS,mismode,seednum); toc(); print("PI done")
# PIanswer.Large = PlugInEstimator(OBS.Large,mismode,seednum)
suppressMessages(source('planid-real-DR.R'))
# tic(); DRanswer = DREstimator(OBS,mismode,seednum); toc(); print("DR done")
tic(); DRanswer = PIanswer; toc(); print("DR done")
# DRanswer.Large = DREstimator(OBS.Large,mismode,seednum)
tic(); WERManswer = WERMEstimator(OBS,mismode,seednum); toc();  print("WERM done")
# WERManswer = rep(0.5,length(PIanswer))

print("answer: naive weighted")
tmp_mat = ResultingPerformanceTable(PIanswer,DRanswer,WERManswer,answer,OBS)
print(tmp_mat)
print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))

print("answer: naive Abs")
tmp_mat = ResultingPerformanceTable_Absolute(PIanswer,DRanswer,WERManswer,answer,OBS)
print(tmp_mat)
print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))


