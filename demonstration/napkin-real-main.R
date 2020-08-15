library(mise)
library(tictoc)
mise()

source('napkin-real-data.R')
source('napkin-real-naive.R')
source('napkin-real-plugin.R')
source('napkin-real-DR.R')
source('napkin-real-WERM.R')  
source('napkin-real-DR-naive.R')
source('napkin-real-asBD.R')


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

args = commandArgs(trailingOnly = TRUE)

mismode = as.numeric(args[1])
N = as.numeric(args[2])

print(N)

N = 10000
mismode = 0

seednum = sample(1:10000000,1)
Nmax = 1000
tmp = dataGen(seednum,N,Nmax)
OBS.Large = tmp[[1]]
OBS = tmp[[2]]
# OBS = OBS
tic(); answer1 = NaiveEstimator(OBS.Large); toc(); print('Done: Answer 1')
# answer1 = NaiveEstimator(OBS)
tic(); answer2 = DRNaiveEstimator(OBS.Large,0); toc(); print('Done: Answer 2')
answer = (answer1 + answer2)/2

tic(); PIanswer = PlugInEstimator(OBS,mismode,seednum); toc(); print('Done: PIanswer')

tic(); DRanswer = DREstimator(OBS,mismode,seednum); toc(); print('Done: DRanswer')

tic(); WERManswer = WERMEstimator(OBS,mismode,seednum); toc(); print('Done: WERManswer')
# asBDanswer = asBDEstimator(OBS,mismode,seednum)

performance_PI = computePerformance(OBS,answer,PIanswer)
performance_DR = computePerformance(OBS,answer,DRanswer)
performance_WERM = computePerformance(OBS,answer,WERManswer)
# performance_asBD = mean(abs(answer-asBDanswer))

tmp_mat = matrix(round(c(performance_PI,performance_DR,performance_WERM),3),ncol=3)
colnames(tmp_mat) = c('Plug-in','DR','WERM')
rownames(tmp_mat) = 'Error'
print(tmp_mat)
print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))
