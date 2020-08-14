library(mise)
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

mismode = 2

seednum = sample(1:10000000,1)
N = 1000; Nmax = 1000
tmp = dataGen(seednum,N,Nmax)
OBS.Large = tmp[[1]]
# OBS = tmp[[2]]
OBS = OBS.Large
answer1 = NaiveEstimator(OBS.Large)
# answer1 = NaiveEstimator(OBS)
# answer2 = DRNaiveEstimator(OBS.Large)
answer = answer1 

PIanswer = PlugInEstimator(OBS,mismode,seednum)
DRanswer = DREstimator(OBS,mismode,seednum)
WERManswer = WERMEstimator(OBS,mismode,seednum)
# asBDanswer = asBDEstimator(OBS,mismode,seednum)

performance_PI = computePerformance(OBS.Large,answer,PIanswer)
performance_DR = computePerformance(OBS.Large,answer,DRanswer)
performance_WERM = computePerformance(OBS.Large,answer,WERManswer)
# performance_asBD = mean(abs(answer-asBDanswer))

tmp_mat = matrix(round(c(performance_PI,performance_DR,performance_WERM),3),ncol=3)
colnames(tmp_mat) = c('Plug-in','DR','WERM')
rownames(tmp_mat) = 'Error'
print(tmp_mat)
print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))
