library(mise)
mise()

source('planid-real-data.R')
source('planid-real-naive.R')
source('planid-real-plugin.R')
source('planid-real-WERM.R')
source('planid-real-DR.R')
source('planid-real-DR-naive.R')

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

mismode = 2

seednum = sample(1:10000000,1)
# seednum = 4536437
# N = 500; Nmax = 1000
OBS = dataGen(seednum)
# OBS.Large = tmp[[1]]
# OBS = tmp[[2]]
# OBS = OBS.Large
answer1 = NaiveEstimator(OBS)
# answer2 = DRNaiveEstimator(OBS,1)
# answer = (answer1+answer2)/2
answer = answer1

PIanswer = PlugInEstimator(OBS,mismode,seednum)
# PIanswer.Large = PlugInEstimator(OBS.Large,mismode,seednum)
DRanswer = DREstimator(OBS,mismode,seednum)
# DRanswer.Large = DREstimator(OBS.Large,mismode,seednum)
WERManswer = WERMEstimator(OBS,mismode,seednum) 
# WERManswer = rep(0.5,length(PIanswer))

performance_PI = computePerformance(OBS,answer,PIanswer)
performance_DR = computePerformance(OBS,answer,DRanswer)
performance_WERM = computePerformance(OBS,answer,WERManswer)
# performance_OBS = mean(abs(answer-obsans))

tmp_mat = matrix(round(c(performance_PI,performance_DR,performance_WERM),3),ncol=3)
colnames(tmp_mat) = c('Plug-in','DR','WERM')
rownames(tmp_mat) = 'Error'
print(paste("Mismode:",mismode))
print(tmp_mat)


print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))
