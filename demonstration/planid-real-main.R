library(mise)
mise()

mismode = 0
distortval = 0.3

source('planid-real-data.R')
source('planid-real-naive.R')
source('planid-real-plugin.R')
source('planid-real-WERM.R')  
source('planid-real-DR.R')

seednum = sample(1:10000000,1)
N = 1000; Nintv = 100000
tmp = dataGen(seednum,N,Nintv)
OBS = tmp[[1]]
X1unique = unique(OBS$X1)[order(unique(OBS$X1))]
X2unique = unique(OBS$X2)[order(unique(OBS$X2))]
answer = NaiveEstimator(OBS,distortval,mismode)

obsans = rep(0,length(answer))
idx = 1 
for (x1val in X1unique){
  for (x2val in X2unique){
    OBSFiltered = subset(OBS, X1 == x1val & X2 == x2val)
    if (nrow(OBSFiltered) == 0){
      obsans[idx] = 0
    }else{
      obsans[idx] = mean(OBSFiltered$Y,na.rm=T)  
    }
    idx = idx + 1 
  }
}

PIanswer = PlugInEstimator(OBS,distortval,mismode)
DRanswer = DREstimator(OBS,distortval,mismode)
WERManswer = WERMEstimator(OBS,distortval,mismode)

performance_PI = mean(abs(answer-PIanswer))
performance_DR = mean(abs(answer-DRanswer))
performance_WERM = mean(abs(answer-WERManswer))
performance_OBS = mean(abs(answer-obsans))

tmp_mat = matrix(round(c(performance_PI,performance_DR,performance_WERM,performance_OBS),3),ncol=4)
colnames(tmp_mat) = c('Plug-in','DR','WERM','OBS')
rownames(tmp_mat) = 'Error'
print(tmp_mat)
print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))
