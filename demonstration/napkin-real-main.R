library(mise)
mise()

simMode = c("Synthetic","Real")
SimSelect = 2

mismode = 0
distortval = 0.3

if (SimSelect == 1){
  source('napkin1d-data.R')
  source('napkin-real-plugin.R')
  source('napkin-real-DR.R')
  source('napkin-real-WERM.R')  
  seednum = sample(1:10000000,1)
  N = 1000; Nintv = 100000
  tmp = dataGen(seednum,N,Nintv)
  OBS = tmp[[1]]
  INTV = tmp [[2]]
  answer = c(mean(INTV[INTV$X.intv==0,'Y.intv']),mean(INTV[INTV$X.intv==1,'Y.intv']))
  obsans = c(mean(OBS[OBS$X==0,'Y']),mean(OBS[OBS$X==1,'Y']))
}else{
  source('napkin-real-data.R')
  source('napkin-real-naive.R')
  source('napkin-real-plugin.R')
  source('napkin-real-DR.R')
  source('napkin-real-WERM.R')  
  seednum = sample(1:10000000,1)
  N = 1000; Nintv = 100000
  tmp = dataGen(seednum,N,Nintv)
  OBS = tmp[[1]]
  answer = NaiveEstimator(OBS,0,0)
  Xunique = unique(OBS$X)[order(unique(OBS$X))]
  obsans = rep(0,length(Xunique))
  idx = 1
  yval = 1 
  for (xval in Xunique){
    obsans[idx] = nrow(subset(OBS,X==xval & Y==yval))/nrow(subset(OBS,X==xval))
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
