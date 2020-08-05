library(mise)
mise()

simMode = c("Synthetic","Real")

mismode = 0
distortval = 0.3

source('napkin-real-data.R')
source('napkin-real-naive.R')
source('napkin-real-plugin.R')
source('napkin-real-DR.R')
source('napkin-real-WERM.R')  
seednum = sample(1:10000000,1)
N = 50; Nintv = 100000
tmp = dataGen(seednum,N,Nintv)
OBS = tmp[[1]]
answer = NaiveEstimator(OBS)
Xunique = unique(OBS$X)[order(unique(OBS$X))]
obsans = rep(0,length(Xunique))
idx = 1
yval = 1 
for (xval in Xunique){
  obsans[idx] = nrow(subset(OBS,X==xval & Y==yval))/nrow(subset(OBS,X==xval))
  idx = idx + 1 
}

PIanswer = PlugInEstimator(OBS,mismode)
DRanswer = DREstimator(OBS,mismode)
WERManswer = WERMEstimator(OBS,mismode)

performance_PI = mean(abs(answer-PIanswer))
performance_DR = mean(abs(answer-DRanswer))
performance_WERM = mean(abs(answer-WERManswer))
# performance_OBS = mean(abs(answer-obsans))

tmp_mat = matrix(round(c(performance_PI,performance_DR,performance_WERM),3),ncol=3)
colnames(tmp_mat) = c('Plug-in','DR','WERM')
rownames(tmp_mat) = 'Error'
print(tmp_mat)
print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))
