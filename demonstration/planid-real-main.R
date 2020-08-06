library(mise)
mise()

mismode = 0
distortval = 0.3

source('planid-real-data.R')
source('planid-real-naive.R')
source('planid-real-plugin.R')
source('planid-real-WERM.R')  
source('planid-real-DR.R')

mismode = 0
# distortval = 0.1

seednum = sample(1:10000000,1)
N = 1000; Nmax = 5000
tmp = dataGen(seednum,N,Nmax)
OBS.Large = tmp[[1]]
OBS = tmp[[2]]
answer = NaiveEstimator(OBS.Large)

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
