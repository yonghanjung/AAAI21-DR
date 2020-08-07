library(mise)
mise()

source('napkin-real-data.R')
source('napkin-real-naive.R')
source('napkin-real-plugin.R')
source('napkin-real-DR.R')
source('napkin-real-WERM.R')  
source('napkin-real-DR-naive.R')

mismode = 0
distortval = 0.1

seednum = sample(1:10000000,1)
N = 5000; Nmax = 5000
tmp = dataGen(seednum,N,Nmax)
OBS.Large = tmp[[1]]
OBS = tmp[[2]]
# OBS = OBS.Large
# answer1 = NaiveEstimator(OBS.Large)
answer1 = NaiveEstimator(OBS.Large)
answer2 = DRNaiveEstimator(OBS.Large)
answer = (answer1 + answer2)/2

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
