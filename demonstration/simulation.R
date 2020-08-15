library(survey)
library(xgboost)
library(cowplot)
library(boot)
library(ipw)
library(Hmisc)
library(ggplot2)
library(foreach)
library(doParallel)
library(R.utils)
library(dplyr)
library(doSNOW)
library(arm)
library(mgcv)
library(mise)
library(tictoc)

# Log Example
## nohup taskset -c 0-15 Rscript simulation.R 'napkin' 100 500 1 30 0 '0814-2330' >log-napkin-0815-0000-mismode-0.txt & 
## nohup taskset -c 16-31 Rscript simulation.R 'planid' 100 500 1 20 0 '0814-2330' >log-planid-0815-0000-mismode-0.txt & 

computePerformance_planid = function(OBS,answer,prediction){
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

computePerformance_napkin = function(OBS,answer,prediction){
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
timeoutLim = 999999

probleminstance = args[1] # napkin
simRound = as.numeric(args[2]) # 100
NumUnit = as.numeric(args[3]) # 500
nidx.start = as.numeric(args[4]) # 1
nidx.end = as.numeric(args[5]) # 20
mismode = as.numeric(args[6]) # 1
filedate = args[7] # 0811-1800

### Example
probleminstance = 'planid'
simRound = 10
NumUnit = 500
nidx.start = 1
nidx.end = 2
mismode = 0
filedate = 'tmp'
####

Nlist = c(nidx.start:nidx.end)*NumUnit
totalNumUnit = (nidx.end-nidx.start)+1

if(probleminstance == 'napkin'){
  source('napkin-real-data.R')
  source('napkin-real-DR.R')
  source('napkin-real-naive.R')
  source('napkin-real-WERM.R')
  source('napkin-real-plugin.R')
}
if (probleminstance == 'planid'){
  source('planid-real-data.R')
  source('planid-real-DR.R')
  source('planid-real-naive.R')
  source('planid-real-WERM.R')
  source('planid-real-plugin.R')
}

if (probleminstance == 'napkin'){
  computePerformance = computePerformance_napkin
}
if (probleminstance == 'planid'){
  computePerformance = computePerformance_planid
}

probleminstance = paste(probleminstance,"mismode",mismode,sep="-")
filetitle = paste(probleminstance,filedate,sep="-")
print(paste("Mismode:",mismode))

timeoutFun = function(Fun, mytime){
  result = withTimeout({
    Fun
  }, timeout = mytime, onTimeout = "silent")
  return(result)
}

RunFunWithTime = function(TimeFUN, EstFUN, OBS, mismode, seednum ,timelim){
  tic()
  estval = TimeFUN(EstFUN(OBS,mismode,seednum),timelim)
  if(is.null(estval)==T){
    estval = NA
    esttime = NA 
  }else{
    esttime = toc()
    esttime = unname(esttime$toc - esttime$tic)  
  }
  return(estval)
}

returnSummary = function(myArray){
  qt.Array = quantile(myArray,probs=c(0.05,0.25,0.5,0.75,0.95),na.rm = T)
  mean.Array = mean(myArray,na.rm = T)
  return(c(as.numeric(qt.Array),mean.Array))
}

print(probleminstance)
# cl = makeSOCKcluster(numCores,outfile='Result/log-parallel.txt')
# registerDoParallel(numCores)  # use multicore, set to the number of our cores


mat.summary.ANSWER = matrix(0,nrow=totalNumUnit,ncol=6) #5th, 25th, 50th, 75th, 95th, mean 
mat.summary.PlugIn = matrix(0,nrow=totalNumUnit,ncol=6)
mat.summary.WERM = matrix(0,nrow=totalNumUnit,ncol=6)
mat.summary.DR = matrix(0,nrow=totalNumUnit,ncol=6)

mat.total.ANSWER = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.PlugIn = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.WERM = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.DR = matrix(0,nrow=totalNumUnit,ncol=simRound)

Nmax = 1000

for (nidx in nidx.start:nidx.end){
  N = Nlist[nidx]
  print(N)
  
  # pb <- txtProgressBar(max = simRound, style = 3)
  # progress <- function(n) setTxtProgressBar(pb, n)
  # opts <- list(progress = progress)
  val.total = foreach(idx= 1:simRound, .combine = 'rbind', 
                      .packages = c('survey', 'boot', 'ipw', 'Hmisc','R.utils','dplyr','arm','xgboost','tictoc','bnlearn')) %do% {
                        seednum = sample(1:10000000,1)
                        tmp = dataGen(seednum,N,Nmax)
                        OBS.Large = tmp[[1]]
                        OBS = tmp[[2]]
                        answer1 = NaiveEstimator(OBS.Large)
                        # answer2 = DRNaiveEstimator(OBS.Large)
                        answer = answer1 
                        
                        PIanswer = RunFunWithTime(timeoutFun,PlugInEstimator, OBS, mismode, seednum, timeoutLim)
                        DRanswer = RunFunWithTime(timeoutFun,DREstimator, OBS, mismode, seednum, timeoutLim)
                        WERManswer = RunFunWithTime(timeoutFun,WERMEstimator, OBS, mismode, seednum, timeoutLim)
                        
                        performance_PI = computePerformance(OBS.Large,answer,PIanswer)
                        performance_DR = computePerformance(OBS.Large,answer,DRanswer)
                        performance_WERM = computePerformance(OBS.Large,answer,WERManswer)
                        
                        iter_result = c(performance_PI, performance_DR, performance_WERM)
                        system(paste("echo 'Progressing:",idx,"'"))
                        return(c(performance_PI, performance_DR, performance_WERM))
                      }
  
  # val.total = c()
  # for(idx in 1:simRound){
  #   seednum = sample(1:10000000,1)
  #   tmp = dataGen(seednum,N,Nmax)
  #   OBS.Large = tmp[[1]]
  #   OBS = tmp[[2]]
  #   answer1 = NaiveEstimator(OBS.Large)
  #   # answer2 = DRNaiveEstimator(OBS.Large)
  #   answer = answer1
  # 
  #   PIanswer = RunFunWithTime(timeoutFun,PlugInEstimator, OBS, mismode, seednum, timeoutLim); print(paste(N,"of",idx,"th Plug-in Done"))
  #   DRanswer = RunFunWithTime(timeoutFun,DREstimator, OBS, mismode, seednum, timeoutLim); print(paste(N,"of",idx,"th DR Done"))
  #   WERManswer = RunFunWithTime(timeoutFun,WERMEstimator, OBS, mismode, seednum, timeoutLim); print(paste(N,"of",idx,"th WERM Done"))
  # 
  #   performance_PI = computePerformance(OBS.Large,answer,PIanswer)
  #   performance_DR = computePerformance(OBS.Large,answer,DRanswer)
  #   performance_WERM = computePerformance(OBS.Large,answer,WERManswer)
  # 
  #   iter_result = c(performance_PI, performance_DR, performance_WERM)
  #   val.total = rbind(val.total,iter_result)
  #   print(paste("Processing",idx))
  # }
  # rownames(val.total) = c(1:simRound)
  # colnames(val.total) = c("PlugIn","DR","WERM")
  
  val.ANSWER = rep(0,simRound)
  val.PlugIn = as.numeric(val.total[,1])
  val.DR = as.numeric(val.total[,2])
  val.WERM = as.numeric(val.total[,3])
  
  mat.total.ANSWER[nidx,] = val.ANSWER
  mat.total.PlugIn[nidx,] = val.PlugIn
  mat.total.DR[nidx,] = val.DR
  mat.total.WERM[nidx,] = val.WERM
  
  mat.summary.ANSWER[nidx,] = returnSummary(val.ANSWER)
  mat.summary.PlugIn[nidx,] = returnSummary(val.PlugIn)
  mat.summary.DR[nidx,] = returnSummary(val.DR)
  mat.summary.WERM[nidx,] = returnSummary(val.WERM)
}

write.csv(mat.total.ANSWER,paste("Result/",filetitle,"-ANSWER.csv",sep=""))
write.csv(mat.total.DR,paste("Result/",filetitle,"-DR.csv",sep=""))
write.csv(mat.total.PlugIn,paste("Result/",filetitle,"-PlugIn.csv",sep=""))
write.csv(mat.total.WERM,paste("Result/",filetitle,"-WERM.csv",sep=""))

conflist = c(5,25,50,75,95,'mean')
for (idx in 1:length(conflist)){
  confval = conflist[idx]
  assign(paste('Answer.',confval,sep=""),mat.summary.ANSWER[,idx])
  assign(paste('DR.',confval,sep=""),mat.summary.DR[,idx])
  assign(paste('PlugIn.',confval,sep=""),mat.summary.PlugIn[,idx])
  assign(paste('WERM.',confval,sep=""),mat.summary.WERM[,idx])
}
df.result = data.frame(Nlist, Answer.5, Answer.25, Answer.50, Answer.75, Answer.95, Answer.mean,
                       DR.5,DR.25,DR.50,DR.75,DR.95,DR.mean, # Global
                       PlugIn.5,PlugIn.25,PlugIn.50,PlugIn.75,PlugIn.95,PlugIn.mean, # Plugin 
                       WERM.5,WERM.25,WERM.50,WERM.75,WERM.95,WERM.mean
)
write.csv(df.result,paste("Result/",filetitle,"-summary.csv",sep=""))


# stopCluster(cl)  




