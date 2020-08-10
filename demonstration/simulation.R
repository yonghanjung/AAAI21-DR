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
## Rscript simulator.R 'napkin' 20 2 30 20 1 20 5 50 'napkin_0730_2300_D20'
## nohup taskset -c 0-29 Rscript simulator.R 'napkin' 20 2 100 20 1 20 25 50 'napkin-0801-0140-D20' >log-napkin-0801-0140-D20.txt & 

args = commandArgs(trailingOnly = TRUE)
cores = detectCores()
timeoutLim = 100

# probleminstance = args[1]
# D = as.numeric(args[2])
# numCate = as.numeric(args[3])
# simRound = as.numeric(args[4])
# totalN = as.numeric(args[5])
# nidx.start = as.numeric(args[6])
# nidx.end = as.numeric(args[7])
# corenum = as.numeric(args[8])
# NumUnit = as.numeric(args[9])
# filetitle = args[10]

# Example
probleminstance = 'napkin'
simRound = 2
NumUnit = 500
totalNumUnit = 2
nidx.start = 1
nidx.end = totalNumUnit
numCores = 4

filetitle = paste(probleminstance,'-0801-0222',sep="")

if(probleminstance == 'napkin'){
  source('napkin-real-data.R')
  source('napkin-real-DR.R')
  source('napkin-real-naive.R')
  source('napkin-real-WERM.R')
  source('napkin-real-plugin.R')
}

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
registerDoParallel(numCore)  # use multicore, set to the number of our cores
Nlist = c(1:totalNumUnit)*NumUnit

mat.summary.ANSWER = matrix(0,nrow=totalNumUnit,ncol=6) #5th, 25th, 50th, 75th, 95th, mean 
mat.summary.PlugIn = matrix(0,nrow=totalNumUnit,ncol=6)
mat.summary.WERM = matrix(0,nrow=totalNumUnit,ncol=6)
mat.summary.DR = matrix(0,nrow=totalNumUnit,ncol=6)

mat.total.ANSWER = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.PlugIn = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.WERM = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.DR = matrix(0,nrow=totalNumUnit,ncol=simRound)

Nmax = 5000
mismode = 0

for (nidx in nidx.start:nidx.end){
  N = Nlist[nidx]
  print(N)
  
  # pb <- txtProgressBar(max = simRound, style = 3)
  # progress <- function(n) setTxtProgressBar(pb, n)
  # opts <- list(progress = progress)
  
  val.total = foreach(idx= 1:simRound, .combine = 'rbind', 
                      .packages = c('survey', 'boot', 'ipw', 'Hmisc','R.utils','dplyr','arm','xgboost','tictoc')) %dopar% {
                        
                        seednum = sample(1:10000000,1)
                        tmp = dataGen(seednum,N,Nmax)
                        OBS.Large = tmp[[1]]
                        OBS = tmp[[2]]
                        answer1 = NaiveEstimator(OBS.Large)
                        # answer2 = DRNaiveEstimator(OBS.Large)
                        answer = answer1 
                        
                        PIanswer = RunFunWithTime(timeoutFun,PlugInEstimator, OBS, seednum, mismode, timeoutLim)
                        DRanswer = RunFunWithTime(timeoutFun,DREstimator, OBS, seednum, mismode, timeoutLim)
                        WERManswer = RunFunWithTime(timeoutFun,WERMEstimator, OBS, seednum, mismode, timeoutLim)
                        
                        performance_PI = mean(abs(answer-PIanswer), na.rm = T)
                        performance_DR = mean(abs(answer-DRanswer), na.rm = T)
                        performance_WERM = mean(abs(answer-WERManswer), na.rm = T)
                        
                        return(c(performance_PI, performance_DR, performance_WERM))
                      }
  
  val.ANSWER = rep(0,simRound)
  val.PlugIn = val.total[,1]
  val.DR = val.total[,2]
  val.WERM = val.total[,3]
  
  mat.total.ANSWER[nidx,] = val.ANSWER
  mat.total.PlugIn[nidx,] = val.PlugIn
  mat.total.DR[nidx,] = val.DR
  mat.total.WERM[nidx,] = val.WERM
  
  mat.summary.ANSWER[nidx,] = returnSummary(val.ANSWER)
  mat.summary.PlugIn[nidx,] = returnSummary(val.PlugIn)
  mat.summary.DR[nidx,] = returnSummary(val.DR)
  mat.summary.WERM[nidx,] = returnSummary(val.WERM)
}

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
write.csv(mat.total.ANSWER,paste("Result/",filetitle,"-ANSWER.csv",sep=""))
write.csv(mat.total.DR,paste("Result/",filetitle,"-DR.csv",sep=""))
write.csv(mat.total.PlugIn,paste("Result/",filetitle,"-PlugIn.csv",sep=""))
write.csv(mat.total.WERM,paste("Result/",filetitle,"-WERM.csv",sep=""))

stopCluster(cl)  




