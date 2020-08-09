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
Nintv = 10^7

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
  estval = TimeFUN(EstFUN(OBS,D,numCate),timelim)
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
registerDoParallel(numCores)  # use multicore, set to the number of our cores
Nlist = c(1:totalNumUnit)*NumUnit

mat.summary.ANSWER = matrix(0,nrow=totalN,ncol=6) #5th, 25th, 50th, 75th, 95th, mean 
mat.summary.PlugIn = matrix(0,nrow=totalN,ncol=6)
mat.summary.WERM = matrix(0,nrow=totalN,ncol=6)
mat.summary.DR = matrix(0,nrow=totalN,ncol=6)

mat.total.ANSWER = matrix(0,nrow=totalN,ncol=simRound)
mat.total.PlugIn = matrix(0,nrow=totalN,ncol=simRound)
mat.total.WERM = matrix(0,nrow=totalN,ncol=simRound)
mat.total.DR = matrix(0,nrow=totalN,ncol=simRound)

Nmax = 10000

for (nidx in nidx.start:nidx.end){
  N = Nlist[nidx]
  print(N)
  
  pb <- txtProgressBar(max = simRound, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  val.total = foreach(idx= 1:simRound, .combine = 'rbind', 
                      .packages = c('survey', 'boot', 'ipw', 'Hmisc','R.utils','dplyr','arm','xgboost','tictoc'),.options.snow = opts) %dopar% {
                        
                        seednum = sample(1:10000000,1)
                        tmp = dataGen(seednum,N,Nmax)
                        OBS.Large = tmp[[1]]
                        OBS = tmp[[2]]
                        answer = NaiveEstimator(OBS.Large)
                        
                        PIanswer = RunFunWithTime(timeoutFun,PlugInEstimator, OBS, seednum, mismode, timeoutLim)
                        DRanswer = RunFunWithTime(timeoutFun,DREstimator, OBS, seednum, mismode, timeoutLim)
                        WERManswer = RunFunWithTime(timeoutFun,WERMEstimator, OBS, seednum, mismode, timeoutLim)
                        
                        performance_PI = mean(abs(answer-PIanswer), na.rm = T)
                        performance_DR = mean(abs(answer-DRanswer), na.rm = T)
                        performance_WERM = mean(abs(answer-WERManswer), na.rm = T)
                        
                        return(performance_PI, performance_DR, performance_WERM)
                      }
  
  val.ANSWER = rep(0,simRound)
  val.PlugIn = val.total[,1]
  val.DR = val.total[,2]
  val.WERM = val.total[,3]
  
  mat.total.ANSWER[nidx,] = val.ANSWER
  mat.total.PlugIn[nidx,] = val.PlugIn
  mat.total.DR[nidx,] = val.DR
  mat.total.WERM[nidx,] = val.WERM
  
  mat.summary.ANSWER[nidx] = returnSummary(val.ANSWER)
  mat.summary.PlugIN[nidx] = returnSummary(val.PlugIn)
  mat.summary.DR[nidx] = returnSummary(val.DR)
  mat.summary.WERM[nidx] = returnSummary(val.WERM)
}




