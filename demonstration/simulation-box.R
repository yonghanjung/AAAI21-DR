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
## nohup taskset -c 0-20 Rscript simulation.R 'napkin' 100 15 '0812-0100' >log-napkin-0810-1300.txt & 

args = commandArgs(trailingOnly = TRUE)
cores = detectCores()
timeoutLim = 999999

probleminstance = args[1]
simRound = as.numeric(args[2]) # 100
numCores = as.numeric(args[3]) # 15
filedate = args[4]

# Example
# probleminstance = 'napkin'
# simRound = 2
# numCores = 4
# filedate = '0812'

if(probleminstance == 'napkin'){
  source('napkin-real-data.R')
  source('napkin-real-DR.R')
  source('napkin-real-naive.R')
  source('napkin-real-WERM.R')
  source('napkin-real-plugin.R')
}

# probleminstance = paste(probleminstance,"mismode",mismode,sep="-")
filetitle = paste(probleminstance,'box',filedate,sep="-")
filetitle = paste(filetitle,'.csv',sep="")

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

computePerformance_PlanID = function(OBS,answer,prediction){
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
  return(sum(abs(answer-prediction)*proportion_X,na.rm=T))
}

computePerformance_Napkin = function(OBS,answer,prediction){
  Xunique = unique(OBS$X)[order(unique(OBS$X))]
  idx = 1 
  proportion_X = rep(0,length(Xunique))
  for (xval in Xunique){
    proportion_X[idx] = nrow(subset(OBS,X==xval))/nrow(OBS)
    idx = idx + 1 
  }
  return(sum(abs(answer-prediction)*proportion_X,na.rm=T))
}

print(probleminstance)

Nmax = 1000
simResult = foreach(idx= 1:simRound, .combine = 'rbind', 
                    .packages = c('survey', 'boot', 'ipw', 'Hmisc','R.utils','dplyr','arm','xgboost','tictoc','bnlearn')) %do% {
                      
                      seednum = sample(1:10000000,1)
                      OBS = dataGen(seednum,N,Nmax)
                      # OBS = tmp[[1]]
                      # OBS = tmp[[2]]
                      # answer1 = NaiveEstimator(OBS.Large)
                      # answer2 = DRNaiveEstimator(OBS.Large)
                      answer = NaiveEstimator(OBS) 
                      
                      performancePlugIn = c(0,0,0)
                      performanceDR = c(0,0,0)
                      performanceWERM = c(0,0,0)
                      for (mismode in c(0,1,2)){
                        PIanswer = RunFunWithTime(TimeFUN = timeoutFun, EstFUN = PlugInEstimator, OBS = OBS, mismode = mismode, seednum = seednum, timelim = timeoutLim)
                        DRanswer = RunFunWithTime(TimeFUN = timeoutFun, EstFUN = DREstimator, OBS = OBS, mismode = mismode, seednum = seednum, timelim = timeoutLim)
                        WERManswer = RunFunWithTime(TimeFUN = timeoutFun, EstFUN = WERMEstimator, OBS = OBS, mismode = mismode, seednum = seednum, timelim = timeoutLim)
                        
                        if (probleminstance == 'napkin'){
                          computePerformance = computePerformance_Napkin
                        }
                        if (probleminstance == 'planid'){
                          computePerformance = computePerformance_PlanID
                        }
                        performance_PI = computePerformance(OBS.Large,answer,PIanswer)
                        performance_DR = computePerformance(OBS.Large,answer,DRanswer)
                        performance_WERM = computePerformance(OBS.Large,answer,WERManswer)
                        
                        performancePlugIn[mismode+1] = performance_PI
                        performanceDR[mismode+1] = performance_DR
                        performanceWERM[mismode+1] = performance_WERM
                      }
                      return(c(0,performancePlugIn,performanceDR,performanceWERM))  
                    }
colnames(simResult) = c('groundTruth',
                        'PI.Mis0','PI.Mis1','PI.Mis2',
                        'DR.Mis0','DR.Mis1','DR.Mis2',
                        'WERM.Mis0','WERM.Mis1','WERM.Mis2')  
simResult = data.frame(simResult)
write.csv(simResult,filetitle)
