# library(survey)
library(xgboost)
# library(cowplot)
library(boot)
library(ipw)
library(Hmisc)
# library(ggplot2)
library(foreach)
# library(doParallel)
library(R.utils)
# library(dplyr)
# library(doSNOW)
library(arm)
library(mgcv)
library(mise)
library(tictoc)

# Log Example
## nohup taskset -c 0-15 Rscript simulation.R 'napkin' 100 500 1 30 0 '0814-2330' >log-napkin-0815-0000-mismode-0.txt & 
## nohup taskset -c 16-31 Rscript simulation.R 'planid' 100 500 1 20 0 '0814-2330' >log-planid-0815-0000-mismode-0.txt &
## nohup taskset -c 0-15 Rscript simulation.R 'napkin' 100 500 1 30 1 '0815-1100' >log-napkin-0815-1100-mismode-1.txt & 
## nohup taskset -c 16-31 Rscript simulation.R 'planid' 100 500 1 20 1 '0815-1100' >log-planid-0815-1100-mismode-1.txt & 

## nohup taskset -c 0-5 Rscript simulation.R 'napkin' 100 500 1 30 0 '0816-0200' >log-napkin-0816-0200-mismode-0.txt & 
## nohup taskset -c 6-10 Rscript simulation.R 'napkin' 100 500 1 30 1 '0816-0200' >log-napkin-0816-0200-mismode-1.txt & 
## nohup taskset -c 11-15 Rscript simulation.R 'napkin' 100 500 1 30 2 '0816-0200' >log-napkin-0816-0200-mismode-2.txt & 

## nohup taskset -c 16-21 Rscript simulation.R 'planid' 100 500 1 20 0 '0816-0200' >log-planid-0816-0200-mismode-0.txt & 
## nohup taskset -c 22-26 Rscript simulation.R 'planid' 100 500 1 20 1 '0816-0200' >log-planid-0816-0200-mismode-1.txt & 
## nohup taskset -c 27-31 Rscript simulation.R 'planid' 100 500 1 20 2 '0816-0200' >log-planid-0816-0200-mismode-2.txt & 

## nohup taskset -c 0-15 Rscript simulation.R 'napkin' 10 500 1 2 '0816-1900' >log-napkinPractice-0816-1900.txt & 
## nohup taskset -c 16-31 Rscript simulation.R 'planid' 10 500 1 2 '0816-1900' >log-planidPractice-0816-1900.txt & 

## nohup taskset -c 0-15 Rscript simulation.R 'napkin' 100 500 1 30 '0816-2300' >log-napkin-0816-2300.txt & 
## nohup taskset -c 16-31 Rscript simulation.R 'planid' 100 500 1 20 '0816-2300' >log-planid-0816-2300.txt &

## nohup taskset -c 0-15 Rscript simulation.R 'napkin' 100 500 1 30 '0817-2200' >log-napkin-0817-2200.txt &
## nohup taskset -c 16-31 Rscript simulation.R 'planid' 100 500 1 20 '0817-1830' >log-planid-0817-1830.txt & 

## nohup taskset -c 0-15 Rscript simulation.R 'napkin' 10 500 1 30 '0818-2200-tmp' >log-napkintmp-0818-2200.txt &
## nohup taskset -c 16-31 Rscript simulation.R 'planid' 10 500 1 25 '0818-2200-tmp' >log-planidtmp-0818-2200.txt & 
## nohup taskset -c 16-31 Rscript simulation.R 'planid' 10 500 1 25 '0819-0200-tmp' >log-planidtmp-0819-0200.txt & 

## nohup taskset -c 0-15 Rscript simulation.R 'napkin' 100 500 1 30 '0819-0230' >log-napkin-0819-0230.txt & # Previous Final

## nohup taskset -c 0-31 Rscript simulation-ver2.R 'planid' 100 500 1 25 '0820-0930-noreweight' >log-planid-ver2-0820-0930.txt &

## nohup taskset -c 16-31 Rscript simulation-noreweight.R 'planid' 100 500 1 20 '0821-0130-noreweight' >log-planid-noreweight-0821-0130.txt &
## nohup taskset -c 16-31 Rscript simulation-noreweight.R 'planid' 100 500 1 20 '0821-0130-noreweight' >log-planid-noreweight-0821-0130.txt &



### DML presentation -- Kennedy Noise 
## nohup taskset -c 0-15 Rscript simulation.R 'napkin' 100 500 1 20 '0908-1630' >log-napkin-0908-1630.txt & 
## nohup taskset -c 0-15 Rscript simulation.R 'napkin' 100 500 1 20 '0908-2000' >log-napkin-0908-2000.txt & 
## nohup taskset -c 0-15 Rscript simulation.R 'napkin' 100 500 1 20 '0909-0000' >log-napkin-0909-0000.txt & 
## nohup taskset -c 16-31 Rscript simulation.R 'planid' 100 500 1 20 '0908-2000' >log-planid-0908-2000.txt & 
## nohup taskset -c 0-31 Rscript simulation.R 'planid' 50 500 1 20 '0909-0800' >log-planid-0909-0800.txt & 
## nohup taskset -c 0-15 Rscript simulation.R 'napkin' 50 500 1 20 '0909-1930' >log-napkin-0909-1930.txt & 
## nohup taskset -c 16-31 Rscript simulation.R 'napkin' 50 500 1 20 '0909-2015' >log-napkin-0909-2015.txt &  // correct-PI vs. mis-DML

reportingPerformance_planid = function(OBS,answer,prediction){
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

reportingPerformance_napkin = function(OBS,answer,prediction){
  Xunique = unique(OBS$X)[order(unique(OBS$X))]
  idx = 1 
  proportion_X = rep(0,length(Xunique))
  for (xval in Xunique){
    proportion_X[idx] = nrow(subset(OBS,X==xval))/nrow(OBS)
    idx = idx + 1 
  }
  return(sum(abs(answer-prediction)*proportion_X))
}

reportingPerformanceAbsolute = function(OBS,answer,prediction){
  return(mean(abs(answer-prediction),na.rm=T))
}

args = commandArgs(trailingOnly = TRUE)
timeoutLim = 999999

probleminstance = args[1] # napkin
simRound = as.numeric(args[2]) # 100
NumUnit = as.numeric(args[3]) # 500
nidx.start = as.numeric(args[4]) # 1
nidx.end = as.numeric(args[5]) # 20
filedate = args[6] # 0811-1800

### Example
# probleminstance = 'napkin'
# simRound = 10
# NumUnit = 500
# nidx.start = 1
# nidx.end = 2
# filedate = "1817tmp"
####

Nlist = c(nidx.start:nidx.end)*NumUnit
totalNumUnit = (nidx.end-nidx.start)+1

if(probleminstance == 'napkin'){
  source('napkin-real-data.R')
  source('napkin-real-DR.R')
  source('napkin-real-WERM.R')
  source('napkin-real-plugin.R')
  source('napkin-real-asBDNaive-groundtruth.R')
}
if (probleminstance == 'planid'){
  source('planid-real-data.R')
  source('planid-real-DR.R')
  source('planid-real-WERM-noreweight.R')
  source('planid-real-plugin.R')
  source('planid-real-asBDNaive-groundtruth.R')
}

if (probleminstance == 'napkin'){
  # reportPerformance = reportingPerformance_napkin
  reportPerformance_Weight = reportingPerformance_napkin
}
if (probleminstance == 'planid'){
  # reportPerformance = reportingPerformance_planid
  reportPerformance_Weight = reportingPerformance_planid
}

# probleminstance = paste(probleminstance,sep="-")
filetitle = paste(probleminstance,filedate,sep="-")
# print(paste("Mismode:",mismode))

timeoutFun = function(Fun, mytime){
  result = withTimeout({
    Fun
  }, timeout = mytime, onTimeout = "silent")
  return(result)
}

RunFunWithTime = function(TimeFUN, EstFUN, OBS, mismode, seednum ,timelim){
  if (identical(EstFUN,PlugInEstimator)){
    printedmsg = paste("PlugIn",mismode,seednum,sep="-")
  }
  if (identical(EstFUN,DREstimator)){
    printedmsg = paste("DR",mismode,seednum,sep="-")
  }
  if (identical(EstFUN,WERMEstimator)){
    printedmsg = paste("WERM",mismode,seednum,sep="-")
  }
  tic(msg=printedmsg)  
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


mat.summary.ANSWER.mis0.weight = matrix(0,nrow=totalNumUnit,ncol=6) #5th, 25th, 50th, 75th, 95th, mean 
mat.summary.PlugIn.mis0.weight = matrix(0,nrow=totalNumUnit,ncol=6)
mat.summary.WERM.mis0.weight = matrix(0,nrow=totalNumUnit,ncol=6)
mat.summary.DR.mis0.weight = matrix(0,nrow=totalNumUnit,ncol=6)

mat.summary.ANSWER.mis1.weight = matrix(0,nrow=totalNumUnit,ncol=6) #5th, 25th, 50th, 75th, 95th, mean 
mat.summary.PlugIn.mis1.weight = matrix(0,nrow=totalNumUnit,ncol=6)
mat.summary.WERM.mis1.weight = matrix(0,nrow=totalNumUnit,ncol=6)
mat.summary.DR.mis1.weight = matrix(0,nrow=totalNumUnit,ncol=6)

mat.summary.ANSWER.mis2.weight = matrix(0,nrow=totalNumUnit,ncol=6) #5th, 25th, 50th, 75th, 95th, mean 
mat.summary.PlugIn.mis2.weight = matrix(0,nrow=totalNumUnit,ncol=6)
mat.summary.WERM.mis2.weight = matrix(0,nrow=totalNumUnit,ncol=6)
mat.summary.DR.mis2.weight = matrix(0,nrow=totalNumUnit,ncol=6)

mat.total.ANSWER.mis0.weight = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.PlugIn.mis0.weight = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.WERM.mis0.weight = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.DR.mis0.weight = matrix(0,nrow=totalNumUnit,ncol=simRound)

mat.total.ANSWER.mis1.weight = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.PlugIn.mis1.weight = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.WERM.mis1.weight = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.DR.mis1.weight = matrix(0,nrow=totalNumUnit,ncol=simRound)

mat.total.ANSWER.mis2.weight = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.PlugIn.mis2.weight = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.WERM.mis2.weight = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.DR.mis2.weight = matrix(0,nrow=totalNumUnit,ncol=simRound)

mat.summary.ANSWER.mis0.abs = matrix(0,nrow=totalNumUnit,ncol=6) #5th, 25th, 50th, 75th, 95th, mean 
mat.summary.PlugIn.mis0.abs = matrix(0,nrow=totalNumUnit,ncol=6)
mat.summary.WERM.mis0.abs = matrix(0,nrow=totalNumUnit,ncol=6)
mat.summary.DR.mis0.abs = matrix(0,nrow=totalNumUnit,ncol=6)

mat.summary.ANSWER.mis1.abs = matrix(0,nrow=totalNumUnit,ncol=6) #5th, 25th, 50th, 75th, 95th, mean 
mat.summary.PlugIn.mis1.abs = matrix(0,nrow=totalNumUnit,ncol=6)
mat.summary.WERM.mis1.abs = matrix(0,nrow=totalNumUnit,ncol=6)
mat.summary.DR.mis1.abs = matrix(0,nrow=totalNumUnit,ncol=6)

mat.summary.ANSWER.mis2.abs = matrix(0,nrow=totalNumUnit,ncol=6) #5th, 25th, 50th, 75th, 95th, mean 
mat.summary.PlugIn.mis2.abs = matrix(0,nrow=totalNumUnit,ncol=6)
mat.summary.WERM.mis2.abs = matrix(0,nrow=totalNumUnit,ncol=6)
mat.summary.DR.mis2.abs = matrix(0,nrow=totalNumUnit,ncol=6)

mat.total.ANSWER.mis0.abs = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.PlugIn.mis0.abs = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.WERM.mis0.abs = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.DR.mis0.abs = matrix(0,nrow=totalNumUnit,ncol=simRound)

mat.total.ANSWER.mis1.abs = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.PlugIn.mis1.abs = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.WERM.mis1.abs = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.DR.mis1.abs = matrix(0,nrow=totalNumUnit,ncol=simRound)

mat.total.ANSWER.mis2.abs = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.PlugIn.mis2.abs = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.WERM.mis2.abs = matrix(0,nrow=totalNumUnit,ncol=simRound)
mat.total.DR.mis2.abs = matrix(0,nrow=totalNumUnit,ncol=simRound)

Nmax = 1000

for (nidx in nidx.start:nidx.end){
  N = Nlist[nidx]
  print(N)
  
  # pb <- txtProgressBar(max = simRound, style = 3)
  # progress <- function(n) setTxtProgressBar(pb, n)
  # opts <- list(progress = progress)
  
  # val.total = c()
  # for (idx in 1:simRound){
  #   seednum = sample(1:10000000,1)
  #   tmp = dataGen(seednum,N,Nmax)
  #   DATA = tmp[[1]]
  #   OBS.Large = tmp[[2]]
  #   OBS = tmp[[3]]
  #   answer = BDNaiveEstimator(DATA)
  #   # answer2 = DRNaiveEstimator(OBS.Large)
  # 
  #   iter_result = c()
  #   print(c(N,idx,seednum))
  #   for (mismode in c(0,1,2)){
  #     PIanswer = RunFunWithTime(timeoutFun,PlugInEstimator, OBS, mismode, seednum, timeoutLim)
  #     DRanswer = RunFunWithTime(timeoutFun,DREstimator, OBS, mismode, seednum, timeoutLim)
  #     WERManswer = RunFunWithTime(timeoutFun,WERMEstimator, OBS, mismode, seednum, timeoutLim)
  # 
  #     performance_PI_Weight = reportPerformance_Weight(OBS,answer,PIanswer)
  #     performance_DR_Weight = reportPerformance_Weight(OBS,answer,DRanswer)
  #     performance_WERM_Weight = reportPerformance_Weight(OBS,answer,WERManswer)
  # 
  #     performance_PI_Abs = reportingPerformanceAbsolute(OBS,answer,PIanswer)
  #     performance_DR_Abs = reportingPerformanceAbsolute(OBS,answer,DRanswer)
  #     performance_WERM_Abs = reportingPerformanceAbsolute(OBS,answer,WERManswer)
  # 
  #     iter_result = c(iter_result,c(performance_PI_Weight, performance_DR_Weight, performance_WERM_Weight, performance_PI_Abs,performance_DR_Abs,performance_WERM_Abs))
  #   }
  #   val.total = rbind(val.total,iter_result)
  #   # iter_result = c(performance_PI, performance_DR, performance_WERM)
  #   system(paste("echo 'Progressing:",idx,"'"))
  #   # return(iter_result)
  # # }
  # }
  
  val.total = foreach(idx= 1:simRound, .combine = 'rbind',
                      .packages = c('survey', 'boot', 'ipw', 'Hmisc','R.utils','dplyr','arm','xgboost','tictoc','bnlearn')) %do% {
                        seednum = sample(1:10000000,1)
                        tmp = dataGen(seednum,N,Nmax)
                        DATA = tmp[[1]]
                        OBS.Large = tmp[[2]]
                        OBS = tmp[[3]]
                        answer = BDNaiveEstimator(DATA)
                        # answer2 = DRNaiveEstimator(OBS.Large)
                        
                        iter_result = c()
                        for (mismode in c(0,1,2)){
                          PIanswer = RunFunWithTime(timeoutFun,PlugInEstimator, OBS, mismode, seednum, timeoutLim)
                          DRanswer = RunFunWithTime(timeoutFun,DREstimator, OBS, mismode, seednum, timeoutLim)
                          # WERManswer = RunFunWithTime(timeoutFun,WERMEstimator, OBS, mismode, seednum, timeoutLim)
                          WERManswer = rep(0,length(DRanswer))
                          
                          performance_PI_Weight = reportPerformance_Weight(OBS,answer,PIanswer)
                          performance_DR_Weight = reportPerformance_Weight(OBS,answer,DRanswer)
                          performance_WERM_Weight = reportPerformance_Weight(OBS,answer,WERManswer)
                          
                          performance_PI_Abs = performance_PI_Weight
                          performance_DR_Abs = performance_DR_Weight
                          performance_WERM_Abs = performance_WERM_Weight
                          
                          # performance_PI_Abs = reportingPerformanceAbsolute(OBS,answer,PIanswer)
                          # performance_DR_Abs = reportingPerformanceAbsolute(OBS,answer,DRanswer)
                          # performance_WERM_Abs = reportingPerformanceAbsolute(OBS,answer,WERManswer)
                          
                          iter_result = c(iter_result,c(performance_PI_Weight, performance_DR_Weight, performance_WERM_Weight, performance_PI_Abs,performance_DR_Abs,performance_WERM_Abs))
                        }
                        # iter_result = c(performance_PI, performance_DR, performance_WERM)
                        system(paste("echo 'Progressing:",idx,"'"))
                        return(iter_result)
                      }
  # 
  ################ Mis0, Weight ################
  val.ANSWER.mis0.weight = rep(0,simRound)
  val.PlugIn.mis0.weight = as.numeric(val.total[,1])
  val.DR.mis0.weight = as.numeric(val.total[,2])
  val.WERM.mis0.weight = as.numeric(val.total[,3])
  
  mat.total.ANSWER.mis0.weight[nidx,] = val.ANSWER.mis0.weight
  mat.total.PlugIn.mis0.weight[nidx,] = val.PlugIn.mis0.weight
  mat.total.DR.mis0.weight[nidx,] = val.DR.mis0.weight
  mat.total.WERM.mis0.weight[nidx,] = val.WERM.mis0.weight
  
  mat.summary.ANSWER.mis0.weight[nidx,] = returnSummary(val.ANSWER.mis0.weight)
  mat.summary.PlugIn.mis0.weight[nidx,] = returnSummary(val.PlugIn.mis0.weight)
  mat.summary.DR.mis0.weight[nidx,] = returnSummary(val.DR.mis0.weight)
  mat.summary.WERM.mis0.weight[nidx,] = returnSummary(val.WERM.mis0.weight)
  
  ################ Mis0, Abs ################
  val.ANSWER.mis0.abs = rep(0,simRound)
  val.PlugIn.mis0.abs = as.numeric(val.total[,4])
  val.DR.mis0.abs = as.numeric(val.total[,5])
  val.WERM.mis0.abs = as.numeric(val.total[,6])
  
  mat.total.ANSWER.mis0.abs[nidx,] = val.ANSWER.mis0.abs
  mat.total.PlugIn.mis0.abs[nidx,] = val.PlugIn.mis0.abs
  mat.total.DR.mis0.abs[nidx,] = val.DR.mis0.abs
  mat.total.WERM.mis0.abs[nidx,] = val.WERM.mis0.abs
  
  mat.summary.ANSWER.mis0.abs[nidx,] = returnSummary(val.ANSWER.mis0.abs)
  mat.summary.PlugIn.mis0.abs[nidx,] = returnSummary(val.PlugIn.mis0.abs)
  mat.summary.DR.mis0.abs[nidx,] = returnSummary(val.DR.mis0.abs)
  mat.summary.WERM.mis0.abs[nidx,] = returnSummary(val.WERM.mis0.abs)
  
  
  ################ Mis1, Weight ################
  val.ANSWER.mis1.weight = rep(0,simRound)
  val.PlugIn.mis1.weight = as.numeric(val.total[,7])
  val.DR.mis1.weight = as.numeric(val.total[,8])
  val.WERM.mis1.weight = as.numeric(val.total[,9])
  
  mat.total.ANSWER.mis1.weight[nidx,] = val.ANSWER.mis1.weight
  mat.total.PlugIn.mis1.weight[nidx,] = val.PlugIn.mis1.weight
  mat.total.DR.mis1.weight[nidx,] = val.DR.mis1.weight
  mat.total.WERM.mis1.weight[nidx,] = val.WERM.mis1.weight
  
  mat.summary.ANSWER.mis1.weight[nidx,] = returnSummary(val.ANSWER.mis1.weight)
  mat.summary.PlugIn.mis1.weight[nidx,] = returnSummary(val.PlugIn.mis1.weight)
  mat.summary.DR.mis1.weight[nidx,] = returnSummary(val.DR.mis1.weight)
  mat.summary.WERM.mis1.weight[nidx,] = returnSummary(val.WERM.mis1.weight)
  
  ################ Mis1, Abs ################
  val.ANSWER.mis1.abs = rep(0,simRound)
  val.PlugIn.mis1.abs = as.numeric(val.total[,10])
  val.DR.mis1.abs = as.numeric(val.total[,11])
  val.WERM.mis1.abs = as.numeric(val.total[,12])
  
  mat.total.ANSWER.mis1.abs[nidx,] = val.ANSWER.mis1.abs
  mat.total.PlugIn.mis1.abs[nidx,] = val.PlugIn.mis1.abs
  mat.total.DR.mis1.abs[nidx,] = val.DR.mis1.abs
  mat.total.WERM.mis1.abs[nidx,] = val.WERM.mis1.abs
  
  mat.summary.ANSWER.mis1.abs[nidx,] = returnSummary(val.ANSWER.mis1.abs)
  mat.summary.PlugIn.mis1.abs[nidx,] = returnSummary(val.PlugIn.mis1.abs)
  mat.summary.DR.mis1.abs[nidx,] = returnSummary(val.DR.mis1.abs)
  mat.summary.WERM.mis1.abs[nidx,] = returnSummary(val.WERM.mis1.abs)
  
  ################ mis2, Weight ################
  val.ANSWER.mis2.weight = rep(0,simRound)
  val.PlugIn.mis2.weight = as.numeric(val.total[,13])
  val.DR.mis2.weight = as.numeric(val.total[,14])
  val.WERM.mis2.weight = as.numeric(val.total[,15])
  
  mat.total.ANSWER.mis2.weight[nidx,] = val.ANSWER.mis2.weight
  mat.total.PlugIn.mis2.weight[nidx,] = val.PlugIn.mis2.weight
  mat.total.DR.mis2.weight[nidx,] = val.DR.mis2.weight
  mat.total.WERM.mis2.weight[nidx,] = val.WERM.mis2.weight
  
  mat.summary.ANSWER.mis2.weight[nidx,] = returnSummary(val.ANSWER.mis2.weight)
  mat.summary.PlugIn.mis2.weight[nidx,] = returnSummary(val.PlugIn.mis2.weight)
  mat.summary.DR.mis2.weight[nidx,] = returnSummary(val.DR.mis2.weight)
  mat.summary.WERM.mis2.weight[nidx,] = returnSummary(val.WERM.mis2.weight)
  
  ################ mis2, Abs ################
  val.ANSWER.mis2.abs = rep(0,simRound)
  val.PlugIn.mis2.abs = as.numeric(val.total[,16])
  val.DR.mis2.abs = as.numeric(val.total[,17])
  val.WERM.mis2.abs = as.numeric(val.total[,18])
  
  mat.total.ANSWER.mis2.abs[nidx,] = val.ANSWER.mis2.abs
  mat.total.PlugIn.mis2.abs[nidx,] = val.PlugIn.mis2.abs
  mat.total.DR.mis2.abs[nidx,] = val.DR.mis2.abs
  mat.total.WERM.mis2.abs[nidx,] = val.WERM.mis2.abs
  
  mat.summary.ANSWER.mis2.abs[nidx,] = returnSummary(val.ANSWER.mis2.abs)
  mat.summary.PlugIn.mis2.abs[nidx,] = returnSummary(val.PlugIn.mis2.abs)
  mat.summary.DR.mis2.abs[nidx,] = returnSummary(val.DR.mis2.abs)
  mat.summary.WERM.mis2.abs[nidx,] = returnSummary(val.WERM.mis2.abs)
}

write.csv(mat.total.ANSWER.mis0.weight,paste("Result/",filetitle,"mis0-ANSWER_weight.csv",sep=""))
write.csv(mat.total.DR.mis0.weight,paste("Result/",filetitle,"mis0-DR_weight.csv",sep=""))
write.csv(mat.total.PlugIn.mis0.weight,paste("Result/",filetitle,"mis0-PlugIn_weight.csv",sep=""))
write.csv(mat.total.WERM.mis0.weight,paste("Result/",filetitle,"mis0-WERM_weight.csv",sep=""))

write.csv(mat.total.ANSWER.mis0.abs,paste("Result/",filetitle,"mis0-ANSWER_abs.csv",sep=""))
write.csv(mat.total.DR.mis0.abs,paste("Result/",filetitle,"mis0-DR_abs.csv",sep=""))
write.csv(mat.total.PlugIn.mis0.abs,paste("Result/",filetitle,"mis0-PlugIn_abs.csv",sep=""))
write.csv(mat.total.WERM.mis0.abs,paste("Result/",filetitle,"mis0-WERM_abs.csv",sep=""))

write.csv(mat.total.ANSWER.mis1.weight,paste("Result/",filetitle,"mis1-ANSWER_weight.csv",sep=""))
write.csv(mat.total.DR.mis1.weight,paste("Result/",filetitle,"mis1-DR_weight.csv",sep=""))
write.csv(mat.total.PlugIn.mis1.weight,paste("Result/",filetitle,"mis1-PlugIn_weight.csv",sep=""))
write.csv(mat.total.WERM.mis1.weight,paste("Result/",filetitle,"mis1-WERM_weight.csv",sep=""))

write.csv(mat.total.ANSWER.mis1.abs,paste("Result/",filetitle,"mis1-ANSWER_abs.csv",sep=""))
write.csv(mat.total.DR.mis1.abs,paste("Result/",filetitle,"mis1-DR_abs.csv",sep=""))
write.csv(mat.total.PlugIn.mis1.abs,paste("Result/",filetitle,"mis1-PlugIn_abs.csv",sep=""))
write.csv(mat.total.WERM.mis1.abs,paste("Result/",filetitle,"mis1-WERM_abs.csv",sep=""))

write.csv(mat.total.ANSWER.mis2.weight,paste("Result/",filetitle,"mis2-ANSWER_weight.csv",sep=""))
write.csv(mat.total.DR.mis2.weight,paste("Result/",filetitle,"mis2-DR_weight.csv",sep=""))
write.csv(mat.total.PlugIn.mis2.weight,paste("Result/",filetitle,"mis2-PlugIn_weight.csv",sep=""))
write.csv(mat.total.WERM.mis2.weight,paste("Result/",filetitle,"mis2-WERM_weight.csv",sep=""))

write.csv(mat.total.ANSWER.mis2.abs,paste("Result/",filetitle,"mis2-ANSWER_abs.csv",sep=""))
write.csv(mat.total.DR.mis2.abs,paste("Result/",filetitle,"mis2-DR_abs.csv",sep=""))
write.csv(mat.total.PlugIn.mis2.abs,paste("Result/",filetitle,"mis2-PlugIn_abs.csv",sep=""))
write.csv(mat.total.WERM.mis2.abs,paste("Result/",filetitle,"mis2-WERM_abs.csv",sep=""))

conflist = c(5,25,50,75,95,'mean')
for (idx in 1:length(conflist)){
  confval = conflist[idx]
  assign(paste('Answer.',confval,sep=""),mat.summary.ANSWER.mis0.weight[,idx])
  assign(paste('DR.',confval,sep=""),mat.summary.DR.mis0.weight[,idx])
  assign(paste('PlugIn.',confval,sep=""),mat.summary.PlugIn.mis0.weight[,idx])
  assign(paste('WERM.',confval,sep=""),mat.summary.WERM.mis0.weight[,idx])
}
df.result.mis0 = data.frame(Nlist, Answer.5, Answer.25, Answer.50, Answer.75, Answer.95, Answer.mean,
                            DR.5,DR.25,DR.50,DR.75,DR.95,DR.mean, # Global
                            PlugIn.5,PlugIn.25,PlugIn.50,PlugIn.75,PlugIn.95,PlugIn.mean, # Plugin 
                            WERM.5,WERM.25,WERM.50,WERM.75,WERM.95,WERM.mean
)
write.csv(df.result.mis0,paste("Result/",filetitle,"mis0-summary_weight.csv",sep=""))

conflist = c(5,25,50,75,95,'mean')
for (idx in 1:length(conflist)){
  confval = conflist[idx]
  assign(paste('Answer.',confval,sep=""),mat.summary.ANSWER.mis0.abs[,idx])
  assign(paste('DR.',confval,sep=""),mat.summary.DR.mis0.abs[,idx])
  assign(paste('PlugIn.',confval,sep=""),mat.summary.PlugIn.mis0.abs[,idx])
  assign(paste('WERM.',confval,sep=""),mat.summary.WERM.mis0.abs[,idx])
}
df.result.mis0 = data.frame(Nlist, Answer.5, Answer.25, Answer.50, Answer.75, Answer.95, Answer.mean,
                            DR.5,DR.25,DR.50,DR.75,DR.95,DR.mean, # Global
                            PlugIn.5,PlugIn.25,PlugIn.50,PlugIn.75,PlugIn.95,PlugIn.mean, # Plugin 
                            WERM.5,WERM.25,WERM.50,WERM.75,WERM.95,WERM.mean
)
write.csv(df.result.mis0,paste("Result/",filetitle,"mis0-summary_abs.csv",sep=""))

conflist = c(5,25,50,75,95,'mean')
for (idx in 1:length(conflist)){
  confval = conflist[idx]
  assign(paste('Answer.',confval,sep=""),mat.summary.ANSWER.mis1.weight[,idx])
  assign(paste('DR.',confval,sep=""),mat.summary.DR.mis1.weight[,idx])
  assign(paste('PlugIn.',confval,sep=""),mat.summary.PlugIn.mis1.weight[,idx])
  assign(paste('WERM.',confval,sep=""),mat.summary.WERM.mis1.weight[,idx])
}
df.result.mis1 = data.frame(Nlist, Answer.5, Answer.25, Answer.50, Answer.75, Answer.95, Answer.mean,
                            DR.5,DR.25,DR.50,DR.75,DR.95,DR.mean, # Global
                            PlugIn.5,PlugIn.25,PlugIn.50,PlugIn.75,PlugIn.95,PlugIn.mean, # Plugin 
                            WERM.5,WERM.25,WERM.50,WERM.75,WERM.95,WERM.mean
)
write.csv(df.result.mis1,paste("Result/",filetitle,"mis1-summary_weight.csv",sep=""))

conflist = c(5,25,50,75,95,'mean')
for (idx in 1:length(conflist)){
  confval = conflist[idx]
  assign(paste('Answer.',confval,sep=""),mat.summary.ANSWER.mis1.abs[,idx])
  assign(paste('DR.',confval,sep=""),mat.summary.DR.mis1.abs[,idx])
  assign(paste('PlugIn.',confval,sep=""),mat.summary.PlugIn.mis1.abs[,idx])
  assign(paste('WERM.',confval,sep=""),mat.summary.WERM.mis1.abs[,idx])
}
df.result.mis1 = data.frame(Nlist, Answer.5, Answer.25, Answer.50, Answer.75, Answer.95, Answer.mean,
                            DR.5,DR.25,DR.50,DR.75,DR.95,DR.mean, # Global
                            PlugIn.5,PlugIn.25,PlugIn.50,PlugIn.75,PlugIn.95,PlugIn.mean, # Plugin 
                            WERM.5,WERM.25,WERM.50,WERM.75,WERM.95,WERM.mean
)
write.csv(df.result.mis1,paste("Result/",filetitle,"mis1-summary_abs.csv",sep=""))


conflist = c(5,25,50,75,95,'mean')
for (idx in 1:length(conflist)){
  confval = conflist[idx]
  assign(paste('Answer.',confval,sep=""),mat.summary.ANSWER.mis2.weight[,idx])
  assign(paste('DR.',confval,sep=""),mat.summary.DR.mis2.weight[,idx])
  assign(paste('PlugIn.',confval,sep=""),mat.summary.PlugIn.mis2.weight[,idx])
  assign(paste('WERM.',confval,sep=""),mat.summary.WERM.mis2.weight[,idx])
}
df.result.mis2 = data.frame(Nlist, Answer.5, Answer.25, Answer.50, Answer.75, Answer.95, Answer.mean,
                            DR.5,DR.25,DR.50,DR.75,DR.95,DR.mean, # Global
                            PlugIn.5,PlugIn.25,PlugIn.50,PlugIn.75,PlugIn.95,PlugIn.mean, # Plugin 
                            WERM.5,WERM.25,WERM.50,WERM.75,WERM.95,WERM.mean
)
write.csv(df.result.mis2,paste("Result/",filetitle,"mis2-summary_weight.csv",sep=""))

conflist = c(5,25,50,75,95,'mean')
for (idx in 1:length(conflist)){
  confval = conflist[idx]
  assign(paste('Answer.',confval,sep=""),mat.summary.ANSWER.mis2.abs[,idx])
  assign(paste('DR.',confval,sep=""),mat.summary.DR.mis2.abs[,idx])
  assign(paste('PlugIn.',confval,sep=""),mat.summary.PlugIn.mis2.abs[,idx])
  assign(paste('WERM.',confval,sep=""),mat.summary.WERM.mis2.abs[,idx])
}
df.result.mis2 = data.frame(Nlist, Answer.5, Answer.25, Answer.50, Answer.75, Answer.95, Answer.mean,
                            DR.5,DR.25,DR.50,DR.75,DR.95,DR.mean, # Global
                            PlugIn.5,PlugIn.25,PlugIn.50,PlugIn.75,PlugIn.95,PlugIn.mean, # Plugin 
                            WERM.5,WERM.25,WERM.50,WERM.75,WERM.95,WERM.mean
)
write.csv(df.result.mis2,paste("Result/",filetitle,"mis2-summary_abs.csv",sep=""))