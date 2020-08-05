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
timeoutLim = 1000000
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
simRound = 100
NumUnit = 50
totalNumUnit = 20
nidx.start = 1
nidx.end = totalNumUnit
corenum = 15

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

RunFunWithTime = function(TimeFUN, EstFUN, OBS,D,numCate,timelim){
  tic()
  estval = TimeFUN(EstFUN(OBS,D,numCate),timelim)
  if(is.null(estval)==T){
    estval = NA
    esttime = NA 
  }else{
    esttime = toc()
    esttime = unname(esttime$toc - esttime$tic)  
  }
  return(c(estval,esttime))
}






