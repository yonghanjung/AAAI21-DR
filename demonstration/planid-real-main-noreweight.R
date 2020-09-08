library(mise)
library(tictoc)
# mise()

suppressMessages(source('planid-real-data.R'))
# suppressMessages(source('planid-real-naive.R'))
suppressMessages(source('planid-real-plugin.R'))
suppressMessages(source('planid-real-WERM-noreweight.R'))
suppressMessages(source('planid-real-DR.R'))
# suppressMessages(source('planid-real-DR-naive.R'))
suppressMessages(source('planid-real-asBD-groundtruth.R'))
suppressMessages(source('planid-real-asBDNaive-groundtruth.R'))

args = commandArgs(trailingOnly = TRUE)

computePerformance = function(OBS,answer,prediction){
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

computePerformanceAbsolute = function(OBS,answer,prediction){
  return(mean(abs(answer-prediction),na.rm=T))
  # X1unique = unique(OBS$X1)[order(unique(OBS$X1))]
  # X2unique = unique(OBS$X2)[order(unique(OBS$X2))]
  # idx = 1 
  # proportion_X = rep(0,length(X1unique)*length(X2unique))
  # for (x1val in X1unique){
  #   for (x2val in X2unique){
  #     proportion_X[idx] = nrow(subset(OBS,X1==x1val & X2==x2val))/nrow(OBS)
  #     idx = idx + 1 
  #   }
  # }
  # return(sum(abs(answer-prediction)*proportion_X))
}

ResultingPerformanceTable_Absolute = function(PIanswer,DRanswer,WERManswer,answer,OBS){
  finalperformance = c()
  for (mismode in c(0,1,2)){
    performance_PI = computePerformanceAbsolute(OBS,answer,get(paste("PIanswer.",mismode,sep="")))
    performance_DR = computePerformanceAbsolute(OBS,answer,get(paste("DRanswer.",mismode,sep="")))
    performance_WERM = computePerformanceAbsolute(OBS,answer,get(paste("WERManswer.",mismode,sep="")))
    myresult = c(performance_PI,performance_DR,performance_WERM)
    finalperformance = rbind(finalperformance,myresult)
  }
  tmp_mat = matrix(round(finalperformance,3),ncol=3)
  colnames(tmp_mat) = c('Plug-in','DR','WERM')
  rownames(tmp_mat) = c('(abs) mis0','(abs) mis1','(abs) mis2')
  
  winner_array = c()
  for (mismode in c(0,1,2)){
    winner_name = colnames(tmp_mat)[which.min(tmp_mat[(mismode+1),])]
    winner_array = c(winner_array,winner_name)
  }
  tmp_mat = cbind(tmp_mat,winner_array)
  colnames(tmp_mat)[ncol(tmp_mat)] = "Winner"
  
  return(tmp_mat)
}

ResultingPerformanceTable_Weight = function(PIanswer,DRanswer,WERManswer,answer,OBS){
  finalperformance = c()
  for (mismode in c(0,1,2)){
    performance_PI = computePerformance(OBS,answer,get(paste("PIanswer.",mismode,sep="")))
    performance_DR = computePerformance(OBS,answer,get(paste("DRanswer.",mismode,sep="")))
    performance_WERM = computePerformance(OBS,answer,get(paste("WERManswer.",mismode,sep="")))
    myresult = c(performance_PI,performance_DR,performance_WERM)
    finalperformance = rbind(finalperformance,myresult)
  }
  tmp_mat = matrix(round(finalperformance,3),ncol=3)
  colnames(tmp_mat) = c('Plug-in','DR','WERM')
  rownames(tmp_mat) = c('(weight) mis0','(weight) mis1','(weight) mis2')
  
  winner_array = c()
  for (mismode in c(0,1,2)){
    winner_name = colnames(tmp_mat)[which.min(tmp_mat[(mismode+1),])]
    winner_array = c(winner_array,winner_name)
  }
  tmp_mat = cbind(tmp_mat,winner_array)
  colnames(tmp_mat)[ncol(tmp_mat)] = "Winner"
  
  return(tmp_mat)
  # print(paste("Mismode:",mismode))
  # print(tmp_mat)
  
  # print(paste("Winner: ",colnames(tmp_mat)[which.min(tmp_mat)],sep=""))
}


# mismode = as.numeric(args[1])

# seednum = sample(1:10000000,1)
# N = as.numeric(args[1])
# set.seed(as.numeric(Sys.time()))

N = 10000
seednum = 8987


# mismode = 0
# N = 10000

Nmax = 1000
tmp = dataGen(seednum,N,Nmax)
DATA = tmp[[1]]
OBS.Large = tmp[[2]]
OBS = tmp[[3]]

print(c(seednum,N,nrow(OBS),nrow(OBS.Large)))

answer = BDNaiveEstimator(DATA)

X1unique = unique(OBS.Large$X1)[order(unique(OBS.Large$X1))]
X2unique = unique(OBS.Large$X2)[order(unique(OBS.Large$X2))]
proportion_X = rep(0,length(X1unique)*length(X2unique)); idx = 1
for (x1val in X1unique){
  for (x2val in X2unique){
    proportion_X[idx] = nrow(subset(OBS.Large,X1==x1val & X2==x2val))/nrow(OBS.Large)
    idx = idx + 1 
  }
}
proportion_X = round(proportion_X,4)

# answer = BDEstimator(DATA)


for (mismode in c(0,1,2)){
  tic(); PIanswer = PlugInEstimator(OBS,mismode,seednum); toc(); print('Done: PIanswer')
  tic(); DRanswer = DREstimator(OBS,mismode,seednum); toc(); print('Done: DRanswer')
  tic(); WERManswer = WERMEstimator(OBS,mismode,seednum); toc(); print('Done: WERManswer')  
  
  assign(paste("PIanswer.",mismode,sep=""),PIanswer)
  assign(paste("DRanswer.",mismode,sep=""),DRanswer)
  assign(paste("WERManswer.",mismode,sep=""),WERManswer)
}
PIanswer = c(PIanswer.0,PIanswer.1,PIanswer.2)
DRanswer = c(DRanswer.0,DRanswer.1,DRanswer.2)
WERManswer = c(WERManswer.0,WERManswer.1,WERManswer.2)

# asBDanswer = asBDEstimator(OBS,mismode,seednum)

print("answer: Weighted")
tmp_mat = ResultingPerformanceTable_Weight(PIanswer,DRanswer,WERManswer,answer,OBS)
print(tmp_mat)

print("answer: Absolute")
tmp_mat = ResultingPerformanceTable_Absolute(PIanswer,DRanswer,WERManswer,answer,OBS)
print(tmp_mat)
