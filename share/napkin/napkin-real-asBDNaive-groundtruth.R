source('RID_functions.R')
source('WERM_Heuristic.R')

returnUniqueFromVector = function(myvector){
  return(unique(myvector)[order(unique(myvector))])
}
myDivisionTruth = function(a,b){
  val = a/b
  if (is.nan(val)){
    val = 0.5
  }
  return(val)
}

BDNaiveEstimator = function(DATA){
  ####################################################
  # Main 
  ####################################################
  K = DATA$K
  R = DATA$R 
  X = DATA$X # CO 
  Y = DATA$Y # BP 
  
  unique.K = returnUniqueFromVector(K)
  unique.R = returnUniqueFromVector(R)
  unique.X = returnUniqueFromVector(X)
  unique.Y = returnUniqueFromVector(Y)
  
  # Setting
  tmp = c()
  tmp = append(tmp,list(unique.K)) # K
  tmp = append(tmp, list(unique.R)) # R
  tmp = append(tmp,list(unique.X)) # X
  tmp = append(tmp,list(unique.Y)) # Y
  allpossible = expand.grid(tmp)
  colnames(allpossible) = c('K','R','X','Y')

  # Compute the P(Y=1 | all)
  for (rowidx in 1:nrow(allpossible)){
    # Joint prob P(Y=1,all)
    val.K = allpossible[rowidx,'K']
    val.R = allpossible[rowidx,'R']
    # val.TPR = allpossible[rowidx,'TPR']
    val.X = allpossible[rowidx,'X']
    val.Y = allpossible[rowidx,'Y']
    
    # Compute the P(Y|rest)
    filterJoint = subset(DATA, K==val.K & R==val.R & X==val.X & Y==val.Y)
    filterCondition = subset(DATA, K==val.K & R==val.R & X==val.X)
    allpossible[rowidx,'prob.Y'] = myDivisionTruth(nrow(filterJoint),nrow(filterCondition))
    
    # Compute the P(PMB,INT,KINK,VTUB)
    filterJoint2 = subset(DATA, K==val.K & R==val.R)
    allpossible[rowidx,'prob.PA'] = myDivisionTruth(nrow(filterJoint2),nrow(DATA))
  }
  allpossible[,'val1'] =  allpossible[,'prob.Y'] * allpossible[,'prob.PA']
  
  YxNaive = rep(0,length(unique.X))
  idx = 1
  for (xval in unique.X){
    YxNaive[idx] = sum(allpossible[allpossible$X == xval & allpossible$Y == 1, 'val1'],na.rm=T)
    idx = idx + 1 
  }
  return(YxNaive)
}




