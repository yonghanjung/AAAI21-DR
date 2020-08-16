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
  Ndata = nrow(DATA)
  
  PMB = DATA$PMB
  INT = DATA$INT
  KINK = DATA$KINK
  VTUB = DATA$VTUB
  SHNT = DATA$SHNT # X1
  VLNG = DATA$VLNG # X2
  Y = DATA$CCHL
  DATA = DATA[,c('PMB','INT','KINK','VTUB','SHNT','VLNG','CCHL')]
  
  PMBUnique = returnUniqueFromVector(PMB)
  INTUnique = returnUniqueFromVector(INT)
  KINKUnique = returnUniqueFromVector(KINK)
  VTUBUnique = returnUniqueFromVector(VTUB)
  SHNTUnique = returnUniqueFromVector(SHNT)
  VLNGUnique = returnUniqueFromVector(VLNG)
  YUnique = returnUniqueFromVector(Y)
  
  X1unique = SHNTUnique
  X2unique = VLNGUnique
  
  # Setting
  tmp = c()
  tmp = append(tmp,list(PMBUnique)) # PMB
  tmp = append(tmp, list(INTUnique)) # INT
  tmp = append(tmp,list(KINKUnique)) # KINK
  tmp = append(tmp,list(VTUBUnique)) # VTUB
  tmp = append(tmp,list(SHNTUnique)) # SHNT, X1
  tmp = append(tmp,list(VLNGUnique)) # VLNG, X2
  tmp = append(tmp,list(YUnique))
  allpossible = expand.grid(tmp)
  colnames(allpossible) = c('PMB','INT','KINK','VTUB','SHNT','VLNG','Y')
  
  # Compute the P(Y=1 | all)
  for (rowidx in 1:nrow(allpossible)){
    # Joint prob P(Y=1,all)
    val.PMB = allpossible[rowidx,'PMB']
    val.INT = allpossible[rowidx,'INT']
    val.KINK = allpossible[rowidx,'KINK']
    val.VTUB = allpossible[rowidx,'VTUB']
    val.SHNT = allpossible[rowidx,'SHNT']
    val.VLNG = allpossible[rowidx,'VLNG']
    val.Y = allpossible[rowidx,'Y']
    # Compute the P(Y|rest)
    filterJoint = subset(DATA,PMB==val.PMB & INT==val.INT & KINK==val.KINK & VTUB==val.VTUB & SHNT==val.SHNT & VLNG==val.VLNG & CCHL==val.Y)
    filterCondition = subset(DATA,PMB==val.PMB & INT==val.INT & KINK==val.KINK & VTUB==val.VTUB & SHNT==val.SHNT & VLNG==val.VLNG )
    allpossible[rowidx,'prob.Y'] = myDivisionTruth(nrow(filterJoint),nrow(filterCondition))
    
    # Compute the P(PMB,INT,KINK,VTUB)
    filterJoint2 = subset(DATA,PMB==val.PMB & INT==val.INT & KINK==val.KINK & VTUB==val.VTUB)
    allpossible[rowidx,'prob.PA'] = myDivisionTruth(nrow(filterJoint2),nrow(DATA))
  }
  allpossible[,'val1'] =  allpossible[,'prob.Y'] * allpossible[,'prob.PA']
  
  YxNaive = rep(0,length(X1unique)*length(X2unique))
  idx = 1
  for (x1val in X1unique){
    for (x2val in X2unique){
      YxNaive[idx] = sum(allpossible[allpossible$SHNT == x1val & allpossible$VLNG == x2val & allpossible$Y == 1, 'val1'],na.rm=T)
      idx = idx + 1   
    }
  }
  return(YxNaive)
}




