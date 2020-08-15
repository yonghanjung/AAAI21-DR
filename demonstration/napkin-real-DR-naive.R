source('RID_functions.R')

DRNaiveEstimator = function(OBS,mismode){
  W = OBS[,1] # High dim surrogate
  R = OBS[,2] # Cofounder 0-numCate
  X = OBS[,3]
  Y = OBS[,4]
  
  Wunique = unique(W)[order(unique(W))]
  Runique = unique(R)[order(unique(R))]
  Xunique = unique(X)[order(unique(X))]
  Yunique = unique(Y)[order(unique(Y))]
  
  # Setting
  DATA = data.frame(cbind(W,R,X,Y))
  DATA = subset(DATA,(is.na(W) == FALSE)&(is.na(R) == FALSE)&(is.na(X) == FALSE)&(is.na(Y) == FALSE))
  Ndata = nrow(DATA)
  # Construct ProbMat 
  # Enumerate all possible values of column
  tmp = c()
  tmp = append(tmp,list(Wunique)) # W
  tmp = append(tmp, list(Runique)) # R
  tmp = append(tmp,list(Xunique)) # X
  tmp = append(tmp,list(Yunique)) # Y
  allpossible = expand.grid(tmp)
  colnames(allpossible) = c('W','R','X','Y')
  
  # Compute the Joint prob 
  for (rowidx in 1:nrow(allpossible)){
    wval = allpossible[rowidx,'W']
    rval = allpossible[rowidx,'R']
    xval = allpossible[rowidx,'X']
    yval = allpossible[rowidx,'Y']
    filterDATA = subset(DATA,W==wval&R==rval&X==xval&Y==yval)
    allpossible[rowidx,'prob'] = nrow(filterDATA)/Ndata
  }
  
  
  ComputeDR = function(xfix){
    # Choose the fixed R 
    RProb = rep(0,0,0)
    idx = 1
    for (rval in Runique){
      filtered_DATA = subset(DATA,R==rval & X==xfix)
      RProb[idx] = nrow(filtered_DATA)/nrow(DATA)
      idx = idx + 1 
    } 
    rfix = Runique[which.max(RProb)]
    # IF_1 
    
    yfix = 1 
    Ir = (R == rfix)*1
    Ix = (X == xfix)*1
    Iy = (Y == yfix)*1
    Ixy = Ix*Iy
    
    ## Compute P(x,y|R,W)
    returnArray = mapply(function(rval,wval){
      # P(x,y,R,W)
      jointprob.xfix.yfix.rval.wval = allpossible[allpossible$X==xfix & allpossible$Y==yfix & allpossible$R==rval & allpossible$W == wval,'prob']
      # P(x,R,W)
      jointprob.xfix.rval.wval = sum(allpossible[allpossible$X==xfix & allpossible$R==rval & allpossible$W == wval,'prob'])
      # P(R,W)
      jointprob.rval.wval = sum(allpossible[allpossible$R==rval & allpossible$W == wval,'prob'])
      # P(W)
      jointprob.wval = sum(allpossible[allpossible$W == wval,'prob'])
      
      # P(x,y,r,W)
      jointprob.xfix.yfix.rfix.wval = allpossible[allpossible$X==xfix & allpossible$Y==yfix & allpossible$R==rfix & allpossible$W == wval,'prob']
      # P(x,r,W)
      jointprob.xfix.rfix.wval = sum(allpossible[allpossible$X==xfix & allpossible$R==rfix & allpossible$W == wval,'prob'])
      # P(r,W)
      jointprob.rfix.wval = sum(allpossible[allpossible$R==rfix & allpossible$W == wval,'prob'])
      
      # P(x,y|R,W)
      val.prob.xfix_yfix.rval_wval = (jointprob.xfix.yfix.rval.wval/jointprob.rval.wval) 
      # P(x,y|r,W)
      val.prob.xfix_yfix.rfix_wval = (jointprob.xfix.yfix.rfix.wval/jointprob.rfix.wval) 
      # P(x|R,W)
      val.prob.xfix.rval_wval = (jointprob.xfix.rval.wval/jointprob.rval.wval) 
      # P(x|r,W)
      val.prob.xfix.rfix_wval = (jointprob.xfix.rfix.wval/jointprob.rfix.wval) 
      # P(R|W)
      val.prob.rval.wval = jointprob.rval.wval/jointprob.wval
      # P(r|W)
      val.prob.rfix.wval = jointprob.rfix.wval/jointprob.wval
      
      returnval = c(val.prob.xfix_yfix.rval_wval, # P(x,y|R,W)
                    val.prob.xfix_yfix.rfix_wval, # P(x,y|r,W)
                    val.prob.xfix.rval_wval, # P(x|R,W)
                    val.prob.xfix.rfix_wval, # P(x|r,W)
                    val.prob.rval.wval, # P(R|W)
                    val.prob.rfix.wval # P(r|W)
                    ) 
      return(returnval)
    },DATA$R,DATA$W)
    returnArray = t(returnArray)
    
    prob.xy.RW = returnArray[,1]
    prob.xy.rW = returnArray[,2]
    prob.x.RW = returnArray[,3]
    prob.x.rW = returnArray[,4]
    prob.R.W = returnArray[,5]
    prob.r.W = returnArray[,6]
    
    if (mismode == 1){
      prob.xy.RW = fix_pred(mis_pred(prob.xy.RW,distortval))
      prob.xy.rW = fix_pred(mis_pred(prob.xy.rW,distortval))
      prob.x.rW = fix_pred(mis_pred(prob.x.rW,distortval))
      prob.x.RW = fix_pred(mis_pred(prob.x.RW,distortval))
    }
    
    if (mismode == 2){
      prob.R.W = fix_pred(mis_pred(prob.R.W,distortval))
      prob.r.W = fix_pred(mis_pred(prob.r.W,distortval))
    }

    # Compute this by IPW 
    smallval = 0
    prob.xy.dor = mean((Ir*(Ixy - prob.xy.RW)/(prob.R.W+smallval)) + prob.xy.rW)
    prob.x.dor = mean((Ir*(Ix - prob.x.RW)/(prob.R.W+smallval)) + (prob.x.rW))
    prob.y.dox = prob.xy.dor/prob.x.dor
    prob.xy.dor.Truth = prob.xy.dor
    prob.x.dor.Truth = prob.x.dor
    prob.y.dox.Truth = prob.y.dox
    
    UIF_M1 = (Ir*(Ixy - prob.xy.RW)/(prob.R.W+smallval)) + prob.xy.rW 
    UIF_M2 = (Ir*(Ix - prob.x.RW)/(prob.R.W+smallval)) + (prob.x.rW)
    EIF_M1 = (Ir*(Ixy - prob.xy.RW)/(prob.R.W+smallval)) + prob.xy.rW - prob.xy.dor
    EIF_M2 = (Ir*(Ix - prob.x.RW)/(prob.R.W+smallval)) + (prob.x.rW - prob.x.dor)
    
    UIF = (1/(prob.x.dor))*(UIF_M1 - EIF_M2*prob.y.dox)
    # UIF = (UIF_M1/prob.x.dor) - (EIF_M2/prob.x.dor)*(prob.xy.dor/prob.x.dor)
    return(mean(UIF))
  }
  YxDR = rep(0,length(Xunique))
  idx = 1 
  for (xval in Xunique){
    YxDR[idx] = ComputeDR(xval)  
    idx = idx + 1 
  }
  return(YxDR)
}




