source('RID_functions.R')
source('WERM_Heuristic.R')

myDivision = function(a,b){
  val = a/b
  if (is.nan(val)){
    val = 0
  }
  return(val)
}



DRNaiveEstimator = function(OBS,mismode){
  X1 = OBS[,1] 
  Z = OBS[,2] 
  R = OBS[,3] 
  X2 = OBS[,4]  
  Y = OBS[,5]
  
  X1unique = unique(X1)[order(unique(X1))]
  Zunique = unique(Z)[order(unique(Z))]
  Runique = unique(R)[order(unique(R))]
  X2unique = unique(X2)[order(unique(X2))]
  Yunique = unique(Y)[order(unique(Y))]
  
  # Setting
  DATA = data.frame(cbind(X1,Z,R,X2,Y))
  DATA = subset(DATA,(is.na(X1) == FALSE)&(is.na(Z) == FALSE)&(is.na(R) == FALSE)&(is.na(X2) == FALSE)&(is.na(Y) == FALSE))
  Ndata = nrow(DATA)
  # Construct ProbMat 
  # Enumerate all possible values of column
  tmp = c()
  tmp = append(tmp,list(X1unique)) # X1
  tmp = append(tmp,list(Zunique)) # Z
  tmp = append(tmp,list(Runique)) # R
  tmp = append(tmp,list(X2unique)) # X2
  tmp = append(tmp,list(Yunique)) # Y
  allpossible = expand.grid(tmp)
  colnames(allpossible) = c('X1','Z','R','X2','Y')
  
  # Compute the Joint prob 
  for (rowidx in 1:nrow(allpossible)){
    x1val = allpossible[rowidx,'X1']
    zval = allpossible[rowidx,'Z']
    rval = allpossible[rowidx,'R']
    x2val = allpossible[rowidx,'X2']
    yval = allpossible[rowidx,'Y']
    filterDATA = subset(DATA,X1==x1val & Z==zval & R==rval & X2==x2val & Y==yval)
    allpossible[rowidx,'prob'] = nrow(filterDATA)/Ndata
  }
  
  
  ############################################
  # Learn Models 
  ############################################
  yval = 1
  yvalfix = yval
  Iy = (Y==yval)*1
  
  distortval = 0.5
  
  compute_UIF_M1_Naive = function(DATA,allpossible,yval,rvalfix,x2valfix){
    # X1 = DATA[,1] 
    # Z = DATA[,2] 
    # R = DATA[,3] 
    # X2 = DATA[,4]  
    # Y = DATA[,5]
    
    yval = 1
    yvalfix = yval
    Iy = (Y==yval)*1
    
    
    ############################################
    # M1 =  M[y | (r,x2);(x1;z)]
    ############################################
    # rvalfix = 0 
    # x2valfix = 0
    Ir = (R==rvalfix)*1
    Ix2 = (X2==x2valfix)*1
    
    returnArray = mapply(function(x1valiter,zvaliter,rvaliter,x2valiter){
      # P(y,r,x2,X1,Z)
      jointprob.yrx2X1Z = sum(allpossible[allpossible$X1==x1valiter & allpossible$Z==zvaliter & allpossible$R==rvalfix & allpossible$X2 == x2valfix & allpossible$Y == yval,'prob'],na.rm=T)
      # P(y,R,X2,X1,Z)
      jointprob.yRX2X1Z = allpossible[allpossible$X1==x1valiter & allpossible$Z==zvaliter & allpossible$R==rvaliter & allpossible$X2 == x2valiter & allpossible$Y == yval,'prob']
      # P(r,x2,X1,Z)
      jointprob.rx2X1Z = sum(allpossible[allpossible$X1==x1valiter & allpossible$Z==zvaliter & allpossible$R==rvalfix & allpossible$X2 == x2valfix,'prob'],na.rm=T)
      # P(R,X2,X1,Z)
      jointprob.RX2X1Z = sum(allpossible[allpossible$X1==x1valiter & allpossible$Z==zvaliter & allpossible$R==rvaliter & allpossible$X2 == x2valiter,'prob'],na.rm=T)
      # P(R,X1)
      jointprob.RX1 = sum(allpossible[allpossible$R==rvaliter & allpossible$X1 == x1valiter,'prob'],na.rm=T)
      # P(X1)
      jointprob.X1 = sum(allpossible[allpossible$X1 == x1valiter,'prob'],na.rm=T)
      # P(X2,X1,Z)
      jointprob.X2ZX1 = sum(allpossible[allpossible$X2==x2valiter & allpossible$Z == zvaliter & allpossible$X1 == x1valiter,'prob'],na.rm=T)
      # P(X1,Z)
      jointprob.ZX1 = sum(allpossible[allpossible$Z == zvaliter & allpossible$X1 == x1valiter,'prob'],na.rm=T)
      
      # P(y|r,x2,X1,Z)
      condprob.y.rx2X1Z = myDivision(jointprob.yrx2X1Z,jointprob.rx2X1Z)
      # P(y|R,X2,X1,Z)
      condprob.y.RX2X1Z = myDivision(jointprob.yRX2X1Z,jointprob.RX2X1Z)
      # P(R|X1)
      condprob.R.X1 = myDivision(jointprob.RX1,jointprob.X1)
      # P(X2|X1,Z)
      condprob.X2.X1Z = myDivision(jointprob.X2ZX1,jointprob.ZX1)
      
      returnval = c(condprob.y.rx2X1Z, # P(y|r,x2,X1,Z)
                    condprob.y.RX2X1Z, # P(y|R,X2,X1,Z)
                    condprob.R.X1, # P(R|X1)
                    condprob.X2.X1Z # P(X2|X1,Z)
                    ) 
      return(returnval)
    },DATA$X1,DATA$Z,DATA$R,DATA$X2)
    returnArray = t(returnArray)
    
    pred.Y.rx2.X1Z = returnArray[,1]
    pred.Y.RX2.X1Z = returnArray[,2]
    prob.R.X1 = returnArray[,3]
    prob.X2.ZX1 = returnArray[,4]
    
    if (mismode == 1){
      pred.Y.rx2.X1Z = fix_pred(mis_pred(pred.Y.rx2.X1Z,distortval))
      pred.Y.RX2.X1Z = fix_pred(mis_pred(pred.Y.RX2.X1Z,distortval))
    }
    if (mismode == 2){
      prob.R.X1 = fix_pred(mis_pred(prob.R.X1,distortval))
    }
    
    UIF_M1 = pred.Y.rx2.X1Z + ((Ir*Ix2)/(prob.X2.ZX1*prob.R.X1))*(Iy - pred.Y.RX2.X1Z)
    # EIF_M1 = UIF_M1 - mean(UIF_M1)
    return(UIF_M1)
  }
  
  compute_UIF_M2_Naive = function(DATA,allpossible,x1valfix,rvalfix){
    X1 = DATA[,1] 
    Z = DATA[,2] 
    R = DATA[,3] 
    X2 = DATA[,4]  
    Y = DATA[,5]
    
    ############################################
    # M2 =  M[r | x1;]
    ############################################
    # x1valfix = 0 
    # rvalfix = 0 
    Ir = (R==rvalfix)*1
    Ix1 = (X1==x1valfix)*1
    
    returnArray = mapply(function(x1valiter){
        # P(X1)
        jointprob.X1 = sum(allpossible[allpossible$X1==x1valiter,'prob'],na.rm=T)
        # P(x1)
        jointprob.x1 = sum(allpossible[allpossible$X1==x1valfix,'prob'],na.rm=T)
        # P(r,X1)
        jointprob.rX1 = sum(allpossible[allpossible$R==rvalfix & allpossible$X1==x1valiter,'prob'],na.rm=T)
        # P(r,x1)
        jointprob.rx1 = sum(allpossible[allpossible$R==rvalfix & allpossible$X1==x1valfix,'prob'],na.rm=T)
        
        # P(r|X1)
        condprob.r.X1 = myDivision(jointprob.rX1,jointprob.X1)
        # P(r|x1)
        condprob.r.x1 = myDivision(jointprob.rx1,jointprob.x1)
        
        returnval = c(jointprob.X1,
                      condprob.r.X1, # P(r|X1)
                      condprob.r.x1 # P(r|x1)
        ) 
    },DATA$X1)
    returnArray = t(returnArray)
    
    prob.X1 = returnArray[,1]
    prob.r.X1 = returnArray[,2]
    prob.r.x1 = returnArray[,3]
    if (mismode == 2){
      prob.r.X1 = fix_pred(mis_pred(prob.r.X1,0.2))
      prob.r.x1 = fix_pred(mis_pred(prob.r.x1,0.2))
    }
    
    UIF_M2 = (Ix1/prob.X1)*(Ir-prob.r.X1)+prob.r.x1
    # EIF_M2 = UIF_M2 - mean(UIF_M2)
    return(UIF_M2)
  }
  
  compute_Yx = function(DATA,allpossible,yvalfix,x1valfix,x2valfix){
    ############################################
    # UIF =  \sum_{r}UIF_M1*mean(UIF_M2) + EIF_M2 * mean(UIF_M1)
    # M1 =  M[y | (r,x2);(x1;z)]
    # M2 =  M[r | x1;]
    ############################################
    UIF = rep(0,nrow(DATA))
    for (rvalfix in Runique){
      UIF_M1 = compute_UIF_M1_Naive(DATA=DATA,allpossible = allpossible, yval=yvalfix, rvalfix=rvalfix, x2valfix=x2valfix)
      EIF_M1 = UIF_M1 - mean(UIF_M1,na.rm=T)
      
      UIF_M2 = compute_UIF_M2_Naive(DATA=DATA,allpossible = allpossible, x1valfix=x1valfix,rvalfix=rvalfix)
      EIF_M2 = UIF_M2 - mean(UIF_M2,na.rm=T)
      
      UIF = UIF + UIF_M1*mean(UIF_M2,na.rm=T) + EIF_M2*mean(UIF_M1,na.rm=T)
    }
    return(mean(UIF,na.rm=T))
  }
  
  Yx = rep(0,length(X1unique)*length(X2unique))
  idx = 1 
  for (x1valfix in X1unique){
    for (x2valfix in X2unique){
      Yx[idx] = compute_Yx(DATA=DATA,allpossible = allpossible, yvalfix=yvalfix,x1valfix=x1valfix,x2valfix=x2valfix)    
      Yx[idx] = max(Yx[idx],0)
      Yx[idx] = min(Yx[idx],1)
      idx = idx + 1 
    }
  }
  return(Yx)
}
