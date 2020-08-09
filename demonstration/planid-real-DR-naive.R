source('RID_functions.R')
source('WERM_Heuristic.R')

returnUnique = function(OBS){
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
  return(list(X1unique,Zunique,Runique,X2unique,Yunique))
}

GoodSplit = function(OBS){
  totalidx = c(1:nrow(OBS))
  
  while(1){
    splitidx_1 = sample(c(1:nrow(OBS)),size=nrow(OBS)/2)
    splitidx_2 = setdiff(totalidx,splitidx_1)
    OBS_1 = OBS[splitidx_1,]
    OBS_2 = OBS[splitidx_2,]
    
    Unique_1 = returnUnique(OBS_1)
    Unique_2 = returnUnique(OBS_2)
    stopSwitch = TRUE 
    for (idx in 1:5){
      if (identical(Unique_1,Unique_2) == FALSE){
        stopSwitch = FALSE 
      }
    }
    if (stopSwitch == TRUE){
      break
    }
  }
  rownames(OBS_1) = c(1:nrow(OBS_1))
  rownames(OBS_2) = c(1:nrow(OBS_2))
  return(list(OBS_1,OBS_2))
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
    
    # Learn P(y | r,x2,X1,Z)
    pred.Y.rx2.X1Z = mapply(function(x1valiter,zvaliter,rvalfix,x2valfix){
      jointprob.yrx2x1z = sum(allpossible[allpossible$X1==x1valiter & allpossible$Z==zvaliter & allpossible$R==rvalfix & allpossible$X2 == x2valfix & allpossible$Y == yval,'prob'])
      jointprob.rx2x1z = sum(allpossible[allpossible$X1==x1valiter & allpossible$Z==zvaliter & allpossible$R==rvalfix & allpossible$X2 == x2valfix,'prob'])
      resultval = jointprob.yrx2x1z/jointprob.rx2x1z
      if(is.nan(resultval)){
        resultval = 0
      }
      return(resultval)
    },DATA$X1,DATA$Z,rep(rvalfix,nrow(DATA)),rep(x2valfix,nrow(DATA)))
    # pred.Y.rx2.X1Z = rep(0,nrow(DATA))
    
    # pred.Y.RX2.X1Z = rep(0,nrow(DATA))
    pred.Y.RX2.X1Z = mapply(function(x1valiter,zvaliter,rvaliter,x2valiter){
      jointprob.yrx2x1z = allpossible[allpossible$X1==x1valiter & allpossible$Z==zvaliter & allpossible$R==rvaliter & allpossible$X2 == x2valiter & allpossible$Y == yval,'prob']
      jointprob.rx2x1z = sum(allpossible[allpossible$X1==x1valiter & allpossible$Z==zvaliter & allpossible$R==rvaliter & allpossible$X2 == x2valiter,'prob'])
      resultval = jointprob.yrx2x1z/jointprob.rx2x1z
      if(is.nan(resultval)){
        resultval = 0
      }
      return(resultval)
    },DATA$X1,DATA$Z,DATA$R,DATA$X2)
    
    if (mismode == 1){
      pred.Y.rx2.X1Z = fix_pred(mis_pred(pred.Y.rx2.X1Z,0.2))
      pred.Y.RX2.X1Z = fix_pred(mis_pred(pred.Y.RX2.X1Z,0.2))
    }
    
    prob.R.X1 = mapply(function(rvaliter,x1valiter){
      jointprob.RX1 = sum(allpossible[allpossible$R==rvaliter & allpossible$X1 == x1valiter,'prob'])
      jointprob.X1 = sum(allpossible[allpossible$X1 == x1valiter,'prob'])
      return(jointprob.RX1/jointprob.X1)
    },DATA$R,DATA$X1)
    if (mismode == 2){
      prob.R.X1 = fix_pred(mis_pred(prob.R.X1,0.2))
    }
    
    prob.X2.ZX1 = mapply(function(x2valiter,zvaliter,x1valiter){
      jointprob.X2ZX1 = sum(allpossible[allpossible$X2==x2valiter & allpossible$Z == zvaliter & allpossible$X1 == x1valiter,'prob'])
      jointprob.ZX1 = sum(allpossible[allpossible$Z == zvaliter & allpossible$X1 == x1valiter,'prob'])
      return(jointprob.X2ZX1/jointprob.ZX1)
    },DATA$X2,DATA$Z,DATA$X1)
    
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
    
    prob.X1 = mapply(function(x1valiter){
      return(sum(allpossible[allpossible$X1==x1valiter,'prob']))
    },DATA$X1)
    
    prob.r.X1 = mapply(function(rvalfix,x1valiter){
      jointprob.rX1 = sum(allpossible[allpossible$R==rvalfix & allpossible$X1==x1valiter,'prob'])
      jointprob.X1 = sum(allpossible[allpossible$X1==x1valiter,'prob'])
      return(jointprob.rX1/jointprob.X1)
    },rep(rvalfix,nrow(DATA)),DATA$X1)
    
    prob.r.x1 = mapply(function(rvalfix,x1valfix){
      jointprob.rX1 = sum(allpossible[allpossible$R==rvalfix & allpossible$X1==x1valfix,'prob'])
      jointprob.X1 = sum(allpossible[allpossible$X1==x1valfix,'prob'])
      return(jointprob.rX1/jointprob.X1)
    },rep(rvalfix,nrow(DATA)),rep(x1valfix,nrow(DATA)))
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
