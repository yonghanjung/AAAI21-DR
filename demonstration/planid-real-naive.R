source('RID_functions.R')
source('WERM_Heuristic.R')
NaiveEstimator = function(OBS){
  X1 = OBS[,1] 
  Z = OBS[,2] 
  R = OBS[,3] 
  X2 = OBS[,4]  
  Y = OBS[,5]
  DATA = data.frame(X1,Z,R,X2,Y)
  
  X1unique = unique(X1)[order(unique(X1))]
  Zunique = unique(Z)[order(unique(Z))]
  Runique = unique(R)[order(unique(R))]
  X2unique = unique(X2)[order(unique(X2))]
  Yunique = unique(Y)[order(unique(Y))]
  
  # Setting
  tmp = c()
  tmp = append(tmp,list(X1unique)) # X1 # RiskAversion: 0 highest 
  tmp = append(tmp, list(Zunique)) # Z
  tmp = append(tmp,list(Runique)) # R
  tmp = append(tmp,list(X2unique)) # X2 # Accident: 3 highest 
  allpossible = expand.grid(tmp)
  colnames(allpossible) = c('X1','Z','R','X2')
  
  zeroval_handle = 0
  # Compute P(Y=1 | x1,z,r,x2)
  Expect.Y = function(myallpossible,DATA,yval){
    newcol = (ncol(myallpossible)+1)
    for (x1val in X1unique){
      for (zval in Zunique){
        for(rval in Runique){
          for(x2val in X2unique){
            filtered_DATA = subset(DATA,X1==x1val & Z==zval & R==rval & X2==x2val)
            filtered_DATA_Y = subset(DATA,X1==x1val & Z==zval & R==rval & X2==x2val & Y==yval)
            if (nrow(filtered_DATA) > 0){
              probY = nrow(filtered_DATA_Y)/nrow(filtered_DATA)   
            }else{
              probY = zeroval_handle 
            }
            myallpossible[myallpossible$X1==x1val & myallpossible$Z==zval & myallpossible$R==rval & myallpossible$X2==x2val,newcol] = probY 
          }
        }
      }
    }
    colnames(myallpossible)[ncol(myallpossible)] = 'prob'
    return(myallpossible)
  }
  
  # Compute P(r|x1)
  ProbR.X1.Param_Real = function(myallpossible,DATA){
    newcol = (ncol(myallpossible)+1)
    for (rval in Runique){
      for (x1val in X1unique){
        filtered_DATA.rx1 = subset(DATA,(X1==x1val & R==rval))
        filtered_DATA.x1 = subset(DATA,(X1==x1val))
        if (nrow(filtered_DATA.rx1) > 0){
          prob_r.x1 = nrow(filtered_DATA.rx1)/nrow(filtered_DATA.x1)   
        }else{
          prob_r.x1 = zeroval_handle
        }
        myallpossible[myallpossible$X1==x1val & myallpossible$R==rval,newcol] = prob_r.x1
      }
    }
    colnames(myallpossible)[ncol(myallpossible)] = 'prob'
    return(myallpossible)
  }
  
  # Compute P(z,x1)
  ProbZXParam_Real = function(myallpossible,DATA){
    newcol = (ncol(myallpossible)+1)
    for (zval in Zunique){
      for (x1val in X1unique){
        filtered_DATA = subset(DATA,(X1==x1val & Z==zval))
        if (nrow(filtered_DATA) > 0){
          prob_x1z = nrow(filtered_DATA)/nrow(DATA)   
        }else{
          prob_x1z = zeroval_handle
        }
        myallpossible[myallpossible$X1==x1val & myallpossible$Z==zval,newcol] = prob_x1z
      }
    }
    colnames(myallpossible)[ncol(myallpossible)] = 'prob'
    return(myallpossible)
  }
  
  Ytable = allpossible
  Ytable = Expect.Y(Ytable,DATA,1)
  
  # Compute P(z,x1)
  Pzx1Table = allpossible 
  Pzx1Table = ProbZXParam_Real(Pzx1Table,DATA)
  
  # Compute P(r|x1)
  Pr.x1Table = allpossible 
  Pr.x1Table = ProbR.X1.Param_Real(Pr.x1Table,DATA)

  allpossibleOrig = allpossible
  
  ComputeVal = allpossible
  ComputeVal$val1 = Ytable$prob * Pzx1Table$prob
  ComputeVal$val2 = Pr.x1Table$prob

  ## Marginalizing over X1,Z
  tmp = c()
  tmp = append(tmp,list(Runique)) # R
  tmp = append(tmp,list(X2unique)) # X2
  allpossible.Marginover.X1Z = expand.grid(tmp)
  colnames(allpossible.Marginover.X1Z) = c('R','X2')
  allpossible.Marginover.X1Z[,(ncol(allpossible.Marginover.X1Z)+1)] = 0
  colnames(allpossible.Marginover.X1Z)[ncol(allpossible.Marginover.X1Z)] = 'val1'
  
  for (rval in Runique){
    for(x2val in X2unique){
      allpossible.Marginover.X1Z[allpossible.Marginover.X1Z$R==rval & allpossible.Marginover.X1Z$X2==x2val,'val1'] = sum(ComputeVal[ComputeVal$X2==x2val & ComputeVal$R==rval,'val1'],na.rm=T)
    }
  }
  
  ## For all R,X1
  tmp = c()
  tmp = append(tmp,list(Runique)) # R
  tmp = append(tmp,list(X1unique)) # X1
  allpossible.X1R = expand.grid(tmp)
  colnames(allpossible.X1R) = c('R','X1')
  allpossible.X1R[,(ncol(allpossible.X1R)+1)] = 0
  colnames(allpossible.X1R)[ncol(allpossible.X1R)] = 'val1'
  
  for (rval in Runique){
    for(x1val in X1unique){
      allpossible.X1R[allpossible.X1R$R==rval & allpossible.X1R$X1==x1val,'val1'] = unique(ComputeVal[ComputeVal$X1==x1val & ComputeVal$R==rval,'val2'])
    }
  }
  
  Yx = rep(0,length(X1unique)*length(X2unique))
  idx = 1 
  for (x1val in X1unique){
    for (x2val in X2unique){
      Yx[idx] = sum(allpossible.Marginover.X1Z[allpossible.Marginover.X1Z$X2==x2val,'val1'] * allpossible.X1R[allpossible.X1R$X1==x1val,'val1'])
      idx = idx + 1 
    }
  }
  return(Yx)
}

