source('WERM_Heuristic.R')
source('RID_functions.R')

WERMEstimator = function(OBS,mismode,seednum){
  X1 = OBS[,1] 
  Z = OBS[,2] 
  R = OBS[,3] 
  X2 = OBS[,4]  
  Y = OBS[,5]
  
  ############################################
  # DATA setup 
  ############################################
  DATA = data.frame(X1,Z,R,X2,Y)
  inVar = data.frame(X1,Z,R,X2)
  
  X1unique = unique(X1)[order(unique(X1))]
  Zunique = unique(Z)[order(unique(Z))]
  Runique = unique(R)[order(unique(R))]
  X2unique = unique(X2)[order(unique(X2))]
  Yunique = unique(Y)[order(unique(Y))]
  
  for (x1val in X1unique){
    for (zval in Zunique){
      variable_name = paste("DATA.X1_",x1val,".Z_",zval,sep="")
      inVar_name = paste("inVar.X1_",x1val,".Z_",zval,sep="")
      myvariable = data.frame(X1=rep(x1val,nrow(DATA)),Z=rep(zval,nrow(DATA)),R,X2,Y)
      myinVar = data.frame(X1=rep(x1val,nrow(DATA)),Z=rep(zval,nrow(DATA)),R,X2)
      assign(variable_name,myvariable)
      assign(inVar_name,myinVar)
    }
  }
  
  ############################################
  # Learn Prob model 
  ############################################
  # P(y | r,x1,x2,z)
  yvalfix = 1 
  Iy.Train = (Y==yvalfix)*1
  X2train = X2 
  Rtrain = R 
  if (mismode == 1){
    Iy.Train = distortVar(Iy.Train,seednum)
  }
  if (mismode == 2){
    Rtrain = distortVar(Rtrain,seednum)
    X2train = distortVar(X2train,seednum)
  }
  mylambda = rep(100/sqrt(nrow(DATA)),nrow(DATA))
  
  model.Y = learnXG(inVar = data.matrix(data.frame(X1,Z,R,X2)),labelval = Iy.Train, regval = mylambda,binommode = 1)
  # Prob.weighted.Y.rx2 = rep(0,nrow(DATA))
  # for (x1val in X1unique){
  #   for (zval in Zunique){
  #     myModelName = paste('model.Y.X1_',x1val,".Z_",zval,sep="")
  #     myInVarName = paste("inVar.X1_",x1val,".Z_",zval,sep="")
  #     myInVar = get(myInVarName)
  #     myModel.Y = learnXG(inVar = data.matrix(myInVar),labelval = Iy.Train, regval = rep(0,nrow(DATA)),binommode = 1) # P(Y|x1,z,R,X2)
  #     myPred.Y = predict(myModel.Y,newdata=data.matrix(myInVar),type='response')
  #     myProb.Y = myPred.Y*Y + (1-myPred.Y)*(1-Y)
  #     prob.X_1.Z = nrow(subset(DATA,X1==x1val & Z == zval))/nrow(DATA)
  #     Prob.weighted.Y.rx2 = Prob.weighted.Y.rx2 + (myProb.Y * prob.X_1.Z)
  #   }
  # }
  pred.Y = predict(model.Y,newdata=data.matrix(data.frame(X1,Z,R,X2)),type='response')
  prob.Y = pred.Y*Y + (1-pred.Y)*(1-Y)
  
  # if (mismode == 1){
  #   prob.Y = fix_pred(mis_pred(prob.Y))
  # }
  # SW = prob.Y/Prob.weighted.Y.rx2
  
  # P(r | x1)
  model.R.X1 = learnXG(inVar = data.matrix(data.frame(X1)),labelval = Rtrain, regval = mylambda, binommode = 0)
  pred.R = t(matrix(predict(model.R.X1, newdata=as.matrix(DATA$X1), type='response'), nrow=length(Runique)))
  prob.R.X1 = mapply(function(rowidx,rval){
    return(pred.R[rowidx,rval+1])
  },c(1:nrow(DATA)),R)
  
  # if (mismode == 2){
  #   prob.R.X1 = fix_pred(mis_pred(prob.R.X1))
  # }
  
  # P(x2 | x1,z)
  model.X2.X1Z = learnXG(inVar = data.matrix(data.frame(X1,Z)),labelval = X2train, regval = mylambda, binommode = 0)
  pred.X2 = t(matrix(predict(model.X2.X1Z, newdata=as.matrix(data.frame(X1,Z)), type='response'), nrow=length(X2unique)))
  prob.X2.X1Z = mapply(function(rowidx,x2val){
    return(pred.X2[rowidx,x2val+1])
  },c(1:nrow(DATA)),X2)
  
  # if (mismode == 2){
  #   prob.X2.X1Z = fix_pred(mis_pred(prob.X2.X1Z))
  # }

  # P(r) 
  prob.R.Array = rep(0,length(Runique))
  for (rval in Runique){
    prob.R.Array[rval+1] =  sum((Rtrain==rval)*1) / nrow(OBS) 
  }
  prob.R = mapply(function(rval){
    return(prob.R.Array[rval+1])
  },R)
  
  # P(x2)
  prob.X2.Array = rep(0,length(X2unique))
  for (x2val in X2unique){
    prob.X2.Array[x2val+1] = sum((X2train==x2val)*1) / nrow(OBS) 
  }
  prob.X2 = mapply(function(x2val){
    return(prob.X2.Array[x2val+1])
  },X2)
  
  # Weight SW_Y  
  SW_Y = (prob.R * prob.X2)/(prob.R.X1 * prob.X2.X1Z)
  # if (mismode == 2){
  #   SW_Y = SW_Y + 0.2  
  # }
  
  
  # Learn P^{SW_Y}(y|r,x2)
  Ybox = rep(0,nrow(DATA))
  bootstrap_iter = 1
  for (idx in 1:bootstrap_iter){
    sampled_df = WERM_Sampler(DATA,SW_Y)
    # Learn Pw(y|r,x2)
    model.weighted.Y.rx2 = learnXG(inVar = data.matrix(data.frame(R = sampled_df$R, X2 = sampled_df$X2)),labelval = Iy.Train, regval = mylambda, binommode = 1)
    pred.weighted.Y.rx2 = predict(model.weighted.Y.rx2, newdata = data.matrix(data.frame(R, X2)),type='response')
    Prob.weighted.Y.rx2 = Y*pred.weighted.Y.rx2 + (1-Y)*(1-pred.weighted.Y.rx2)
    Ybox = Ybox + Prob.weighted.Y.rx2
  }
  Prob.weighted.Y.rx2 = Ybox/bootstrap_iter
  # if (mismode == 1){
  #   Prob.weighted.Y.rx2 = fix_pred(mis_pred(Prob.weighted.Y.rx2))
  # }
  # Weight SW
  SW = prob.Y/Prob.weighted.Y.rx2
  # if (mismode == 2){
  #   SW = SW + 0.2   
  # }
  
  
  ################################################################
  # Learn h and W.
  ################################################################
  # regvallist = seq(0,10,by=0.2)
  # lambda_W = learnHyperParam(regvallist,data.matrix(DATA),SW,0,TFcontinuous = 0)
  lambda_W = mylambda
  learned_W = learnWdash(SW,data.matrix(DATA),lambda_W)
  # lambda_h = learnHyperParam(regvallist,data.matrix(data.frame(X1,Z,X2)),Y,1,TFcontinuous = 0)
  lambda_h = rep(0,nrow(DATA))
  
  Yx = rep(0,length(X1unique)*length(X2unique))
  idx = 1
  for (x1val in X1unique){
    for (x2val in X2unique){
      Yx[idx] = WERM_Heuristic(inVar_train = data.matrix(data.frame(X1,X2,Z)), 
                     inVar_eval = data.matrix(data.frame(X1 = rep(x1val,length(X1)),X2 = rep(x2val,length(X2)),Z)), 
                     Y = Iy.Train, Ybinary = 1, lambda_h = lambda_h, learned_W = learned_W 
                     )
      idx = idx + 1 
    }
  }
  
  return(Yx)
}





