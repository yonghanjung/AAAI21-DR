source('WERM_Heuristic.R')

WERMEstimator = function(OBS,distortval,mismode){
  X1 = OBS[,1] 
  Z = OBS[,2] 
  R = OBS[,3] 
  X2 = OBS[,4]  
  Y = OBS[,5]
  
  ############################################
  # DATA setup 
  ############################################
  DATA = data.frame(X1,Z,R,X2,Y)
  X1unique = unique(X1)[order(unique(X1))]
  Zunique = unique(Z)[order(unique(Z))]
  Runique = unique(R)[order(unique(R))]
  X2unique = unique(X2)[order(unique(X2))]
  Yunique = unique(Y)[order(unique(Y))]
  
  ############################################
  # Learn Prob model 
  ############################################
  # P(y | r,x1,x2,z)
  model.Y = learnXG(inVar = data.matrix(data.frame(X1,Z,R,X2)),labelval = Y, regval = rep(0,nrow(DATA)),binommode = 1)
  pred.Y = predict(model.Y,newdata=data.matrix(data.frame(X1,Z,R,X2)),type='response')
  prob.Y = pred.Y*Y + (1-pred.Y)*(1-Y)
  if (mismode == 1){
    prob.Y = fix_pred(mis_pred(prob.Y,distortval))
  }
  
  # P(r | x1)
  model.R.X1 = learnXG(inVar = data.matrix(data.frame(X1)),labelval = R, regval = rep(0,nrow(DATA)),binommode = 0)
  pred.R = t(matrix(predict(model.R.X1, newdata=as.matrix(DATA$X1), type='response'), nrow=length(Runique)))
  prob.R.X1 = mapply(function(rowidx,rval){
    return(pred.R[rowidx,rval+1])
  },c(1:nrow(DATA)),R)
  
  # P(x2 | x1,z)
  model.X2.X1Z = learnXG(inVar = data.matrix(data.frame(X1,Z)),labelval = X2, regval = rep(0,nrow(DATA)),binommode = 0)
  pred.X2 = t(matrix(predict(model.X2.X1Z, newdata=as.matrix(data.frame(X1,Z)), type='response'), nrow=length(X2unique)))
  prob.X2.X1Z = mapply(function(rowidx,x2val){
    return(pred.X2[rowidx,x2val+1])
  },c(1:nrow(DATA)),X2)
  
  # P(r) 
  prob.R.Array = rep(0,length(Runique))
  for (rval in Runique){
    prob.R.Array[rval+1] = nrow(subset(DATA,R==rval))/nrow(DATA) 
  }
  prob.R = mapply(function(rval){
    return(prob.R.Array[rval+1])
  },R)
  
  # P(x2)
  prob.X2.Array = rep(0,length(X2unique))
  for (x2val in X2unique){
    prob.X2.Array[x2val+1] = nrow(subset(DATA,X2==x2val))/nrow(DATA) 
  }
  prob.X2 = mapply(function(x2val){
    return(prob.X2.Array[x2val+1])
  },X2)
  
  # Weight SW_Y  
  SW_Y = (prob.R * prob.X2)/(prob.R.X1 * prob.X2.X1Z)
  
  # Learn P^{SW_Y}(y|r,x2)
  Ybox = rep(0,nrow(DATA))
  bootstrap_iter = 10
  for (idx in 1:bootstrap_iter){
    sampled_df = WERM_Sampler(DATA,SW_Y)
    # Learn Pw(y|r,x2)
    model.weighted.Y.rx2 = learnXG(inVar = data.matrix(data.frame(R = sampled_df$R, X2 = sampled_df$X2)),labelval = Y, regval = rep(0,length(Y)), binommode = 1)
    pred.weighted.Y.rx2 = predict(model.weighted.Y.rx2, newdata = data.matrix(data.frame(R, X2)),type='response')
    Prob.weighted.Y.rx2 = Y*pred.weighted.Y.rx2 + (1-Y)*(1-pred.weighted.Y.rx2)
    Ybox = Ybox + Prob.weighted.Y.rx2
  }
  Prob.weighted.Y.rx2 = Ybox/bootstrap_iter
  # Weight SW
  SW = prob.Y/Prob.weighted.Y.rx2
  
  ################################################################
  # Learn h and W.
  ################################################################
  regvallist = seq(0,10,by=0.2)
  lambda_W = learnHyperParam(regvallist,data.matrix(DATA),SW,0)
  learned_W = learnWdash(SW,data.matrix(DATA),lambda_W)
  lambda_h = learnHyperParam(regvallist,data.matrix(data.frame(X1,Z,X2)),Y,1)
  
  Yx = rep(0,length(X1unique)*length(X2unique))
  idx = 1
  for (x1val in X1unique){
    for (x2val in X2unique){
      Yx[idx] = WERM_Heuristic(inVar_train = data.matrix(data.frame(X1,X2,Z)), 
                     inVar_eval = data.matrix(data.frame(X1 = rep(x1val,length(X1)),X2 = rep(x2val,length(X2)),Z)), 
                     Y = Y, Ybinary = 1, lambda_h = lambda_h, learned_W = learned_W, mismode = mismode, distortval = distortval 
                     )
      idx = idx + 1 
    }
  }
  
  return(Yx)
}

# 
# source('planid-data.R')
# source('planid-param.R')
# # source('napkin-est.R')
# # library('mise')
# #
# N = 1000
# Nintv = 1000000
# D = 1
# numCate = 2
# C = numCate - 1
# mismode = 2
# distortval = 0.2
# 
# seednum = sample(1:10000000,1)
# # seednum = 1234
# mytmp = dataGen(seednum,N,Nintv,D)
# OBS = mytmp[[1]]
# INTV = mytmp[[2]]
# answer = c(mean(INTV[INTV$X1intv==0 & INTV$X2intv==0,'Yintv']),
#            mean(INTV[INTV$X1intv==0 & INTV$X2intv==1,'Yintv']),
#            mean(INTV[INTV$X1intv==1 & INTV$X2intv==0,'Yintv']),
#            mean(INTV[INTV$X1intv==1 & INTV$X2intv==1,'Yintv'])
# )
# 
# 
# # naiveanswer = timeoutFun(naiveAdj(OBS,D,numCate),timelim)
# dlanswer = dlAdj(OBS,D,distortval,mismode)
# paramanswer = paramAdj(OBS,D,distortval,mismode)
# 
# dlperformance = mean(abs(answer-dlanswer),na.rm=T)
# paramperformance = mean(abs(answer-paramanswer),na.rm=T)
# print(round(dlperformance,4))
# print(round(paramperformance,4))
# 
# # print(c(paste("Naive: ",round(naiveperformance,4),sep=""),paste("Multi: ",round(multiperformance,4),sep="")))
# # if (naiveperformance < multiperformance){
# #   print("Win: Naive")
# # }else{
# #   print("Win: Multi")
# # }





