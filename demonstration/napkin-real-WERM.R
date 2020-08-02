
source('WERM_Heuristic.R')


WERMEstimator = function(OBS,distortval,mismode){
  W = OBS[,1] # High dim surrogate
  R = OBS[,2] # Cofounder 0-numCate
  X = OBS[,3]
  Y = OBS[,4]
  Wunique = unique(W)[order(unique(W))]
  Runique = unique(R)[order(unique(R))]
  Xunique = unique(X)[order(unique(X))]
  Yunique = unique(Y)[order(unique(Y))]
  
  # Setting
  DATA = data.frame(W,R,X,Y)
  DATA = subset(DATA,(is.na(W) == FALSE)&(is.na(R) == FALSE)&(is.na(X) == FALSE)&(is.na(Y) == FALSE))
  Ndata = nrow(DATA)
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
  
  ### Compute P(R|W)
  model.R.W = learnXG(inVar = data.matrix(data.frame(W)),labelval = R, regval = rep(0,nrow(DATA)),binommode = 0)
  pred.R.W = predict(model.R.W,newdata=data.matrix(data.frame(W)),type='response')
  pred.R.W  = t(matrix(pred.R.W,nrow=3))
  prob.R.W = rep(0,nrow(DATA))
  for (idx in 1:nrow(DATA)){
    prob.R.W[idx] = pred.R.W[idx,(R[idx]+1)]
  }
  if (mismode == 2){
    prob.R.W = fix_pred(mis_pred(prob.R.W,distortval))
  }
  
  prob.R = mapply(function(rval){
    jointprob.R = sum(allpossible[allpossible$R==rval,'prob'])
    return(jointprob.R)
  },DATA$R)
  
  SW_importance_sampling = prob.R/prob.R.W
  
  Iy = (Y == 1)*1
  regvallist = seq(0,10,by=0.2)
  lambda_W = learnHyperParam(regvallist,data.matrix(data.frame(W=W,R=R)),SW_importance_sampling,0)
  learned_W = learnWdash(SW_importance_sampling,data.matrix(data.frame(W=W,R=R)),lambda_W)
  lambda_h = learnHyperParam(regvallist,data.matrix(data.frame(X=X)),Iy,1)
  
  YxWERM = rep(0,length(Xunique))
  idx = 1 
  for (xval in Xunique){
    # Choose the fixed R 
    RProb = rep(0,length(Xunique))
    rval_idx = 1
    for (rval in Runique){
      filtered_DATA = subset(DATA,R==rval & X==xval)
      RProb[rval_idx] = nrow(filtered_DATA)/nrow(DATA)
      rval_idx = rval_idx + 1 
    } 
    rfix = Runique[which.max(RProb)]
    YxWERM[idx] = WERM_Heuristic(inVar_train=data.frame(X=X,R=R),inVar_eval=data.frame(X=rep(xval,nrow(OBS)),R=R),Y = Iy, Ybinary = 1, lambda_h = lambda_h, learned_W=learned_W,mismode,distortval) 
    idx = idx + 1 
  }
  return(YxWERM)
}







