myXOR = function(W){
  myval = xor(W[,2],W[,1])
  for(idx in 3:(ncol(W))){
    myval = xor(myval,W[,idx])
  }
  return(myval)
}

dataGen = function(seednum,N,Nintv,D,C){
  # W should be high dim 
  
  set.seed(seednum)
  varval = 2
  
  c1 = rnorm(D,1,1)
  c2 = rnorm(D,-2,1)
  cz = rnorm(D,2,1)
  
  u1mean = -1 
  u2mean = 2 
  u3mean = -2 
  
  U1 = rnorm(N,u1mean,varval) # X1,Z 
  U2 = rnorm(N,u2mean,varval) # X1,Y 
  U3 = rnorm(N,u3mean,varval) # Z,Y 
  
  U1.intv = rnorm(Nintv,u1mean,varval)
  U2.intv = rnorm(Nintv,u2mean,varval)
  U3.intv = rnorm(Nintv,u3mean,varval)
  
  fX1 = function(N,U1,U2){
    Ux = rnorm(N,0,1)
    # X = rbinom(N,size=1,inv.logit(log(abs(1*(U1*Z))+1)* 2*(2*Z-1)*Ux - 4*Ux*exp(U1-2) -3))
    X = rbinom(N,size=1,inv.logit(1*U1 - 2*U2 + Ux   ))
    return(X)
  }
  
  fZ = function(N,U1,U3){
    Uz = rnorm(N,0,1)
    Z = matrix(0,ncol=D,nrow=N)
    for (idx in 1:D){
      Z[,idx] = rbinom(N,size=1,prob=inv.logit(c1[idx]*U1+c2[idx]*U3))
    }
    Z = data.frame(Z)
    colnames(Z) = paste('Z',1:D,sep="")
    return(Z)
  }
  
  fR = function(N,X1){
    Ur = rnorm(N,0,1)
    R = rbinom(N,size=1,inv.logit(-1*X1 - Ur*(2*X1-1)   ))
    return(R)
  }
  
  fX2 = function(N,X1,Z){
    Ux2 = rnorm(N,0,0.5)
    Zmat = as.matrix(2*Z-1)
    czmat = as.matrix(cz)
    X2val = inv.logit(Zmat %*% czmat)
    X2 = rbinom(N,size=1,inv.logit(-1*X2val + Ux2-1 ))
    return(X2)
  }
  
  fY = function(N,R,X2,U2,U3){
    Uy = rnorm(N,-2,1)
    Y = rbinom(N,size=1,inv.logit(0.5*R - 2*(2*X2-1) + 2*U2 - 0.5*U3 ))
    # ind.X = 2*X - 1 
    # Y =1*(2*U2 + ind.X- Uy )
    # Y = inv.logit(Y)
    
    
    # Y = rbinom(N,size=1,inv.logit(-1*ind.X* U2*log(abs(U2*ind.X)+1) + 0.1*U2-Uy))
    return(Y)
  }
  
  # OBS construction 
  X1 = fX1(N,U1,U2)
  Z = fZ(N,U1,U3)
  R = fR(N,X1)
  X2 = fX2(N,X1,Z)
  Y = fY(N,R,X2,U2,U3)
  
  OBS = data.frame(Z,R,X1,X2,Y)
  
  # OBS construction 
  X1.intv = c(rep(0,Nintv/4),rep(0,Nintv/4),rep(1,Nintv/4),rep(1,Nintv/4))
  Z.intv = fZ(Nintv,U1.intv,U3.intv)
  R.intv = fR(Nintv,X1.intv)
  X2.intv = c(rep(0,Nintv/4),rep(1,Nintv/4),rep(0,Nintv/4),rep(1,Nintv/4))
  Y.intv = fY(Nintv,R.intv,X2.intv,U2.intv,U3.intv)
  
  INTV = data.frame(Z.intv,R.intv,X1.intv,X2.intv,Y.intv)
  
  return(list(OBS,INTV))
}

numStrata = function(D,numCate){
  return(2 * (numCate^D) * 2)
}