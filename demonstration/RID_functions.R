fix_pred = function(pred){
  smallval = 1e-1
  pred[pred<0] = smallval; pred[pred>1] = 1-smallval
  return(pred)
}

mis_pred = function(predval){
  distortval = 0.5
  sign_rv = 2*rbinom(n=length(predval),size=1,prob=0.5)-1
  predval = predval + sign_rv*rnorm(length(predval),distortval,0.1)
  predval = fix_pred(predval)
  return(predval)
}

# distortVar = function(myvar,seednum){
#   set.seed(seednum)
#   while(1){
#     numSample =  length(unique(myvar))
#     numSize = length(myvar)
#     myprob = c(runif(1,min=0,max=1))
#     sumprob = myprob[1]
#     for (idx in 1:(numSample-1)){
#       lastprob = runif(1,min=0,max=(1-sumprob))
#       myprob = c(myprob,lastprob)
#       sumprob = sumprob + lastprob
#     }
#     myprob[length(myprob)] = 1-sum(myprob[1:(length(myprob)-1)])
#     mySample = sample(x=c(0:(numSample-1)),size=numSize,replace=T,prob=myprob)
#     if (length(unique(mySample)) == length(unique(myvar)) ){
#       return(mySample)
#     }
#   }
# }


distortVar = function(myvar,seednum){
  set.seed(seednum)
  while(1){
    mysample = sample(c(0:(length(unique(myvar))-1)),size=length(myvar),replace=T)
    if (length(unique(mysample)) == length(unique(myvar)) ){
      return(mysample)
    }
  }
}

