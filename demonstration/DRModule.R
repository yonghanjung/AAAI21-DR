returnUnique = function(OBS){
  mylist = c()
  for (colidx in 1:ncol(OBS)){
    myVar = OBS[,colidx]
    myVarUnique = unique(myVar)[order(unique(myVar))]
    mylist = append(mylist,list(myVarUnique))
  }
  return(mylist)
}

recursiveFilter = function(OBS,allpossible,rowidx){
  myFilter = OBS 
  for (colidx in c(1:ncol(OBS))){
    myFilter = myFilter[myFilter[,colidx] == allpossible[rowidx,colidx],]
  }
  return(myFilter)
}

nonRandomSampleSplit = function(OBS,mysize){
  # tic(msg="nonRandomSample"); toc()
  # print("nonRandomSample")
  OBSuniqueList = returnUnique(OBS)
  # Enumerate all possible values of column
  tmp = c()
  for (idx in 1:length(OBSuniqueList)){
    tmp = append(tmp,list(OBSuniqueList[[idx]]))
  }
  allpossible = expand.grid(tmp)
  colnames(allpossible) = colnames(OBS)
  mycollect = c()
  for (rowidx in 1:nrow(allpossible)){
    # filtered_OBS = subset(OBS, W == allpossible[rowidx,'W'] & R == allpossible[rowidx,'R'] & X == allpossible[rowidx,'X'] & Y == allpossible[rowidx,'Y'] )
    filtered_OBS = recursiveFilter(OBS,allpossible,rowidx)
    if (nrow(filtered_OBS) > 0){
      mycollect = rbind(mycollect,filtered_OBS[1,]) 
    }
  }
  if (nrow(mycollect) < mysize){
    for (rowidx in 1:nrow(OBS)){
      if (rowidx %in% as.numeric(rownames(mycollect)) == FALSE){
        mycollect = rbind(mycollect,OBS[rowidx,])  
      }
      if (nrow(mycollect) == mysize){
        break 
      }
    }
  }
  OBS_1 = mycollect
  totalidx = c(1:nrow(OBS))
  mycollect2 = setdiff(totalidx,as.numeric(rownames(mycollect)))
  
  OBS_2 = OBS[mycollect2,]
  
  # OBS_2_unique = returnUnique(OBS_2)
  # OBS_1_unique = returnUnique(OBS_1)
  # if (identical( OBS_2_unique[[length(OBS_2_unique)]], OBS_1_unique[[length(OBS_1_unique)]] ) == FALSE){
  #   while(1){
  #     splitidx = sample(c(1:nrow(OBS)),size=nrow(OBS)/2)
  #     OBS_2 = OBS[splitidx,]
  #     if (identical(OBS_2_unique[[length(OBS_2_unique)]], OBS_1_unique[[length(OBS_1_unique)]] )){
  #       break
  #     }
  #   }  
  # }
  # OBS_1 = mycollect
  rownames(OBS_1) = c(1:nrow(OBS_1))
  rownames(OBS_2) = c(1:nrow(OBS_2))
  return(list(OBS_1,OBS_2))
}



GoodSplit = function(OBS){
  totalidx = c(1:nrow(OBS))
  
  iteridx = 0
  iterMax = 10
  
  while(1){
    # print(iteridx)
    iteridx = iteridx + 1 
    
    splitidx_1 = sample(c(1:nrow(OBS)),size=nrow(OBS)/2)
    splitidx_2 = setdiff(totalidx,splitidx_1)
    OBS_1 = OBS[splitidx_1,]
    OBS_2 = OBS[splitidx_2,]
    
    Unique_1 = returnUnique(OBS_1)
    Unique_2 = returnUnique(OBS_2)
    stopSwitch = TRUE 
    if (identical(Unique_1,Unique_2) == FALSE){
      stopSwitch = FALSE 
    }
    if (stopSwitch == TRUE){
      break
    }
    if (iteridx > iterMax){
      mytmp = nonRandomSampleSplit(OBS,mysize=nrow(OBS)/2)
      OBS_1 = mytmp[[1]]
      OBS_2 = mytmp[[2]]
      break 
    }
  }
  rownames(OBS_1) = c(1:nrow(OBS_1))
  rownames(OBS_2) = c(1:nrow(OBS_2))
  if (nrow(OBS_1) < nrow(OBS_2)){
    OBS_2 = OBS_2[c(1:nrow(OBS_1)),]
  }
  if (nrow(OBS_2) < nrow(OBS_1)){
    OBS_1 = OBS_2[c(1:nrow(OBS_2)),]
  }
  return(list(OBS_1,OBS_2))
}

RunTryCatchProb_plugin = function(FUN, allpossible, DATA_Train, Data_Eval, mylambda){
  tryCatch(
    expr = {
      myresult = FUN(allpossible,Data_Eval,mylambda)  
    },
    error = function(e){
      print("Error in the Function")
    },
    finally = {
      myresult = FUN(allpossible,DATA_Train,mylambda)  
    }
  )
  return(myresult)
}

RunTryCatchProb_WERM = function(FUN, DATA_Train, Data_Eval, mylambda){
  tryCatch(
    expr = {
      myresult = FUN(Data_Eval,DATA_Train,mylambda)  
    },
    error = function(e){
      print("Error in the Function")
    },
    finally = {
      myresult = FUN(DATA_Train,Data_Eval,mylambda)  
    }
  )
  return(myresult)
}
