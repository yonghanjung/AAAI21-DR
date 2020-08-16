
# a = P(y|x1,x2,z,r)
# b = P(z,x1) or P(x1)
# c = P(r|x1)
# d = P(x2|x1,z)

myFun = function(a,b,c,d){
  return(((a & b) | (c & d)) & (c | b))
}

myresult = c()
for (a in c(T,F)){
  for (b in c(T,F)){
    for (c in c(T,F)){
      for (d in c(T,F)){
        myresult = rbind(myresult, c(a,b,c,d, myFun(a,b,c,d)))
      }
    }
  }
}
myresult = data.frame(myresult)
colnames(myresult) = c("a","b","c","d","result")
rownames(myresult) = c(1:nrow(myresult))
print(myresult[myresult$c==F & myresult$d==F,])
# print(myresult[myresult$a==F & myresult$d==F,])

colnames(myresult) = c("P(y|x1,x2,z,r)","P(z,x1)","P(r|x1)","P(x2|x1,z)","result")