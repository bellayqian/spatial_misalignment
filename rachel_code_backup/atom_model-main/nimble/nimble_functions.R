
createProbVec<-nimbleFunction(
  run=function(vec=double(1),a=integer(0),M=integer(0),mat=double(2),K=integer(0)){
    ans<-numeric(K)
    for (m in 1:M){
      ans[mat[m,1]:mat[m,2]]<-vec[a+(mat[m,1]:mat[m,2])]/sum(vec[a+(mat[m,1]:mat[m,2])])
    }
    return(ans[])
    returnType(double(1))
  }
)

CcreateProbVec<-compileNimble(createProbVec)