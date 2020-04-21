# This section of code computes the daily doses of aluminum for safety and other purposes
safedose5 = 0.0078*(2.65899- (1.86774/sqrt(1 + dailytimes/30)) + 1.59926*sqrt(1 + dailytimes/30))
safedose50 = 0.0078*(3.35319+ 1.74026*sqrt(1 + dailytimes/30) + 0.618471*log(0.1+ 0.1*dailytimes/30))
breastmilkdose = 0.0078*c(rep(0.03,181),rep(0.7,220))
formulamilkdose = 0.0078*c(rep(0.15,181),rep(0.7,220))


# Here is a function which will compute bb under immediate absorption from daily doses
bbfunc<-function(doses,injectiontimes){
  q=length(injectiontimes)
  mm = length(times)
  tempmat = matrix(0,nrow=mm,ncol=q)
  offsetmat = matrix(0,nrow=mm,ncol=q)
  for(j in 1:q){
    tempmat[,j]=doses[j]*bbfast
    index = which(times==injectiontimes[j])
    offsetmat[index:mm,j] = tempmat[1:(mm-index+1),j]
  }
  out=apply(offsetmat,1,sum)
  out
}


# Here is a function which will compute bb under slow absorption from daily doses
bbslowfunc<-function(doses,injectiontimes){
  q=length(injectiontimes)
  mm = length(times)
  tempmat = matrix(0,nrow=mm,ncol=q)
  offsetmat = matrix(0,nrow=mm,ncol=q)
  for(j in 1:q){
    tempmat[,j]=doses[j]*bbslow
    index = which(times==injectiontimes[j])
    offsetmat[index:mm,j] = tempmat[1:(mm-index+1),j]
  }
  out=apply(offsetmat,1,sum)
  out
}


# Here is a function which will compute bb for ALOH under slow absorption
ALOHfunc<-function(doses,injectiontimes){
  q=length(injectiontimes)
  mm = length(times)
  tempmat = matrix(0,nrow=mm,ncol=q)
  offsetmat = matrix(0,nrow=mm,ncol=q)
  for(j in 1:q){
    tempmat[,j]=doses[j]*bbslowAlOH
    index = which(times==injectiontimes[j])
    offsetmat[index:mm,j] = tempmat[1:(mm-index+1),j]
  }
  out=apply(offsetmat,1,sum)
  out
}