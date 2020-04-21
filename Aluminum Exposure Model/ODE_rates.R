library("deSolve")

dosedat <- cbind(c(0, 60, 120, 180, 365), c(0.250, 1.20, 1.20, 0.975, 0.600))
injectiontimes = dosedat[,1]
doses=dosedat[,2]
m = length(doses)
endtime = 400
dt =0.1
times <- seq(0,endtime,by = dt)
mm = length(times)
dailytimes = seq(0,endtime,1)


#ODE for 100% instantaneous absorption in children
childrenmod <- function(time, state, pars) {
  with(as.list(c(state, pars)), {
    dX1 <- -k10*( 0.361 + 0.639*(time/(time+231.462)))*X1 + k21*X2 + k31*X3 -k12*X1 -k13*X1
    dX2 <- k12*X1 - k21*X2
    dX3 <- k13*X1 - k31*X3
    return(list(c(dX1,dX2,dX3)))
  })
}


# ODE for 100% instantaneous absorption in adults
adultmod <- function(time, state, pars) {
  with(as.list(c(state, pars)), {
    dX1 <- -k10*X1 + k21*X2 + k31*X3 -k12*X1 -k13*X1
    dX2 <- k12*X1 - k21*X2
    dX3 <- k13*X1 - k31*X3
    return(list(c(dX1,dX2,dX3)))
  })
}


#ODE for slow release in children
vaccinechildrenmod <- function(time, state, pars) {
  with(as.list(c(state, pars)), {
    dM <- -ka* M 
    dX1 <- ka*M -k10*( 0.361 + 0.639*(time/(time+231.462)))*X1 + k21*X2 + k31*X3 -k12*X1 -k13*X1
    dX2 <- k12*X1 - k21*X2
    dX3 <- k13*X1 - k31*X3
    return(list(c(dM,dX1,dX2,dX3)))
  })
}


# ODE for slow release in adults
vaccineadultmod <- function(time, state, pars) {
  with(as.list(c(state, pars)), {
    dM <- -ka*M
    dX1 <- ka*M -k10*X1 + k21*X2 + k31*X3 -k12*X1 -k13*X1
    dX2 <- k12*X1 - k21*X2
    dX3 <- k13*X1 - k31*X3
    return(list(c(dM,dX1,dX2,dX3)))
  })
}