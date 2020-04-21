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
Due to the nonlinear nature of the compartmental pharmacokinetic model, the solutions to the ODEs specified above are computed using the “ODE” function from the package “deSolve”, and the results are stored as a list of values for each time period.

#compute the solution to the ode for slow release in children
pars=c(ka=(-log(1-0.51)/28),k10=0.411,k12=0.076,k21=0.215,k31=0.0005,k13=0.065)
stateini = c(M=1,X1=0,X2=0,X3=0)
mat <- ode( func = vaccinechildrenmod, y = stateini, parms = pars, times = times)
bbslow = apply(mat[,3:5],1,sum) # total bb for children


#compute the solution to the ode for slow release in children (AlOH)
pars=c(ka=(-log(1-0.17)/28),k10=0.411,k12=0.076,k21=0.215,k31=0.0005,k13=0.065)
stateini = c(M=1,X1=0,X2=0,X3=0)
mat <- ode( func = vaccinechildrenmod, y = stateini, parms = pars, times = times)
bbslowAlOH = apply(mat[,3:5],1,sum) # total bb for children


#compute the solution to the ode for fast release in children
pars=c(k10=0.411,k12=0.076,k21=0.215,k31=0.0005,k13=0.065)
stateini = c(X1=1,X2=0,X3=0)
mat <- ode( func = childrenmod, y = stateini, parms = pars, times = times)
bbfast = apply(mat[,2:4],1,sum) # total bb for children


#compute the solution to the ode for background mass in children
backgroundstate = c(X1=0,X2=0,X3=0.4)
background <- ode( func = childrenmod, y = backgroundstate, parms = pars, times = times)
backgroundbb=apply(background[,2:4],1,sum)
The minimal risk values (body burdens) of aluminum for the 5th and 50th percentile of infants body weight was calculated based on the safety threshold set by the Agency for Toxic Substances and Disease Registry (ATSDR). Exposure to aluminum from breastmilk and infant formula are also estimated as detailed in Mitkus et al. 2011. Functions were written to compute the body burden of aluminum following instantaneous (100% absorption) as well as the slow absorption from intramuscular injection of vaccines.

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
Combining all the values and functions that are set above, we finally generate the the values reflecting the body burden (for the 5th and 50th body weight percentile of infants).

bbdose = bbslowfunc(doses,injectiontimes)
bbdoseAlOH = ALOHfunc(doses,injectiontimes)
bolusbb = bbfunc(doses,injectiontimes)
safebb5 = bbfunc(safedose5,dailytimes)
safebb50 = bbfunc(safedose50,dailytimes)
bmbb= bbfunc(breastmilkdose,dailytimes)
fmbb= bbfunc(formulamilkdose,dailytimes)


# adjust everything for background
netbodyburden = backgroundbb+bbdose
netbodyburdenAlOH = backgroundbb+bbdoseAlOH
bolusbb = bolusbb + backgroundbb
safebb5 = safebb5 + backgroundbb
safebb50 = safebb50 + backgroundbb
bmbb= bmbb + backgroundbb
fmbb= fmbb + backgroundbb
Producing the figures.

# produce the plots

plot(times,netbodyburden,xlab="Days of age",ylim=c(0,7),ylab="Body burden (mg)",type="l")
lines(times,safebb50,col="blue")
lines(times,bmbb,col="orange")
lines(times,fmbb,col="green4")
legend(0,7,c("MRL 50","MRL 5","Vaccines (AlP04)","Formula","Breastmilk"), col = c("blue","red","black","green4","orange"),
       text.col = "black", lwd=2,merge = TRUE)
polygon(c(times,rev(times)),c(safebb5,rev(netbodyburden)),col=cm.colors(16)[9])
lines(times,safebb5,col="red")


plot(times,netbodyburdenAlOH,xlab="Days of age",ylim=c(0,7),ylab="Body burden (mg)",type="l")
lines(times,safebb50,col="blue")
lines(times,bmbb,col="orange")
lines(times,fmbb,col="green4")
legend(0,7,c("MRL 50","MRL 5","Vaccines (AlOH)","Formula","Breastmilk"), col = c("blue","red","black","green4","orange"),
       text.col = "black", lwd=2,merge = TRUE)
polygon(c(times,rev(times)),c(safebb5,rev(netbodyburdenAlOH)),col=cm.colors(16)[9])
lines(times,safebb5,col="red")


plot(times,bolusbb,xlab="Days of age",ylim=c(0,7),ylab="Body burden (mg)",type="l")
lines(times,safebb50,col="blue")
lines(times,bmbb,col="orange")
lines(times,fmbb,col="green4")
legend(0,7,c("MRL 50","MRL 5","Vaccines (Bolus)","Formula","Breastmilk"), col = c("blue","red","black","green4","orange"),
       text.col = "black", lwd=2,merge = TRUE)
# polygon(c(times,rev(times)),c(safebb5,rev(bolusbb)),col=cm.colors(16)[9])
lines(times,safebb5,col="red")


dev.off()
## null device 
##           1
Creating an “output” dataset for export and further analysis if needed. This datset is essentially the data that generated the figures above.

output = cbind(times,backgroundbb,bbdose,netbodyburden,netbodyburdenAlOH,safebb5,safebb50,bmbb,fmbb)
colnames(output)<-c("time","Background BB","BB due to Doses","Net BB","Net BB AlOH","Safe BB for 5%","Safe BB for 50%","Brest milk BB","Formula milk BB")