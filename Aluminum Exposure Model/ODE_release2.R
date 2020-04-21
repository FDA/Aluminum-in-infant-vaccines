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