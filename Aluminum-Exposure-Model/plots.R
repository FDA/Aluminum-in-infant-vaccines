# produce the plots

#body burden v age in days [Vaccines (AlP04)]
png("Body_burden_v_age_in_days_[Vaccines_(AlP04)].png")
plot(times,netbodyburden,xlab="Days of age",ylim=c(0,7),ylab="Body burden (mg)",type="l")
lines(times,safebb50,col="blue")
lines(times,bmbb,col="orange")
lines(times,fmbb,col="green4")
legend(0,7,c("MRL 50","MRL 5","Vaccines (AlP04)","Formula","Breastmilk"), col = c("blue","red","black","green4","orange"),
       text.col = "black", lwd=2,merge = TRUE)
polygon(c(times,rev(times)),c(safebb5,rev(netbodyburden)),col=cm.colors(16)[9])
lines(times,safebb5,col="red")
dev.off()

#body burden v age in days [Vaccines (AlOH)]
png("Body_burden_v_age_in_days_[Vaccines_(AlOH)].png")
plot(times,netbodyburdenAlOH,xlab="Days of age",ylim=c(0,7),ylab="Body burden (mg)",type="l")
lines(times,safebb50,col="blue")
lines(times,bmbb,col="orange")
lines(times,fmbb,col="green4")
legend(0,7,c("MRL 50","MRL 5","Vaccines (AlOH)","Formula","Breastmilk"), col = c("blue","red","black","green4","orange"),
       text.col = "black", lwd=2,merge = TRUE)
polygon(c(times,rev(times)),c(safebb5,rev(netbodyburdenAlOH)),col=cm.colors(16)[9])
lines(times,safebb5,col="red")
dev.off()

#Body burden v age in days [Vaccines (Bolus)]
png("Body_burden_v_age_in_days_[Vaccines_(Bolus)].png")
plot(times,bolusbb,xlab="Days of age",ylim=c(0,7),ylab="Body burden (mg)",type="l")
lines(times,safebb50,col="blue")
lines(times,bmbb,col="orange")
lines(times,fmbb,col="green4")
legend(0,7,c("MRL 50","MRL 5","Vaccines (Bolus)","Formula","Breastmilk"), col = c("blue","red","black","green4","orange"),
       text.col = "black", lwd=2,merge = TRUE)
# polygon(c(times,rev(times)),c(safebb5,rev(bolusbb)),col=cm.colors(16)[9])
lines(times,safebb5,col="red")
dev.off()

