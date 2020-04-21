#Creating an “output” dataset for export and further analysis if needed. This datset is essentially the data that generated the figures.

output = cbind(times,backgroundbb,bbdose,netbodyburden,netbodyburdenAlOH,safebb5,safebb50,bmbb,fmbb)
colnames(output)<-c("time","Background BB","BB due to Doses","Net BB","Net BB AlOH","Safe BB for 5%","Safe BB for 50%","Brest milk BB","Formula milk BB")
write.csv(output, 'outputs.csv')