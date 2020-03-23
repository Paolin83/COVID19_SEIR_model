
SEIR_betavar<-function(theta, B, status,fw){
#pars=parameters
  #beta=expit(est$par[3])-0.5
  #status=initials
  #forecast=14
seir1<-NULL
if(length(B)==1) B=rep(B,fw)
for(i in 1:fw){
parameters <- c(mu=theta[1], beta=B[i], sigma=theta[2], gamma=theta[3])
if( i==1) initials <- c(S=status[1], E=status[2], I=status[3], R=status[4])
if( i>1) initials <- c(S = seir1_temp$results$S[2], E = seir1_temp$results$E[2], I =seir1_temp$results$I[2], R = seir1_temp$results$R[2])
seir1_temp <- SEIR(pars = parameters, init = initials, time = 0:1)
seir1 <- rbind(seir1,SEIR(pars = parameters, init = initials, time = 0:1)$results[2,])
}
seir1$time<-1:fw
seir1
}
