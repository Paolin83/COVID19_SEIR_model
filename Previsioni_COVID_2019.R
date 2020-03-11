rm(list=ls())
###import italian dataset updated 10 March 2020
dat_csv<-read.csv("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv",header=T)
dat_csv$t<-1:dim(dat_csv)[1]
##### SEIR model

#stimo r0
#vedi https://kingaa.github.io/clim-dis/parest/parest.html

fit1 <- lm(log(totale_casi)~t,data=dat_csv)

# stimo 

plot(dat_csv$t,log(dat_csv$totale_casi))
abline(coef(summary(fit1))[,1])
summary(fit1)
# R0
slope <-coef(summary(fit1))[2,1]; slope
slope.se <- coef(summary(fit1))[2,2]; slope.se
###R0 stima + o meno co durata contagio 14 giorni
R_0=slope*14+1;R_0
(slope+c(-1,1)*1.96*slope.se)*14+1

I0<-max(dat_csv$totale_casi)
R0<-max(dat_csv$dimessi_guariti)
N<-60480000
library(EpiDynamics)
#beta
beta0<-R_0/(14)
# durata contagio 14 giorni
duratac<-14
#ipotiziamo un sigma (coronavirus transmission rate) del 5%  (la metà dell'influenza)
mu0<-1/(82*365.25) # 1/durata della vita
parameters <- c(mu = mu0, beta = 1/14, sigma = 0.05, gamma = 1/duratac)
# numero medio di persone esposte per ogni infetto
f1<-10
initials <- c(S = 0.95, E = (f1*I0/N), I = I0/N, R = R0/N)
seir1 <- SEIR(pars = parameters, init = initials, time = 0:14)

f2<-5
initials <- c(S = 0.95, E = (f2*I0/N), I = I0/N, R = R0/N)
seir2 <- SEIR(pars = parameters, init = initials, time = 0:14)

f3<-3
initials <- c(S = 0.95, E = (f3*I0/N), I = I0/N, R = R0/N)
seir3 <- SEIR(pars = parameters, init = initials, time = 0:14)



plot(c(dat_csv$totale_casi,seir1$results$I[-1]*60480000),type="l",ylab="Number of infectus",xlab="time")
lines(c(dat_csv$totale_casi,seir2$results$I[-1]*60480000),col=2)
lines(c(dat_csv$totale_casi,seir3$results$I[-1]*60480000),col=3)
legend("topleft",c("first scenario","second scenario","third scenario"),lty=1,col=1:3)

