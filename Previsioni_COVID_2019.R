rm(list=ls())
#base R code for
###import italian dataset updated 11 March 2020
dat_csv<-read.csv("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv",header=T)
days<-dim(dat_csv)[1]
dat_csv$t<-1:days

##### SEIR model

#estimate r0
#see https://kingaa.github.io/clim-dis/parest/parest.html
#calculate r0 based in the last XX observation 
beta_vec<-NULL
sd_vec<-NULL
for (i in 1:(days-6)){
fit <- lm(log(totale_attualmente_positivi)~t,data=dat_csv[i:days,])
beta_vec<-c(beta_vec,coef(fit)[2])
sd_vec<-c(sd_vec,coef(summary(fit))[2,2])
}

label <-days:7
mean  <- (beta_vec*14+1)
lower <- ((beta_vec-1.96*sd_vec)*14+1)
upper <- ((beta_vec+1.96*sd_vec)*14+1)

df <- data.frame(label, mean, lower, upper)

# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=rev(df$label))

library(ggplot2)
fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Last days used for the calculation") + ylab("R0 Mean (95% CI)") +
  theme_bw()  # use a white background
print(fp)




fit1 <- lm(log(totale_attualmente_positivi)~t,data=dat_csv)

# stimo 

plot(dat_csv$t,log(dat_csv$totale_attualmente_positivi))
abline(coef(summary(fit1))[,1])
summary(fit1)
# R0
slope <-coef(summary(fit1))[2,1]; slope
slope.se <- coef(summary(fit1))[2,2]; slope.se
###R0 incubation time of 14 days
R_0=slope*14+1;R_0
(slope+c(-1,1)*1.96*slope.se)*14+1

I0<-dat_csv$totale_attualmente_positivi[dim(dat_csv)[1]]
R0<-dat_csv$dimessi_guariti[dim(dat_csv)[1]] 
N<-60480000
library(EpiDynamics)
#beta
beta0<-R_0/(14)
# durata contagio 14 giorni
duratac<-14
#ipotiziamo un sigma (coronavirus transmission rate) del 5%  (la metà dell'influenza)
mu0<-1/(82*365.25) # 1/durata della vita
parameters <- c(mu = mu0, beta = beta0, sigma = 0.05, gamma = 1/duratac)
# numero medio di persone esposte per ogni infetto
f1<-10
initials <- c(S = 0.95, E = (f1*I0/N), I = I0/N, R = R0/N)
seir1 <- SEIR(pars = parameters, init = initials, time = 0:14)

parameters <- c(mu = mu0, beta = beta0, sigma = 0.05, gamma = 1/duratac)
f2<-5
initials <- c(S = 0.95, E = (f2*I0/N), I = I0/N, R = R0/N)
seir2 <- SEIR(pars = parameters, init = initials, time = 0:14)

parameters <- c(mu = mu0, beta = beta0*1/2, sigma = 0.05, gamma = 1/duratac)
f3<-5
initials <- c(S = 0.95, E = (f3*I0/N), I = I0/N, R = R0/N)
seir3 <- SEIR(pars = parameters, init = initials, time = 0:14)



plot(c(dat_csv$totale_attualmente_positivi,seir1$results$I[-1]*N),type="l",ylab="Number of infectus",xlab="time",main="Infected")
lines(c(dat_csv$totale_attualmente_positivi,seir2$results$I[-1]*N),col=2)
lines(c(dat_csv$totale_attualmente_positivi,seir3$results$I[-1]*N),col=3)
legend("topleft",c("first scenario","second scenario","third scenario"),lty=1,col=1:3)

date<-seq(as.Date("2020-02-24"),as.Date("2020-02-24")+14-1+dim(dat_csv)[1],1)
plot(date,c(dat_csv$totale_casi,(seir1$results$I[-1]+seir1$results$R[-1])*N),type="l",ylab="Cases",xlab="time",main="Cumulative Infected")
lines(date,c(dat_csv$totale_casi,(seir2$results$I[-1]+seir2$results$R[-1])*N),col=2)
lines(date,c(dat_csv$totale_casi,(seir3$results$I[-1]+seir3$results$R[-1])*N),col=3)
lines(date[1:dim(dat_csv)[1]],dat_csv$totale_casi,lwd=2)
legend("topleft",c("first scenario","second scenario","third scenario"),lty=1,col=1:3)
PlotMods(seir2)
