#Load libraries
library(spdep)
library(maptools)
library(INLA)
library(rgdal)
library(foreign)
library("dygraphs")
library("xts")
#define work directory
################################
rm(list=ls())
wd<-"~/Dropbox/COVID19"
setwd(wd)
################################ IMPORT DATASETS
# ISTAT Province dataset
db_istat<-read.csv("Elenco-comuni-italiani2.csv",sep=";",header=T,dec=",")
pop_provincia<-data.frame(codice_provincia=as.numeric(names(table(db_istat[,3]))),pop=tapply(db_istat[,20],db_istat[,3],sum))
###import international dataset
dat_int<-read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv",header=T)
hubei<-as.numeric(dat_int[dat_int$Province.State=="Hubei",5:dim(dat_int)[2]])

###import updated italian dataset  for each province
dat_csv<-read.csv("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-province/dpc-covid19-ita-province.csv",header=T)
dat_csv<-merge(dat_csv,pop_provincia,by=c("codice_provincia"),all.y=TRUE)
dat_csv<-dat_csv[dat_csv$codice_provincia<112,]
#select a Region, in this example "Lombardia" Code region 3
#Region<-3
#dat_csv<-dat_csv[dat_csv$codice_regione==Region,]
dat_csv_n$denominazione_provincia<-droplevels(dat_csv$denominazione_provincia)
#### number of province
nprov<-length(table(dat_csv$denominazione_provincia)[table(dat_csv$denominazione_provincia)>0])

#order the obsercation for date and space, useful for INLA later....
dat_csv<-dat_csv[order(dat_csv$data),]
dat_csv<-dat_csv[order(dat_csv$codice_provincia),]
###space ID
dat_csv$ID<-as.factor(dat_csv$codice_provincia)
levels(dat_csv$ID)<-1:length(levels(dat_csv$ID))
dat_csv$ID<-as.numeric(dat_csv$ID)
#calculate new_cases
dat_csv$cumulative_cases<-dat_csv$totale_casi
dat_csv$totale_casi[dat_csv$data!="2020-02-24 18:00:00"]<-unlist(tapply(dat_csv$totale_casi,dat_csv$ID,diff))
dat_csv$totale_casi[dat_csv$totale_casi<0]<-0

#import dataset shapefile Italian province
nc.province <- readShapePoly(paste(wd,"/ProvCM01012019_g/ProvCM01012019_g_WGS84.shp",sep=""))
nc.province<-nc.province[nc.province$COD_REG==Region,]
nc.province<-nc.province[order(nc.province$COD_PROV),]
nc.province$ID<-as.factor(nc.province$COD_PROV)
levels(nc.province$ID)<-1:length(levels(nc.province$ID))
nc.province$ID<-as.numeric(nc.province$ID)
#Create adjacency matrix
nc.nb <- poly2nb(nc.province,snap=0.01)
nb2INLA(paste(wd,"nc.adj",sep=""), nc.nb)
file.adj <- paste(wd,"nc.adj",sep="")
################################

################################
#prio distribution for iid and besag components
prior.iid = c(1,0.01)
prior.besag = c(1,0.001)
initial.iid = 4
initial.besag = 3
####COVID19 MAPS BYM with INLA
##########################indexing, for 
#### number of days
dat_csv$t<-as.numeric(dat_csv$data)
days<-max(dat_csv$t)
dat_csv$t3<-dat_csv$t2<-dat_csv$t
#### province
dat_csv$ID2<-dat_csv$ID
##########################
formula.bym= totale_casi ~ 1+f(ID, model = "bym", graph = file.adj, param = c(prior.iid, prior.besag), initial = c(initial.iid, initial.besag))
fit_1<-inla(formula.bym, family="poisson", data=dat_csv, E=pop,control.compute = list(dic=T))
summary(fit_1)
##incidence rate ratio for each province
IRR_mean <- exp(fit_1$summary.random$ID$mean[1:nprov])
nc.province$IRR_mean<-IRR_mean
spplot(nc.province, c( "IRR_mean"))

########Pure temporal model (RW2)
#RW2 temporal model
Date<-seq(as.Date("2020-02-24"),as.Date("2020-02-24")+days-1,1)
formula_t = totale_casi ~ 1+f(t,model="rw2", constr = FALSE)+f(t2,model="iid")
fit_2<-inla(formula_t, family="poisson", data=dat_csv, E=pop,control.compute = list(dic=T))
summary(fit_2)
plot(Date,fit_2$summary.random$t$mean,ylab="Temp Coef",xlab="Date")
#estimating a linear model without the first obs
lm_fit<-lm(fit_2$summary.random$t$mean~I(1:days),subset=-1)
coef(lm_fit)[2]
plot(1:days,fit_2$summary.random$t$mean,ylab="Temp Coef",xlab="Date")
abline(lm_fit)
##Incidence rate ratios, the rate of increase od cases respect to time 1
plot(Date,exp(fit_2$summary.random$t$mean),ylab="Indicence of cases",xlab="Date")

####interaction model type 4 (complete time trend interaction to estimate specific curves for each province)
####specification
formula.intIV<- totale_casi ~f(ID,model="bym",graph=file.adj) +
  f(t,model="rw2",constr = FALSE) +
  f(t2,model="iid") +
  f(ID2,model="besag", graph=file.adj,
    group=t3,
    control.group=list(model="rw2"))
#estimation
fit_st4<-inla(formula.intIV, family="poisson", data=dat_csv, E=pop,control.compute = list(dic=T))
summary(fit_st4)

#### components
t_1<-fit_st4$summary.random$t$mean #overall trend
t_2<-fit_st4$summary.random$t2$mean #iid trend
s_1<-fit_st4$summary.random$ID$mean[1:nprov] # spatial CAR component
s_2<-fit_st4$summary.random$ID$mean[(nprov+1):(nprov*2)] # spatial idd
st<-fit_st4$summary.random$ID2$mean # spatio-temporal model with car and ar1 specification 
## Mean IRR for each province
nc.province$IRR_st<-exp(s_1+apply(t(matrix(st,nrow=nprov)),2,mean))
spplot(nc.province, c( "IRR_st"),main="IRR")
## Trends for province

trends<-t_1+t(matrix(st,nrow=nprov))
# trends+ specific
matplot(trends)
#specific for each province
matplot(t(matrix(st,nrow=nprov)))


##################predictions 3 days forward
Forecast=3
dat_csv2<-dat_csv 
dat_csv2$t<-rep(1:days+days,nprov)
dat_csv2$totale_casi<-NA
dat_csv2<-dat_csv2[dat_csv2$t<=(days+Forecast),]
dat_csv_n<-rbind(dat_csv,dat_csv2)
dat_csv_n$t3<-dat_csv_n$t2<-dat_csv_n$t
#introducting hubei starting on 22/01, quarantine started on 25/01
# applying a smoothing , data has errors
hubei_new_cases<-filter(diff(hubei), rep(1, 5))

plot(hubei/58.5,type="l",ylab="COVID19 cases for milion of inhabitant",ylim=c(0,1530),xlab="Days since over 100 cases")
lines(c(rep(NA,2),cumsum(tapply(dat_csv_n$totale_casi,dat_csv_n$t,sum))/10),col=2)
abline(v=2,col=1,lty=2)
abline(v=4,col=2,lty=2)
abline(v=10,col=3,lty=3)
legend("bottomright",c("Hubei","Lombardy"),lty=1,col=1:2)
text(10,1500,"Hubei quarantine")
text(12,1350,"Lodi home lock",col=2)
text(20,1200,"Lombardy home lock",col=3)
#province of lodi restrictive measure started on 24/02 lag 28 days, at the beginning of the time series
#other provinces on 11/03 lagged started their restrictions
hubei_lodi<-hubei_new_cases[1:(days+Forecast)]
hubei_others<-c(rep(NA,8),hubei_new_cases[1:(days+Forecast-8)])
dat_csv_n<-dat_csv_n[order(dat_csv_n$codice_provincia,dat_csv_n$t),]
dat_csv_n$hubei<-rep(hubei_others,12)
dat_csv_n$hubei[dat_csv_n$denominazione_provincia=="Lodi"]<-hubei_lodi
#overdispersion
dat_csv_n$over<-1:dim(dat_csv_n)[1]

## we model the number of cases
## in the rw2 I set constr = FALSE is set to FALSE and that, for this reason, the intercept is not included in the linear predictor.
## include overdispersion parameter
formula.intIVn<- totale_casi ~ log(hubei+1)+log(pop)+
  f(t,model="rw2", constr = FALSE) +
  f(t2,model="iid") +f(over,model="iid") +
f(ID,model="bym",graph=file.adj) +
  f(ID2,model="besag", graph=file.adj,
    group=t3,
    control.group=list(model="rw2"))+ f(over,model="iid") 

fit_st4n<-inla(formula.intIVn, family="poisson", data=dat_csv_n,control.compute = list(dic=T),control.predictor = list(link = 1))
summary(fit_st4n)


#### components
t_1<-fit_st4n$summary.random$t$mean #overall trend
t_2<-fit_st4n$summary.random$t2$mean #iid trend
s_1<-fit_st4n$summary.random$ID$mean[1:nprov] # spatial CAR component
s_2<-fit_st4n$summary.random$ID$mean[(nprov+1):(nprov*2)] # spatial idd
st<-fit_st4n$summary.random$ID2$mean # spatio-temporal coefficients 
## Trends for province
Date_n<-seq(as.Date("2020-02-24"),as.Date("2020-02-24")+days+Forecast-1,1)
trends<-t_1+t(matrix(st,nrow=nprov))
plot(Date_n,t_1,ylab="RW2 Time coefficient",xlab="Date")
matplot(t(matrix(st,nrow=nprov)),ylab="RW2 Time coefficient for Province",xlab="Date")
#extract number of cases for lodi
est_lodi<-fit_st4n$summary.fitted.values[dat_csv_n$denominazione_provincia=="Lodi",]


days.before<-Date_n[1:days]
days.ahead<-Date_n[(days+1):(days+Forecast)]
mu.lower<-est_lodi$mean-1.96*est_lodi$sd
mu.lower[mu.lower<0]<-0
mu.upper<-est_lodi$mean+1.96*est_lodi$sd
mu.med<-xts(est_lodi$mean,order.by = c(days.before,days.ahead),frequency = 7)
counts<-mu.med
step.ahead<-Forecast
mu<-xts(x = as.matrix(cbind(counts,mu.lower,mu.upper)) , order.by = c(days.before,days.ahead))
p <- dygraph(mu,main=paste("Lodi (Credible Interval ",100*0.95,"%)",sep = ""),ylab=" Infected",xlab="Day",height=400,width=800) %>%  dySeries(c("mu.lower", "counts", "mu.upper"),label="counts")
p<-p %>% dyLegend(show = "always", hideOnMouseOut = FALSE) %>%  dyShading(from = days.ahead[1], to = days.ahead[step.ahead], color = "#CCEBD6")%>% dyEvent(days.ahead[1], "Prediction", labelLoc = "bottom")
p


