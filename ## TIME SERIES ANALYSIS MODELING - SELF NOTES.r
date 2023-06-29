## TIME SERIES ANALYSIS - SELF NOTES - PART - I 


library()   #list all available packages#

install.packages("TSA")     #Install packages- TSA#
library("TSA")     #Load packages - TSA#
library(help=TSA)

setwd("  ")

#1#
win.graph(width=4.875, height=2.5,pointsize=8)
data(larain)
larain
plot(larain,ylab='Inches',xlab='Year',type='o')



#2#
mel=ts(read.csv("melanoma.csv",header=F),start=1936)
x=plot(mel, ylab='per',xlab='Year',type='o')
savePlot(filename=”x”,type=c("jpg"))



#3#
temp=ts(read.csv("temperature.csv",header=F),start=1880,frequency=12)
plot(temp, ylab='temperature',xlab='month',type='o')


#4#
data(oilfilters)
plot(oilfilters,type='l',ylab='Sales')
points(y=oilfilters,x=time(oilfilters),pch=as.vector(season(oilfilters)))

#5#
temp=ts(read.table("17130.txt",header=F),start=1950,frequency=12)
plot(temp, ylab='temperature',xlab='month',type='l')
plot(ts(window(temp,c(1980,1),c(1981,12),12)), ylab='temperature',xlab='month',type='l')
tempmatrix=matrix(temp,57,12,byrow = TRUE)
ts.plot(ts(tempmatrix,start=1950), ylab='temperature',xlab='month',col=1:12)
ts.plot(ts(tempmatrix[,2],start=1950), ylab='temperature',xlab='month',col=1:1)
ts.plot(ts(tempmatrix[,2:5],start=1950), ylab='temperature',xlab='month',col=1:4)
par(mfrow=c(3,2))
plot(ts(tempmatrix[,1],start=1950), main="january", ylab="temperature")
plot(ts(tempmatrix[,2],start=1950), main="February", ylab="temperature")
plot(ts(tempmatrix[,3],start=1950), main="March", ylab="temperature")
plot(ts(tempmatrix[,4],start=1950), main="April", ylab="temperature")
plot(ts(tempmatrix[,5],start=1950), main="May", ylab="temperature")
plot(ts(tempmatrix[,6],start=1950), main="june", ylab="temperature")

#6#
data(AirPassengers)
AP=AirPassengers
layout(1:2)
plot(aggregate(AP))
boxplot(AP ~ cycle(AP))

#PART 2

#model identification example from acf pacf (correlogram)

#1 Simulating AR(1) process with phi = 0.5 
ar1sim=arima.sim(list(order = c(1,0,0), ar = 0.5), n = 800)
plot(ar1sim,ylab=expression(Y[t]),type='o');
acf(ar1sim);
pacf(ar1sim);

#2 Simulating non-stationary AR(1) or ARIMA(1,1,0) process with phi = 0.5 
ar1sim=arima.sim(list(order = c(1,1,0), ar = 0.5), n = 800)
plot(ar1sim,ylab=expression(Y[t]),type='o');
acf(ar1sim);
pacf(ar1sim);

# after differencing
plot(diff(ar1sim))
acf(diff(ar1sim))
pacf(diff(ar1sim))


#3 Simulating AR(2) process with phi(1) = 1.5 phi(2) = -0.75
ar2sim=arima.sim(n = 200, list(order = c(2,0,0), ar = c(1.5,-0.75)))
acf(ar2sim)
pacf(ar2sim)
plot(ar2sim, ylab=expression(Y[t]),type='o')

#4 Simulating MA(1) process with theta = 0.7 
ma1sim=arima.sim(list(order = c(0,0,1), ma = 0.7), n = 200)
plot(ma1sim,ylab=expression(Y[t]),type='o')
acf(ma1sim)
pacf(ma1sim)


#5 Simulating ARMA(1,1) process with theta(1) = 0.7 and phi=0.5 
armasim=arima.sim(n=200,list(order = c(1,0,1),ar=c(0.5),ma=c(0.7)))
acf(armasim);
pacf(armasim);



#6 Non-stationary model ARIMA
arma2sim=arima.sim(n=800,list(order = c(2,1,2),ar=c(1.1,-0.44),ma=c(1.7,-0.5)))
plot(arma2sim,ylab=expression(Y[t]),type='o');
acf(arma2sim);
pacf(arma2sim);

# after differencing

plot(diff(arma2sim))
acf(diff(arma2sim))
pacf(diff(arma2sim))


#7 acf and pacf of mel series

mel=ts(read.csv("melanoma.csv",header=F),start=1936)
plot(mel)
acf(mel,lag.max=36) 
pacf(mel,lag.max=36)
acf(diff(mel),lag.max=36) 
pacf(diff(mel),lag.max=36) 

#PART 3

# EACF example
library("TSA")
w=ts(read.table("w.txt",header=F), start=1818)
plot(w, ylab='per',xlab='Year',type='o')
acf(w)
pacf(w)
eacf(w)

# We need caschrono package.
# Install caschrono at first, 
library("caschrono")
armaselect(w, nbmod=5)


# Simulation of an ARMA(1,2) model
yc = arima.sim(n = 300, list( ar = -0.8, ma= c(-0.3, 0.6)), sd = 1)                                                                           
armaselect(yc, nbmod=5)
eacf(yc)



# Real data example
rec=ts(read.table("recruit.dat",header=F))
plot(rec)
acf(rec)
pacf(rec)
armaselect(rec, nbmod=5)
eacf(rec)




# Describing Stationarity and Trend Characteristics via KPSS test

library("tseries")

# a) Deterministic Trend
DT=ts(read.table("DT.txt",header=F)) 
plot(DT)
kpss.test(DT,null=c("Level"))
kpss.test(DT,null=c("Trend"))
model1=lm(DT~time(DT))
summary(model1)
plot(y=rstudent(model1),x=as.vector(time(DT)), ylab='Standardized Residuals',xlab='Time',type='o')
z=rstudent(model1)
# Trend removed from the series. Residuals seems stationary.


# b) Stochastic Trend
ST=ts(read.table("Stoch.txt",header=F)) 
plot(ST)
kpss.test(ST,null=c("Level"))
kpss.test(ST,null=c("Trend"))
model2=lm(ST~time(ST))
summary(model2)
plot(y=rstudent(model2),x=as.vector(time(ST)), ylab='Standardized Residuals',xlab='Time',type='o')
# Trend can not removed via detrending, residuals seems nonstationary. Differencing is needed.
dif=diff(ST)
plot(dif)
# Differenced seeries seems stationary.



# Random walk with drift simulation  "x[t] = x[t - 1] + rho + w[t]"
x=w=rnorm(1000)
for (t in 2:1000) x[t] = x[t - 1] + 0.15 + w[t]
plot(x, type = "l")
acf(x)
pacf(x)



#example on transformation. when you need stabilize the variance.
data(airpass)
plot(airpass)
BoxCox.ar(airpass,method = c("yule-walker"))
BoxCox.ar(airpass,method = c("ols"))
plot(log(airpass))
plot(sqrt(airpass))

# when you want to transform the series for specific value of lambda use following codes
library(forecast)
lambda <- BoxCox.lambda(airpass)
lambda
plot(BoxCox(airpass,lambda))

#PART 4

library(TSA)

#1 Augmented Dickey Fuller (ADF) unit root test

ST=ts(read.table("Stoch.txt",header=F)) 
acf(ST)
pacf(ST)
library(fUnitRoots)
adfTest(ST, lags=1, type="c")
adfTest(ST, lags=1, type="ct")

std=diff(ST)
plot(std)
adfTest(std, type="nc")


# philip-perron unit root test
pp.test(ST)
pp.test(std)



#2
data(nporg)  # Data set contains the fourteen U.S. economic time series used by Nelson Plosser #
gnp = ts(na.omit(nporg[, "gnp.r"]))
plot(gnp)

adfTest(gnp, type="c")
adfTest(gnp, type="ct")

pp.test(gnp)

kpss.test(gnp,null=c("Level"))
kpss.test(gnp,null=c("Trend"))

dfgnp=diff(gnp)
plot(dfgnp)

adfTest(dfgnp, type="nc")

pp.test(dfgnp)

kpss.test(dfgnp,null=c("Level"))

df2gnp=diff(dfgnp)
plot(df2gnp)
kpss.test(df2gnp,null=c("Level"))


#3
ar=ts(read.table("ar.txt",header=T))

plot(ar,ylab=expression(Y[t]),type='o');
acf(ar);
pacf(ar);

kpss.test(ar,null=c("Level"))
kpss.test(ar,null=c("Trend"))


adfTest(ar, type="c")
adfTest(ar, type="ct")

pp.test(ar)


ardif=diff(ar)
plot(ardif)
kpss.test(ardif,null=c("Level"))
adfTest(ardif, type="nc")
pp.test(ardif)

acf(ardif)
pacf(ardif)
eacf(ardif)
library(caschrono)
armaselect(ardif)

##PART 5

library(TSA)
library(forecast)

# Seasonal  Time Series – Seasonal Unit root 

# HEGY test

install.packages("pdR")
library(pdR)

?HEGY.test
# HEGY.test (wts, itsd, regvar=0, selectlags=list(mode="signf", Pmax=NULL))

#1#
data(SP)
plot(SP)
lambda <- BoxCox.lambda(SP)
lambda
plot(BoxCox(SP,lambda))
SP=(BoxCox(SP,lambda))
acf(as.vector(SP))
pacf(as.vector(SP))
outsp=HEGY.test(wts=SP, itsd=c(1,1,c(1:3)),regvar=0, selectlags=list(mode="aic", Pmax=8))
outsp
outsp$stats
outsp$regvarcoefs
outsp$hegycoefs
summary(outsp)


#2#
data(tempdub)
plot(tempdub)
lambda <- BoxCox.lambda(tempdub)
lambda
plot(BoxCox(tempdub,lambda))
tempdub=(BoxCox(tempdub,lambda))
acf(as.vector(tempdub))
pacf(as.vector(tempdub))
outtempdub=HEGY.test (tempdub, itsd= c(1,0,c(1:11)), regvar=0, selectlags=list(mode="aic", Pmax=12))
outtempdub
outtempdub$stats
outtempdub$regvarcoefs
outtempdub$hegycoefs



# ndiffs(x, alpha=0.05, test=c("kpss","adf", "pp"), max.d=2)
# nsdiffs(x, m=frequency(x), test=c("ch"), max.D=1)
# nsdiffs(x, m=frequency(x), test=c("ocsb"), max.D=1)
help(nsdiffs)

#3#
prod=ts(read.table("prod.dat",header=F),start=1950,frequency=12)
plot(prod)
nsdiffs(prod, m=12, test=c("ocsb"), max.D=1)
nsdiffs(prod, m=12, test=c("ch"), max.D=1)

ndiffs(prod, test=c("kpss"), max.d=4)
ndiffs(diff(prod), test=c("kpss"), max.d=4)
nsdiffs(diff(prod), m=12, test=c("ocsb"), max.D=1)
nsdiffs(diff(prod), m=12, test=c("ch"), max.D=1)

##PART 6

library(forecast)
library(fpp)

# Simple Exponential Smoothing (SES)

data(oil)
oildata <- window(oil,start=1996,end=2007)
plot(oildata, ylab="Oil (millions of tonnes)",xlab="Year")

fit1 <- ses(oildata, alpha=0.2, initial="simple", h=3)
fit2 <- ses(oildata, alpha=0.6, initial="simple", h=3)
fit3 <- ses(oildata, h=3)


plot(window(oil,start=1996,end=2010))
lines(fitted(fit1), col="blue", type="o")
lines(fitted(fit2), col="red", type="o")
lines(fitted(fit3), col="green", type="o")
lines(fit1$mean, col="blue", type="o")
lines(fit2$mean, col="red", type="o")
lines(fit3$mean, col="green", type="o")
legend("topleft",lty=1, col=c(1,"blue","red","green"), 
       c("data", expression(alpha == 0.2), expression(alpha == 0.6),
         expression(alpha == 0.89)),pch=1)




# DES with exponential trend or dampened trend (additive or multiplicative)
livestock2 <- window(livestock,start=1970,end=2000)
fit1 <- ses(livestock2)
fit2 <- holt(livestock2)
fit3 <- holt(livestock2,exponential=TRUE)
fit4 <- holt(livestock2,damped=TRUE)


accuracy(fit1) 
accuracy(fit2)
accuracy(fit3)
accuracy(fit4)


plot(fit2$model$state)
plot(fit3$model$state)
plot(fit4$model$state)


plot(fit3, type="o", ylab="Livestock, sheep in Asia (millions)", flwd=1, plot.conf=FALSE)
lines(window(livestock,start=2001),type="o")
lines(fit1$mean,col=2)
lines(fit2$mean,col=3)
lines(fit4$mean,col=5)

legend("topleft", lty=1, pch=1, col=1:5, c("Data","SES","Holt's","Exponential", "Additive Damped"))



# State space representation of Exponential Smoothing Methods
# We need this for obtaining forecast error variance to construct CI.

# ets(y, model="ZZZ", damped=NULL, alpha=NULL, beta=NULL,
#     gamma=NULL, phi=NULL, additive.only=FALSE, lambda=NULL,
#     lower=c(rep(0.0001,3), 0.8), upper=c(rep(0.9999,3),0.98),
#     opt.crit=c("lik","amse","mse","sigma","mae"), nmse=3,
#     bounds=c("both","usual","admissible"),
#     ic=c("aicc","aic","bic"), restrict=TRUE)

## Example quarterly data
x=ts(read.table("w.txt",header=F),start=1950,frequency=4)
plot(x)

fit1=ets(ts(x[1:180],frequency=4,start=1950), model="ZZZ")
fr1=forecast(fit1)
plot(forecast(fit1))

fit2=ets(ts((x[1:180]),frequency=4,start=1950),model='MMM')
fr2=forecast(fit2)
plot(forecast(fit2))

fit3=auto.arima(ts(x[1:180],frequency=4,start=1950))
fr3=forecast(fit3)
plot(forecast(fit3))

accuracy(fr1)
accuracy(fr2)
accuracy(fr3)


## Example monthly data
prcp=ts(read.table("prc.txt",header=F),start=1950,frequency=12)
tsdisplay(prcp)
plot(prcp)


fit1=ets(ts(prcp[1:708],frequency=12,start=1950), model="ZZZ")
plot(fit1)
fr1=forecast(fit1)
plot(forecast(fit1))

fit2=ets(ts((prcp[1:708]+0.001),frequency=12,start=1950),model='MMM')
plot(fit2)
fr2=forecast(fit2)
plot(forecast(fit2))

fit3=auto.arima(ts(prcp[1:708],frequency=12,start=1950))
fr3=forecast(fit3)
plot(forecast(fit3))

accuracy(fr1,ts(prcp[709:732],frequency=12, start=2009))
accuracy(fr2,ts(prcp[709:732],frequency=12, start=2009))
accuracy(fr3,ts(prcp[709:732],frequency=12, start=2009))