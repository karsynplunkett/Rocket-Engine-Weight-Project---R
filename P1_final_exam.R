#Karsyn Plunkett
#Final Exam Project 
#STAT 4610
#PROBLEM 1
setwd("C:/Users/Karsyn/Documents/STAT 4610")
engine <- read.csv("LiquidRocketEngineData.csv", header=TRUE)
attach(engine)
tWeight = log(Weight)
lm1<-lm(Weight~Thrust)
par(mfrow=c(1,1))
plot(engine$Thrust,engine$Weight,xlab="Thrust", ylab="Weight")

pairs(Weight~Thrust+StagedCombustion+GasGenerator+PressureFed+ExpanderCycle+ISP+US+Russia+France+Germany, data=engine,gap=0.4,cex.labels=1.5)
m1 <- lm(Weight~Thrust+StagedCombustion+GasGenerator+PressureFed+ExpanderCycle+ISP+US+Russia+France+Germany)
summary(m1)

par(mfrow=c(2,2))
plot(m1)

StanRes1 <- rstandard(m1)
par(mfrow=c(3,3))
plot(Thrust,StanRes1, ylab="Standardized Residuals")
plot(Cycle,StanRes1, ylab="Standardized Residuals")
plot(ISP,StanRes1, ylab="Standardized Residuals")
plot(Country,StanRes1, ylab="Standardized Residuals")
plot(Length,StanRes1, ylab="Standardized Residuals")
plot(Diameter,StanRes1, ylab="Standardized Residuals")
plot(CP,StanRes1, ylab="Standardized Residuals")
plot(Nozzle,StanRes1, ylab="Standardized Residuals")

par(mfrow=c(1,1))
fit1 <- m1$fitted.values
m2 <- lm(Weight~fit1 + I(fit1^2))
plot(fit1,Weight,xlab="Fitted Values")
fitnew <- seq(-15,60,len=76)
lines(fitnew,predict(m2,newdata=data.frame(fit1=fitnew)))
abline(lsfit(m1$fitted.values,Weight),lty=2)

library(alr3)
inverseResponsePlot(m1,key=TRUE)

install.packages("car")
install.packages("alr3")
library(alr3)
summary(tranxy <- powerTransform(Thrust + ISP))

avPlots(m1)

par(mfrow=c(1,1))
lnslr=lm(log(engine$Weight)~log(Thrust)+StagedCombustion+GasGenerator+PressureFed+ExpanderCycle+log(ISP)+US+Russia+France+Germany)
plot(log(engine$Thrust), log(engine$Weight), xlab = 'LN Thrust', ylab= 'LN Weight')
summary(lnslr)
abline(lsfit(log(engine$Thrust)+log(engine$ISP), log(engine$Weight)))

StanRes2 <- rstandard(lnslr)
par(mfrow=c(3,3))
plot(Thrust,StanRes2, ylab="Standardized Residuals")
plot(ISP,StanRes2, ylab="Standardized Residuals")

avPlots(lnslr)

par(mfrow=c(2,2))
plot(lnslr)
summary(lnslr)

library(leaps)
b <- regsubsets(as.matrix(X),log(Weight))
rs <- summary(b)
par(mfrow=c(1,2))
plot(1:7,rs$adjr2,xlab="Subset Size",ylab="Adjusted R-squared")
library(car)
subsets(b,statistic=c("adjr2"))
rs$adjr2


n <- length(lnslr$residuals)
#MODEL BASED ON BACKWARD SELECTION PARTB
backwardAIC <- step(lnslr,direction="backward", data=engine)
backwardBIC <- step(lnslr,direction="backward", data=engine, k=log(n))

#MODEL BASED ON FORWARD SELECTION PARTC
mint <- lm(log(Weight)~1,data=engine)
forwardAIC <- step(mint,scope=list(lower=~1, upper=~log(Thrust)+StagedCombustion+GasGenerator+PressureFed+ExpanderCycle+log(ISP)+US+Russia+France+Germany),direction="forward", data=engine)
forwardBIC <- step(mint,scope=list(lower=~1, upper=~log(Thrust)+StagedCombustion+GasGenerator+PressureFed+ExpanderCycle+log(ISP)+US+Russia+France+Germany),direction="forward", data=engine,k=log(n))


newmod =lm(log(engine$Weight)~log(Thrust)+log(ISP)+ExpanderCycle+Russia)
plot(log(engine$Thrust)+log(ISP), log(engine$Weight), xlab = 'LN Thrust', ylab= 'LN Weight')
summary(newmod)
newmod2 = lm(log(engine$Weight) ~ log(Thrust) + StagedCombustion + GasGenerator + PressureFed + log(ISP))
summary(newmod2)
par(mfrow=c(2,2))
plot(newmod)
par(mfrow=c(2,2))
plot(newmod2)
detach(engine)
