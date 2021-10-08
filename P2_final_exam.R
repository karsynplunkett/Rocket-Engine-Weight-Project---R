#Karsyn Plunkett
#Final Exam Project 
#STAT 4610
#PROBLEM 2
setwd("C:/Users/Karsyn/Documents/STAT 4610")
motor <- read.csv("SolidMotor.csv", header=TRUE)
attach(motor)
pairs(y2~x1+x2+x3+x4+x5+x6+x7+x8, data=motor,gap=0.4,cex.labels=1.5)
m1 <- lm(y2~x1+x2+x3+x4+x5+x6+x7+x8)
summary(m1)

library(alr3)
summary(tranxy <- powerTransform(x1+x2+x3+x4+x5+x6+x7+x8))

summary(tranxy <- powerTransform(Thrust + ISP))

StanRes1 <- rstandard(m1)
par(mfrow=c(3,3))
plot(x1,StanRes1, ylab="Standardized Residuals")
plot(x2,StanRes1, ylab="Standardized Residuals")
plot(x3,StanRes1, ylab="Standardized Residuals")
plot(x4,StanRes1, ylab="Standardized Residuals")
plot(x5,StanRes1, ylab="Standardized Residuals")
plot(x6,StanRes1, ylab="Standardized Residuals")
plot(x7,StanRes1, ylab="Standardized Residuals")
plot(x8,StanRes1, ylab="Standardized Residuals")


par(mfrow=c(1,1))
fit1 <- m1$fitted.values
m12 <- lm(y2~fit1 + I(fit1^2))
plot(fit1,y2,xlab="Fitted Values")
fitnew <- seq(-15,60,len=76)
lines(fitnew,predict(m12,newdata=data.frame(fit1=fitnew)))
abline(lsfit(m1$fitted.values,y2),lty=2)

par(mfrow=c(2,2))
plot(m1)

m2<- lm(y2~x1+x3+x4+x5+x6+x7+x8)
summary(m2)
library(alr3)
inverseResponsePlot(m1,key=TRUE)

library(alr3)
summary(tranxy <- powerTransform(x1 +x2+x3+x4+x5+x6+x7+x8 ))


model1 <-lm(y2~(x1+x2+x3+x4+x5+x6+x7+x8)^2+I(x1^2))
summary(model1)

model2<-lm(y2~(x1+x2+x3+x4+x5+x6+x7+x8)^3+I(x3^2))
summary(model2)

par(mfrow=c(2,2))
plot(model2)

logmod<-lm(log(y2)~(x1+x2+x3+x4+x5+x6+x7+x8)^2+I(x1^2))
summary(logmod)

par(mfrow=c(2,2))
plot(logmod)

StanRes2 <- rstandard(model2)
par(mfrow=c(3,3))
plot(x1,StanRes2, ylab="Standardized Residuals")
plot(x2,StanRes2, ylab="Standardized Residuals")
plot(x3,StanRes2, ylab="Standardized Residuals")
plot(x4,StanRes2, ylab="Standardized Residuals")
plot(x5,StanRes2, ylab="Standardized Residuals")
plot(x6,StanRes2, ylab="Standardized Residuals")
plot(x7,StanRes2, ylab="Standardized Residuals")
plot(x8,StanRes2, ylab="Standardized Residuals")



n <- length(logmod$residuals)
#MODEL BASED ON BACKWARD SELECTION PARTB
backwardAIC <- step(logmod,direction="backward", data=motor)
backwardBIC <- step(logmod,direction="backward", data=motor, k=log(n))

#MODEL BASED ON FORWARD SELECTION PARTC
mint <- lm(log(y2)~1,data=motor)
forwardAIC <- step(mint,scope=list(lower=~1, upper=~(x1+x2+x3+x4+x5+x6+x7+x8)^2),direction="forward", data=motor)
forwardBIC <- step(mint,scope=list(lower=~1, upper=~x1+x2+x3+x4+x5+x6+x7+x8)^2+I(x1^2),direction="forward", data=motor,k=log(n))
bmod = lm(log(y2) ~ x7 + x8 + x6 + x5 + x3 + x4 + x1 + x7:x8 + x7:x3 + x7:x1 + x7:x4 + x6:x3 + x7:x5 + x8:x6 + x7:x6 + x8:x5 + x8:x3)
summary(bmod)

par(mfrow=c(2,2))
plot(bmod)

###################################################
#PART B
pairs(y5~x1+x2+x3+x4+x5+x6+x7+x8, data=motor,gap=0.4,cex.labels=1.5)
mod <- lm(y5~x1+x2+x3+x4+x5+x6+x7+x8)
summary(mod)

par(mfrow=c(2,2))
plot(mod)

par(mfrow=c(1,1))
fit1 <- mod$fitted.values
mod1 <- lm(y5~fit1 + I(fit1^2))
plot(fit1,y5,xlab="Fitted Values")
fitnew <- seq(-15,60,len=76)
lines(fitnew,predict(mod1,newdata=data.frame(fit1=fitnew)))
abline(lsfit(mod$fitted.values,y5),lty=2)

library(alr3)
summary(tranxy <- powerTransform(x1 +x2+x3+x4+x5+x6+x7+x8 ))

model11 <-lm(y5~(x1+x2+x3+x4+x5+x6+x7+x8)^2)
summary(model11)

model12 <-lm(y5~(x1+x2+x3+x4+x5+x6+x7+x8)^3)
summary(model12)

par(mfrow = c(2,2))
plot(model12)

n <- length(model11$residuals)
#MODEL BASED ON BACKWARD SELECTION PARTB
backwardAIC <- step(model11,direction="backward", data=motor)
backwardBIC <- step(model11,direction="backward", data=motor, k=log(n))

#MODEL BASED ON FORWARD SELECTION PARTC
mint <- lm(y5~1,data=motor)
forwardAIC <- step(mint,scope=list(lower=~1, upper=~(x1+x2+x3+x4+x5+x6+x7+x8)^2),direction="forward", data=motor)
forwardBIC <- step(mint,scope=list(lower=~1, upper=~x1+x2+x3+x4+x5+x6+x7+x8)^2+I(x1^2),direction="forward", data=motor,k=log(n))

newmod = lm(y5 ~ x8 + x7 + x3 + x4 + x6 + x5 + x7:x3 + x8:x7 + x7:x4 + x7:x6 + x8:x3 + x3:x4)
summary(newmod)
newmod2 = lm(y5 ~ x3 + x4 + x5 + x6 + x7 + x8 + x3:x4 + x3:x7 + x3:x8 + x4:x7 + x6:x7 + x7:x8)
summary(newmod2)

par(mfrow =c(2,2))
plot(newmod)
