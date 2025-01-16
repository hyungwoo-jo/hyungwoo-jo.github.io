library(survival);library(dplyr);library(ggplot2);library(purrr);library(cmprsk)

# makedata --------------------

set.seed(10)
ftime <- rexp(200)
fstatus <- sample(0:2,200,replace=TRUE)
cov <- matrix(runif(200),nrow=200)
dimnames(cov)[[2]] <- c('x1')
cov
print(z <- crr(ftime,fstatus,cov))
summary(z)
z.p <- predict(z,rbind(c(.1,.5,.8),c(.1,.5,.2)))
plot(z.p,lty=1,color=2:3)
crr(ftime,fstatus,cov,failcode=2)
# quadratic in time for first cov
crr(ftime,fstatus,cov,cbind(cov[,1],cov[,1]),function(Uft) cbind(Uft,Uft^2))

### survival package
pdata <- finegray(Surv(ftime, fstatus) ~ cov, data=cbind(ftime, fstatus, cov))
pdata


fgfit <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 1,
               weight=fgwt, data=pdata, model = T)


