#gglasso using a unified BMD algorithm
library(gglasso)

#############first-step:grouup lasso for change-points#########
##序列生成
Y <- c()
set.seed(123)
e<-rnorm(1024)
Y[1] <- e[1]
for (i in 2:512) Y[i] <- 0.9*Y[i-1]+e[i]
for (i in 513:768) Y[i] <- 1.69*Y[i-1]-0.81*Y[i-2]+e[i]
for (i in 769:1024) Y[i] <- 1.32*Y[i-1]-0.81*Y[i-2]+e[i]
n <- length(Y)
plot(1:1024, Y, type='l')
abline(v=512, col="red")
abline(v=769, col="red")

#design matrix and respond vector
order <- 5  #maximum order of AR model
YY <- Y[seq(-1,-order,by=-1)]
x <- cbind(Y[c(-1,-2,-3,-4,-1024)],
           Y[c(-1,-2,-3,-1023,-1024)],
           Y[c(-1,-2,-1022,-1023,-1024)],
           Y[c(-1,-1021,-1022,-1023,-1024)],
           Y[c(-1020,-1021,-1022,-1023,-1024)])
XX<-x
system.time(
  for (i in 1:(n-order-1)){
    x[seq(1,i), ] <- rep(0, time=order)
    XX <- cbind(XX, x)
  }
)
XX <- as.matrix(XX)

#group lasso
grp <- rep(1:(n-order), each=order)
system.time(fit <- gglasso(x=XX, y=YY, group=grp, loss='ls', nlambda=100, intercept=FALSE)) 
system.time(fit.cv <- cv.gglasso(x=XX, y=YY, group=grp, nfolds=10)) 
plot(fit.cv)
lambda <- fit.cv$lambda.1se 
coef <- coef.gglasso(object=fit, s=lambda)[-1]
length(which(coef!=0))/order-1 #the number of change-points in first-step

#coefficient matrix/box
cname <- c()
for (i in 1:order) cname[i] <- paste("c", i, sep="")
rname <- c()
for (i in 1:(n-order)) rname[i] <- paste("r", i, sep="")
matrix.coef<-matrix(coef, ncol=order, byrow=TRUE, dimnames=list(rname, cname))
location.cp <- which(matrix.coef[ ,1]!=0)[-1] #the location of change-points in first-step

#fitting value VS real value
Y.pre <- predict(fit, XX, s=lambda, type="link")
plot(1:(n-order), YY, type="l")
lines(1:(n-order), Y.pre, type="l", col=2)
legend("topleft", c("real", "prediction"), horiz=TRUE, pch=c(21,21), col=c(1,2))
abline(v=507)
abline(v=764)
for (i in 1:length(location.cp)) abline(v=location.cp[i], col=2) #估计变点

#group.coef
group.coef <- matrix.coef[which(matrix.coef[ ,1]!=0), ]
for (i in 2:dim(group.coef)[1]) group.coef[i, ] <- group.coef[i, ]+group.coef[i-1, ]


#############second-step:IC for optimal change-points################
###IC=n*ln(SSR/n)+ln(n)*K
ncoef <- coef
ic <- c()
init.ic <- (n-order)*log((sum((Y.pre-YY)^2))/(n-order))+length(which(coef!=0))*log(n-order)
ic[1] <- init.ic
m <- ic[1]
j <- 1
while (m<=ic[j]) {
  j <- j+1
  ic[j] <- m
  ma.noid <- matrix(which(ncoef!=0), ncol=order, byrow=TRUE)[-1, ]
  icc <- c()
  for (i in 1:dim(ma.noid)[1]) {
    nncoef <- ncoef
    nncoef[ma.noid[i, ]] <- rep(0, time=order)
    icc[i] <- (n-order)*log(sum((XX%*%nncoef-YY)^2)/(n-order))+log(n-order)*length(which(nncoef!=0))
  }
  m <- min(icc)
  print(m)
  last.coef <- ncoef[ma.noid[which(icc==m), ]]
  ncoef[ma.noid[which(icc==m), ]] <- rep(0, time=order)
}
ncoef[ma.noid[which(icc==m), ]] <- last.coef
e.coef <- matrix(ncoef, ncol=order, byrow=TRUE)
location.ecp <- which(e.coef[ ,1]!=0)[-1] #the location of change-points in second-step
end.coef <- matrix(ncoef[which(ncoef!=0)], ncol=order, byrow=TRUE)
for (i in 2:dim(end.coef)[1]) end.coef[i, ] <- end.coef[i, ]+end.coef[i-1, ]

###results comparision
plot(1:(n-order), YY, type="l", col=1, lwd=2, xlab="time", ylab="amount")
lines(1:(n-order), Y.pre, type="l", col=2, lwd=2)
lines(1:(n-order), XX%*%ncoef, type="l", col=3, lwd=2)
legend("bottomright", cex=0.75, c("original","first-step","second-step"), 
       lwd=c(2,2), col=c(1,2,3), bty="n")
for (i in 1:length(location.ecp)) abline(v=location.ecp[i], col=2)
