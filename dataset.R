#use: cf function.R, Algorithm1.R ,SARS.R
##cf_fit() result:A1_opt(change delta), SARS_opt(change lm for no -1)

#1 dim
library(splines)
library(data.table)
library(foreach)

set.seed(1)
n <- 100
order <- 3
x <- seq(0, 1, length.out=n)
e <- rnorm(n, 0, 0.1)
e1 <- 0.1
e2 <- 0.01
n_data <- 100

raw_data <- function(x, knot, order, alpha, e) {
  basis <- bs(x=x, knots=knot, degree=order, Boundary.knots=c(0,1),intercept=TRUE)
  y.spline <- basis%*%alpha + e
  dt <- data.table(x, y.spline=as.vector(y.spline))
  return(dt=dt)
}

#1
knot <- c(0.3, 0.45, 0.9)
alpha <- scale(c(-10, -5, 5, 10, 15, 20, 1))
cf_fit(knot, alpha, e, "f1", 5) #0.007389765(0.007387377) 0.007503943(0.007170740)
legend("topleft",ncol=1, col=c("black","red","blue"), lwd=c(2,2,2) ,lty=c(1,1,1), legend=c("real func.","Algorithm1","Algorithm2"), cex=0.8, bty="n")

cf_location(knot, alpha, e, "f1", 5) 
cf_delta(knot, alpha, e, 3)
cf_time(knot, alpha, e, 3)

cf_dataset(knot, alpha, e1, 5)
cf_dataset(knot, alpha, e2, 5)

#2
knot <- c(0.16, 0.22, 0.5, 0.7)
alpha <- c(0.8, -0.16, -0.22, 0.35, -0.53, 0.71, -0.22, 0.8)
cf_fit(knot, alpha, e, "f2", 3) #0.007083954(0.007288528) 0.007279812(0.007084091)
legend("topleft",ncol=3, col=c("black","red","blue"), lwd=c(2,2,2) ,lty=c(1,1,1), legend=c("real func.","Algorithm1","Algorithm2"), cex=0.8, bty="n")

cf_location(knot, alpha, e, "f2", 3) 
cf_delta(knot, alpha, e, 3)
cf_time(knot, alpha, e, 3)

cf_dataset(knot, alpha, e1, 3)
cf_dataset(knot, alpha, e2, 5)

#3
data3_e0 <- data.frame(x, y.spline=sin(20*x))
knot <- c(0.16 ,0.25, 0.51, 0.83)
alpha <- scale(lm(y.spline~bs(x, knots=knot, degree=order),data3_e0)$coe)
cf_fit(knot, alpha, e, "f3", 5) #0.006944444(0.006977317) 0.020681079(0.020681079)
legend("topleft",ncol=3, col=c("black","red","blue"), lwd=c(2,2,2) ,lty=c(1,1,1), legend=c("real func.","Algorithm1","Algorithm2"), cex=0.8, bty="n")

cf_location(knot, alpha, e, "f3", 3) 
cf_delta(knot, alpha, e, 3)
cf_time(knot, alpha, e, 3)

cf_dataset(knot, alpha, e1, 5)
cf_dataset(knot, alpha, e2, 5)

#4
knot <- c(0.2, 0.5, 0.5, 0.5, 0.5, 0.8)
alpha <- scale(c(-5,-3,-2,5,25,25,5,-2,-3,-5))
cf_fit(knot, alpha, e, "f4", 5) #0.007288702(0.007291319) 0.007646512(0.007642966)
legend("topleft",ncol=1, col=c("black","red","blue"), lwd=c(2,2,2) ,lty=c(1,1,1), legend=c("real func.","Algorithm1","Algorithm2"), cex=0.8, bty="n")

cf_location(knot, alpha, e, "f4", 5) 
cf_delta(knot, alpha, e, 3)
cf_time(knot, alpha, e, 5)

cf_dataset(knot, alpha, e1, 4)
cf_dataset(knot, alpha, e2, 5)

#5
knot <- c(0.1, 0.5, 0.5, 0.5, 0.5, 0.8)
alpha <- scale(c(-10,-5,-2,20,1,20,25,10,1,-1))
cf_fit(knot, alpha, e, "f5", 5) #0.008295713(0.008336688) 0.008362527(0.008295930)
legend("topleft",ncol=1, col=c("black","red","blue"), lwd=c(2,2,2) ,lty=c(1,1,1), legend=c("real func.","Algorithm1","Algorithm2"), cex=0.8, bty="n")

cf_location(knot, alpha, e, "f5", 5) 
cf_delta(knot, alpha, e, 3)
cf_time(knot, alpha, e, 5)

cf_dataset(knot, alpha, e1, 5)
cf_dataset(knot, alpha, e2, 5)

#6
knot <- c(0.3, 0.3, 0.3, 0.4, 0.6, 0.7)
alpha <- scale(c(3, 5, 4, 10, 4, 7, 8, 8, 6, 4))
cf_fit(knot, alpha, e, "f6", 5) #0.007193820(0.006945821) 0.006980891(0.006620028)
legend("topleft",ncol=1, col=c("black","red","blue"), lwd=c(2,2,2) ,lty=c(1,1,1), legend=c("real func.","Algorithm1","Algorithm2"), cex=0.8, bty="n")

cf_location(knot, alpha, e, "f6", 5) 
cf_delta(knot, alpha, e, 3)
cf_time(knot, alpha, e, 5)

cf_dataset(knot, alpha, e1, 5)
cf_dataset(knot, alpha, e2, 5)

###diss. special case.
knot <- c(0.4, 0.4, 0.4, 0.4, 0.7)#normal
alpha <- scale(c(-1.78, -0.33, 1.62, -5, 5, 1.85, -3, -0.35, 1.39))
cf_fit(knot, alpha, e, "f7", 5) #0.008466138 0.123289523(0.008452876)

knot <- c(0.3, 0.3, 0.3, 0.3, 0.7, 0.8)#normal, but no -1 SARS all better
alpha <- scale(c(2, 5, 4, 20, 4, 16, 10, 8, 5, 4))
cf_fit(knot, alpha, e, "f7", 5) #0.00718391 0.06184805(0.007058172)
cf_delta(knot, alpha, e, 3)
cf_location(knot, alpha, e, "f7", 5)

knot <- c(0.1, 0.5, 0.5, 0.5, 0.5, 0.7)#normal, but no -1 SARS all better
alpha <- scale(c(-10, -5, -2, 20, 1, 20, 25, 10, 1, -1))
cf_fit(knot, alpha, e, "f7", 5) #0.007489104 0.048184573(0.009312252)
cf_delta(knot, alpha, e, 3)
cf_location(knot, alpha, e, "f7", 5)


#2 dim
library(scatterplot3d)
library(rgl)
n <- 100
order <- 3
x <- seq(0, 1, length.out = n)
X <- expand.grid(x, x)

coe_basis2 <- function(x, knot){
  basis <- (x-knot)^order
  basis[x<knot] <- 0
  return(basis)
}

uni_basis <- function(x, knot){
  basis1 <- outer(x, 0:order, "^")
  basis2 <- matrix(NA, length(x), length(knot))
  for(i in 1:length(knot)){
    basis2[,i] <- coe_basis2(x, knot[i])
  }
  bs <- cbind(basis1, basis2)
  return(bs)
}

n_basis <- function(X, knot1, knot2){
  b1 <- uni_basis(X[,1], knot1)
  b2 <- uni_basis(X[,2], knot2)
  
  cross <- expand.grid(1:dim(b1)[2],1:dim(b2)[2])
  
  basis <- matrix(NA, dim(b1)[1], nrow(cross))
  for(i in 1:nrow(cross)){
    basis[,i] <- b1[,cross$Var1[i]]*b2[,cross$Var2[i]]
  }
  return(basis)
}
cf <- function(k, name){
  knot1 <- k[[2]]
  knot2 <- k[[3]]
  
  # plot(r1, r2, main="fitted vs raw knot")
  # points(knot1, knot2, col="red")
  n_basis <- function(X, knot1, knot2){
    b1 <- uni_basis(X[,1], knot1)
    b2 <- uni_basis(X[,2], knot2)
    
    cross <- expand.grid(1:dim(b1)[2],1:dim(b2)[2])
    
    basis <- matrix(NA, dim(b1)[1], nrow(cross))
    for(i in 1:nrow(cross)){
      basis[,i] <- b1[,cross$Var1[i]]*b2[,cross$Var2[i]]
    }
    return(basis)
  }
  fitted_y <- function(alpha_basis, fitted_basis){
    #summary(lm(y~B-1)$coefficients)
    alpha_fit <- matrix(lm(y~alpha_basis-1)$coefficients)
    alpha_fit[is.na(alpha_fit)] <- 0
    y_hat <- fitted_basis%*%alpha_fit
    
    return(y_hat)
  }
  
  cf_plot <- function(y1=y, y2, X1, X2, name){
    pmat <- persp(x, x, y0_mat, theta=45)
    points(trans3d(X1[,1],X1[,2],y1,pmat),col="red",pch=20)
    points(trans3d(X2[,1],X2[,2],y2,pmat),col="blue",pch=20)
    title(main=name)
    
    if(length(y1)==length(y2)){
      mse <- sum((y1-y2)^2)/n^2
      #title(sub=paste0("mse=",mse))
      return(mse)
    }
  }
  
  B_fit <- n_basis(X, knot1, knot2)
  y_hat <- fitted_y(B_fit, B_fit)
  cf_plot(y, y_hat, X, X, name)
  #cf_plot3d(y, y_hat,X,X)
}

#1
x <- seq(0, 1, length.out = n)
X <- expand.grid(x, x)
y0 <- sin(20*X[,1])+sin(20*X[,2])
e <- rnorm(n^2, 0, 0.1)
y <- y0 + e
y0_mat <- matrix(y0, n)
pmat <- persp(x, x, y0_mat, theta=45)
points(trans3d(X[,1], X[,2], y, pmat),col="red", pch=20)
title(main="f1")
scatterplot3d(X[,1], X[,2], y0)
plot3d(X[,1], X[,2], y0)

data1 <- data.frame(y.spline=y, x1=X[,1], x2=X[,2])
a1 <- Algorithm1_2dim(data1, 7, 1)
cf(a1, "f1")
legend("bottomleft",ncol=1, col=c("red","blue"),pch = c(19, 19), legend=c("real func.","Algorithm1"), cex=0.8, bty="n")

s1 <- SARS_2dim(data1)
cf(s1, "f1")
legend("bottomleft",ncol=1, col=c("red","blue"), pch = c(19, 19), legend=c("real func.","Algorithm2"), cex=0.8, bty="n")

cf_2time(data1)

#2
x <- seq(-1, 1, length.out = n)
X <- expand.grid(x, x)
e <- rnorm(n^2, 0, 0.1)
y0 <- exp(X[,1]*sin(pi*X[,2]))
y <- y0+e
y0_mat <- matrix(y0, n)
pmat <- persp(x, x, y0_mat, theta=45)
points(trans3d(X[,1], X[,2], y, pmat),col="red", pch=20)
title(main="f2")
scatterplot3d(X[,1], X[,2], y0)
plot3d(X[,1], X[,2], y0)


data2 <- data.frame(y.spline=y, x1=X[,1], x2=X[,2])
a2 <- Algorithm1_2dim(data2, 7, 1)
cf(a2, "f2")
legend("bottomleft",ncol=1, col=c("red","blue"),pch = c(19, 19), legend=c("real func.","Algorithm1"), cex=0.8, bty="n")

s2 <- SARS_2dim(data2)
cf(s2, "f2")
legend("bottomleft",ncol=1, col=c("red","blue"), pch = c(19, 19), legend=c("real func.","Algorithm2"), cex=0.8, bty="n")

cf_2time(data2)

#3
x <- seq(-2, 2, length.out = n)
X <- expand.grid(x, x)
e <- rnorm(n^2, 0, 0.1)
y0 <- 3*sin(X[,1]*X[,2])
y <- y0+e
y0_mat <- matrix(y0, n)
pmat <- persp(x, x, y0_mat, theta=45)
points(trans3d(X[,1], X[,2], y, pmat),col="red", pch=20)
title(main="f3")
scatterplot3d(X[,1], X[,2], y0)
plot3d(X[,1], X[,2], y0)

data3 <- data.frame(y.spline=y, x1=X[,1], x2=X[,2])
k3 <- Algorithm1_2dim(data3, 7, 2)
cf(k3, "f3")
legend("bottomleft",ncol=1, col=c("red","blue"),pch = c(19, 19), legend=c("real func.","Algorithm1"), cex=0.8, bty="n")

s3 <- SARS_2dim(data3)
cf(s3, "f3")
legend("bottomleft",ncol=1, col=c("red","blue"), pch = c(19, 19), legend=c("real func.","Algorithm2"), cex=0.8, bty="n")

cf_2time(data3)

#4
x <- seq(-1, 1, length.out = n)
X <- expand.grid(x, x)
e <- rnorm(n^2, 0, 0.1)
y0 <- sin(2*pi*X[,1])*X[,2]
y <- y0+e
y0_mat <- matrix(y0, n)
pmat <- persp(x, x, y0_mat, theta=45)
points(trans3d(X[,1], X[,2], y, pmat),col="red", pch=20)
title(main="f4")
scatterplot3d(X[,1], X[,2], y0)
plot3d(X[,1], X[,2], y0)

data4 <- data.frame(y.spline=y, x1=X[,1], x2=X[,2])
k4 <- Algorithm1_2dim(data4, 7, 1)
cf(k4, "f4")
legend("bottomleft",ncol=1, col=c("red","blue"),pch = c(19, 19), legend=c("real func.","Algorithm1"), cex=0.8, bty="n")

s4 <- SARS_2dim(data4)
cf(s4, "f4")
legend("bottomleft",ncol=1, col=c("red","blue"), pch = c(19, 19), legend=c("real func.","Algorithm2"), cex=0.8, bty="n")

cf_2time(data4)

#5
x <- seq(0, 1, length.out = n)
X <- expand.grid(x, x)
e <- rnorm(n^2, 0, 0.1)
y0 <- 40*exp(8*(X[,1]-0.5)^2+(X[,2]-0.5)^2)/(exp(8*(X[,1]-0.2)^2+(X[,2]-0.7)^2)+exp(8*(X[,1]-0.7)^2+(X[,2]-0.2)^2))
y <- y0+e
y0_mat <- matrix(y0, n)
pmat <- persp(x, x, y0_mat, theta=45)
points(trans3d(X[,1], X[,2], y, pmat),col="red", pch=20)
title(main="f5")
data5 <- data.frame(y.spline=y, x1=X[,1], x2=X[,2])
k5 <- Algorithm1_2dim(data5, 7, 1)
cf(k5)
s5 <- SARS_2dim(data5)
cf(s5)
