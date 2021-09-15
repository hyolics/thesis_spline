#8/26: 先確定是不是delta大小問題
#改成選到節點後先去配適，計算出的殘差當作新的Y，再重新計算p value

library(splines)
library(data.table)
library(foreach)

#Algorithm 1.
Algorithm1_test <- function(data,J){
  coe_alpha <- function(data, order=2){
    coe <- lm(y.spline~x, data)$coe
    return(coe)
  }
  
  coe_sigma <- function(data1, data2, order=2){
    cov1 <- summary(lm(y.spline~x, data1))$cov 
    cov2 <- summary(lm(y.spline~x, data2))$cov
    if((nrow(cov1)!=2)|(nrow(cov2)!=2)){
      sigma <- matrix(0, 2, 2)
    }else{
      sigma <- matrix(cov1+cov2, 2, 2)
    }
    return(sigma)
    
    sigma <- matrix(cov1+cov2, order, order)
    return(sigma)
  }
  
  phi_test <- function(data, knot, delta, order=2){
    data1 <- subset(data,(x>=(knot-delta))&(x<=knot))
    data2 <- subset(data,(x<=(knot+delta))&(x>=knot))
    
    if ((nrow(data1)>1)&(nrow(data2)>1)){
      D <- coe_alpha(data1, order)-coe_alpha(data2, order)
      S <- coe_sigma(data1, data2, order)
      sigma_hat <- median(abs(data$y.spline[-1]-data$y.spline[-length(data$y.spline)]))/(0.6745*sqrt(2))
      W <- (1/sigma_hat^2)*t(D)%*%solve(S)%*%D
      
      p_value <- pchisq(W, order, lower.tail = F)
      
    }else{
      p_value <- NA
    }
    
    return(p_value)
    #return(list(W,p_value))
  }
  
  Sj <- NULL
  bic <- NULL
  for(j in 1:J){
    delta <- 1/2^j
    set_T <- foreach (i = 1:n, .combine='rbind') %do% {
      phi_test(data, data$x[i], delta, order=2)
    }
    
    tmp <- data.frame(x, set_T)
    potential <- subset(tmp, set_T<0.05)
    
    S <- NULL
    SS <- NULL
    while(nrow(potential)!=0){
      xi <- which.min(potential$set_T)
      S <- potential[xi,]$x
      SS <- c(S, SS)
      potential <- subset(potential, (x<=(potential[xi,]$x-delta))|(x>=(potential[xi,]$x+delta)))
      y.spline <-lm(y.spline~bs(x, knots=SS, degree=3),data)$residual
      
      set_T <- foreach (i = 1:nrow(potential), .combine='rbind') %do% {
        phi_test(potential, potential$x[i], delta, order=2)
      }
      potential$set_T <- set_T
      
      
      #print(S)
    }
    Sj[[j]] <- list(c(j, delta), SS)
    #print(Sj)
    #print(Sj[[j]][[2]])
    if (is.null(Sj[[j]][[2]])==FALSE){
      fit <- lm(y.spline~bs(x, df = NULL, degree=3, knots=Sj[[j]][[2]]), data)
      k <- length(Sj[[j]][[2]])
      rss <- sum(fit$residuals^2)
      bic[[j]] <- n*log(rss/n)+k*log(n)
      
    }else{
      bic[[j]] <- NA
    }
  }
  #Sj
  #bic
  return(Sj[[which.min(bic)]])
}

#fit
fit_Algorithm1 <- function(data_e0,data){
  plot(data,lwd=2,main="Algorithm 1")
  lines(data_e0, lwd=2)
  abline(v=knot, lty=3)
  
  result <- Algorithm1_test(data, 7)[[2]]
  abline(v=result, col="red", lty=3)
  fit <- lm(y.spline~bs(x, knots=result, degree=order),data)
  lines(x, predict(fit), type="l", col="red" ,lwd=2)
  
  fit_0knot <- lm(y.spline~bs(x, knots=NULL, degree=order),data)
  lines(x, predict(fit_0knot), type="l", col="green" ,lwd=2)
  legend("topleft",ncol=2, col=c("black","red","black","red"), lwd=c(2,2,3,3) ,lty=c(1,1,3,3), legend=c("func.","fit func.","raw knot","result"), cex=0.8, bty="n")
  
  #summary knot
  #list(raw=sort(knot),result=sort(result))
  
  #MSE
  mse <- sum((predict(fit)-data$y.spline)^2)/n
  return(mse)
}

#init.
set.seed(1)
n <- 100
order <- 3
x <- seq(0, 1, length.out=n)
e <- rnorm(n ,0, 0.05)

#generate data
raw_data <- function(x, knot, order, alpha, e) {
  basis <- bs(x=x, knots=knot, degree=order, Boundary.knots=c(0,1),intercept=TRUE)
  y.spline <- basis%*%alpha + e
  dt <- data.table(x, y.spline=as.vector(y.spline))
  return(dt)
}

knot <- c(0.02, 0.3, 0.65, 0.96)
alpha <- c(0.6, 0.3, 0.5, 0.2, 0.8, 0.4, 0.2, 0.7)
data_e0 <- raw_data(x, knot, order, alpha, e=0)
data <- raw_data(x, knot, order, alpha, e)

plot(data,lwd=2,main="Algorithm 1")
lines(data_e0, lwd=2)
abline(v=knot, lty=3)

fit_Algorithm1(data_e0, data)

