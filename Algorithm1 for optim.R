library(splines)
library(data.table)
library(foreach)
#Algorithm 1. opt
Algorithm1_opt <- function(data, J){
  coe_alpha <- function(data,order=2){
    coe <- lm(y.spline~x, data)$coe
    return(coe)
  }
  
  coe_sigma <- function(data1, data2, order=2){
    cov1 <- summary(lm(y.spline~x, data1))$cov 
    cov2 <- summary(lm(y.spline~x, data2))$cov
    sigma <- matrix(cov1+cov2, order, order)
    return(sigma)
  }
  
  rss <- function(init_knot){
    k <- 3
    order <- 2
    #sigma_hat <- median(abs(data$y.spline[-1]-data$y.spline[-length(data$y.spline)]))/(0.6745*sqrt(2))
    f_hat <- lm(y.spline~bs(x, knots=init_knot, degree=k),data)
    R <- sum((data$y.spline-predict(f_hat))^2)/n #+ 2*(3+order+1)*sigma_hat/n
    return(R)
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
  
  
  J <- 7
  Sj <- NULL
  bic <- NULL
  for(j in 1:J){
    delta <- 1/2^j
    #delta <- (max(data$x)-min(data$x))/n*J*1.2
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
      # S <- last(optim(potential[xi,]$x, fn=gradient, method = "CG")$par)
      # S <- potential$x[which.min(abs(S-potential$x))]
      # potential <- subset(potential, (x<=(xi-delta))|(x>=(xi+delta)))
      potential <- subset(potential, (x<=(potential[xi,]$x-delta))|(x>=(potential[xi,]$x+delta)))
      SS <- c(S, SS)
      #print(S)
    }
    Sj[[j]] <- list(j, SS)
    #print(Sj)
    #print(Sj[[j]][[2]])
    if (is.null(Sj[[j]][[2]])==FALSE){
      fit <- lm(y.spline~bs(x, degree=3, knots=Sj[[j]][[2]]), data)
      k <- length(Sj[[j]][[2]])
      r <- sum(fit$residuals^2)
      bic[[j]] <- n*log(r/n)+k*log(n)
      
    }else{
      bic[[j]] <- NA
    }
    
  }
  #Sj
  #bic
  
  #Sj[[which.min(bic)]]
  
  result_op <- optim(Sj[[which.min(bic)]][[2]], fn=rss, method = "CG")$par
  return(result_op)
}

#fit
fit_Algorithm1_opt <- function(data_e0,data){
  plot(data,lwd=2,main="Algorithm 1 with opt")
  lines(data_e0, lwd=2)
  abline(v=knot, lty=3)
  
  result <- Algorithm1_opt(data, 7)
  abline(v=result, col="red", lty=3)
  fit <- lm(y.spline~bs(x, knots=result, degree=3),data)
  lines(predict(fit), type="l", col="red" ,lwd=2)
  legend("topleft",ncol=2, col=c("black","red","black","red"), lwd=c(2,2,3,3) ,lty=c(1,1,3,3), legend=c("func.","fit func.","raw knot","result"), cex=0.8, bty="n")
  
  #summary knot
  #list(raw=sort(knot),result=sort(result_a1_op))
  
  #MSE
  mse <- sum((predict(fit)-data$y.spline)^2)/n
  return(mse)
}

#init.
set.seed(1)
n <- 100
order <- 3
x <- seq(0, 1, length.out=n)

#generate data
raw_data <- function(x, knot, order, alpha, e) {
  basis <- bs(x=x, knots=knot, degree=order, Boundary.knots=c(0,1),intercept=TRUE)
  y.spline <- basis%*%alpha + e
  dt <- data.table(x, y.spline=as.vector(y.spline))
  return(dt)
}

knot <- c(0.02, 0.3, 0.65, 0.96)
knot <- runif(4, 0, 1)
alpha <- c(0.6, 0.3, 0.5, 0.2, 0.8, 0.4, 0.2, 0.7)
alpha <- rnorm(8, 0, 1)
data_e0 <- raw_data(x, knot, order, alpha, e=0)
e <- rnorm(n ,0, 0.05)
data <- raw_data(x, knot, order, alpha, e)

plot(data,lwd=2,main="Algorithm 1 with opt")
lines(data_e0, lwd=2)
abline(v=knot, lty=3)

fit_Algorithm1_opt(data_e0, data)
