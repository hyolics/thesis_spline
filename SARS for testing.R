SARS_opt <- function(data, C){
  basis_X <-function(data, xi, l, order=3){
    X <- NULL
    x1 <- ifelse((data$x-xi)>0,(data$x-xi)^l, 0)
    x2 <- outer(data$x-xi,0:(order-1),"^")
    X <- cbind(x1, x2)
    return(X)
  }
  
  coe_gamma <- function(data, data_j){
    alpha_hat <- lm(y.spline~.-1, data_j)$coe[1]
    sigma_hat <- median(abs(data$y.spline[-1]-data$y.spline[-length(data$y.spline)]))/(0.6745*sqrt(2))
    cov_matrix <- summary(lm(y.spline~.-1,data_j))$cov[1]
    var_hat <- cov_matrix*sigma_hat^2
    gamma <- abs(alpha_hat/sqrt(var_hat))
    return(gamma)
  }
  
  rss <- function(init_knot){
    sigma_hat <- median(abs(data$y.spline[-1]-data$y.spline[-length(data$y.spline)]))/(0.6745*sqrt(2))
    f_hat <- lm(y.spline~bs(x, knots=init_knot, degree=3),data)
    R <- sum((data$y.spline-predict(f_hat))^2)/n + 2*(3+order+1)*sigma_hat/n
    return(R)
  }
  
  l <- 1
  
  delta_0 <- (max(data$x)-min(data$x))/5
  delta <- ((max(data$x)-min(data$x))/n)*1.2
  data_selection <- data
  Sp <- NULL
  p <- 0
  potential <- NULL
  while(delta<delta_0){
    delta <- ((max(data$x)-min(data$x))/n)*1.2^p

    for(i in 1:nrow(data_selection)){
      J <- subset(data_selection ,(x>=(x[i]-delta))&(x<=(x[i]+delta)))
      
      if(nrow(J)>5){
        X <- basis_X(J, x[i], l, 2)
        data_j <- data.frame(y.spline=J$y.spline, X)
        gamma <- coe_gamma(data, data_j)
        
        if((!is.na(gamma)) & (gamma>C)){
          tmp <- data_selection$x[i]
          potential <- c(tmp, potential)
          data_selection <- subset(data_selection, (x>=(tmp+delta))|(x<=(tmp-delta)))
        }
      }
    }
    p <- p+1
  }
  return(potential)
}

SARS_f(data, 5)

#fit
fit_SARS_f <- function(data_e0, data, C){
  plot(data, lwd=2, main="SARS. with f")
  lines(data_e0, lwd=2)
  abline(v=knot, lty=3)
  
  result <- SARS_f(data, C)
  abline(v=result , col="red", lty=3)
  
  fit <- lm(y.spline~bs(x, knots=result, degree=3),data)
  lines(x, predict(fit), type="l", col="red" ,lwd=2)
  legend("topleft",ncol=2, col=c("black","red","black","red"), lwd=c(2,2,3,3) ,lty=c(1,1,3,3), legend=c("func.","fit func.","raw knot","result"), cex=0.8, bty="n")
  
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
alpha <- c(0.6, 0.3, 0.5, 0.2, 0.8, 0.4, 0.2, 0.7)
data_e0 <- raw_data(x, knot, order, alpha, e=0)
e <- rnorm(n ,0, 0.05)
data <- raw_data(x, knot, order, alpha, e)

plot(data, lwd=2, main="SARS. with f")
lines(data_e0, lwd=2)
abline(v=knot, lty=3)

fit_SARS_f(data_e0, data, 5)

