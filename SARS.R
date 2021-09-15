library(splines)
library(data.table)
#SARS.
SARS <- function(data, C){
    basis_X <-function(data, xi, l, order=2){
    X <- NULL
    x1 <- ifelse((data$x-xi)>0,(data$x-xi)^l, 0)
    x2 <- outer(data$x-xi,0:(order-1),"^")
    X <- cbind(x1, x2)
    return(X)
  }
  
  #X <- basis_X(data, 0.3, 0, 2)
  
  coe_gamma <- function(data, data_j){
    alpha_hat <- lm(y.spline~., data_j)$coe[1]
    sigma_hat <- median(abs(data$y.spline[-1]-data$y.spline[-length(data$y.spline)]))/(0.6745*sqrt(2))
    cov_matrix <- summary(lm(y.spline~.,data_j))$cov[1]
    var_hat <- cov_matrix*sigma_hat^2
    gamma <- abs(alpha_hat/sqrt(var_hat))
    #return(list(alpha_hat, sigma_hat, var_hat))
    return(gamma)
  }
  #coe_gamma(data, data_j)
  
  l <- 1
  delta_0 <- (max(data$x)-min(data$x))/5
  delta <- ((max(data$x)-min(data$x))/n)*1.2
  data_selection <- data
  Sp <- NULL
  bic <- NULL
  p <- 0
  while((nrow(data_selection)>1)&(delta < delta_0)){
    potential <- NULL
    data_selection <- data
    p <- p + 1
    delta <- ((max(data$x)-min(data$x))/n)*1.2^p
    
    for(i in 1:nrow(data_selection)){
      J <- subset(data_selection ,(x>=(x[i]-delta))&(x<=(x[i]+delta)))
  
      if(nrow(J)>4){
        X <- basis_X(J, x[i], l, 2)
        data_j <- data.frame(y.spline=J$y.spline, X)
        gamma <- coe_gamma(data, data_j)
      
        if(gamma>C){
          tmp <- data_selection$x[i]
          potential <- c(tmp, potential)
          
          data_selection <- subset(data_selection, (x>=(tmp+delta))|(x<=(tmp-delta)))
          #print(nrow(data_selection))
        }
      }
    }
    
    Sp[[p]] <- list(c(p, delta), potential)
    if (is.null(Sp[[p]][[2]])==FALSE){
      fit <- lm(y.spline~bs(x, df = NULL, degree=3, knots=Sp[[p]][[2]]), data)
      k <- length(Sp[[p]][[2]])
      sigma_hat <- median(abs(data$y.spline[-1]-data$y.spline[-length(data$y.spline)]))/(0.6745*sqrt(2))
      rss <- sum(fit$residuals^2)+ 2*(k+4)*sigma_hat/n
      bic[[p]] <- n*log(rss/n)+k*log(n)
      
    }else{
      bic[[p]] <- NA
    }
  }
  #Sp
  #bic
  return(Sp[[which.min(bic)]])
}

#fit
fit_SARS <- function(data_e0, data, C){
  plot(data, lwd=2, main="SARS.")
  lines(data_e0, lwd=2)
  abline(v=knot, lty=3)
  
  result <- SARS(data, C)[[2]]
  abline(v=result, col="red", lty=3)
  fit <- lm(y.spline~bs(x, knots=result, degree=(length(result)-1)),data)
  lines(x, predict(fit), type="l", col="red" ,lwd=2)
  legend("topleft",ncol=2, col=c("black","red","black","red"), lwd=c(2,2,3,3) ,lty=c(1,1,3,3), legend=c("func.","fit func.","raw knot","result"), cex=0.8, bty="n")
  
  #summary knot
  #list(raw=sort(knot),result=sort(result_s1))
  
  #MSE:0.002177784
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

plot(data, lwd=2, main="SARS.")
lines(data_e0, lwd=2)
abline(v=knot, lty=3)

fit_SARS(data_e0, data, 5)
