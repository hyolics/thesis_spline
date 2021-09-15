#1 dim
library(splines)
library(data.table)
library(foreach)

raw_data <- function(x, knot, order, alpha, e) {
  basis <- bs(x=x, knots=knot, degree=order, Boundary.knots=c(0,1),intercept=TRUE)
  y.spline <- basis%*%alpha + e
  dt <- data.table(x, y.spline=as.vector(y.spline))
  return(dt=dt)
}

t1 <- function(data){
  result_a1 <- Algorithm1_opt(data, 7)
  fit_a1 <- lm(y.spline~bs(x, knots=result_a1, degree=order),data)
  lines(x, predict(fit_a1), type="l", col="red" ,lwd=2)
  mse <- sum((predict(fit_a1)-data$y.spline)^2)/n
  return(mse)
}

t2 <- function(data, C){
  result_s1 <- SARS_opt(data, C)
  fit_s1 <- lm(y.spline~bs(x, knots=result_s1, degree=order),data)
  #abline(v=result_s1, lty=2, col="red")
  lines(x, predict(fit_s1), type="l", col="blue" ,lwd=2)
  mse <- sum((predict(fit_s1)-data$y.spline)^2)/n
  return(mse)
}

t2_org <- function(data){
  result_s1 <- SARS_org(data)
  fit_s1 <- lm(y.spline~bs(x, knots=result_s1, degree=order),data)
  #abline(v=result_s1, lty=2, col="red")
  lines(x, predict(fit_s1), type="l", col="green" ,lwd=2)
  mse <- sum((predict(fit_s1)-data$y.spline)^2)/n
  return(mse)
}

l1 <- function(data){
  result_a1 <- Algorithm1_opt(data, 7)
  abline(v=result_a1, col="red", lty=3)
  abline(v=knot,lty=2)
  
  fit_a1 <- lm(y.spline~bs(x, knots=result_a1, degree=order),data)
  lines(x, predict(fit_a1), type="l", col="red" ,lwd=2)
  #fit_a1_0knot <- lm(y.spline~bs(x, knots=NULL, degree=order),data)
  #lines(x, predict(fit_a1_0knot), type="l", col="red" ,lwd=2)
  #legend("topleft",ncol=2, col=c("black","red","black","red"), lwd=c(2,2,3,3) ,lty=c(1,1,3,3), legend=c("func.","fit func.","raw knot","result"), cex=0.8, bty="n")
  return(length(result_a1))
}

l2 <- function(data, C){
  result_a1 <- SARS_opt(data, C)
  abline(v=result_a1, col="blue", lty=3)
  
  fit_a1 <- lm(y.spline~bs(x, knots=result_a1, degree=order),data)
  lines(x, predict(fit_a1), type="l", col="blue" ,lwd=2)
  # fit_a1_0knot <- lm(y.spline~bs(x, knots=NULL, degree=order),data)
  # lines(x, predict(fit_a1_0knot), type="l", col="red" ,lwd=2)
  #legend("topleft",ncol=2, col=c("black","red","black","red"), lwd=c(2,2,3,3) ,lty=c(1,1,3,3), legend=c("func.","fit func.","raw knot","result"), cex=0.8, bty="n")
  return(length(result_a1))
}

l2_org <- function(data){
  result_a1 <- SARS_org(data)
  abline(v=result_a1, col="blue", lty=3)
  
  fit_a1 <- lm(y.spline~bs(x, knots=result_a1, degree=order),data)
  lines(x, predict(fit_a1), type="l", col="green" ,lwd=2)
  # fit_a1_0knot <- lm(y.spline~bs(x, knots=NULL, degree=order),data)
  # lines(x, predict(fit_a1_0knot), type="l", col="red" ,lwd=2)
  #legend("topleft",ncol=2, col=c("black","red","black","red"), lwd=c(2,2,3,3) ,lty=c(1,1,3,3), legend=c("func.","fit func.","raw knot","result"), cex=0.8, bty="n")
  return(length(result_a1))
}

cf_time <- function(knot, alpha, e, C){
  data_e0 <- raw_data(x, knot, order, alpha, e=0)
  data <- raw_data(x, knot, order, alpha, e)
  t <- Sys.time()
  a1 <- Algorithm1_opt(data, 7)
  t_a1 <- Sys.time() - t
  t <- Sys.time()
  s1 <- SARS_opt(data, C)
  t_s1 <- Sys.time() - t
  #return(list(t_a1, t_s1))
  return(list(100*(t_a1)/60, 100*(t_s1)/60))
}

cf_location <- function(knot, alpha, e, title, C){
  data_e0 <- raw_data(x, knot, order, alpha, e=0)
  data <- raw_data(x, knot, order, alpha, e)
  plot(data,ylab = "y",main = title)
  lines(data_e0 ,lwd=3)
  return(c(l1(data),l2(data, C)))
}

cf_delta <- function(knot, alpha, e, C){
  data <- raw_data(x, knot, order, alpha, e)
  return(c(Algorithm1(data, 7)[[1]],SARS(data, C)[[1]]))
}

cf_fit <- function(knot, alpha, e, title, C){
  data_e0 <- raw_data(x, knot, order, alpha, e=0)
  data <- raw_data(x, knot, order, alpha, e)
  plot(data_e0,type = "l",ylab = "y",main = title)
  abline(v=knot,lty=2)
  
  plot(data ,ylab = "y",main = title)
  lines(data_e0 ,lwd=3)
  return(c(t1(data),t2(data, C)))
}

cf_dataset <- function(knot, alpha, e, C){
  generate_error <- function(seed,error){
    set.seed(seed)
    e <- matrix(rnorm(n*n_data ,0, error),n_data,n)
    return(e)
  }
  
  MSE <- function(data, result_knot_1, result_knot_2){
    fit1 <- lm(y.spline~bs(x, knots=result_knot_1, df=3+length(result_knot_1)),data)
    fit2 <- lm(y.spline~bs(x, knots=result_knot_2, df=3+length(result_knot_2)),data)
    mse1 <- sum((predict(fit1)-data$y.spline)^2)/n
    mse2 <- sum((predict(fit2)-data$y.spline)^2)/n
    return(c(mse1, mse2))
  }
  
  cf <- matrix(NA, n_data, 2)
  num <- NULL
  for(i in 1:n_data){
    data <- data.frame(raw_data(x, knot, 3, alpha, generate_error(1, e)[i,]))
    
    num[[i]] <-c(length(Algorithm1_opt(data, 10)),length(SARS_opt(data, C)))
    cf[i,] <- MSE(data, Algorithm1_opt(data, 10), SARS_opt(data, C))
  }
  r1 <- apply(matrix(unlist(num), 2, n_data), 1, mean)
  
  result <- data.frame(cf, better=apply(cf,1,which.min))
  colnames(result) <- c("Algorithm1","SARS","better")
  r2 <- table(result$better)
  r3 <- apply(result,2,mean)
  r4 <- apply(result,2,median)
  r5 <- apply(result,2,var)
  return(list(r1, r2, r3, r4, r5))
}


#2 dim
n <- 100
order <- 3
x <- seq(0, 1, length.out = n)
X <- expand.grid(x, x)
n_2data <- 10

cf_2dim <- function(k){
  knot1 <- k[[2]]
  knot2 <- k[[3]]
  
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
  
  fitted_y <- function(alpha_basis, fitted_basis){
    alpha_fit <- matrix(lm(y~alpha_basis-1)$coefficients)
    alpha_fit[is.na(alpha_fit)] <- 0
    y_hat <- fitted_basis%*%alpha_fit
    
    return(y_hat)
  }
  
  
  MSE <- function(y1=y, y2, X1, X2){
    if(length(y1)==length(y2)){
      mse <- sum((y1-y2)^2)/n^2
      return(mse)
    }
  }
  
  B_fit <- n_basis(X, knot1, knot2)
  y_hat <- fitted_y(B_fit, B_fit)
  return (MSE(y, y_hat, X, X))
}

cf_2plot <- function(y1=y, y2, X1, X2, name){
  pmat <- persp(x, x, y0_mat, theta=45)
  points(trans3d(X1[,1],X1[,2],y1,pmat),col="red",pch=20)
  points(trans3d(X2[,1],X2[,2],y2,pmat),col="blue",pch=20)
  title(main=name)
}

cf_2time <- function(data){
  t <- Sys.time()
  a1 <- Algorithm1_2dim(data, 7, 1)
  t_a1 <- Sys.time() - t
  t <- Sys.time()
  s1 <- SARS_2dim(data)
  t_s1 <- Sys.time() - t
  return(list(10*(t_a1)/60, 10*(t_s1)/60))
}

cf_2dataset <- function(data, error){
  generate_error <- function(seed,error){
    set.seed(seed)
    e <- matrix(rnorm(n^2*n_2data ,0, error),n_2data,n)
    return(e)
  }
  
  cf <- matrix(NA, n_2data, 2)
  num <- NULL
  for(i in 1:n_2data){
    data <- data.frame(y.spline=y0+generate_error(1, error)[i,], x1=X[,1], x2=X[,2])
    
    num[[i]] <-c(length(Algorithm1_2dim(data, 7, 1)[[2]]),length(SARS_2dim(data)[[2]]))
    cf[i,] <- c(cf_2dim(Algorithm1_2dim(data, 7, 1)), cf_2dim(SARS_2dim(data)))
  }
  
  r1 <- apply(matrix(unlist(num), 2, n_2data), 1, mean)
  
  result <- data.frame(cf, better=apply(cf,1,which.min))
  colnames(result) <- c("Algorithm1","SARS","better")
  r2 <- table(result$better)
  r3 <- apply(result,2,mean)
  r4 <- apply(result,2,median)
  r5 <- apply(result,2,var)
  return(list(r1, r2, r3, r4, r5))
}

