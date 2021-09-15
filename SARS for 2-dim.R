#init.
library(rgl)

n <- 100
order <- 3
x <- seq(0, 1, length.out = n)
plot(x,sin(20*x))
data <- data.frame(x, y.spline=sin(20*x))
knot1 <- r1 <- Algorithm1_opt(data,7)[[2]]
knot2 <- r2 <- Algorithm1_opt(data,7)[[2]]
fit <- lm(y.spline~bs(x, knots=knot1, degree=order),data)
lines(x, predict(fit), type="l", col="red" ,lwd=2)

n_knot <- length(knot1)
n_alpha <- ((order+1)+n_knot)

e <- rnorm(n^2, 0, 0.05)
X <- expand.grid(x, x)
y <- sin(20*X[,1])+sin(20*X[,2])

#generate data
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

B <- n_basis(X, knot1, knot2)
alpha <- lm(y~B-1)$coefficients
y0 <- B%*%alpha
y <- y0 + e

data <- data.frame(y.spline=y, x1=X[,1], x2=X[,2])

y0_mat <- matrix(y0, n)
pmat <- persp(x, x, y0_mat, theta=45)
points(trans3d(X[,1], X[,2], y, pmat),col="red", pch=20)
title(main="real data")

# plot3d(X[,1], X[,2], y0)
# points3d(X[,1], X[,2], y, col="red")

##check infulence of error: real knot wiht fitted
fitted_y <- function(alpha_basis, fitted_basis){
  #summary(lm(y~B-1)$coefficients)
  alpha_fit <- matrix(lm(y~alpha_basis-1)$coefficients)
  alpha_fit[is.na(alpha_fit)] <- 0
  y_hat <- fitted_basis%*%alpha_fit
  
  return(y_hat)
}

SARS_2dim <- function(data){
  basis_X <-function(x, xi, l, order=2){
    X <- NULL
    x1 <- ifelse((x-xi)>0,(x-xi)^l, 0)
    x2 <- outer(x-xi,0:(order-1),"^")
    X <- cbind(x1, x2)
    return(X)
  }
  
  fitted_y <- function(alpha_basis, fitted_basis){
    #summary(lm(y~B-1)$coefficients)
    alpha_fit <- matrix(lm(y~alpha_basis-1)$coefficients)
    alpha_fit[is.na(alpha_fit)] <- 0
    y_hat <- fitted_basis%*%alpha_fit
    
    return(y_hat)
  }
  
  fit_model <- function(data_j){
    m1 <- lm(y.spline~.-1, data_j)
    m2 <- lm(y.spline~(X1+X2+X3+X4)-1, data_j)
    m3 <- lm(y.spline~(X1+X2+X3)-1, data_j)
    m4 <- lm(y.spline~(X1+X2)-1, data_j)
    m5 <- lm(y.spline~(X1)-1, data_j)
    model <- list(m5, m4, m3, m2, m1)
    f <- anova(m5, m4, m3, m2, m1)$`Pr(>F)`
    f[c(2:5)][is.na(f[c(2:5)])] <- 1
    if(na.omit(all(f[c(2:5)]>=0.05))){
      return(m1)
    }else{
      return(model[[which.min(f)]])
    }
  }
  
  coe_gamma <- function(data, data_j){
    alpha_hat <- na.omit(fit_model(data_j)$coe)
    sigma_hat <- median(abs(data$y.spline[-1]-data$y.spline[-length(data$y.spline)]))/(0.6745*sqrt(2))
    cov_matrix <- na.omit(diag(summary(lm(y.spline~.-1,data_j))$cov))
    var_hat <- cov_matrix*sigma_hat^2
    gamma <- abs(alpha_hat/sqrt(var_hat))
    #return(list(alpha_hat, sigma_hat, var_hat))
    return(gamma)
  }
  
  circle_region <- function(data, knot_x1, knot_x2, delta){
    d_circle <- (data$x1-knot_x1)^2 + (data$x2-knot_x2)^2
    data <- data.frame(data, d_circle)
    data_selection <- subset(data, d_circle<=delta^2)
    return(data_selection)
  }
  
  l <- 1
  
  d <- dist(X)
  delta_0 <- (max(d)-min(d))/5
  delta <- ((max(d)-min(d))/n)*1.5
  
  data_selection <- data
  Sp <- NULL
  bic <- NULL
  p <- 0
  while((nrow(data_selection)>1)&(delta<delta_0)){
    data_selection <- data
    S_x1 <- NULL
    S_x2 <- NULL
    SS_x1 <- NULL
    SS_x2 <- NULL
    p <- p + 1
    delta <- ((max(d)-min(d))/n)*1.5^p
    
    for(i in 1:nrow(data_selection)){
      J <- circle_region(data_selection, data_selection[i,]$x1, data_selection[i,]$x2, delta)
      
      if(nrow(J)>5){
        X1 <- basis_X(J$x1, data_selection[i,]$x1, l, 2)
        X2 <- basis_X(J$x2, data_selection[i,]$x2, l, 2)
        
        XX <- matrix(c(X1[,1]*X2[,1], X1[,1]*X2[,c(2,3)], X1[,c(2,3)]*X2[,1]), nrow=dim(X1)[1])
        data_j <- data.frame(y.spline=J$y.spline, XX)
        
        gamma <- coe_gamma(data, data_j)
        
        if(all(gamma>5)){
          S_x1 <- data_selection$x1[i]
          S_x2 <- data_selection$x2[i]
          data_selection <- subset(data_selection, !(x1 %in% J$x1)&!(x2 %in% J$x2))
          SS_x1 <- c(S_x1, SS_x1)
          SS_x2 <- c(S_x2, SS_x2)
        }
      }
    }
    
    Sp[[p]] <- list(c(p, delta), SS_x1, SS_x2)
    if (is.null(Sp[[p]][[2]])==FALSE){
      knot1 <- Sp[[p]][[2]]
      knot2 <- Sp[[p]][[3]]
      K <- length(knot1)
      B_fit <- n_basis(X, knot1, knot2)
      y_hat <- fitted_y(B_fit, B_fit)
      sigma_hat <- median(abs(data$y.spline[-1]-data$y.spline[-length(data$y.spline)]))/(0.6745*sqrt(2))
      rss <- sum((y_hat-y)^2) #+ 2*(2*K+4)*sigma_hat/n^2
      bic[[p]] <- n^2*log(rss/n^2)+(K*2)*log(n^2)
      
    }else{
      bic[[p]] <- NA
    }
  }
  #return(c(Sp, bic))
  return(Sp[[which.min(bic)]])
}

k <- SARS_2dim(data)
k

cf <- function(k){
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
      title(sub=paste0("mse=",mse))
      return(mse)
    }
  }
  
  B_fit <- n_basis(X, knot1, knot2)
  y_hat <- fitted_y(B_fit, B_fit)
  cf_plot(y, y_hat, X, X, "fitted vs raw data")
  #cf_plot3d(y, y_hat,X,X)
}
cf(k)
#####################################################
