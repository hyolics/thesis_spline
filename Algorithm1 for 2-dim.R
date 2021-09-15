#二維想法: 10000個點
#一個圓內的分群隨機畫一條直線，確定為節點後，變成剔除一個圓內的點 
#不要用bs function，改成自己寫，才能將比較重要的基底挑出來，讓效果比較明顯
#演算法部分，先用線性，挑出來結果是一組(x1,x2)

#4/29測試:knot是否真的需要，把所有節點拿掉fit
#5/13修正:
#產生資料直接放X，產出10000*2矩陣
#fitted是用lm產生估計係數
#畫圖要畫1.真實資料(y0+e) 2.用真實點fitted 3.估計點fitted 4.沒用點fitted
#5/20:原來方程式效果不佳，改成推廣sin(20*x)至二維
#5/27:因現在分群是隨機抽，為使結果穩定，改多抽幾次，並取其中最小p作為最終結果再比較。

#init.
library(rgl)

n <- 100
order <- 3
x <- seq(0, 1, length.out = n)
plot(x,sin(20*x))
data <- data.frame(x, y.spline=sin(20*x))
knot1 <- r1 <- c(0.08579784, 0.86429079, 0.54973985, 0.23467282, 0.39295810, 0.70674221)
knot2 <- r2 <- c(0.08579784, 0.86429079, 0.54973985, 0.23467282, 0.39295810, 0.70674221)
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

cf_plot3d <- function(y1=y, y2, X1, X2, name){
  plot3d(X1[,1], X1[,2], y1,col="red")
  points3d(X2[,1], X2[,2], y2, col="blue")
}

y_hat <- fitted_y(B, B)
cf_plot(y, y_hat,X,X, "influence of error(real knot fitted)")
#cf_plot3d(y, y_hat,X,X)

##check infulence of knot: no knot with fitted 
uni_basis_0knot <- function(x){
  basis1 <- outer(x, 0:order, "^")
  bs <- cbind(basis1)
  return(bs)
}

n_basis_0knot <- function(X) {
  b1 <- uni_basis_0knot(X[,1])
  b2 <- uni_basis_0knot(X[,2])
  
  cross <- expand.grid(1:dim(b1)[2],1:dim(b2)[2])
  basis <- matrix(NA, dim(b1)[1], nrow(cross))
  for(i in 1:nrow(cross)){
    basis[,i] <- b1[,cross$Var1[i]]*b2[,cross$Var2[i]]
  }
  return(basis)
}

B_0knot <- n_basis_0knot(X)#10000*16
y_0knot <- fitted_y(B_0knot, B_0knot)
cf_plot(y, y_0knot, X, X, "influence of knot(no knot fitted)")
#cf_plot3d(y, y_0knot, X, X)

cf_plot(y_hat, y_0knot, X, X, "real knot vs. no knot fitted")
#cf_plot3d(y_hat, y_0knot, X, X)

# ###test data
# x0 <- seq(0, 1, length.out = 20)
# X0 <- expand.grid(x0, x0)#400*2
# B_fit <- n_basis(X0, knot1, knot2)
# 
# y_hat <- fitted_y(B, B_fit)
# cf_plot(y, y_hat, X, X0, "influence of error(real knot fitted)")
# #cf_plot3d(y, y_hat, X, X0)
# 
# B_fit_0knot <- n_basis_0knot(X0)
# y_0knot <- fitted_y(B_0knot, B_fit_0knot )
# cf_plot(y, y_0knot,X,X0, "influence of knot(no knot fitted)")
# #cf_plot3d(y, y_0knot,X,X0)
# 
# cf_plot(y_hat, y_0knot, X0, X0, "real knot vs. no knot fitted")
# #cf_plot3d(y_hat, y_0knot, X0, X0)

#Algorithm1
Algorithm1_2dim <- function(data,J,n_sample){
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
  
  fitted_y <- function(alpha_basis, fitted_basis){
    #summary(lm(y~B-1)$coefficients)
    alpha_fit <- matrix(lm(y~alpha_basis-1)$coefficients)
    alpha_fit[is.na(alpha_fit)] <- 0
    y_hat <- fitted_basis%*%alpha_fit
    
    return(y_hat)
  }
  
  coe_alpha <- function(data, order=2){
    coe <- lm(y.spline~x1+x2, data)$coe
    return(coe)
  }
  
  coe_sigma <- function(data1, data2){
    cov1 <- summary(lm(y.spline~x1+x2, data1))$cov 
    cov2 <- summary(lm(y.spline~x1+x2, data2))$cov
    if((nrow(cov1)!=3)|(nrow(cov2)!=3)){
      sigma <- matrix(0, 3, 3)
    }else{
      sigma <- matrix(cov1+cov2, 3, 3)
    }
    return(sigma)
  }
  
  circle_region <- function(data, knot_x1, knot_x2, delta){
    d_circle <- (data$x1-knot_x1)^2 + (data$x2-knot_x2)^2
    data <- data.frame(data, d_circle)
    data_selection <- subset(data, d_circle<=delta^2)
    return(data_selection)
  }
  
  split_region <- function(data, knot_x1, knot_x2, a){
    criterion <- a*data$x1+data$x2-knot_x1+knot_x2 #ax+y-x1+x2
    data_selection <- data.frame(data, criterion)
    return(data_selection)
  }
  
  phi_test <- function(data, knot_x1, knot_x2, delta, order=2){
    set.seed(NULL)
    a <- rnorm(1)
    data <- split_region(circle_region(data, knot_x1, knot_x2, delta), knot_x1, knot_x2, a)
    data1 <- subset(data, criterion<0)
    data2 <- subset(data, criterion>0)
    
    if ((nrow(data1)>2)&(nrow(data2)>2)){
      D <- coe_alpha(data1, order=2)-coe_alpha(data2, order=2)
      S <- coe_sigma(data1, data2)
      sigma_hat <- median(abs(data$y.spline[-1]-data$y.spline[-length(data$y.spline)]))/(0.6745*sqrt(2))
      if(det(S)!=0){
        W <- (1/sigma_hat^2)%*%t(D)%*%solve(S)%*%D
      
        p_value <- pchisq(W, order, lower.tail = F)
        }else{
          p_value <- 1
          }
      }else{
        p_value <- 1
      }
    
    return(p_value)
  }
  

  Sj <- NULL
  bic <- NULL
  for(j in 1:J){
    delta <- 1/2^j
    set_T <- matrix(NA, n_sample, nrow(data))
    for (time in 1:n_sample){
      for(i in 1:nrow(data)){
        set_T[time,i]  <- phi_test(data, data[i,]$x1, data[i,]$x2, delta, order=2)
      }
    }
    set_T <- apply(set_T,2,min)
    
    tmp <- data.frame(x1=X[,1], x2=X[,2], set_T)
    potential <- subset(tmp, set_T<0.05)
    
    S_x1 <- NULL
    S_x2 <- NULL
    SS_x1 <- NULL
    SS_x2 <- NULL
    while(nrow(potential)!=0){
      xi <- which.min(potential$set_T)
      S_x1 <- potential[xi,]$x1
      S_x2 <- potential[xi,]$x2
      region <- circle_region(data, S_x1, S_x2, delta)
      potential <- subset(potential, !(x1 %in% region$x1)&!(x2 %in% region$x2))
      SS_x1 <- c(S_x1, SS_x1)
      SS_x2 <- c(S_x2, SS_x2)
    }
    Sj[[j]] <- list(j, SS_x1, SS_x2)
    
    if (length(Sj[[j]][[2]])>1){
      knot1 <- Sj[[j]][[2]]
      knot2 <- Sj[[j]][[3]]
      K <- length(knot1)
      B_fit <- n_basis(X, knot1, knot2)
      y_hat <- fitted_y(B_fit, B_fit)
      rss <- sum((y_hat-y)^2)
      bic[[j]] <- n^2*log(rss/n^2)+(K*2)*log(n^2)
      
    }else{
      bic[[j]] <- NA
    }
  }
  #Sj
  #bic
  #Sj[[which.min(bic)]]
  return(Sj[[which.min(bic)]])
}

#MSE
##with fitted
# k0 <- Algorithm1_2dim(data, 7, 1)
# k1 <- Algorithm1_2dim(data, 7, 2)
# k2 <- Algorithm1_2dim(data, 7, 3)
# k3 <- Algorithm1_2dim(data, 7, 4)


cf <- function(k,time){
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
  cf_plot(y, y_hat, X, X, paste("fitted vs raw data",time))
  #cf_plot3d(y, y_hat,X,X)
}

cf(k0,"1")
cf(k1,"2")
cf(k2,"3")
cf(k3,"4")


