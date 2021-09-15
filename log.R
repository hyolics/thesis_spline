##係數要調
x1 <- outer(x,1:order-1,"^")
knot <- c(0.25, 0.5, 0.75)
x2 <- matrix(rep(0,3*n),c(3,1))
for(i in 1:3){
  x2[i,] <- ifelse(x-knot[i]>=0,(x-knot[i])^2,0)
}
###v2
# x2<- outer(x,knot,">")*outer(x,knot,"-")^order
# X <- cbind(x1,x2)

X <- cbind(x1, t(x2))
e <- rnorm(n,0,1^2)
alpha <- c(0.6, 1, 0.5, 0.2, 0.8, 0.3)
y <- rowSums(X*alpha)+e
y <- rowSums(X*alpha)+e
plot(x, y)
# curve()
# abline()

#----------------------------------
basis <- function(x, degree, i, knots) {
  if(degree == 0){
    B <- ifelse((x >= knots[i]) & (x < knots[i+1]), 1, 0)
  } else {
    if((knots[degree+i] - knots[i]) == 0) {
      alpha1 <- 0
    } else {
      alpha1 <- (x - knots[i])/(knots[degree+i] - knots[i])
    }
    if((knots[i+degree+1] - knots[i+1]) == 0) {
      alpha2 <- 0
    } else {
      alpha2 <- (knots[i+degree+1] - x)/(knots[i+degree+1] - knots[i+1])
    }
    B <- alpha1*basis(x, (degree-1), i, knots) + alpha2*basis(x, (degree-1), (i+1), knots)
  }
  return(B)
}
bs <- function(x, degree=3, interior.knots=NULL, intercept=FALSE, Boundary.knots = c(0,1)) {
  if(missing(x)) stop("You must provide x")
  if(degree < 1) stop("The spline degree must be at least 1")
  Boundary.knots <- sort(Boundary.knots)
  interior.knots.sorted <- NULL
  if(!is.null(interior.knots)) interior.knots.sorted <- sort(interior.knots)
  knots <- c(rep(Boundary.knots[1], (degree+1)), interior.knots.sorted, rep(Boundary.knots[2], (degree+1)))
  K <- length(interior.knots) + degree + 1
  B.mat <- matrix(0,length(x),K)
  for(j in 1:K) B.mat[,j] <- basis(x, degree, j, knots)
  if(any(x == Boundary.knots[2])) B.mat[x == Boundary.knots[2], K] <- 1
  if(intercept == FALSE) {
    return(B.mat[,-1])
  } else {
    return(B.mat)
  }
}

#---------------------------------------------------
plot.spline <- function(basisdata, points = FALSE) {
  p <- ggplot(data = basisdata$dt)
  if (points) p <- p + geom_point(aes(x=x, y = y), color = "grey75")  
  p <- p + 
    geom_line(aes(x = x, y = y.spline), color = "red", size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1), breaks = knots) +
    theme(panel.grid.minor = element_blank())
  return(p)
}

####################################

coe_p <- function(data, n=100){
  p <- 0
  delta <- 1.2^p*((max(data$x)-min(data$x))/n)
  delta_0 <- (max(data$x)-min(data$x))/5
  while(delta < delta_0){
    p <- p + 1
    delta <- 1.2^p*((max(data$x)-min(data$x))/n)
  }
  return(p-1)
}

l <- 1
gamma <- NULL
Sj <- NULL
bic <- NULL
for(p in 1:16){
  delta <- ((max(data$x)-min(data$x))/n)*1.2^p
  for(i in 1:n){
    J <- subset(data,(x>=(x[i]-delta))&(x<=(x[i]+delta)))
    
    X <- basis_X(J, x[i], l, 2)
    data_j <- data.frame(y.spline=J$y.spline, X)
    gamma[i] <- coe_gamma(data, data_j)
    
  }
  tmp <- data.frame(x=x, gamma=gamma)
  potential <- subset(tmp, gamma>5)
  #gamma
  S <- NULL
  SS <- NULL
  while(nrow(potential)!=0){
    xi <- which.max(potential$gamma)
    S <- potential[xi,]$x
    potential <- subset(potential, (x<=(potential[xi,]$x-delta))|(x>=(potential[xi,]$x+delta)))
    SS <- c(S, SS)
    #print(S)
  }
  
  Sj[[p]] <- list(p, SS)
  #print(Sj)
  #print(Sj[[p]][[2]])
  if (is.null(Sj[[p]][[2]])==FALSE){
    fit <- lm(y.spline~bs(x, df = NULL, degree=3, knots=Sj[[p]][[2]]), data)
    k <- length(Sj[[p]][[2]])
    rss <- sum(fit$residuals^2)
    bic[[p]] <- n*log(rss/n)+k*log(n)
    
  }else{
    bic[[p]] <- NA
  }
}
###################################################
#2-dim
library(splines)
x <- seq(0, 1,length.out = 100)
X <- expand.grid(x,x)
knot <- seq(0, 1,length.out = 10)

uni_basis <- function(x) {
  basis <- splineDesign(knot,x,outer.ok = T)
  return(basis)
}
uni_basis(x)
dim(t(uni_basis(X[,1]))%*%uni_basis(X[,2]))#6*6

n_basis <- function(x) {
  basis <- t(uni_basis(x[1]))%*%uni_basis(x[2])
  return(basis)
}
n_basis(X[1,]) #6*6
B <- apply(X, 1, n_basis) #36*10000
y0 <- t(B)%*%alpha
y <- y0 + e

bs_uni <- function(x, knot, order){
  basis <- bs(x=x, knots=knot, degree=order-1, Boundary.knots=c(0,1) ,intercept=TRUE)
  return(basis)
}

bs_tensor <- function(X, order, bs_uni){
  if (is.null(dim(X))){ 
    return(bs_uni(X, order)) 
  }
  n <- dim(X)[1]
  d <- dim(X)[2]
  
  basis_m<- NULL
  for(i in 1:d){
    basis_m[[i]] <- bs_uni(X[,i], knot[,i], order)
  }
  
  total <- (dim(basis_m[[1]])[2])^d
  t <- NULL
  for(i in 1:d){
    t[[i]] <- 1:total
  }
  
  ind_matrix <- as.matrix(expand.grid(t))
  bt <- matrix(0, n, total)
  for(j in 1:total){
    ind <- ind_matrix[j,]
    b_prod <- rep(1, n)
    for(i in 1:d){
      b_prod <- b_prod*basis_m[[i]][,ind[i]]
    }
    bt[,j] <- b_prod
    return(bt)
  }
}

generate_f_hat <- function(data_x, data_y, order, bs_uni){
  bs_data <- bs_tensor(data_x, order, bs_uni)
  coef_hat <- lm(data_y~bs_data-1)$coef
  f_hat <- function(x){
    bx_x <- bs_tensor(x, order, bs_uni)
    return(as.numeric(bx_x %*% coef_hat))
  }
  return(f_hat)
}
f_hat <- generate_f_hat(X, y, order-1, bs_uni)
mean((f_hat(X)-f)^2)

#random basis
n_basis <- function(x) {
  set.seed(NULL)
  s <- sample(1:n_alpha, n_alpha, replace = T)
  b1 <- uni_basis(x[1], knot1)
  b1_s <- matrix(0, n_alpha, n_alpha)
  for(i in length(s)){
    b1_s[,i] <- b1[,s[i]]
  }
  
  basis <- t(b1)%*%uni_basis(x[2], knot2)
  return(basis)
}


## check same length every knot
mse_knot <-NULL
mse_knot[1] <- sum((y0_0knot-y0)^2)/n^2
for(i in 1:length(Sj[[k]][[2]])){
  knot1 <- Sj[[k]][[2]][1:i]
  knot2 <- Sj[[k]][[3]][1:i]
  
  mse_knot[i+1] <- cf_func(knot1, knot2, paste0("knot num=",i+1))
  
}
mse_knot

plot(mse_knot, main="mse with knot num.")
lines(mse_knot)
