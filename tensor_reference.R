set.seed(1)
f <- function(x){ dnorm(x[,1]-0.5, sd=0.2)*dnorm(x[,2]-0.5, sd=0.2) }
n <- 5000
X <- matrix(runif(n*2), n,2)
y <- f(X) + rnorm(n,sd=0.4)
data <- data.frame(y.spline=y, x1=X[,1], x2=X[,2])

bx.trigo <- function(x, m){
  n <- length(x)
  b.trigo <- matrix(0, n, 2*m)
  for (k in 1:m){
    b.trigo[ ,k] <- cos(2*pi*k*x)
    b.trigo[ ,k+m] <- sin(2*pi*k*x)
  }
  return(cbind(rep(1,n), b.trigo))
}
### bx.tensor: compute tensor basis functions at x, x can be a matrix
bx.tensor <- function(x, m, bx.uni){
  if ( is.null(dim(x)) ) { return( bx.uni(x,m) ) }
  n <- dim(x)[1]
  d <- dim(x)[2]
  mat.list <- vector("list", d)
  for (i in 1:d){ mat.list[[i]] <- bx.uni(x[,i],m) }
  n.basis1 <- dim(mat.list[[1]])[2]
  n.basisd <- n.basis1^d
  v <- vector("list", d)
  for (i in 1:d){ v[[i]] <- 1:n.basis1 }
  ind.mat <- as.matrix(expand.grid(v))
  bx <- matrix(0, n, n.basisd)
  for (j in 1:n.basisd) {
    ind <- ind.mat[j,]
    b.prod <- rep(1, n)
    for (k in 1:d) { b.prod <- b.prod*mat.list[[k]][,ind[k]] }
    bx[, j] <- b.prod
  }
  return(bx)
}

get.fhat <- function(data.x, data.y, m, bx.uni){
  bx.data <- bx.tensor(data.x, m, bx.uni)
  coef.hat <- lm(data.y~bx.data-1)$coef
  fhat <- function(x){
    bx.x <- bx.tensor(x, m, bx.uni)
    return( as.numeric( bx.x %*% coef.hat ) )
  }
  return(fhat)
}
### generate data and compute fhat
set.seed(1)
f <- function(x){ dnorm(x[,1]-0.5, sd=0.2)*dnorm(x[,2]-0.5, sd=0.2) }
n <- 1000
X <- matrix(runif(n*2), n,2)

y <- f(X) + rnorm(n,sd=0.4)
fhat <- get.fhat(X, y, 3, bx.trigo)
### compute ISE using monte carlo integration based on 10000 monte carlo samples
set.seed(2)
n.mc <- 10000
u <- matrix(runif(n.mc*2),n.mc, 2)
mean((fhat(u) - f(u))^2 )
