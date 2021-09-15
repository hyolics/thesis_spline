library(splines)
Knots   = seq(-1,1,length=10)
# number of splines
M = (length(Knots)-4)^2 
bspline = function(x) {splineDesign(Knots,x,outer.ok = T)}
btens   = function(x) {t(bspline(x[1]))%*%bspline(x[2])}
# numebr of points to plot 
ng      = 51
# create vectors for plotting 
xgr     = seq(-1,1,length=ng)
xgr2    = expand.grid(xgr,xgr)
alpha   = rnorm(M) #alpha
Bx      = apply(xgr2,1,btens)
Bmat    = matrix(t(Bx)%*%alpha,ng)
# contour(xgr,xgr,Bmat)
# persp(xgr,xgr,Bmat,theta=15)
library(rgl)
open3d()
plot3d(xgr,xgr,Bmat, col = rainbow(n), size=2)

###################################
library(splines)
n <- 50
n_knot <- 10
M <- (n_knot-4)^2 # number of splines
e <- rnorm(n, 0, 0.05)
x<-seq(0, 1, length.out=n)
X <- expand.grid(x,x)

knot <- seq(0,1,length=n_knot)
alpha <- rnorm(M) 

# knot <- matrix(runif(n_knot), n_knot/2, 2)
# alpha <- runif(M, 0, 1)

bspline <- function(x) {
  splineDesign(knot,x,outer.ok = T)
}
btens <- function(x) {
  t(bspline(x[1]))%*%bspline(x[2])
}



Bx <- apply(X_grid, 1, btens)
y0 <- t(Bx)%*%alpha
y <- y0 + e

data <- data.frame(y.spline=y, X)
colnames(data) <- c("y.spline", "x1", "x2")
nrow(data)
head(data)

library(rgl)
open3d()
plot3d(y, X, X, col = rainbow(n), size=2)
plot3d(data$y.spline, data$x1, data$x2, col = rainbow(n), size=2)
Bmat    = matrix(t(Bx)%*%alpha,n)
contour(X, X, Bmat)
persp(X, X, Bmat,theta=15)


data <- data.frame(y.spline=y, x1= X, x2=X)
nrow(data)
head(data)