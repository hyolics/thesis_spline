library(rgl)
open3d()
x <- sort(rnorm(1000))
y <- rnorm(1000)
z <- rnorm(1000) + atan2(x,y)
plot3d( x, y, z, 
        col = rainbow(1000), size = 2)
M <- par3d("userMatrix")
play3d( par3dinterp( 
  userMatrix = list( M, rotate3d( M, pi/2, 1, 0, 0), 
                     rotate3d( M, pi/2, 0, 1, 0) ) ), duration = 4)

swissroll = function( n, sigma = 0.05){
  angle = (3*pi/2)*(1+2*runif(n));
  height = runif(n);
  xdata = cbind( angle*cos(angle), 
                 height, angle*sin(angle))
  xdata = scale(xdata) + matrix(rnorm(n*3, 0, sigma), n, 3)
  order.angle = order(angle)
  sort.angle = sort(order.angle, index.return=TRUE)
  col.id = rainbow(n)
  my.color = col.id[sort.angle$ix]
  colnames(xdata) = paste( "x", 1:3, sep = "")
  return( list( xdata = xdata, 
                angle = angle, 
                color = my.color))
}
swissdata = swissroll(500)
xdata = swissdata$xdata
x.color = swissdata$color
open3d()
plot3d( xdata[,1], xdata[,2], xdata[,3], 
        col=x.color, size=3, xlab="",
        ylab= "", zlab="", axes = FALSE)
