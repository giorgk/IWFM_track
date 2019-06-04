# http://www.corrmap.com/features/homography_transformation.php

Xo <- c(5.0115,  0.3571, 6.0253, 8.2834)
Yo <- c(0.7143,  2.6676, 9.3440, 3.8047)
Zo <- c(1,1,1,1)
Xomat <- matrix(data = c(Xo,Yo,Zo), nrow=4)

Xd <- c(0,  0, 1, 1)
Yd <- c(0,  1, 1, 0)
Zd <- c(1,1,1,1)
Xdmat <- matrix(data = c(Xd,Yd,Zd), nrow=4)

XL <- matrix(data = 0, nrow = 8, ncol = 8)

XL[1:4,1:3] <- Xomat
XL[5:8,4:6] <- Xomat
XL[1:4,7] <- -Xomat[,1]*Xdmat[,1]
XL[5:8,7] <- -Xomat[,1]*Xdmat[,2]
XL[1:4,8] <- -Xomat[,2]*Xdmat[,1]
XL[5:8,8] <- -Xomat[,2]*Xdmat[,2]


cc <- solve(XL,c(Xd,Yd))

htransf <- function(x,y,cf){
  xp <- (cf[1]*x+cf[2]*y+cf[3])/(cf[7]*x+cf[8]*y+1)
  yp <- (cf[4]*x+cf[5]*y+cf[6])/(cf[7]*x+cf[8]*y+1)
  return(c(xp,yp))
}

tst <- htransf(0.9643,    2.7843,cc)
