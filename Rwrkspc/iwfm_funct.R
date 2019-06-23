library(h5)
library(pracma)
library(lubridate)


iwfm.findElemId <- function(x, y, BC, MSH, XY){
  dst_bc <- order(sqrt((BC[,1] - x)^2 + (BC[,2] - y)^2))
  #return(dst_bc)
  for (i in 1:10) {
    if (MSH[dst_bc[i],4] == 0)
      np <- 3
    else
      np <- 4

    xp <- XY[, 1][MSH[dst_bc[i],1:np]]
    yp <- XY[, 2][MSH[dst_bc[i],1:np]]

    if (inpolygon(x, y, xp, yp))
      return (dst_bc[i])
  }
  return(-9)
}

iwfm.findLayer <- function(xp, yp, zp, elid, MSH, XY, STRAT){
  xv <- XY[,1][MSH[elid,1:4]]
  yv <- XY[,2][MSH[elid,1:4]]
  zv <- STRAT[,][MSH[elid,1:4][1:4],]

  zlay <- iwfm.idw_interp(xp,yp,xv,yv,zv)
  #return(length(zlay)-1)

  if (zp > zlay[1]){
    out <- list("Lay" = 0)
    return(out)
  }
  if (zp < zlay[length(zlay)]){
    out <- list("Lay" = -9)
    return (-out)
  }

  for (i in 1:length(zlay)-1) {
    if (zp < zlay[i] && zp > zlay[i+1]){
      u = (zp - zlay[i+1])/(zlay[i] - zlay[i+1])
      out <- list("Lay" = i, "t" = u)
      return(out)
    }
  }
}

iwfm.loadData <- function(file){
  hf <- h5file(name =  file, mode = "r")
  hpartGroups <- list.groups(hf)
  hpartDataSets <- list.datasets(hf)
  XY <- hf["/geodata/XY"]
  MSH <- hf["/geodata/MSH"]
  BC <- hf["/geodata/BCEL"]
  STRAT <- hf["/geodata/STRAT"]
  FI <- hf["/geodata/FI"]
  VFLOW <- hf["flowdata/VFLOW"]
  HFLOW <- hf["flowdata/HFLOW"]
  HTCF <- hf["geodata/HTCF"]
  NRML <- hf["geodata/NRML"]
  FACEZ <- hf["geodata/FACEZ"]
  BCZ <- hf["geodata/bcZ"]

  # convert the flow data to arrays for faster indexing
  tmpH <- array(dim = dim(HFLOW))
  for (i in 1:dim(tmpH)[3]) {
    tmpH[,,i] <- HFLOW[,,i]
  }
  tmpV <- array(dim = dim(VFLOW))
  for (i in 1:dim(tmpV)[3]) {
    tmpV[,,i] <- VFLOW[,,i]
  }

  output <- list("XY" = XY, "MSH" = MSH, "BC" = BC, "STRAT" = STRAT, "FI" = FI, 
                 "HFLOW" = tmpH, "VFLOW" = tmpV, "HTCF" = HTCF, "NRML" = NRML, 
                 "FACEZ" = FACEZ, "BCZ" = BCZ )
  h5close(hf)
  return(output)
}

iwfm.idw_interp <- function(x,y, xv, yv, v, thres = 0.001){
  # given a set of points (xv, yv) with some values v associated
  # calculates the values on the point (x,y)
  # For each column of v returns a value
  dst <- sqrt((xv - x)^2 + (yv - y)^2)
  id <- which(dst <= thres)
  w <- vector(mode = "numeric", length = length(xv))
  if (length(id) == 1){
    w[id] <- 1
  }else{
    w <- dst/sum(dst)
  }

  if (is.null(dim(v))){
    vRet <- vector(mode = "numeric", length = 1)
    vRet[1] <- sum(w*v)
  }
  else{
    vRet <- vector(mode = "numeric", length = dim(v)[2])
    for (i in 1:dim(v)[2]) {
      vRet[i] <- sum(w*v[,i])
    }
  }

  return(vRet)
}

iwfm.calcVzRBF <- function(p, TM, ZV, BC, BCZ, S){
  XYdif <- repmat((BC[,1] - p[1])^2 + (BC[,2] - p[2])^2, dim(BCZ)[2], 1)
  Zdif <- array( (BCZ[,] - p[3])^2, dim = c(dim(BCZ)[1]*dim(BCZ)[2],1) )
  
  w <- exp(- ( XYdif/(2*S$XY^2) + Zdif/(2*S$Z^2)) )
  sumw <- sum(w)

  if (length(dim(ZV)) == 3){
    va <- ZV[, TM$T1, ]
    vb <- ZV[, TM$T2, ]
    vab <- va*TM$t + vb*(1-TM$t)
    vab <- array(vab, dim = c(dim(vab)[1]*dim(vab)[2],1) )
  }else if (length(dim(ZV)) == 2){
    vab <- array( ZV[,], dim = c(dim(ZV)[1]*dim(ZV)[2],1) )
  }
  vab <- c(vab, matrix(data = 0, nrow = dim(BCZ)[1], ncol = 1))

  return( sum(w*vab)/sumw )
}


iwfm.calcVxyRBF <- function(p, TM, HV, NRML, FCZ, S){
  XYdif <- repmat((NRML[,1] - p[1])^2 + (NRML[,2] - p[2])^2, dim(FCZ)[2], 1)
  Zdif <- array( (FCZ[,] - p[3])^2, dim = c(dim(FCZ)[1]*dim(FCZ)[2],1) )
  
  w <- exp(- ( XYdif/(2*S$XY^2) + Zdif/(2*S$Z^2)) ) 
  
  sumw <- sum(w)
  
  if (length(dim(HV)) == 3){
    va <- HV[, TM$T1, ]
    vb <- HV[, TM$T2, ]
    vab <- va*TM$t + vb*(1-TM$t)
    vab <- array(vab, dim = c(dim(vab)[1]*dim(vab)[2],1) )
  }else if (length(dim(HV)) == 2){
    vab <- array( HV[,], dim = c(dim(HV)[1]*dim(HV)[2],1) )
  }

  nx <- sum(w*repmat(NRML[,3], dim(FCZ)[2], 1 ) )/sumw
  ny <- sum(w*repmat(NRML[,4], dim(FCZ)[2], 1 ) )/sumw

  nxy <- c(nx,ny)/sqrt((nx^2+ny^2))

  return( (sum(w*vab)/sumw)*nxy )
}

iwfm.calcVz <- function(pz, el, L, TM, ZV){
  if (length(dim(ZV)) == 3){
    vTOPa <- ZV[el,TM$T1,L$Lay]
    vTOPb <- ZV[el,TM$T2,L$Lay]
    if (L$Lay < 4){
      vBOTa <- ZV[el, TM$T1, L$Lay+1]
      vBOTb <- ZV[el, TM$T2, L$Lay+1]
    }else{
      vBOTa <- 0
      vBOTb <- 0
    }
    # Calculate the vertical velocity between the top and bottom
    Vza <- vTOPa*L$t + vBOTa*(1-L$t)
    Vzb <- vTOPb*L$t + vBOTb*(1-L$t)
    
    vz <- Vza*TM$t + Vzb*(1-TM$t)
  }else if (length(dim(ZV)) == 2){
    vt <- ZV[el,L$Lay]
    if (L$Lay < 4){
      vb <- ZV[el,L$Lay + 1]
    }
    else{
      vb <- 0
    }
    vz <- vt*L$t + vb*(1-L$t)
  }
  return(vz)
}

iwfm.calcVel <- function(p, TM, el, L, HV, ZV, FI, cf){
  
  pu <- iwfm.UnitCoords(p[1], p[2], cf)
  # Check if we are inside a quad or triangle
  np <- 4
  if (FI[el,4] == 0){
    np <- 3
  }
  
  if (length(dim(HV)) == 3){
    # Horizontal Velocity at first time step
    va <- sign(FI[el,1:np]) * HV[abs(FI[el,1:np]), TM$T1, L$Lay]
    # Horizontal Velocity at second time step
    vb <- sign(FI[el,1:np]) * HV[abs(FI[el,1:np]), TM$T2, L$Lay]
    # interpolate between the time steps
    vab <- va*TM$t + vb*(1-TM$t)
  }else if (length(dim(HV)) == 2){
    vab <- sign(FI[el,1:np]) * HV[abs(FI[el,1:np]), L$Lay]
  }
  

  # Using the quad shape functions interpolate the xy velocities
  vxy <- iwfm.Quadinterp(pu[1], pu[2], vab)
  # take a step on the unit space
  pu2 <- pu + vxy
  # convert the point to real space
  p2 <- iwfm.CoordFromUnit(pu2[1], pu2[2], cf)
  
  # This is the normalized direction of velocity on the real space 
  vnrm <- (p2[1:2]-p[1:2])/sqrt(sum((p2[1:2]-p[1:2])^2))
  
  vxy <- sqrt(sum(vxy^2))*vnrm
  
  
  if (length(dim(ZV)) == 3){
    # Vertical velocity at the two time steps
    vTOPa <- ZV[el,TM$T1,L$Lay]
    vTOPb <- ZV[el,TM$T2,L$Lay]
    if (L$Lay < 4){
      vBOTa <- ZV[el, TM$T1, L$Lay+1]
      vBOTb <- ZV[el, TM$T2, L$Lay+1]
    }else{
      vBOTa <- 0
      vBOTb <- 0
    }
    # Calculate the vertical velocity between the top and bottom
    Vza <- vTOPa*L$t + vBOTa*(1-L$t)
    Vzb <- vTOPb*L$t + vBOTb*(1-L$t)
  
    vz <- Vza*TM$t + Vzb*(1-TM$t)
  }else if (length(dim(ZV)) == 2){
    vt <- ZV[el,L$Lay]
    if (L$Lay < 4){
      vb <- ZV[el,L$Lay + 1]
    }
    else{
      vb <- 0
    }
    vz <- vt*L$t + vb*(1-L$t)
  }
  

  return(c(vxy,vz)/30)
}

iwfm.pntVelRBF <- function(p, TM, D, simTime, po){
 
  v1xy <- iwfm.calcVxyRBF(p, TM, D$HFLOW, D$NRML, D$FACEZ, po$RBFSTD)
  v2z <- iwfm.calcVzRBF(p, TM, D$VFLOW, D$BC, D$BCZ, po$RBFSTD)
  return(c(v1xy,v2z))

}

iwfm.pntVel <- function(p, t, HV, ZV, FI, BC, MSH, XY, STRAT, HTCF, simTime){
  Ttp <- iwfm.interpTime(simTime, t)
  if ((Ttp$T1 == -10000 || Ttp$T1 == 10000) && length(dim(HV)) == 3){
    return(list("v" = c(0,0,0), "reas" = "outTime"))
  }
  else{
    elid <- iwfm.findElemId(p[1], p[2], BC, MSH, XY)
    if (elid < 0){
      return(list("v" = c(0,0,0), "reas" = "outDomain"))
    }
    else{
      lay <- iwfm.findLayer(p[1], p[2], p[3] ,elid, MSH, XY, STRAT)
      if (lay$Lay == 0){
        return(list("v" = c(0,0,0), "reas" = "outTop"))
      } else if (lay$Lay == -9){
        return(list("v" = c(0,0,0), "reas" = "outBot"))
      }
      else{
        
        #punit <- iwfm.UnitCoords(p[1], p[2], HTCF[elid,])
        v <- iwfm.calcVel(p, Ttp, elid, lay, HV, ZV, FI, HTCF[elid,])
        return(list("v" = v, "reas" = "outNO"))
      }
    }
  }
}

iwfm.SimTime <- function(sy = 1973, sm = 10, ey = 2015, em = 9){
  nmonths <- length(seq.Date(from = as.Date(paste0(sy,"/",sm,"/1")),to = as.Date(paste0(ey,"/",em,"/1")),by = "month"))

  tmp = vector(mode = "list", length = nmonths)
  tmp <-seq(ISOdate(sy,sm,30), by="month", length.out = nmonths)

  i <- 1
  y <- sy
  m <- sm

  while(TRUE){
    ndays <- days_in_month(as.Date(paste0(y,"/",m,"/",1)))
    tmp[i] <- ISOdate(y,m,ndays) # as.POSIXlt(paste0(y,"-",m,"-",ndays, " 00:00:00"),format="%Y-%m-%d %H:%M:%S", tz = "GMT")
    i <- i + 1
    m <- m+1
    if (m == 13){
      m <- 1
      y <- y +1
    }

    if (y == ey && m > em)
      break
  }
  return(tmp)
}

iwfm.interpTime <- function(simTS, tm){
  tmp <- which(simTS <= tm)
  if (length(tmp) == 0){
    return(list("T1" = -10000, "T2" = -10000, "t" = 0))
  }
  else{
    ta <- max(tmp)
  }

  tmp <- which(simTS >= tm)
  if (length(tmp) == 0){
    return(list("T1" = 10000, "T2" = 10000, "t" = 0))
  }
  else{
    tb <- min(tmp)
  }

  if (ta == tb){
    if (ta == 1){
      tb <- 2
    }
    else{
      ta <- ta - 1
    }
  }

  t <- as.double(tm - simTS[ta])/as.double(simTS[tb] - simTS[ta])
  return(list("T1" = ta, "T2" = tb, "t" = t))
}

# Calculate Homographic transformation Coefficients
iwfm.calc_MSH_HTC <- function(XY, MSH){
  HTCF <- matrix(data = NA, nrow = dim(MSH)[1], ncol = 8)

  for (i in 1:dim(MSH)[1]) {
    x <- vector(mode= "numeric", length = 4L)
    y <- vector(mode= "numeric", length = 4L)
    x[1:3] <- XY[,1][MSH[i,1:3]]
    y[1:3] <- XY[,2][MSH[i,1:3]]
    if (MSH[i,4] == 0){
      # This is a triangle and we make it to dummy quad by projecting the first point
      # across the line of 2dn and 3rd points
      xm <- (x[2]+x[3])/2
      ym <- (y[2]+y[3])/2
      x[4] <- x[3]
      y[4] <- y[3]
      x[3] <- xm*2+x[1]*(1-2)
      y[3] <- ym*2+y[1]*(1-2)
    }else{
      x[4] <- XY[,1][MSH[i,4]]
      y[4] <- XY[,2][MSH[i,4]]
    }
    HTCF[i,] <- iwfm.calcHTC(x,y)
  }

  return(HTCF)
}

iwfm.calcHTC <- function(x,y){
  xt <- x/100000
  yt <- y/1000000
  XL <- matrix(data = 0, nrow = 8, ncol = 8)
  xd <- c(0,  1, 1, 0)
  yd <- c(0,  0, 1, 1)
  z <- c(1,1,1,1)
  Xmat <- matrix(data = c(xt,yt,z), nrow=4)
  Xdmat <- matrix(data = c(xd,yd,z), nrow=4)

  XL[1:4,1:3] <- Xmat
  XL[5:8,4:6] <- Xmat
  XL[1:4,7] <- -Xmat[,1]*Xdmat[,1]
  XL[5:8,7] <- -Xmat[,1]*Xdmat[,2]
  XL[1:4,8] <- -Xmat[,2]*Xdmat[,1]
  XL[5:8,8] <- -Xmat[,2]*Xdmat[,2]

  cc <- solve(XL,c(xd,yd))
  return(cc)
}

iwfm.UnitCoords <- function(x,y,cf){
  x <- x/100000
  y <- y/1000000
  xp <- (cf[1]*x+cf[2]*y+cf[3])/(cf[7]*x+cf[8]*y+1)
  yp <- (cf[4]*x+cf[5]*y+cf[6])/(cf[7]*x+cf[8]*y+1)
  return(c(xp,yp))
}

iwfm.CoordFromUnit <- function(xu, yu, cf){
  #      1 2 3 4 5 6 7 8
  # syms a b c d e f g h X Y x y
  # A=[X*g-a X*h-b;Y*g-d Y*h -e]
  # B=[c-X;f-Y]
  # XX = linsolve(A,B)
  # -(Y*b - X*e - b*f + c*e - Y*c*h + X*f*h)/(a*e - b*d - Y*a*h + Y*b*g + X*d*h - X*e*g)
  # (Y*a - X*d - a*f + c*d - Y*c*g + X*f*g)/(a*e - b*d - Y*a*h + Y*b*g + X*d*h - X*e*g)
  
  x <- -(yu*cf[2] - xu*cf[5] - cf[2]*cf[6] + cf[3]*cf[5] - yu*cf[3]*cf[8] + xu*cf[6]*cf[8])/(cf[1]*cf[5] - cf[2]*cf[4] - yu*cf[1]*cf[8] + yu*cf[2]*cf[7] + xu*cf[4]*cf[8] - xu*cf[5]*cf[7])
  y <- (yu*cf[1] - xu*cf[4] - cf[1]*cf[6] + cf[3]*cf[4] - yu*cf[3]*cf[7] + xu*cf[6]*cf[7])/(cf[1]*cf[5] - cf[2]*cf[4] - yu*cf[1]*cf[8] + yu*cf[2]*cf[7] + xu*cf[4]*cf[8] - xu*cf[5]*cf[7])
  x <- 100000*x
  y <- 1000000*y
  return(c(x, y))
}

iwfm.Quadinterp <- function(xu, yu, v){
  # Have to check the interpolation that the shape functions correspond to the right sides.
  # Especially for the triangles
  if (length(v) == 4){
    vx <- xu*v[2] + (1-xu)*v[4]
    vy <- (1-yu)*v[1] + yu*v[3]
    return(c(vx,vy))
  } else if (length(v) == 3){
    vx <- xu*v[1] + sqrt(2)*xu*v[2] + (xu-1)*v[3]
    vy <- (yu-1)*v[1] + sqrt(2)*yu*v[2] + yu*v[3]
    return(c(vx,vy))
  }
}

iwfm.findNextpointRBF <- function(p,v, t, D, po, dist_step, simTime){
  
  tm_step <- dist_step/sqrt(sum(v^2))
  
  d <- sign(po$Dir)
  
  k1 <- d*v
  
  # Calculate p2 by taking 1/4 step using the velocity K1
  p2 <- p + (1/4) * tm_step * k1
  t2 <- t +  d*seconds(24*(1/4)*tm_step*3600)
  k2 <- d*iwfm.pntVelRBF(p2,t2,D,simTime,po)
  
  # Calculate p3 by taking 3/8 step using a linear combination of velocities k1 and k2
  p3 <- p + (3/8)*tm_step * ((3/32)*k1 + (9/32)*k2)
  t3 <- t + d*seconds(24*(3/8)*tm_step*3600)
  k3 <- d*iwfm.pntVelRBF(p3,t3,D,simTime,po)
  
  # Calculate p4 by taking 12/13 step using a linear combination of k1 - k3
  p4 <- p + (12/13)*tm_step * ((1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3)
  t4 <- t + d*seconds(24*(12/13)*tm_step*3600)
  k4 <- d*iwfm.pntVelRBF(p4,t4,D,simTime,po)
  
  # Calculate p5 by taking a full step using a linear combination of k1 - k4
  p5 <- p + (1)*tm_step * ((439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4)
  t5 <- t + d*seconds(24*(1)*tm_step*3600)
  k5 <- d*iwfm.pntVelRBF(p5,t5,D,simTime,po)
  
  # Calculate p6 by taking a half step using a linear combination of k1 - k4
  p6 <- p + (1/2)*tm_step * ( -(8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5)
  t6 <- t + d*seconds(24*(1/2)*tm_step*3600)
  k6 <- d*iwfm.pntVelRBF(p6,t6,D,simTime,po)
  
  # Calculate two approximations of the next point using different combinations of velocities
  yn <- p + tm_step * ( (25/216)*k1 + (1408/2565)*k3 + (2197/4101)*k4 - (1/5)*k5)
  zn <- p + tm_step * ( (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6)
  
  
  
  R <- (1/tm_step)*sqrt(sum((zn-yn)^2))
  
  if (is.nan(R)){
    return(list("p" = NA, "Status" = "Repeat", "sh" = 0.5 ))
  }
  
  delta <- 0.84*(po$Tolerance/R)^(1/4)
  if (R < po$Tolerance){
    return(list("p" = zn, "Status" = "OK", "sh" = delta*po$Tolerance ))
  }
  else{
    return(list("p" = NA, "Status" = "Repeat", "sh" = delta*tol ))
  }
}


iwfm.findNextpoint <- function(p, v, t, HV, ZV, FI, BC, MSH, XY, STRAT, HTCF, simTime, tm_step, tol, drct){
  d <- sign(drct)
  HV <- d*HV
  ZV <- d*ZV

  k1 <- v

  # Calculate p2 by taking 1/4 step using the velocity K1
  p2 <- p + (1/4) * tm_step * k1
  t2 <- t +  d*seconds(24*(1/4)*tm_step*3600)
  VV <- iwfm.pntVel(p2, t2, HV, ZV, FI, BC, MSH, XY, STRAT, HTCF, simTime)
  if (VV$reas == "outNO"){
    k2 <- VV$v
    # Calculate p3 by taking 3/8 step using a linear combination of velocities k1 and k2
    p3 <- p + (3/8)*tm_step * ((3/32)*k1 + (9/32)*k2)
    t3 <- t + d*seconds(24*(3/8)*tm_step*3600)
    VV <- iwfm.pntVel(p3, t3, HV, ZV, FI, BC, MSH, XY, STRAT, HTCF, simTime)
    if (VV$reas == "outNO"){
      k3 <- VV$v
      # Calculate p4 by taking 12/13 step using a linear combination of k1 - k3
      p4 <- p + (12/13)*tm_step * ((1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3)
      t4 <- t + d*seconds(24*(12/13)*tm_step*3600)
      VV <- iwfm.pntVel(p4, t4, HV, ZV, FI, BC, MSH, XY, STRAT, HTCF, simTime)
      if (VV$reas == "outNO"){
        k4 <- VV$v
        # Calculate p5 by taking a full step using a linear combination of k1 - k4
        p5 <- p + (1)*tm_step * ((439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4)
        t5 <- t + d*seconds(24*(1)*tm_step*3600)
        VV <- iwfm.pntVel(p5, t5, HV, ZV, FI, BC, MSH, XY, STRAT, HTCF, simTime)
        if (VV$reas == "outNO"){
          k5 <- VV$v
          # Calculate p6 by taking a half step using a linear combination of k1 - k4
          p6 <- p + (1/2)*tm_step * ( -(8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5)
          t6 <- t + d*seconds(24*(1/2)*tm_step*3600)
          VV <- iwfm.pntVel(p6, t6, HV, ZV, FI, BC, MSH, XY, STRAT, HTCF, simTime)
          if (VV$reas == "outNO"){
            k6 <- VV$v
            yn <- p + tm_step * ( (25/216)*k1 + (1408/2565)*k3 + (2197/4101)*k4 - (1/5)*k5)
            zn <- p + tm_step * ( (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6)

            R <- (1/tm_step)*sqrt(sum((zn-yn)^2))
            delta <- 0.84*(tol/R)^(1/4)
            if (R < tol){
              return(list("p" = zn, "reas" = VV$reas, "sh" = delta*tol ))
            }
            else{
              return(list("p" = NA, "reas" = "Repeat", "sh" = delta*tol ))
            }
          }
          else{
            return(list("p" = NA, "reas" = VV$reas))
          }
        }
        else{
          return(list("p" = NA, "reas" = VV$reas))
        }
      }
      else{
        return(list("p" = NA, "reas" = VV$reas))
      }
    }
    else{
      return(list("p" = NA, "reas" = VV$reas))
    }
  }else{
    return(list("p" = NA, "reas" = VV$reas))
  }
}
 

iwfm.wellParticles <- function(wells, npart, rd, MSH, XY, STRAT, BC){

  t <- seq(from = 1, to = npart, by = 1)
  particles <- list("Eid" = NA, "Sid" = NA, "X" = NA, "Y" = NA, "Z" = NA, "TM" = ISOdate(1900,10,1))

  sind <- 1
  eind <- npart

  for (iw in 1:dim(wells)[1]) {
    elw <- iwfm.findElemId(wells$X[iw], wells$Y[iw], BC, MSH, XY)
    xv <- XY[,1][MSH[elw,1:4]]
    yv <- XY[,2][MSH[elw,1:4]]
    zv <- STRAT[,1][MSH[elw,1:4][1:4],]
    gse <- iwfm.idw_interp(wells$X[iw], wells$Y[iw], xv, yv, zv)
    p1 <- c(wells$X[iw], wells$Y[iw], gse - wells$DWT[iw])
    p2 <- c(wells$X[iw], wells$Y[iw], p1[3] - wells$SL[iw])

    xu <- rd*cos(t)
    yu <- rd*sin(t)
    zu <- t*(1/npart)

    xu <- xu + wells$X[iw]
    yu <- yu + wells$Y[iw]
    zu <- zu*(p1[3] - p2[3]) + p2[3]

    particles$Eid[sind:eind] <- iw
    particles$Sid[sind:eind] <- t
    particles$X[sind:eind] <- xu
    particles$Y[sind:eind] <- yu
    particles$Z[sind:eind] <- zu
    particles$TM[sind:eind] <- ISOdate(wells$YY[iw], wells$MM[iw], wells$DD[iw])

    sind <- sind + npart
    eind <- eind + npart
  }

  return(data.frame(particles))
}

iwfm.traceParticles <- function(prt, D, po = iwfm.options(), simTime = iwfm.SimTime()){

  #streamline <- iwfm.calcStreamline(particles[1,3:6], D, po, simTime)
  #streamline <- iwfm.calcStreamlineRBF(particles[1,3:6], D, po, simTime)
  
  cnt <- 1
  for (i in 1:dim(prt)[1]) {
    print(i)
    streamline <- iwfm.calcStreamlineRBF(prt[i,3:6], D, po, simTime, cnt)
    
    if (i == 1){
      strmlns <- data.frame(prt$Eid[i], prt$Sid[i], list("STRM" = streamline))
      names(strmlns) <- c("Eid", "Sid", "CNT", "X", "Y", "Z", "TM", "Exit")
    }else{
      df <- data.frame(prt$Eid[i], prt$Sid[i], list("STRM" = streamline))
      names(df) <- c("Eid", "Sid", "CNT", "X", "Y", "Z", "TM", "Exit")
      strmlns <- rbind(strmlns,df)
    }
    cnt <- strmlns$CNT[length(strmlns$CNT)] + 1
  }
  
  return(strmlns)
}


iwfm.checkPoint <- function(p, BC, MSH, XY, STRAT){
  elid <- iwfm.findElemId(p[1], p[2], BC, MSH, XY)
  if (elid < 0){
    return(list("Status" = "outDomain"))
  }
  else{
    lay <- iwfm.findLayer(p[1], p[2], p[3] ,elid, MSH, XY, STRAT)
    if (lay$Lay == 0){
      return(list("Status" = "outTop"))
    } else if (lay$Lay == -9){
      return(list("Status" = "outBot"))
    }
    else{
      return( list("Status" = "IN", "ElId" = elid, "Lay" = lay) )
    }
  }
}

iwfm.calcStreamlineRBF <- function(prt, D, po, simTime, cnt ){
  iter <- 1
  pos <- 1
  timeStep <- po$Step
  
  stepIncrease <- TRUE
  prt <- data.frame(cnt, as.numeric(prt[1,1]), as.numeric(prt[1,2]), as.numeric(prt[1,3]), prt$TM[pos])
  names(prt) <- c("CNT", "X", "Y", "Z", "TM")
  
  chck <- iwfm.checkPoint(c(as.numeric(prt[1,2]), as.numeric(prt[1,3]), as.numeric(prt[1,4])), D$BC, D$MSH, D$XY, D$STRAT)
  if (chck$Status != "IN")
    return(list("STRMLN" = prt, "Exit" = chck$Status))
  
  
  while (iter < po$MaxIter){
    #print(iter)
    p <- c(as.numeric(prt[pos,2]), as.numeric(prt[pos,3]), as.numeric(prt[pos,4]))
    TM <- iwfm.interpTime(simTime, prt[pos,4])
    if ((TM$T1 == -10000 || TM$T1 == 10000) && length(dim(D$HFLOW)) == 3){
      return(list("STRMLN" = prt, "Exit" = "outTime"))
    }

    v <- iwfm.pntVelRBF(p, TM, D, simTime, po)
    PN <- iwfm.findNextpointRBF(p,v, prt[pos,4], D, po, timeStep, simTime)
    
    if (PN$Status == "Repeat"){
      timeStep <- timeStep*PN$sh
    }else{
      chck <- iwfm.checkPoint(PN$p, D$BC, D$MSH, D$XY, D$STRAT)
      if (chck$Status == "IN" ){
        cnt <- cnt + 1
        df <- data.frame(cnt, PN$p[1], PN$p[2], PN$p[3], prt$TM[pos] + sign(po$Dir)*seconds(24*timeStep*3600))
        names(df) <- c("CNT", "X", "Y", "Z", "TM")
        prt <- rbind(prt,df)
        pos <- pos +1
        if (stepIncrease)
          timeStep <- min(po$MaxStep, timeStep*PN$sh)
      }
      else{
        if (chck$Status == "outTop"){
          if (timeStep < po$MinStep){
            cnt <- cnt + 1
            df <- data.frame(cnt, p[1], p[2], p[3], prt$TM[pos] + sign(po$Dir)*seconds(24*timeStep*3600))
            names(df) <- c("CNT", "X", "Y", "Z", "TM")
            prt <- rbind(prt,df)
            return(list("STRMLN" = prt, "Exit" = chck$Status))
          }
          else{
            stepIncrease <- FALSE
            timeStep <- timeStep/2
          }
        }
      }
    }
    iter <- iter + 1
  }
  return(list("STRMLN" = prt, "Exit" = "outIter"))
}




iwfm.calcStreamline <- function(prt, D, po, simTime ){
  iter <- 1
  pos <- 1
  timeStep <- po$Step


  while (iter < po$MaxIter){
    p <- c(as.numeric(prt[pos,1]), as.numeric(prt[pos,2]), as.numeric(prt[pos,3]))
    VV <- iwfm.pntVel(p, prt[pos,4], D$HFLOW, D$VFLOW,D$FI, D$BC, D$MSH, D$XY, D$STRAT, D$HTCF, simTime)
    if (VV$reas == "outNO"){
      PN <- iwfm.findNextpoint(p, VV$v, prt[pos,4], D$HFLOW, D$VFLOW, D$FI, D$BC, D$MSH, D$XY, D$STRAT, D$HTCF, simTime, timeStep, po$Tolerance, po$Dir)
      if(PN$reas == "Repeat"){

        timeStep <- timeStep*PN$sh

      } else if (PN$reas == "outNO"){

        df <- data.frame(PN$p[1], PN$p[2], PN$p[3], prt$TM[pos] + sign(po$Dir)*seconds(24*timeStep*3600))
        names(df) <- c("X", "Y", "Z", "TM")
        prt <- rbind(prt,df)
        pos <- pos +1

        timeStep <- min(po$MaxStep, timeStep*PN$sh)

      } else{
        if (PN$reas != "outTime"){
          timeStep <- timeStep/2
          if (timeStep < po$MinStep){
            return(list("STRMLN" = prt, "reas" = PN$reas))
          }
        }
        return(list("STRMLN" = prt, "reas" = PN$reas))
      }
    }
    else{
       return(list("STRMLN" = prt, "reas" = VV$reas))
    }
    iter <- iter + 1
  }
  return(list("STRMLN" = prt, "reas" = "outIter"))
}

iwfm.options <- function(){
  return(
    list(
         "Dir" = -1,
         "MaxIter" = 1000,
         "Step" = 10,
         "MaxStep" = 500,
         "MinStep" = 1,
         "Tolerance" = 1,
         "RBFSTD" = list("XY" = 5000, "Z" = 100)
         )
    )
}

iwfm.AverageVelField <- function(startDate, endDate, D, simTime){
  st <- iwfm.interpTime(simTime, startDate)
  et <- iwfm.interpTime(simTime, endDate)
  HF <- matrix(data = NA, nrow = dim(D$HFLOW)[1], ncol = 4)
  VF <- matrix(data = NA, nrow = dim(D$VFLOW)[1], ncol = 4)
  for (i in 1:4) {
    HF[,i] <- apply(D$HFLOW[,st$T2:et$T2,i],1,mean)
    VF[,i] <- apply(D$VFLOW[,st$T2:et$T2,i],1,mean)
  }
  D$HFLOW <- HF
  D$VFLOW <- VF
  return(D)
}

iwfm.sourceOutline <- function(wellid, strmln){
  sid <- unique(strmln$Sid[ which(strmln$Eid == wellid)])
  df <- data.frame()
  for (i in 1:length(sid)){
    iend <- max(which(strmln$Eid == wellid & strmln$Sid == sid[i]))
    tmp <- data.frame(strmln[iend,4:6])
    df <- rbind(df, tmp)
  }
  return(df)
  
}


iwfm.writeData<-function(filename, A){
  if (is.null(dim(A))){
    nr <- length(A)
    nc <- 1
  }
  else{
    nr <- dim(A)[1]
    nc <- dim(A)[2]
  }
  
  fprintf("%d %d\n", nr, nc, file = filename, append = FALSE)
  
  write(formatC(t(A), format = "f", digits = 11, drop0trailing = TRUE), #format(t(A), nsmall = nsml, digits = 15, justify = "left", scientific = FALSE), 
        file = filename, ncolumns = nc, append = TRUE)
}