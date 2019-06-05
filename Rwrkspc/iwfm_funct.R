library(h5)
library(pracma)
library(lubridate)

#iwfm.calc_velocity(x,y,z)

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
  
  # convert the flow data to arrays for faster indexing
  tmpH <- array(dim = dim(HFLOW))
  for (i in 1:dim(tmpH)[3]) {
    tmpH[,,i] <- HFLOW[,,i]
  }
  tmpV <- array(dim = dim(VFLOW))
  for (i in 1:dim(tmpV)[3]) {
    tmpV[,,i] <- VFLOW[,,i]
  }
  
  output <- list("XY" = XY, "MSH" = MSH, "BC" = BC, "STRAT" = STRAT, "FI" = FI, "HFLOW" = tmpH, "VFLOW" = tmpV, "HTCF" = HTCF )
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

iwfm.calcVel <- function(pu, TM, el, L, HV, ZV, FI){
  # Check if we are inside a quad or triangle
  np <- 4
  if (FI[el,4] == 0){
    np <- 3
  }
  
  # Horizontal Velocity at first time step
  va <- sign(FI[el,1:np]) * HV[abs(FI[el,1:np]), TM$T1, L$Lay]
  # Horizontal Velocity at second time step
  vb <- sign(FI[el,1:np]) * HV[abs(FI[el,1:np]), TM$T2, L$Lay]
  # interpolate between the time steps
  vab <- va*TM$t + vb*(1-TM$t)
  
  # Using the quad shape functions interpolate the xy velocities
  vxy <- iwfm.Quadinterp(pu[1], pu[2], vab)
  
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
  
  return(c(vxy,vz)/30)  
}

iwfm.pntVel <- function(p, t, HV, ZV, FI, BC, MSH, XY, STRAT, HTCF, simTime){
  Ttp <- iwfm.interpTime(simTime, t)
  if (Ttp$T1 == -10000 || Ttp$T1 == 10000){
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
        punit <- iwfm.UnitCoords(p[1], p[2], HTCF[elid,])
        v <- iwfm.calcVel(punit, Ttp, elid, lay, HV, ZV, FI)
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

iwfm.findNextpoint <- function(p, v, t, HV, ZV, FI, BC, MSH, XY, STRAT, HTCF, simTime, tm_step, tol){
  k1 <- v
  
  # Calculate p2 by taking 1/4 step using the velocity K1
  p2 <- p + (1/4) * tm_step * k1
  VV <- iwfm.pntVel(p2, t+seconds(24*(1/4)*tm_step*3600), HV, ZV, FI, BC, MSH, XY, STRAT, HTCF, simTime)
  if (VV$reas == "outNO"){
    k2 <- VV$v
    # Calculate p3 by taking 3/8 step using a linear combination of velocities k1 and k2
    p3 <- p + (3/8)*tm_step * ((3/32)*k1 + (9/32)*k2)
    VV <- iwfm.pntVel(p3, t+seconds(24*(3/8)*tm_step*3600), HV, ZV, FI, BC, MSH, XY, STRAT, HTCF, simTime)
    if (VV$reas == "outNO"){
      k3 <- VV$v
      # Calculate p4 by taking 12/13 step using a linear combination of k1 - k3
      p4 <- p + (12/13)*tm_step * ((1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3)
      VV <- iwfm.pntVel(p4, t+seconds(24*(12/13)*tm_step*3600), HV, ZV, FI, BC, MSH, XY, STRAT, HTCF, simTime)
      if (VV$reas == "outNO"){
        k4 <- VV$v
        # Calculate p5 by taking a full step using a linear combination of k1 - k4
        p5 <- p + (1)*tm_step * ((439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4)
        VV <- iwfm.pntVel(p5, t+seconds(24*(1)*tm_step*3600), HV, ZV, FI, BC, MSH, XY, STRAT, HTCF, simTime)
        if (VV$reas == "outNO"){
          k5 <- VV$v
          # Calculate p6 by taking a half step using a linear combination of k1 - k4
          p6 <- p + (1/2)*tm_step * ( -(8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5)
          VV <- iwfm.pntVel(p6, t+seconds(24*(1/2)*tm_step*3600), HV, ZV, FI, BC, MSH, XY, STRAT, HTCF, simTime)
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
  particles <- list("X" = NA, "Y" = NA, "Z" = NA, "TM" = ISOdate(1900,10,1))
  
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
    
    particles$X[sind:eind] <- xu
    particles$Y[sind:eind] <- yu
    particles$Z[sind:eind] <- zu
    particles$TM[sind:eind] <- ISOdate(wells$YY[iw], wells$MM[iw], wells$DD[iw])
    
    sind <- sind + npart
    eind <- eind + npart
  }
  
  return(particles)
}


iwfm.trackParticle <- function(prt, D, partOpt ){
  iter <- 1
  
  while (iter < partOpt$MaxIter){
    
  }
}

iwfm.options <- function(){
  return(
    list("Step" = 0.2, 
         "Dir" = -1,
         "MaxIter" = 1000
         )
    )
}

  