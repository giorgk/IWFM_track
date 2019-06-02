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
  
  if (zp > zlay[1])
    return(0)
  if (zp < zlay[length(zlay)])
    return (-9)
  
  for (i in 1:length(zlay)-1) {
    if (zp < zlay[i] && zp > zlay[i+1])
      return(i)
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
  
  # convert the flow data to arrays for faster indexing
  tmpH <- array(dim = dim(HFLOW))
  for (i in 1:dim(tmpH)[3]) {
    tmpH[,,i] <- HFLOW[,,i]
  }
  tmpV <- array(dim = dim(VFLOW))
  for (i in 1:dim(tmpV)[3]) {
    tmpV[,,i] <- VFLOW[,,i]
  }
  
  output <- list("XY" = XY, "MSH" = MSH, "BC" = BC, "STRAT" = STRAT, "FI" = FI, "HFLOW" = tmpH, "VFLOW" = tmpV)
  h5close(hf)
  return(output)
}

iwfm.idw_interp <- function(x,y, xv, yv, v, thres = 0.001){
  dst <- sqrt((xv - x)^2 + (yv - y)^2)
  id <- which(dst <= thres)
  w <- vector(mode = "numeric", length = length(xv))
  if (length(id) == 1){
    w[id] <- 1
  }else{
    w <- dst/sum(dst)
  }
  
  vRet <- vector(mode = "numeric", length = dim(v)[2])
  for (i in 1:dim(v)[2]) {
    vRet[i] <- sum(w*v[,i])
  }
  return(vRet)
}

iwfm.calcVel <- function(x, y, z, t, elid, lay, hVel, vVel){
  
  return(0)  
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

iwfm.interpTime <- function(simTS, t){
  ta <- max(which(simTS < t))
  tb <- min(which(simTS > t))
  t <- as.double(t - simTS[ta])/as.double(simTS[tb] - simTS[ta])
  output <- list("T1" = ta, "T2" = tb, "t" = t)
}
