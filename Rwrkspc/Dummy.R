source("iwfm_funct.R")
simTime <- iwfm.SimTime()
D <- iwfm.loadData("PartTrackData.h5")

xp <- 829282
yp <- 3992499
zp <- -100
tp <-ISOdate(2010,5,5)

Ttp <- iwfm.interpTime(simTime, tp) 
elid <- iwfm.findElemId(xp, yp, D$BC, D$MSH, D$XY)
lay <- iwfm.findLayer(xp, yp, zp ,elid, D$MSH, D$XY, D$STRAT)

tmp <- iwfm.calcVel(xp, yp, zp, Ttp, elid, lay, D$HFLOW, D$VFLOW)

D$HFLOW[abs(D$FI[elid][1,]),Ttp$T1,lay]

