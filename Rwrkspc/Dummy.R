setwd("f:/UCDAVIS/IWFM_track/Rwrkspc")
source("iwfm_funct.R")
simTime <- iwfm.SimTime()
D <- iwfm.loadData("PartTrackData.h5")
#     #in triangle# in quad  in first element
xp <- 750633 #835242 #553884.057971
yp <- 4077362 #3958831 #4497546.534653
zp <- -30
tp <- ISOdate(2010,5,5)

p <- c(xp,yp,zp)

Ttp <- iwfm.interpTime(simTime, tp)
VV <- iwfm.calcVxyzRBF(p, Ttp, D$HFLOW, D$VFLOW, D$NRML, D$FACEZ)


VV <- iwfm.pntVel(p, tp, D$HFLOW, D$VFLOW, D$FI, D$BC, D$MSH, D$XY, D$STRAT, D$HTCF, simTime)

PN <- iwfm.findNextpoint(p, VV$v, tp, D$HFLOW, D$VFLOW, D$FI, D$BC, D$MSH, D$XY, D$STRAT, D$HTCF, simTime, 1*0.4352258*0.5, 1)
Ttp <- iwfm.interpTime(simTime, tp)
elid <- iwfm.findElemId(xp, yp, D$BC, D$MSH, D$XY)
lay <- iwfm.findLayer(xp, yp, zp ,elid, D$MSH, D$XY, D$STRAT)

# Transform the point to the reference quad 
punit <- iwfm.UnitCoords(xp,yp,D$HTCF[elid,])



vtmp <- iwfm.calcVel(punit, Ttp, elid, lay, D$HFLOW, D$VFLOW, D$FI)/30
vmag <- sqrt(sum(vtmp^2))

#Runge-Kutta-Fehlberg Method (RKF45)
# http://mathfaculty.fullerton.edu/mathews/n2003/RungeKuttaFehlbergMod.html
# https://ece.uwaterloo.ca/~dwharder/NumericalAnalysis/14IVPs/rkf45/complete.html
# http://maths.cnam.fr/IMG/pdf/RungeKuttaFehlbergProof.pdf
# https://math.okstate.edu/people/yqwang/teaching/math4513_fall11/Notes/rungekutta.pdf
step_m <- 50 # Step size in meters

# based on the velocity calculate the time needed to travel the step size meters
t_m <- step_m/vmag
pinit <- c(xp,yp,zp)

k1 <- vtmp # The velocity at the previous point

# Calculate p2 by taking 1/4 step using the velocity K1
p2 <- pinit + (1/4)*t_m*k1

Ttp2 <- iwfm.interpTime(simTime, tp+seconds(24*(1/4)*t_m*3600))
elid2 <- iwfm.findElemId(p2[1], p2[2], D$BC, D$MSH, D$XY)
lay2 <- iwfm.findLayer(p2[1], p2[2], p2[3] ,elid2, D$MSH, D$XY, D$STRAT)
punit2 <- iwfm.UnitCoords(p2[1], p2[2], D$HTCF[elid2,])
k2 <- iwfm.calcVel(punit2, Ttp2, elid2, lay2, D$HFLOW, D$VFLOW, D$FI)/30

# Calculate p3 by taking 3/8 step using a linear combination of velocities k1 and k2
p3 <- pinit + (3/8)*t_m * ((3/32)*k1 + (9/32)*k2)

Ttp3 <- iwfm.interpTime(simTime, tp+seconds(24*(3/8)*t_m*3600))
elid3 <- iwfm.findElemId(p3[1], p3[2], D$BC, D$MSH, D$XY)
lay3 <- iwfm.findLayer(p3[1], p3[2], p3[3] ,elid3, D$MSH, D$XY, D$STRAT)
punit3 <- iwfm.UnitCoords(p3[1], p3[2], D$HTCF[elid3,])
k3 <- iwfm.calcVel(punit3, Ttp3, elid3, lay3, D$HFLOW, D$VFLOW, D$FI)/30

# Calculate p4 by taking 12/13 step using a linear combination of k1 - k3
p4 <- pinit + (12/13)*t_m * ((1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3)

Ttp4 <- iwfm.interpTime(simTime, tp+seconds(24*(12/13)*t_m*3600))
elid4 <- iwfm.findElemId(p4[1], p4[2], D$BC, D$MSH, D$XY)
lay4 <- iwfm.findLayer(p4[1], p4[2], p4[3] ,elid4, D$MSH, D$XY, D$STRAT)
punit4 <- iwfm.UnitCoords(p4[1], p4[2], D$HTCF[elid4,])
k4 <- iwfm.calcVel(punit4, Ttp4, elid4, lay4, D$HFLOW, D$VFLOW, D$FI)/30

# Calculate p5 by taking a full step using a linear combination of k1 - k4
p5 <- pinit + (1)*t_m * ((439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4)

Ttp5 <- iwfm.interpTime(simTime, tp+seconds(24*(1)*t_m*3600))
elid5 <- iwfm.findElemId(p5[1], p5[2], D$BC, D$MSH, D$XY)
lay5 <- iwfm.findLayer(p5[1], p5[2], p5[3] ,elid5, D$MSH, D$XY, D$STRAT)
punit5 <- iwfm.UnitCoords(p5[1], p5[2], D$HTCF[elid5,])
k5 <- iwfm.calcVel(punit5, Ttp5, elid5, lay5, D$HFLOW, D$VFLOW, D$FI)/30

# Calculate p6 by taking a half step using a linear combination of k1 - k4
p6 <- pinit + (1/2)*t_m * ( -(8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5)

Ttp6 <- iwfm.interpTime(simTime, tp+seconds(24*(1/2)*t_m*3600))
elid6 <- iwfm.findElemId(p6[1], p6[2], D$BC, D$MSH, D$XY)
lay6 <- iwfm.findLayer(p6[1], p6[2], p6[3] ,elid6, D$MSH, D$XY, D$STRAT)
punit6 <- iwfm.UnitCoords(p6[1], p6[2], D$HTCF[elid6,])
k6 <- iwfm.calcVel(punit6, Ttp6, elid6, lay6, D$HFLOW, D$VFLOW, D$FI)/30

# Calculate one approximation
yn <- pinit + t_m * ( (25/216)*k1 + (1408/2565)*k3 + (2197/4101)*k4 - (1/5)*k5)
zn <- pinit + t_m * ( (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6)

iwfm.findNextpoint(p, v, t, HV, ZV, FI, BC, MSH, XY, STRAT, HTCF, simTime, tm_step, tol)


# Commands to run particle tracking
wells <- read.table(file = "welldata.dat", header = TRUE)
particles <- iwfm.wellParticles(wells, 80, 5, D$MSH, D$XY, D$STRAT, D$BC)
po <- iwfm.options()

# Use average velocity
D <- iwfm.AverageVelField(ISOdate(2005,10,31), ISOdate(2015,9,30), D, simTime)

iwfm.traceParticles(particles, D, po)


p1 <- c(566695.358649, 4413576.717052, 46.243613)
V1 <- iwfm.pntVel(p1, tp, D$HFLOW, D$VFLOW, D$FI, D$BC, D$MSH, D$XY, D$STRAT, D$HTCF, simTime)
p2 <- c(566695.808515, 4413575.957847, 46.240942)
V2 <- iwfm.pntVel(p2, tp, D$HFLOW, D$VFLOW, D$FI, D$BC, D$MSH, D$XY, D$STRAT, D$HTCF, simTime)
