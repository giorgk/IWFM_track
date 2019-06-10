library(h5)
library(pracma)

setwd("f:/UCDAVIS/IWFM_track/Rwrkspc")
source("iwfm_funct.R")

# ----------USER input Section---------------
# c2vsim_path is the main folder where the subfolders such as Preproces, Simulation exist
c2vsim_path <-"f:/UCDAVIS/C2VSIM_FG_OR/C2Vsim_FG_v2/C2VSimFG-BETA_PublicRelease/"

# This is the path where the *hdf output files of the C2Vsim simulation are.
# By default is under the c2vsim_path/Results but that's not always the case so set it here explicitly
results_path <- paste0(c2vsim_path, "Results/")

# Output file name with path. Do not add any extension to this 
output_file <- "f:/UCDAVIS/IWFM_track/Rwrkspc/PartTrackData" 

# That's all as far as inputs is concerned. 
# However there are some values that their default values may not work:
ND_linesSkip <- 90 # Number of comment lines in the Preprocessor/C2VSimFG_Nodes.dat file.
ND_nRead <-30179 # Number of lines to read after ND_linesSkip. This is equal to the number of nodes.
Strat_linesSkip <- 105 # Number of comment lines in the Preprocessor/C2VSimFG_Stratigraphy.dat file.
MSH_linesSkip <- 142 # Number of comment lines in the Preprocessor/C2VSimFG_Elements.dat file
MSH_nRead <-32537 # Number of lines to read after MSH_linesSkip This is equal to the number of elements.
# Convertion factor. 
# According to Can the Flow units in the hdf files are in ft^3. we have to multiply by
Flow_CNVRT = 0.0283168 # to convert them to m^3

### Read Ascii files
## Read coordinate file
print("Reading coordinates...")
XY <- read.table(file = paste0(c2vsim_path, "Preprocessor/C2VSimFG_Nodes.dat"),
                 header = FALSE, sep = "", skip = ND_linesSkip, nrows = ND_nRead,
                 quote = "",fill = TRUE,
                 col.names = c("ID", "X", "Y"))

## Read stratigraphy file
print("Reading Stratigraphy...")
strat <- read.table(file = paste0(c2vsim_path, "Preprocessor/C2VSimFG_Stratigraphy.dat"),
                    header = FALSE, sep = "", skip = Strat_linesSkip, nrows = ND_nRead,
                    quote = "",fill = TRUE,
                    col.names = c("ID", "GSE", "A1", "L1", "A2", "L2", "A3", "L3", "A4", "L4"))
# Convert feet to meter
strat[-1] <- strat[-1] * 0.3048
# For each of the 4 layers add the thickness of aquiclude to the layer
for(i in 1:4){
  strat[[ paste0("L", as.character(i))]] <- strat[[ paste0("L", as.character(i))]] + strat[[ paste0("A", as.character(i))]]
}
#Delete the aquicludes and Keep only the layer thicknesses
strat$A1 <- NULL
strat$A2 <- NULL
strat$A3 <- NULL
strat$A4 <- NULL
# Subtract the thickness from the surface to obtain the elevation of the bottom of each layer
strat$L1 = strat$GSE - strat$L1 
strat$L2 = strat$L1 - strat$L2 
strat$L3 = strat$L2 - strat$L3
strat$L4 = strat$L3 - strat$L4

## Read Mesh file
print("Reading Mesh elements...")
MSH <- read.table(file = paste0(c2vsim_path, "Preprocessor/C2VSimFG_Elements.dat"),
                  header = FALSE, sep = "", skip = MSH_linesSkip, nrows = MSH_nRead,
                  quote = "",fill = TRUE,
                  col.names = c("ID", "nd1", "nd2", "nd3", "nd4", "S"))

print("Calculate element areas...")
elemArea <- vector(mode = "numeric", length = length(MSH[[1]]))
bcElem <- matrix(data = 0, nrow = length(MSH[[1]]), ncol = 2)
bcZ <- matrix(data = 0, nrow = length(MSH[[1]]), ncol = 5)
for (i in 1:length(MSH[[1]])){
  x <- XY$X[as.integer(c(MSH[i,2:5]))]
  y <- XY$Y[as.integer(c(MSH[i,2:5]))]
  bcElem[i,] <- c(mean(x), mean(y))
  elemArea[i] <- polyarea(x,y)
  bcZ[i,] <- as.vector(apply(strat[as.integer(c(MSH[i,2:5])),2:6],2,mean))
}

# Calculate the face areas

### Read the HDF5 budget file
GW_BDGinfo <- h5file(name =  paste0(results_path, "C2VSimFG_GW_ZBudget.hdf"), mode = "r")
hfileGroups <- list.groups(GW_BDGinfo)
hfileDataSets <- list.datasets(GW_BDGinfo)
faceElem <- GW_BDGinfo[hfileDataSets[17]]

# initialize a list of matrices to hold the face areas for each layer
faceArea = matrix(data = 0, nrow = dim(faceElem)[1], ncol = 4)
faceIndex = matrix(data = 0L, nrow = dim(MSH)[1], ncol = 4)
faceZ = matrix(data = 0, nrow = dim(faceElem)[1], ncol = 4)

## Calculate the face areas and indices for the inner faces of the Mesh
print("Calculate inner element face areas and indices...")
for (i in 1:dim(faceElem)[1]){
  ela <- faceElem[i][1] # This is the element index of the first column
  if (ela != 0 ){ #If the element index is not zero then find out how many sides this element has
    na <- 4 # number of element sides 
    if (MSH[ela,5] == 0)
      na <- 3
  }
  
  elb <- faceElem[i][2] # Repeat for the element of the second column
  if (elb != 0 ){
    nb <- 4
    if (MSH[elb,5] == 0)
      nb <- 3
  }
  
  if (ela != 0 && elb != 0){
    # find the vertex indices for the common face between the two elements
    breakThis <- FALSE
    for (ii in 1:na){
      if (faceIndex[ela,ii] != 0)
        next
      # find the node indices of this face of the first element
      ia <- MSH[ela,ii+1]
      iii <- ii + 1
      if (ii == na)
        iii = 1
      ib <- MSH[ela, iii+1]
      
      # Loop through the faces of the second element
      for (jj in 1:nb){
        if (faceIndex[elb,jj] != 0)
          next
        ja <- MSH[elb,jj+1]
        jjj <- jj +1
        if (jj == nb)
          jjj <- 1
        jb <- MSH[elb,jjj+1]
        
        if ((ia == ja & ib == jb) | (ia == jb & ib == ja)){
          faceIndex[ela,ii] = i
          faceIndex[elb,jj] = -i
          # calculate the area of the face for each layer
          L <- sqrt((XY$X[ia] - XY$X[ib])^2 + (XY$Y[ia] - XY$Y[ib])^2)
          xv <- c(0, L, L, 0)
          for (k in 1:4){
            yv <- c(strat[ia,k+2], strat[ib,k+2], strat[ib,k+1], strat[ia,k+1])
            faceArea[i,k] <- polyarea(xv,yv)
            faceZ[i,k] <- mean(yv)
          }
          breakThis <- TRUE
          break
        }
      }
      if (breakThis)
        break
    }
  }
}

FcLm = matrix(nrow = dim(faceElem)[1], ncol = dim(faceElem)[2])
for (i in 1:dim(FcLm)[1]){
  FcLm[i,1] <- faceElem[i][1]
  FcLm[i,2] <- faceElem[i][2]
}

print("Calculate outer element face areas and indices...")
## Calculate the face areas and indices for the boundary faces of the Mesh
for (i in 1:dim(faceIndex)[1]){
  # find how many faces do not have index
  zrID <- which(faceIndex[i,] == 0)
  if (MSH[i, 5] == 0){
    nonValid <- which(zrID == 4)
    zrID <- zrID[-nonValid]
  }
  if (!isempty(zrID)){
    # Find the faces that are associated with the element
    elemFCid <- which(FcLm == i,arr.ind = TRUE)
    # remove the faces that they are connected with another element
    keep <- vector(mode = "integer")
    for (j in 1:dim(elemFCid)[1]){
      oth_el <- 2
      if (elemFCid[j,2] == 2)
        oth_el <- 1
      if (FcLm[elemFCid[j,1], oth_el] == 0)
        keep <- append(keep, j)
    }
    elemFCid <- elemFCid[keep,]
    
    if (is.null(dim(elemFCid)[1])){
      if (!(length(elemFCid) == 2 & length(zrID) == 1)){
        print(paste0("Mismatch between faces for element: ",i) )
      }
      else{
        if (elemFCid[2] == 1){
          faceIndex[i,zrID[1]] = elemFCid[1]
        }else{
          faceIndex[i,zrID[1]] = -elemFCid[1]
        }
        ia = MSH[i,zrID[1]+1]
        ib = MSH[i,zrID[1]+2]
        if (ib == 0 | zrID[1]+1 == 5){
          ib = MSH[i,2]
        }
        L <- sqrt((XY$X[ia] - XY$X[ib])^2 + (XY$Y[ia] - XY$Y[ib])^2)
        xv <- c(0, L, L, 0)
        for (k in 1:4){
          yv <- c(strat[ia,k+2], strat[ib,k+2], strat[ib,k+1], strat[ia,k+1])
          faceArea[elemFCid[1], k] <- polyarea(xv,yv)
          faceZ[elemFCid[1], k] <- mean(yv)
        }
      }
    } else{
      if (dim(elemFCid)[1] != length(zrID)){
        print(paste0("Mismatch between faces for element: ",i) )
      }else{
        for (j in 1:length(zrID)){
          if (elemFCid[j,2] == 1){
            faceIndex[i,zrID[j]] = elemFCid[j,1]
          }else{
            faceIndex[i,zrID[j]] = -elemFCid[j,1]
          }
          ia = MSH[i,zrID[j]+1]
          ib = MSH[i,zrID[j]+2]
          if (ib == 0 | zrID[j]+1 == 5){
            ib = MSH[i,2]
          }
          L <- sqrt((XY$X[ia] - XY$X[ib])^2 + (XY$Y[ia] - XY$Y[ib])^2)
          xv <- c(0, L, L, 0)
          for (k in 1:4){
            yv <- c(strat[ia,k+2], strat[ib,k+2], strat[ib,k+1], strat[ia,k+1])
            faceArea[elemFCid[j,1], k] <- polyarea(xv,yv)
            faceZ[elemFCid[j,1], k] <- mean(yv)
          }
        }
      }
    }
  }
}

## Read Horizontal flow
hflow <- array(dim = c(dim(FcLm)[1], 504, 4))
ids <- c(28, 55, 82, 109) 
for (i in 1:length(ids)) {
  print(paste("reading Horizontal flows from: ", hfileDataSets[ids[i]]))
  tempfflow <- GW_BDGinfo[hfileDataSets[ids[i]]]
  for (j in 1:dim(tempfflow)[1]) {
    hflow[,j,i] <- tempfflow[j,]*Flow_CNVRT
  }
}

## Read deep percolation for the layer 1 vertical flow
# Positive values means downward movement therefore we negate the values as
# in the particle tracking negative z ccorresponds to downward movement
print(paste("reading Vertical flow for layer 1 from: ", hfileDataSets[24]))
DeepPerc <- GW_BDGinfo[hfileDataSets[24]]
Vflow <- array(data = 0, dim = c(dim(DeepPerc)[2], dim(DeepPerc)[1], 5))
for (i in 1:dim(DeepPerc)[1]) {
  Vflow[,i,1] <- -DeepPerc[i,]*Flow_CNVRT
}

# The vertical flows are written per node.
# Read the vertical flows
ids <- c(46, 73, 100)
vertflowNodes <- array(dim = c(dim(XY)[1], 504, length(ids)))
for (i in 1:length(ids)) {
  print(paste("reading Vertical flows from: ", hfileDataSets[ids[i]]))
  Vertflow <- GW_BDGinfo[hfileDataSets[ids[i]]]
  for (j in 1:dim(Vertflow)[1]) {
    vertflowNodes[,j,i] <- Vertflow[j,]
  }
}

# First we have to find out how many element share each node
NsharedElem <- vector(mode = "integer", length = dim(XY)[1])
for(i in 1:dim(XY)[1]){
  elemlist <- which(MSH[c(-1,-6)] == i, arr.ind = TRUE)
  NsharedElem[i] <- dim(elemlist)[1]
}

# For each element we will add the vertical flows of the nodes divided by the number of elements each node is connected
for (i in 1:dim(MSH)[1]){
  velemflow <- matrix(data = 0, nrow = 504, ncol = 3)
  for(j in 2:5){
    if (MSH[i,j] == 0)
      break
    velemflow <- velemflow + vertflowNodes[MSH[i,j],,]/as.numeric(NsharedElem[MSH[i,j]])
  }
  for(j in 1:3){
    Vflow[i,,j+1] <- -velemflow[,j]*Flow_CNVRT
  }
}

# Divide the vertical flows with the element areas
for (i in 1:dim(Vflow)[2]) {
  for (j in 1:dim(Vflow)[3]) {
    Vflow[,i,j] <- Vflow[,i,j]/elemArea
  }
}

# Divide the horizontal flows with the face areas
for (i in 1:dim(hflow)[2]) {
  for (j in 1:dim(hflow)[3]) {
    hflow[,i,j] <- hflow[,i,j]/faceArea[,j]
  }
}

# find the point where the velocity is defined and the direction
print("Calculate flow normals ...")
NRML <- matrix(nrow = dim(FcLm)[1], ncol = 4)
for( i in 1:dim(FcLm)[1]){
  el1 <- FcLm[i,1] # The first element of the face
  el2 <- FcLm[i,2] # The second element of the face
  if (el1 == 0){ # if the first element is zero switch the indices
    el1 <- el2
    el2 <- 0
  }
  
  fcind <- which(abs(faceIndex[el1,])==i)
  nd1 <- fcind
  nd2 <- fcind + 1
  if ( nd2 == 5 || (nd2 == 4 && MSH$nd4[el1] == 0))
    nd2 <- 1
  
  p1 <- as.numeric((XY[MSH[el1,nd1+1],2:3]))
  p2 <- as.numeric((XY[MSH[el1,nd2+1],2:3]))
  pm <- (p1 + p2)/2
  
  # Find a 100 m offset from the line
  off <- (100/sqrt(sum( (p1[2] - p2[2])^2 + (p2[1] - p1[1])^2 )))*c(p1[2] - p2[2], p2[1] - p1[1])
  
  if (faceIndex[el1,fcind] > 0){
    pn <- pm - off
  }
  else{
    pn <- pm + off
  }
  # the point pn points to the face flow direction and the vector pm - pn is perpendicular to line p1-p2
  nrm <- (pn - pm)/100
  NRML[i,] <- c(pm,nrm)
}


## ------ Write outputs --------
## Convert dataframes to matrices
XYmat <- matrix(nrow = dim(XY)[1], ncol = 2)
XYmat[,1] <- XY$X
XYmat[,2] <- XY$Y

MSHmat <- matrix(nrow = dim(MSH)[1], ncol = 5)
MSHmat[,1] <- MSH$nd1
MSHmat[,2] <- MSH$nd2
MSHmat[,3] <- MSH$nd3
MSHmat[,4] <- MSH$nd4
MSHmat[,5] <- MSH$S

STRATmat <- matrix(nrow = dim(strat)[1], ncol = 5)
STRATmat[,1] <- strat$GSE
STRATmat[,2] <- strat$L1
STRATmat[,3] <- strat$L2
STRATmat[,4] <- strat$L3
STRATmat[,5] <- strat$L4

print("Calculate Homographic transformation coefficients...")
HTCF <- iwfm.calc_MSH_HTC(XYmat, MSHmat)

#if (1 == 0){
print(paste("Writing data to ", paste0(output_file, ".h5")))
hf <- h5file(name = paste0(output_file, ".h5"), mode = "w")
hf["flowdata/VFLOW"] <- Vflow
hf["flowdata/HFLOW"] <- hflow
hf["geodata/XY"] <- XYmat
hf["geodata/MSH"] <- MSHmat
hf["geodata/STRAT"] <- STRATmat
hf["geodata/FI"] <- faceIndex
hf["geodata/FCEL"] <- FcLm
hf["geodata/BCEL"] <- bcElem
hf["geodata/HTCF"] <- HTCF
hf["geodata/NRML"] <- NRML
hf["geodata/FACEZ"] <- faceZ
hf["geodata/bcZ"] <- bcZ

h5close(hf)
#}

