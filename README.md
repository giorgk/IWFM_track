# IWFM_track
## Introduction
A particle tracking tool for IWFM 

The code considers particle tracking only in groundwater.

It is written in c++ using Visual Studio 2017.
At the moment there is a compiled executable that runs under windows x64.

While there is a working version, the repository is still under development and a bit messy.
The correct source code for the working version is [here](https://github.com/giorgk/IWFM_track/tree/master/VSversion/iwfmTrack/iwfmTrack)

## Building iwfmTrack
The code has very little dependencies.
In fact it requires only boost geometry and boost tree libraries that are header only.
The only library that needs compile is the boost-options library.
The easiest way to prepare the environment is via [vcpck](https://github.com/microsoft/vcpkg)

Follow their instruction to install vcpkg.

Then install the following and you are all set.
```
vcpkg install boost
vcpkg install boost-program-options
```
Do not forget to execute
```
vcpkg integrate install
```
so that VS knows all the paths.

## Prepare Data for particle tracking
The IWFM output data require alot of preprocessing.
For the preprocessing we use R.

The [Preproc.R](https://github.com/giorgk/IWFM_track/blob/master/Rwrkspc/PreprocTrack.R)
generates an h5 file that is needed in the following step.

In the first 16 lines the user needs to define 3 paths. There are explanations in the script.
Once those paths are set, simply run the script to generate the _PartTrackData.h5_ file.

Then run the following commands to write the required data
```
setwd("Path/of/Rwrkspc")
source("iwfm_funct.R")
D <- iwfm.loadData("PartTrackData.h5")
iwfm.writeData("XYfile.dat",D$XY[,])
iwfm.writeData("MSHfile.dat",D$MSH[,])
iwfm.writeData("STRATfile.dat",D$STRAT[,])
iwfm.writeData("NRMLfile.dat",D$NRML[,])
iwfm.writeData("FACEZfile.dat",D$FACEZ[,])
iwfm.writeData("BCfile.dat",D$BC[,])
iwfm.writeData("BCZfile.dat",D$BCZ[,])
```

The current version of particle tracking operates on a steady state field (at least the c++ implementation).

Therefore we need to average the flows. To caclulate the average flows from 2005 to 2015 do the following:
```
simTime <- iwfm.SimTime()
D <- iwfm.AverageVelField(ISOdate(2005,10,31), ISOdate(2015,9,30), D, simTime)
iwfm.writeData("HFLOWfile.dat",D$HFLOW[,])
iwfm.writeData("VFLOWfile.dat",D$VFLOW[,])
```

Another possibility is to use selected months instead of a range. 
First examine the output of the `iwfm.SimTime()` and decide which ids you should use.

For example, to average the values of all April months from 2005 to 2015 do the following
```
IDS <- seq(from=379, to=499, by=12)
D <- iwfm.AverageVelField_IDS(IDS, D)
``` 

## Run particle tracking
To run the particle tracking one need to prepare a configuration file with the following options:
```
XY=path\to\XYfile.dat
MSH=path\to\MSHfile.dat
STRAT=path\to\STRATfile.dat
NORMALS=path\to\NRMLfile.dat
FACEZ=path\to\FACEZfile.dat
HFLOW=path\to\HFLOWfile.dat
VFLOW=path\to\VFLOWfile.dat
BC=path\to\BCfile.dat
BCZ=path\to\BCZfile.dat
WELLS=path\to\welldataVS.dat
OUTFILE=path\to\streamlines.dat
```
Most of the above files are the ones generated from the preprocessing step.

The OUTFILE is the file where the streamlines will be printed
The WELLS file is a file with the following format:
```
4
554756 4497000 20 50 2015 8 1
566343 4414078 10 35 2015 8 1
697362 4156498 5 30 2015 8 1
752641 4024654 10 40 2015 8 1
```
where the first line is the number of wells.
The following lines are
X Y Depth ScreenLength Year month day

The depth is the distance from ground surface to the top of the screen

The year month day are not used in the c++ implementation

To run the code execute under a powershell 
```
.\iwfmTrack.exe -c config.dat
```

The code has many more options to configure the behaviour of particle tracking.
To see all options run the following
```
.\iwfmTrack.exe
```
## OUTPUT
The output of particle tracking is a file with the following format
```
100
1 0 1 17 
554761.403 4497008.415 128.371 0.10783 -0.04751 -0.02639 0.00000
554728.267 4497022.762 137.730 0.10379 -0.04548 -0.02319 0.00000
554690.439 4497038.991 148.292 0.10002 -0.04385 -0.01101 0.00000
.
.
.
554397.764 4497235.767 198.738 0.00000 0.00000 0.00000 0.00000
1 1 1 20 
554751.839 4497009.093 129.371 0.10789 -0.04782 -0.02641 0.00000
554739.843 4497014.398 132.338 0.10636 -0.04703 -0.02458 0.00000
554721.921 4497022.282 136.783 0.10425 -0.04599 -0.02322 0.00000
.
.
.
```
where :

the first line is the number of pathlines.

The next line is
```
Eid Sid exitflag Nsteps
```
where 

- Eid is the entity id (e.g. the well id)
- Sid is the streamline id
- exitflag is the reason why the particle tracking terminated for this particular streamline
- Nsteps Number of steps for this streamline.

Next the flowing line is repeated `Nsteps` times
```
X Y Z VX VY VZ AGE
```
where `X Y Z VX VY VZ` are the pathlines coordinates and velocity vector and 
`AGE` is the time since the start of particle tracking.


-----------------------


**OBSOLETE**

The following is an obsolete attempt to write the particle tracking under linux.
Because it had alot of dependencies many of which were hard to compile I stop the development.
However I keep it because it has some usefull information for other projects.

## Dependencies
The IWFM tracking is built on top of the following libraries:

* **[HDF5](https://www.hdfgroup.org/)**
[Download](https://www.hdfgroup.org/downloads/hdf5/source-code/) the CMake version and build it by executing the `build-unix.sh`. This is going to create a build directory under the main folder.

* **[CGAL](https://www.cgal.org/)**
[Download](https://github.com/CGAL/cgal/archive/releases/CGAL-4.11.3.tar.gz) and extract the source code. We use the version [4.11.3](https://github.com/CGAL/cgal/releases/tag/releases%2FCGAL-4.11.3) in this project, however newer or earlier versions might work as well, but have not tested.
Build the CGAL as follows:
```
cd /path/to/cgal-releases-CGAL-4.11.3
mkdir -p build/release
cd build/release
cmake -DCMAKE_BUILD_TYPE=Release ../..
make
```
Repeat the same steps by replacing the release with with debug to compile a debug version of the library
```
cd /path/to/cgal-releases-CGAL-4.11.3
mkdir -p build/debug
cd build/debug
cmake -DCMAKE_BUILD_TYPE=Debug ../..
make
```
* **[BOOST](https://www.boost.org/)**
Although boost is a header only library we use [Boost.Program_options](https://www.boost.org/doc/libs/1_69_0/doc/html/program_options.html)  library that requires linking with the library.

* **[EIGEN](http://eigen.tuxfamily.org)**
Eigen is a header only library and essensially the cmake commands just copy the files to the install directory. Download the latest stable version. At the moment we use version 3.3.7 and extract it to a folder eg. eigen_3_3_7. This is how we "_build_" Eigen:
```
cd eigen_3_3_7
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/Eigen_install ..
make install
```

* **[OpenSceneGraph](http://www.openscenegraph.org/)**
At the moment this is required. However we may remove or make this optional in the future. The easiest way to build this is to follow [these steps](https://vicrucann.github.io/tutorials/osg-linux-quick-install/). In short, just in case this link won't work, you need to do the following:

```
git clone https://github.com/openscenegraph/osg
cd osg
git checkout tags/OpenSceneGraph-x.y.z
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=/path/to/install/OSG ..
make
make install
```
 
## Building IWFM track
Once you have the above dependencies met you can build IWFMtrack as
```
cmake -G "Unix Makefiles" 
	    -DCMAKE_BUILD_TYPE=Debug 
	    -DHDF5_DIR=/path/to/CMake-hdf5-1.10.5/build/_CPack_Packages/Linux/TGZ/HDF5-1.10.5-Linux/HDF_Group/HDF5/1.10.5/share/cmake/hdf5 
	    -DCGAL_DIR:PATH=/path/to/cgal-releases-CGAL-4.11.3/build/debug 
	    -DOSG_DIR:PATH=/path/to/OSG/install/folder 
	    -DEigen3_DIR=/path/to/Eigen_install/share/eigen3/cmake .
```
 *Note that the above command has to be one line*
 
 
