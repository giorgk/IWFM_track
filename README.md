# IWFM_track
An interactive particle tracking tool for IWFM 

(Work in progress)

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
git osg
git checkout tags/OpenSceneGraph-x.y.z
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=/path/to/install/OSG
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
 
 