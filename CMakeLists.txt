# Set the name of the project and target
set(TARGET "iwfm_track")

# Declare the source files the target consists of
SET(TARGET_SRC
	${TARGET}.cpp
	)

cmake_minimum_required (VERSION 3.10.1)
project( ${TARGET} C CXX )

set (CMAKE_CXX_STANDARD 14)

set (LIB_TYPE STATIC) # or SHARED
string(TOLOWER ${LIB_TYPE} SEARCH_TYPE)

#find_package (HDF5 NAMES hdf5 COMPONENTS C CXX ${SEARCH_TYPE})
#find_package(HDF5 COMPONENTS CXX HL REQUIRED)
find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
find_package( CGAL QUIET COMPONENTS Core )
find_package( OpenThreads )
find_package( osg )
find_package( osgDB )
find_package( osgUtil )
find_package( osgGA )
find_package( osgViewer )

include( ${CGAL_USE_FILE} )

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

find_package(Boost COMPONENTS graph program_options REQUIRED)


message("CGAL_USE_FILE: ${CGAL_USE_FILE}")

link_directories( iwfmHeaders ${HDF5_LIBRARY_DIRS} ${Boost_LIBRARY_DIR} ${OSG_LIB_DIR})

include_directories (${HDF5_INCLUDE_DIR} ${Boost_INCLUDE_DIR} ${OPENTHREADS_INCLUDE_DIR} ${OSG_INCLUDE_DIR})
set (LINK_LIBS ${LINK_LIBS} ${HDF5_C_${LIB_TYPE}_LIBRARY})


message("boost lib: ${Boost_LIBRARIES}")
message("boost inc: ${Boost_INCLUDE_DIR}")

add_executable (${TARGET} ${TARGET_SRC})

#add_to_cached_list( CGAL_EXECUTABLE_TARGETS test01 )

target_link_libraries (${TARGET} ${HDF5_CXX_LIBRARIES} ${Boost_LIBRARIES} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${OSG_LIBRARIES} ${OSGVIEWER_LIBRARIES} ${OSGUTIL_LIBRARIES} ${OSGDB_LIBRARIES} ${OSGGA_LIBRARIES} ${OPENTHREADS_LIBRARIES} Eigen3::Eigen)

