##
#  CMake script
##
#Code structure and usage
##---------------------------------------------------------------------------------------------------------------------
#The openIGA code is basded on deal.ii, trinilios, SuperLU and hdf5
# this code contains 3 parts in general
# 1: IGA base: it is self complete. Basically it read .h5 mesh and provides nurbs basis function. 
#
# 2: FEM generic functions and superLU solver: It includes IGA base, and provides generic FEM functions and superLU solver
#
# 3: supplementary: separately classes providing mesh generator, special data structure and math operators etc.
#
# How to install
# in build folder where you found this CMakeList.txt, change it following guides below
# $cmake CMakeList.txt
# $make release install
# please check "libopenIGA.a" is installed in lib(or other) folder
# delete this build folder if necessary
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------


CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)


#Include deal.II
##---------------------------------------------------------------------------------------------------------------------
#usually DEAL_II_DIR for explicitly deal.II path; deal.II_DIR and DEAL_II_DIR could be preloaded enviromental variable

#SET(ENV{DEAL_II_DIR} "/nfs/mcfs_comp/home/rudraa/software/dealii-8.4.1/installWithTrilinos")
FIND_PACKAGE(deal.II 8.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()
##---------------------------------------------------------------------------------------------------------------------




##Include all source files, add library openIGA 
##---------------------------------------------------------------------------------------------------------------------
##Usually you do NOT need to change anything here
# Set the name of the project and target

SET(TARGET "main")
#this openIGA is NOT installed as a shared lib
#SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

FILE(GLOB IGAbase_SRC "../src/IGA/IGAbase/*.cc" "../src/IGA/solver/superLU_solver.cc"  "../src/supplementary/*.cc"  "../src/IGA/*.cc" "../src/IGA/function/*.cc" "../src/model/*.cc")
ADD_LIBRARY(openIGA ${IGAbase_SRC})
DEAL_II_SETUP_TARGET(openIGA)
FILE (GLOB TARGET_SRC  "main.cc")

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()
##---------------------------------------------------------------------------------------------------------------------



##Set dirdctory for openIGA, please give a full path and make sure the directory is valid 
##---------------------------------------------------------------------------------------------------------------------
install(TARGETS openIGA DESTINATION /home/wzhenlin/software/openIGA/lib)

##---------------------------------------------------------------------------------------------------------------------



#Include external lib as necessary (if they cannot be found automatically);
##---------------------------------------------------------------------------------------------------------------------
SET (hdfLibrary "/sw/arcts/centos7/hdf5/1.8.16-intel-16.0.3/lib/libhdf5_cpp.so")
SET (hdfLibrary2 "/sw/arcts/centos7/hdf5/1.8.16-intel-16.0.3/lib/libhdf5.a")
SET (hdfInclude "/sw/arcts/centos7/hdf5/1.8.16-intel-16.0.3/include")

SET(szip_LIB "/sw/arcts/centos7/szip/2.1/lib/libsz.so")
SET(szip_INC "/sw/arcts/centos7/szip/2.1/include")

SET (SUPERLU_INC "/home/wzhenlin/software/SuperLU_4.3/SRC/")
SET (SUPERLU_LIB "/home/wzhenlin/software/SuperLU_4.3/lib/libsuperlu_4.3.a")
##---------------------------------------------------------------------------------------------------------------------




#Include all pathes for head files and link all libraries, modify when you include external libs differently
##---------------------------------------------------------------------------------------------------------------------
INCLUDE_DIRECTORIES(${hdfInclude} ${SUPERLU_INC} ${szip_INC})
TARGET_LINK_LIBRARIES(main openIGA ${hdfLibrary} ${hdfLibrary2} ${SUPERLU_LIB} ${szip_LIB})













