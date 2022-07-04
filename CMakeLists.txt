# ------------------------------------------------------------------------------
# CMakeLists.txt file Milos Ivanovic, Cranfield 2012.
# ------------------------------------------------------------------------------
#
# -- Define minimum CMake version ----------------------------------------------
cmake_minimum_required (VERSION 2.6)
#
# -- Define project name -------------------------------------------------------
#      (comment: this will become the name of the workspace in Visual Studio)
project (sph07)
#
# -- print the windows status flag ---------------------------------------------
message( STATUS "WIN32: " ${WIN32} )
message( STATUS "UNIX:  " ${UNIX} )
message( STATUS "HOME:  " $ENV{HOME} )
#
# -- Define parameters ---------------------------------------------------------
#
if( WIN32 )
  # -- for Windows platforms (includes 64bit despite the name) -----------------
  set (MPI_INC_DIR "F:/topalovic/PAK/SOLVER/PAKFS/MPICH2x86/include")
  set (MPI_LIBRARY_PATH  "F:/topalovic/PAK/SOLVER/PAKFS/MPICH2x86/lib")
  set (MPI_LIBRARY mpi.lib)
  # -- print the Library Path variable
  message( STATUS "MPI_LIBRARY_PATH: " ${MPI_LIBRARY_PATH} )
  message( STATUS "CMAKE_SYSTEM:     " ${CMAKE_SYSTEM} )
elseif (UNIX)  
  # -- other platforms ---------------------------------------------------------
  set (CMAKE_CXX_COMPILER mpicxx)
  set (CMAKE_LINKER mpicxx)
  # -- print the Library Path variable
  message( STATUS "CMAKE_SYSTEM:  " ${CMAKE_SYSTEM} )
endif()
# ------------------------------------------------------------------------------
#
# -- define compile options ---------------------------------------------------- 
if (WIN32)
  # -- for Windows patch paraview bug ------------------------------------------
  #      MICROSOFT_PATCH: fix number formating
  add_definitions( /DMICROSOFT_PATCH)
  #
  # -- for Windows patch paraview bug ------------------------------------------
  #      PARAVIEW_PATCH: write geometry info each time step rather than once
  #      OMPI_IMPORTS:   linking with OpenMPI
  add_definitions( /DPARAVIEW_PATCH /DOMPI_IMPORTS)
  #
  # -- define path to OpenMPI library ------------------------------------------
  #      Not sure how this is reflected in the Visual Studio Project Properties
  add_library(${MPI_LIBRARY} STATIC IMPORTED)
  set_property(TARGET ${MPI_LIBRARY} PROPERTY IMPORTED_LOCATION ${MPI_LIBRARY_PATH}/${MPI_LIBRARY})
elseif (UNIX)
  add_definitions(-DPARAVIEW_PATCH)
endif()

# ------------------------------------------------------------------------------
#
#set (BOOST_ROOT "E:/boost_1_50_0")
FIND_PACKAGE(Boost 1.46.0)
if (NOT Boost_FOUND)
  message(FATAL_ERROR "Unable to find correct Boost version. Did you set BOOST_ROOT?")
endif()
#
#  -- Specify that all *.cpp files are used for build --------------------------
file(GLOB SRC_FILES "*.cpp" "*.h")
add_executable(sph07 ${SRC_FILES})
#
# -- Define the version number -------------------------------------------------
set (SPH07_VERSION_MAJOR 0)
set (SPH07_VERSION_MINOR 2)
#
# Configure a header file to pass some of the CMake settings to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/SPH07Config.h.in"
  "${PROJECT_BINARY_DIR}/SPH07Config.h"
  )
#
# -- Define include directories ------------------------------------------------
include_directories(${PROJECT_BINARY_DIR} ${MPI_INC_DIR} ${Boost_INCLUDE_DIRS})
#
#  -- Define Libraries that should be linked to --------------------------------
#target_link_libraries(sph07 ${MPI_LIBRARY})
#
# ------------------------------------------------------------------------------
#   END
# ------------------------------------------------------------------------------

