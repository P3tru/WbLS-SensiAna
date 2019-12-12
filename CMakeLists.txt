cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(EventWrapper)

set(CMAKE_CXX_STANDARD 11)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)
#find_package(ROOT CONFIG REQUIRED)
if(NOT ROOT_FOUND)
    message(STATUS "ROOT package not found.")
    return()
endif()

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
#include(${ROOT_USE_FILE})

find_package(Boost 1.57.0 COMPONENTS system filesystem REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# ## NOT NEEDED FOR DICT
include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${ROOT_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
        $ENV{RATROOT}/include
)

#----------------------------------------------------------------------------
# Locate RATROOT Libs
#
LINK_DIRECTORIES($ENV{RATROOT}/lib)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


#----------------------------------------------------------------------------
# Add the executable
#
add_executable(TemplateAnalysis TemplateAnalysis.cc ${sources} ${headers})
target_link_libraries(TemplateAnalysis ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent)

add_executable(EventWrapper pyevent.cpp ${sources} ${headers})
target_link_libraries(EventWrapper ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent)

add_executable(PlotCollectedPE PlotCollectedPE.cc ${sources} ${headers})
target_link_libraries(PlotCollectedPE ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent)

add_executable(CreateEResMatrix CreateEResMatrix.cc ${sources} ${headers})
target_link_libraries(CreateEResMatrix ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent)

add_executable(CreatePosMatrix CreatePosMatrix.cc ${sources} ${headers})
target_link_libraries(CreatePosMatrix ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent)
