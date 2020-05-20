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

find_package(PythonLibs REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# ## NOT NEEDED FOR DICT
include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${ROOT_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
        $ENV{RATROOT}/include
        ${PYTHON_INCLUDE_DIRS}
)

#----------------------------------------------------------------------------
# Locate RATROOT Libs
#
LINK_DIRECTORIES($ENV{RATROOT}/lib)
LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/lib)

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

add_executable(TemplateAnalysisMT TemplateAnalysisMT.cc ${sources} ${headers})
target_link_libraries(TemplateAnalysisMT ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent)

add_executable(EventWrapper pyevent.cpp ${sources} ${headers})
target_link_libraries(EventWrapper ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent)

add_executable(PlotCollectedPE PlotCollectedPE.cc ${sources} ${headers})
target_link_libraries(PlotCollectedPE ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent)

add_executable(CreateEResMatrix CreateEResMatrix.cc ${sources} ${headers})
target_link_libraries(CreateEResMatrix ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent)

add_executable(CreatePosMatrix CreatePosMatrix.cc ${sources} ${headers})
target_link_libraries(CreatePosMatrix ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent)

add_executable(PlotPEVSNHits PlotPEVSNHits.cc ${sources} ${headers})
target_link_libraries(PlotPEVSNHits ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent)

add_executable(VtxRecon VtxRecon.cc ${sources} ${headers})
target_link_libraries(VtxRecon ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent)

add_executable(PECounting PECounting.cc ${sources} ${headers})
target_link_libraries(PECounting ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent)

add_executable(SelectIBD SelectIBD.cc ${sources} ${headers})
target_link_libraries(SelectIBD ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent)

add_executable(CreatePDF CreatePDF.cc ${sources} ${headers})
target_link_libraries(CreatePDF ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent)

add_executable(FlattenHits FlattenHits.cc ${sources} ${headers})
target_link_libraries(FlattenHits ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent)

add_executable(ReconHits ReconHits.cc ${sources} ${headers})
target_link_libraries(ReconHits ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent ${PYTHON_LIBRARIES})



add_executable(FitMacro FitMacro.C ${sources} ${headers})
target_link_libraries(FitMacro ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent -lEVFunctions -lHitFunctions -lLL)

add_executable(TResVSCosTheta TResVSCosTheta.C ${sources} ${headers})
target_link_libraries(TResVSCosTheta ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent -lEVFunctions -lHitFunctions -lLL)

add_executable(NMatrix NMatrix.C ${sources} ${headers})
target_link_libraries(NMatrix ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent -lEVFunctions -lHitFunctions -lLL)

add_executable(PMatrix PMatrix.C ${sources} ${headers})
target_link_libraries(PMatrix ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent -lEVFunctions -lHitFunctions -lLL)

add_executable(ProfileTResCosT ProfileTResCosT.C ${sources} ${headers})
target_link_libraries(ProfileTResCosT ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent -lEVFunctions -lHitFunctions -lLL)

add_executable(dE dE.C ${sources} ${headers})
target_link_libraries(dE ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent -lEVFunctions -lHitFunctions -lLL)

add_executable(PlotQVSNHits.C macros/PlotQVSNHits.C ${sources} ${headers})
target_link_libraries(PlotQVSNHits.C ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent -lEVFunctions -lHitFunctions -lLL)

add_executable(PlotTResCos.C macros/PlotTResCos.C ${sources} ${headers})
target_link_libraries(PlotTResCos.C ${ROOT_LIBRARIES} ${Boost_LIBRARIES} -lRATEvent -lEVFunctions -lHitFunctions -lLL -lMCFunctions -lCalibFunctions -lAnalyzerFunctions)
