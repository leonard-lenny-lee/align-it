cmake_minimum_required(VERSION 3.14)
project(align_it)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
            FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

option(BUILD_RDKIT_SUPPORT "Build against the RDKit instead of OpenBabel3" ON)
option(BUILD_PYTHON_SUPPORT "Build Python wrappers (requires RDKit support)" ON)

if(BUILD_PYTHON_SUPPORT AND (NOT BUILD_RDKIT_SUPPORT))
    message(WARNING "Disabling BUILD_PYTHON_SUPPORT because BUILD_RDKIT_SUPPORT is not enabled")
    set(BUILD_PYTHON_SUPPORT OFF)
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
set(ALIGNIT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(DEFINED ENV{CONDA_PREFIX})
    set(CMAKE_INSTALL_PREFIX $ENV{CONDA_PREFIX})
    set(CONDA_INCLUDE_DIR $ENV{CONDA_PREFIX}/include)
    set(CONDA_LIB_DIR $ENV{CONDA_PREFIX}/lib)
    include_directories(${CONDA_INCLUDE_DIR})
    link_directories(${CONDA_LIB_DIR})
endif()

set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 4)
file(WRITE  ${ALIGNIT_INCLUDE_DIR}/config.h "#define ALIGNIT_VERSION ${VERSION_MAJOR}\n")
file(APPEND ${ALIGNIT_INCLUDE_DIR}/config.h "#define ALIGNIT_RELEASE ${VERSION_MINOR}\n")
file(APPEND ${ALIGNIT_INCLUDE_DIR}/config.h "#define ALIGNIT_SUBRELEASE ${VERSION_PATCH}\n")

set(CMAKE_THREAD_LIBS_INIT "-lpthread")
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(CMAKE_USE_WIN32_THREADS_INIT 0)
set(CMAKE_USE_PTHREADS_INIT 1)
set(THREADS_PREFER_PTHREAD_FLAG ON)

include_directories(${ALIGNIT_INCLUDE_DIR})
if(NOT BUILD_RDKIT_SUPPORT)
    message(STATUS "Building Align-it against OpenBabel3")
    find_package(OpenBabel3 REQUIRED)
    include_directories(${OPENBABEL3_INCLUDE_DIR})
else(NOT BUILD_RDKIT_SUPPORT)  # RDKit
    message(STATUS "Building Align-it against RDKit")
    add_definitions(-DUSE_RDKIT)
    find_package(RDKit REQUIRED)
    find_package(Threads)
    set(RDKit_THREAD_LIBS Threads::Threads)
    find_package(Boost REQUIRED)
endif(NOT BUILD_RDKIT_SUPPORT)

add_library(${PROJECT_NAME} SHARED src/siMath.cpp src/coordinate.cpp src/pharmacophore.cpp
    src/stringTokenizer.cpp src/chargeFuncCalc.cpp src/solutionInfo.cpp src/utilities.cpp
    src/pharMerger.cpp src/hybridCalc.cpp src/alignment.cpp src/functionMapping.cpp
    src/result.cpp src/compScore.cpp
    src/aromFuncCalc.cpp src/hAccFuncCalc.cpp src/hDonFuncCalc.cpp src/lipoFuncCalc.cpp
    src/alignLib.cpp src/mainWar.cpp src/mainErr.cpp src/attaFuncCalc.cpp
)
target_include_directories(${PROJECT_NAME} PUBLIC ${OPENBABEL3_INCLUDE_DIR} ${RDKIT_INCLUDE_DIR} ${Boost_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${OPENBABEL3_LIBRARIES} ${RDKIT_LIBRARIES} ${RDKit_THREAD_LIBS})
install(TARGETS ${PROJECT_NAME} DESTINATION lib)
message("Installation directory of lib${PROJECT_NAME}: ${CMAKE_INSTALL_PREFIX}/lib")

if(BUILD_PYTHON_SUPPORT)
    if(NOT Python3_EXECUTABLE)  # Install Python into virtual environments, instead of globally
        execute_process(
            COMMAND which python
            OUTPUT_VARIABLE Python3_EXECUTABLE
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    endif()
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    add_subdirectory(src/Wrap)
    if(SKBUILD)
        message(STATUS "Building pyalignit using skbuild")
        install(TARGETS cpyalignit DESTINATION .)
    else()
        execute_process(
            COMMAND "${Python3_EXECUTABLE}" -c "if True:
              import sysconfig
              print(sysconfig.get_path('platlib'))"
            OUTPUT_VARIABLE PYTHON_SITE
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "Install location ${PYTHON_SITE}/pyalignit")
        install(DIRECTORY pyalignit DESTINATION ${PYTHON_SITE})
    endif(SKBUILD)
endif(BUILD_PYTHON_SUPPORT)
