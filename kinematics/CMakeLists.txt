cmake_minimum_required(VERSION 3.12)
project(Kinematics)

MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)
set(CMAKE_CXX_STANDARD 14)
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#include_directories(.)
#set example to be build
SET(examples  Kinematics)

# Specify path to boost header files
SET(MC_INC_DIR
        /Users/xuejiaoyang/MC++/include)



# Specify path to SUNDIALS header files
SET(SUNDIALS_INC_DIR
        /Users/xuejiaoyang/sundial/instdir/include
        CACHE STRING
        "Location of SUNDIALS header files")


# Add path to SUNDIALS header files and boost header files
INCLUDE_DIRECTORIES(${SUNDIALS_INC_DIR} )
INCLUDE_DIRECTORIES(${MC_INC_DIR})

#Set libraries
# Set search path for SUNDIALS libraries
SET(SUNDIALS_LIB_DIR /Users/xuejiaoyang/sundial/instdir/lib)

# Find the SUNDIALS solver's library
FIND_LIBRARY(SUNDIALS_SOLVER_LIB
        libsundials_cvode.a ${SUNDIALS_LIB_DIR}
        DOC "CVODE library")

# Find the NVECTOR library
FIND_LIBRARY(SUNDIALS_NVEC_LIB
        libsundials_nvecserial.a ${SUNDIALS_LIB_DIR}
        DOC "NVECTOR library")

# Set additional libraries
SET(SUNDIALS_EXTRA_LIB  -lm CACHE STRING "Additional libraries")

# List of all libraries
SET(SUNDIALS_LIBS ${SUNDIALS_SOLVER_LIB} ${SUNDIALS_NVEC_LIB} ${SUNDIALS_EXTRA_LIB})

#put all sourc files in this folder in a single variable, (they must be in the same root as main function)
aux_source_directory(. DIR_SRCS)

FOREACH(example ${examples})
    add_executable(${example} main.cpp ${DIR_SRCS})
    TARGET_LINK_LIBRARIES(${example} ${SUNDIALS_LIBS} )
ENDFOREACH(example ${examples})
