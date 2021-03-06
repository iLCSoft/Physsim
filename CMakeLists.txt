########################################################
# cmake file for building KalFit libraries
# @author Jan Engels, DESY
# @author F. Gaede, DESY
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 FATAL_ERROR )
########################################################


# project name
PROJECT( Physsim )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 4 )
SET( ${PROJECT_NAME}_VERSION_PATCH 1 )



# project version

cmake_policy(SET CMP0008 NEW)  

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( ROOT REQUIRED )


# macros for generating root dict sources with rootcint
SET( ROOT_DICT_CINT_DEFINITIONS "-DHANDLE_DICT_EXCEPTIONS=IGNORED_FOR_CINT" )
INCLUDE( "${ILCUTIL_ROOT}/cmakemodules/MacroRootDict.cmake" )



### LIBRARY AND TOOLS #######################################################


# require proper c++
#ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
#----- need long long for int64 for now ------
#ADD_DEFINITIONS( "-Wno-long-long -fno-strict-aliasing" )


# Physsim library
ADD_SUBDIRECTORY( ./src )


# --------- kaltest examples -------------------------------------------------
# if this option is set to ON examples will be built with default target 'make'
# if set to OFF examples will only be built with target 'make examples'
#OPTION( BUILD_EXAMPLES "Set to ON to build examples" OFF )
#MESSAGE( STATUS "BUILD_EXAMPLES = ${BUILD_EXAMPLES}" )
#ADD_SUBDIRECTORY( ./examples/kaltest )
# ----------------------------------------------------------------------------



# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( PhyssimConfig.cmake PhyssimConfigVersion.cmake PhyssimLibDeps.cmake )

