########################################################
# cmake file for building KalFit libraries
# @author Jan Engels, DESY
# @author F. Gaede, DESY
CMAKE_MINIMUM_REQUIRED( VERSION 3.15 FATAL_ERROR )
########################################################

# project name
PROJECT( Physsim )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 4 )
SET( ${PROJECT_NAME}_VERSION_PATCH 2 )

SET( ${PROJECT_NAME}_VERSION  "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}" )

#--- Define basic build settings -----------------------------------------------
# Provides install directory variables as defined for GNU software
include(GNUInstallDirs)

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( ROOT 6.14 REQUIRED COMPONENTS Physics )

# macros for generating root dict sources with rootcint
SET( ROOT_DICT_CINT_DEFINITIONS "-DHANDLE_DICT_EXCEPTIONS=IGNORED_FOR_CINT" )
INCLUDE( "${ILCUTIL_ROOT}/cmakemodules/MacroRootDict.cmake" )

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

include(cmake/PhyssimCreateConfig.cmake)
