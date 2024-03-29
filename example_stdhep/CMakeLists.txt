########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 3.15 FATAL_ERROR)
########################################################


# project name
PROJECT( MEMExample )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

SET( ${PROJECT_NAME}_VERSION  "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}" )

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )
FIND_PACKAGE( ROOT 6.14 REQUIRED )
FIND_PACKAGE( Physsim REQUIRED)

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version

# optional package
#FIND_PACKAGE( AIDA )
#IF( AIDA_FOUND )
#    INCLUDE_DIRECTORIES( ${AIDA_INCLUDE_DIRS} )
#    LINK_LIBRARIES( ${AIDA_LIBRARIES} )
#    ADD_DEFINITIONS( "-DMARLIN_USE_AIDA" )
#    MESSAGE( STATUS "AIDA -- found" )
#ELSE()
#    MESSAGE( STATUS "AIDA -- not found" )
#ENDIF()



### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



### LIBRARY #################################################################

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( ./include )
#INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} PUBLIC Physsim::Physsim ${Marlin_LIBRARIES} ROOT::Tree)
TARGET_INCLUDE_DIRECTORIES( ${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    ${Marlin_INCLUDE_DIRS} )


INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()
