# build Physsim library

LIST(APPEND ROOT_DICT_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)

SET(lib_sources
    anlib/ANL4DVector.cxx
    genlib/GENLib.cxx
    hellib/HELLib.cxx
    LCME/LCMEBase.cxx
    LCME/LCMEEEH.cxx
    LCME/LCMEEEZ.cxx
    LCME/LCMENNH.cxx
    LCME/LCMENNHH.cxx
    LCME/LCMEWW.cxx
    LCME/LCMEZH.cxx
    LCME/LCMEZHH.cxx
    LCME/LCMEZZ.cxx
    LCME/LCMEZZH.cxx
    LCME/LCMEZZZ.cxx
    utils/TAttLockable.cxx
    )

SET(headers
    ${CMAKE_SOURCE_DIR}/include/physsim/ANL4DVector.h
    ${CMAKE_SOURCE_DIR}/include/physsim/GENLib.h
    ${CMAKE_SOURCE_DIR}/include/physsim/GENNumCon.h
    ${CMAKE_SOURCE_DIR}/include/physsim/HELLib.h
    ${CMAKE_SOURCE_DIR}/include/physsim/LCMEBase.h
    ${CMAKE_SOURCE_DIR}/include/physsim/LCMEEEH.h
    ${CMAKE_SOURCE_DIR}/include/physsim/LCMEEEZ.h
    ${CMAKE_SOURCE_DIR}/include/physsim/LCMENNH.h
    ${CMAKE_SOURCE_DIR}/include/physsim/LCMENNHH.h
    ${CMAKE_SOURCE_DIR}/include/physsim/LCMEWW.h
    ${CMAKE_SOURCE_DIR}/include/physsim/LCMEZH.h
    ${CMAKE_SOURCE_DIR}/include/physsim/LCMEZHH.h
    ${CMAKE_SOURCE_DIR}/include/physsim/LCMEZZ.h
    ${CMAKE_SOURCE_DIR}/include/physsim/LCMEZZH.h
    ${CMAKE_SOURCE_DIR}/include/physsim/LCMEZZZ.h
    ${CMAKE_SOURCE_DIR}/include/physsim/TAttLockable.h
    ${CMAKE_SOURCE_DIR}/include/physsim/ZHHBases.h
)

ADD_SHARED_LIBRARY(Physsim ${lib_sources})
TARGET_INCLUDE_DIRECTORIES(Physsim
    PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
TARGET_LINK_LIBRARIES(Physsim PUBLIC ROOT::Core ROOT::Tree ROOT::Physics ROOT::RIO)
ROOT_GENERATE_DICTIONARY(Physsim_dict ${headers} MODULE Physsim LINKDEF LinkDef.h)

INSTALL(TARGETS
    Physsim
    EXPORT PhyssimTargets
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )

INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/libPhyssim_rdict.pcm
    ${CMAKE_CURRENT_BINARY_DIR}/libPhyssim.rootmap
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )

INSTALL(DIRECTORY
    ${CMAKE_SOURCE_DIR}/include/physsim
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

IF( APPLE ) #---- need special linker flags for ROOT dictionary on MacOS
    SET_TARGET_PROPERTIES( Physsim PROPERTIES
        LINK_FLAGS "-single_module -undefined dynamic_lookup -bind_at_load"
    )
ENDIF( APPLE )


