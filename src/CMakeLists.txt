# build Physsim library

SET( lib_input_dirs genlib hellib anlib utils LCME)

FOREACH( lib_input_dir ${lib_input_dirs} )
    LIST( APPEND ROOT_DICT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/${lib_input_dir} )
ENDFOREACH()

#MESSAGE( STATUS "ROOT_DICT_INCLUDE_DIRS: ${ROOT_DICT_INCLUDE_DIRS}" )


FOREACH( lib_input_dir ${lib_input_dirs} )

    AUX_SOURCE_DIRECTORY( ${lib_input_dir} lib_sources )

    PREPARE_ROOT_DICT_HEADERS( ${lib_input_dir} )

    INSTALL_DIRECTORY( ${lib_input_dir}/ DESTINATION "include/physsim"
        FILES_MATCHING PATTERN "*.h" PATTERN "LinkDef.h" EXCLUDE
    )

    GEN_ROOT_DICT_SOURCES( ${lib_input_dir}Dict.cxx )

    INSTALL(FILES ${ROOT_DICT_OUTPUT_DIR}/${lib_input_dir}Dict_rdict.pcm DESTINATION lib)
    
    LIST( APPEND lib_sources ${ROOT_DICT_OUTPUT_SOURCES} )

ENDFOREACH()

INCLUDE_DIRECTORIES( ${ROOT_DICT_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )

#MESSAGE( STATUS "Physsim lib sources: ${lib_sources}" )

ADD_SHARED_LIBRARY( Physsim ${lib_sources} )
INSTALL_SHARED_LIBRARY( Physsim DESTINATION lib )
TARGET_LINK_LIBRARIES( Physsim ${ROOT_LIBRARIES} )

IF( APPLE ) #---- need special linker flags for ROOT dictionary on MacOS
    SET_TARGET_PROPERTIES( Physsim PROPERTIES
        LINK_FLAGS "-single_module -undefined dynamic_lookup -bind_at_load"
    )
ENDIF( APPLE )


