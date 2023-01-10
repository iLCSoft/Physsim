set(CMAKE_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

include(CMakePackageConfigHelpers)

# Version file is same wherever we are
write_basic_package_version_file(${PROJECT_BINARY_DIR}/PhyssimConfigVersion.cmake
                                 VERSION ${Physsim_VERSION}
                                 COMPATIBILITY SameMajorVersion )

# Build tree config
# export(EXPORT PhyssimTargets NAMESPACE Physsim:: FILE ${PROJECT_BINARY_DIR}/PhyssimTargets.cmake)

# Install tree config
configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/PhyssimConfig.cmake.in
  ${PROJECT_BINARY_DIR}/PhyssimConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
  )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/PhyssimConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/PhyssimConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_CMAKEDIR})

install(EXPORT PhyssimTargets
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
  NAMESPACE Physsim::)
