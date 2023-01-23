set (CURRENT_PACKAGE_NAME "calqual_exe")
set (CURRENT_EXE_NAME "calqual")

include (${CMAKE_SOURCE_DIR}/cmake/version.cmake)
include (${CMAKE_SOURCE_DIR}/cmake/common.cmake)
include (${CMAKE_SOURCE_DIR}/cmake/workarounds.cmake)

add_executable (calqual qualite2.cpp)

target_link_libraries (calqual PUBLIC qualif)
# INSTALL_RPATH modifie le rpath pour les libs internes au projet :
set_target_properties (calqual PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1	INSTALL_RPATH ${CMAKE_PACKAGE_RPATH_DIR})

# INSTALLATION :
include(CMakePackageConfigHelpers)
# ConfigPackageLocation : c'est plus ou moins standardisé, le défaut étant lib/cmake.
# On nos recommande très vivement d'utiliser ce répertoire afin de limiter les soucis par la suite,
# notamment au niveau des vues.
set (ConfigPackageLocation ${CMAKE_CMAKE_DIR})

install(TARGETS calqual EXPORT calqual DESTINATION ${CMAKE_INSTALL_BINDIR}
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Utilisation par d'autres modules de calqual :
set (RUNTIME_INSTALL_DIR bin/)		# Calqual_RUNTIME_DIR avec calqualConfig.cmake.in
configure_package_config_file(cmake/${CURRENT_EXE_NAME}Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_PACKAGE_NAME}/${CURRENT_EXE_NAME}Config.cmake
        INSTALL_DESTINATION ${CMAKE_CMAKE_DIR}
        PATH_VARS RUNTIME_INSTALL_DIR
    )
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_PACKAGE_NAME}/${CURRENT_EXE_NAME}ConfigVersion.cmake
        VERSION ${QUALIF_VERSION}
        COMPATIBILITY SameMinorVersion
    )
install (
   FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_PACKAGE_NAME}/${CURRENT_EXE_NAME}Config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${CURRENT_PACKAGE_NAME}/${CURRENT_EXE_NAME}ConfigVersion.cmake
   DESTINATION ${ConfigPackageLocation} PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ GROUP_WRITE WORLD_READ
   )


