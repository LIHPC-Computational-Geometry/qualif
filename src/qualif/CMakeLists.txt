set (CURRENT_PACKAGE_NAME "qualif")
set (CURRENT_PACKAGE_MACRO_NAME "QUALIF")

include (${CMAKE_SOURCE_DIR}/cmake/version.cmake)
include (${CMAKE_SOURCE_DIR}/cmake/common.cmake)
include (${CMAKE_SOURCE_DIR}/cmake/workarounds.cmake)

file (GLOB HEADERS public/*.h)
file (GLOB CPP_SOURCES *.cpp)

# On ajoute les en-têtes aux sources. C'est utile pour cmake dans certains cas,
# par exemple lorsqu'ils doivent être pré-processés (moc, ...).
set (ALL_TARGETS qualif)

add_library (qualif ${CPP_SOURCES} ${HEADERS})
foreach (lib ${ALL_TARGETS})
   set_property (TARGET ${lib} PROPERTY VERSION ${PYTHON_UTIL_VERSION})
   set_property (TARGET ${lib} PROPERTY SOVERSION ${PYTHON_UTIL_MAJOR_VERSION})
# Etre capable une fois installée de retrouver les libs utilisées, ... :
   set_target_properties (${lib} PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1)
endforeach (lib)
target_link_libraries (qualif PUBLIC qualifinternal)
target_include_directories (qualif PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/public>$<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)

# INSTALLATION :
include(CMakePackageConfigHelpers)
# Génération du fichier contenant les informations sur la distribution :
configure_file (cmake/QualifDistribution.h.in  ${CMAKE_CURRENT_BINARY_DIR}/QualifDistribution.tmp)
configure_file (${CMAKE_CURRENT_BINARY_DIR}/QualifDistribution.tmp  ${CMAKE_CURRENT_BINARY_DIR}/QualifDistribution.h)


# ConfigPackageLocation : c'est plus ou moins standardisé, le défaut étant lib/cmake.
# On nos recommande très vivement d'utiliser ce répertoire afin de limiter les soucis par la suite,
# notamment au niveau des vues.
set (ConfigPackageLocation ${CMAKE_CMAKE_DIR})

# La génération des règles d'installation :
install (TARGETS ${ALL_TARGETS} EXPORT ${CMAKE_PACKAGE_TARGETS_NAME}
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
         INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Recopie des en-têtes dans le répertoire d'installation :
install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
# Dont QualifDistribution.h généré plus haut :
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/QualifDistribution.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Le fichier contenant les informations de version de ce module :
write_basic_package_version_file(${CMAKE_PACKAGE_VERSION_FILE} VERSION ${QUALIF_VERSION} COMPATIBILITY SameMajorVersion)
export (EXPORT ${CMAKE_PACKAGE_TARGETS_NAME} FILE ${CMAKE_PACKAGE_TARGETS_FILE} NAMESPACE ${CMAKE_MODULE_NAMESPACE})

# Utilisation par d'autres modules de qualif :
set (INCLUDE_INSTALL_DIR include/)		# qualif_INCLUDE_DIR avec qualifConfig.cmake.in
set (LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})	# qualif_LIB_DIR avec qualifConfig.cmake.in
# Création du fichier qualifConfig.cmake qui définira qualif_INCLUDE_DIR, qualif_LIB_DIR et qualif_BIN_DIR :
configure_package_config_file (cmake/${CMAKE_PACKAGE_CONFIG_NAME}.cmake.in ${CMAKE_PACKAGE_CONFIG_FILE}
	INSTALL_DESTINATION ${CMAKE_CMAKE_DIR}
	PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR)

# Copie du fichier qualifTargets.cmake dans prefix/lib/qualif/cmake :
install (EXPORT ${CMAKE_PACKAGE_TARGETS_NAME}
	FILE ${CMAKE_PACKAGE_TARGETS_BASENAME} NAMESPACE ${CMAKE_MODULE_NAMESPACE}
	DESTINATION ${ConfigPackageLocation} PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ GROUP_WRITE WORLD_READ)
# Copie des fichiers qualifConfig.cmake et qualifConfigVersion.cmake dans prefix/lib/qualif/cmake :
install (FILES ${CMAKE_PACKAGE_CONFIG_FILE} ${CMAKE_PACKAGE_VERSION_FILE} DESTINATION ${ConfigPackageLocation} PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ GROUP_WRITE WORLD_READ)

