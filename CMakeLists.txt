cmake_minimum_required (VERSION 3.14)

project (Qualif CXX C)

include (GNUInstallDirs)
include (cmake/version.cmake)
include (cmake/common.cmake)


set (Qualif_VERSION ${QUALIF_VERSION})

set (CMAKE_VERBOSE_MAKEFILE ON)
option (BUILD_SHARED_LIBS "Creation de bibliotheques dynamiques (defaut : ON)" ON)

add_subdirectory (src)


