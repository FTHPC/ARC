cmake_minimum_required(VERSION 3.18)
project(ARC VERSION 0.0.2 LANGUAGES C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
option(BUILD_SHARED_LIBS "prefer shared libraries" ON)

#correct was to set a default build type
# https://blog.kitware.com/cmake-and-the-default-build-type/
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type was set. Setting build type to ${default_build_type}.")
  set(CMAKE_BUILD_TYPE ${default_build_type} CACHE 
    STRING "Choose the type to build" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()


include(GNUInstallDirs)
include(CTest)

find_package(OpenMP REQUIRED)


add_library(arc
  src/arc.c
  src/galois.c
  src/jerasure.c
  src/reed_sol.c
  )
target_link_libraries(
  arc
  PUBLIC
  OpenMP::OpenMP_C
  m
  )
target_include_directories(
  arc
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/ARC>
  )

install(TARGETS arc EXPORT ARCConfig
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
install(EXPORT ARCConfig NAMESPACE ARC:: DESTINATION ${CMAKE_INSTALL_DATADIR}/ARC/cmake)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ARC)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
option(BUILD_EXAMPLES "build example codes" OFF)
if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
