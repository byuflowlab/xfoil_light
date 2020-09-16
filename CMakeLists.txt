cmake_minimum_required(VERSION 3.5)

project(libxfoil)
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "-O3 -fdefault-real-8 -fPIC -shared")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "-O3 -r8 -fPIC -shared")
endif()

file(GLOB_RECURSE sources src/*.f src/*.INC)
file(GLOB_RECURSE sources_cs src_cs/*.f src_cs/*.f90 src_cs/*.INC)
add_library(xfoil SHARED ${sources})
add_library(xfoil_cs SHARED ${sources_cs})
install(TARGETS xfoil xfoil_cs DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
