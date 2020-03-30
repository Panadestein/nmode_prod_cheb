cmake_minimum_required(VERSION 3.5)

project(sop_frb)
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check")
endif()

list(APPEND CMAKE_Fortran_FLAGS_DEBUG ${bounds})
list(APPEND CMAKE_Fortran_FLAGS ${dialect})

#
# Compile.
#
file(GLOB_RECURSE sources  src/*.f90 src/*.h)
add_executable(sop_frb ${sources})