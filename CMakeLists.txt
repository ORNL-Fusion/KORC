cmake_minimum_required (VERSION 3.12)

project (korc_proj C CXX Fortran)

if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
    message (FATAL_ERROR "In-source builds are not allowed. Please create a "
                         "directory and and run cmake from there, passing the "
                         "the path to this source directory as the last "
                         "argument. CMake created a CMakeCache.txt file and "
                         "a CMakeFiles Directory. It is safe to remove these.")
endif ()

#  FIXME: Remove once link issues are resolved with the new macos linker.
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   SET (CMAKE_EXE_LINKER_FLAGS -ld_classic)
endif ()

#-------------------------------------------------------------------------------
#  Setup build types.
#-------------------------------------------------------------------------------
# Set a default build type if none was specified
set (default_build_type "Release")
if (EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set(default_build_type "Debug")
endif ()

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
# Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                 "Debug"
                 "Release"
                 "MinSizeRel"
                 "RelWithDebInfo"
    )
endif ()

#-------------------------------------------------------------------------------
#  Build Options
#-------------------------------------------------------------------------------
option (USE_ACC "Use OpenACC" OFF)
option (USE_PSPLINE "Use PSPLINE interpolation" OFF)
option (USE_FIO "Use FIO interpolation support" OFF)

#-------------------------------------------------------------------------------
#  Sanitizer options
#-------------------------------------------------------------------------------
add_library (sanitizer INTERFACE)

macro (register_sanitizer_option name)
    string (TOUPPER ${name} upper_name)
    
    option (SANITIZE_${upper_name} "Enable the ${name} sanitizer" OFF)

    target_compile_options (sanitizer
                            INTERFACE
                            $<$<BOOL:${SANITIZE_${upper_name}}>:-fsanitize=${name}>
    )
    target_link_options (sanitizer
                         INTERFACE
                         $<$<BOOL:${SANITIZE_${upper_name}}>:-fsanitize=${name}>
    )
endmacro ()

register_sanitizer_option (address)
register_sanitizer_option (leak)
register_sanitizer_option (memory)
register_sanitizer_option (thread)
register_sanitizer_option (undefined)
register_sanitizer_option (float-divide-by-zero)

#-------------------------------------------------------------------------------
#  Build Locations
#-------------------------------------------------------------------------------
#  FIXME: Leave this as is for now but check if we want to use the cmake install
#         functions better.

#  Set directories to build all binary files.
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

option(KORC_BUILD_TESTS "Build tests for KORC" ON)

#  Configure OpenMP
#-------------------------------------------------------------------------------
#  Dependencies
#-------------------------------------------------------------------------------
find_package (OpenMP)

set (MPI_ROOT "" CACHE PATH "Search paths for MPI Library")
find_package (MPI REQUIRED)

set (HDF5_ROOT "" CACHE PATH "Search paths for HDF5 Library")
find_package (HDF5 REQUIRED COMPONENTS Fortran)

add_library (korc_depends INTERFACE)

target_link_libraries (korc_depends

                       INTERFACE

                       sanitizer
		       $<$<BOOL:${OpenMP_Fortran_FOUND}>:OpenMP::OpenMP_Fortran>
                       $<$<BOOL:${MPI_Fortran_FOUND}>:MPI::MPI_Fortran>
                       $<$<BOOL:${HDF5_FOUND}>:hdf5::hdf5_fortran>
                       $<$<BOOL:${USE_FIO}>:${FIO_LIBRARIES}>
                       $<$<BOOL:${USE_FIO}>:${M3DC1_LIBRARIES}>
                       $<$<BOOL:${USE_FIO}>:${NIMFIO_LIBRARIES}>
                       $<$<BOOL:${USE_FIO}>:${NIMBLK_LIBRARIES}>
                       $<$<BOOL:${USE_FIO}>:${NIMBTY_LIBRARIES}>
                       $<$<BOOL:${USE_FIO}>:${NIMLOC_LIBRARIES}>
                       $<$<BOOL:${USE_FIO}>:${NIMITER_LIBRARIES}>
                       $<$<BOOL:${USE_FIO}>:${NIMMAT_LIBRARIES}>
                       $<$<BOOL:${USE_FIO}>:${NIMPAR_LIBRARIES}>
                       $<$<BOOL:${USE_FIO}>:${NIMLIB_LIBRARIES}>
                       $<$<BOOL:${USE_FIO}>:${NIMMPI_LIBRARIES}>
                       $<$<BOOL:${USE_FIO}>:${NIMLAP_LIBRARIES}>
                       $<$<BOOL:${USE_FIO}>:${NIMSLU_LIBRARIES}>
                       $<$<BOOL:${USE_FIO}>:${NIMLUD_LIBRARIES}>
                       
)
target_compile_options (korc_depends

                        INTERFACE

                        $<$<COMPILE_LANGUAGE:Fortran>:-cpp>
)
#  FIXME: Consider making FIO a cmake module and using find_package(FIO) to set
#         the dependencies.
target_include_directories (korc_depends

                            INTERFACE

                            $<$<BOOL:${USE_FIO}>:${FIO_INCLUDE_PATH}>
                            $<$<BOOL:${USE_FIO}>:${NIMROD_INCLUDE_PATH}>
                            $<$<BOOL:${USE_FIO}>:${NIMFIO_INCLUDE_PATH}>
                            $<$<BOOL:${USE_FIO}>:${NIMCORE_INCLUDE_PATH}>
                            $<$<BOOL:${USE_FIO}>:${NIMLIB_INCLUDE_PATH}>
                            $<$<BOOL:${USE_FIO}>:${NIMEXT_INCLUDE_PATH}>
)
target_compile_definitions (korc_depends

                            INTERFACE

                            DOUBLE_PRECISION
                            PARALLEL_RANDOM
#                            $<$<CONFIG:Debug>:DBG_CHECK>
                            $<$<BOOL:${USE_ACC}>:ACC>
                            $<$<BOOL:${USE_PSPLINE}>:PSPLINE>
                            $<$<BOOL:${USE_FIO}>:FIO>
)
target_link_options (korc_depends

                     INTERFACE

                     $<$<BOOL:USE_FIO>:${FIO_LINK_FLAGS}>
                     $<$<BOOL:USE_FIO>:${NIMFIO_LINK_FLAGS}>
                     $<$<BOOL:USE_FIO>:${NIMCORE_LINK_FLAGS}>
                     $<$<BOOL:USE_FIO>:${NIMLIB_LINK_FLAGS}>
                     $<$<BOOL:USE_FIO>:${NIMEXT_LINK_FLAGS}>
)

#-------------------------------------------------------------------------------
#  Source files.
#-------------------------------------------------------------------------------
add_subdirectory (src)

#-------------------------------------------------------------------------------
#  Testing
#-------------------------------------------------------------------------------
if (${KORC_BUILD_TESTS})
    enable_testing ()
    add_subdirectory (test)
endif ()
