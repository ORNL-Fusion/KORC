add_library(fruit ${CMAKE_CURRENT_SOURCE_DIR}/fruit/fruit.f90)

add_executable (xtest)
target_sources(xtest

               PRIVATE

               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_test.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test_io.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test_hpc.f90>
)
target_link_libraries (xtest

                       PUBLIC

                       korc
                       fruit
)

set_property(TARGET xtest PROPERTY LINKER_LANGUAGE Fortran)

add_test (NAME egyro_test
          COMMAND ${CMAKE_SOURCE_DIR}/test/egyro/korc_egyro.sh ${CMAKE_SOURCE_DIR}/test/egyro
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if(USE_PSPLINE)
  add_test (NAME mars_test
            COMMAND ${CMAKE_SOURCE_DIR}/test/mars/korc_mars.sh ${CMAKE_SOURCE_DIR}/test/mars
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif()

foreach(RANK RANGE 1 ${MPIEXEC_MAX_NUMPROCS})
  add_test (NAME unit_testing_${RANK}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${RANK} ./xtest TEST_OUT
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endforeach()

