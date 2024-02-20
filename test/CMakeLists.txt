#  Append all source files to variable korc_sources. As new files are added,
#  this must be updated.

add_library (fruit fruit/fruit.f90)

set (korc_tests "")
list (APPEND korc_tests
     korc_test.f90
     test_io.f90
     test_hpc.f90
)

add_executable (xtest ${korc_tests})
set_property(TARGET xtest PROPERTY LINKER_LANGUAGE Fortran)

target_link_libraries (xtest fruit korc ${MPI_Fortran_LIBRARIES})
target_link_libraries (xtest ${HDF5_LIBRARIES})
target_link_libraries (xtest ${KORCSRC_LIBRARIES})

add_test (NAME egyro_test
          COMMAND ${CMAKE_SOURCE_DIR}/test/egyro/korc_egyro.sh ${CMAKE_SOURCE_DIR}/test/egyro
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/build/bin)

add_test (NAME unit_testing
          COMMAND xtest TEST_OUT
          WORKING_DIRECTORY build/bin)

