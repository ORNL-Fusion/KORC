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

target_link_libraries (xtest fruit korc_lib ${MPI_Fortran_LIBRARIES})
target_link_libraries (xtest ${HDF5_LIBRARIES})
target_link_libraries (xtest ${KORCSRC_LIBRARIES})

