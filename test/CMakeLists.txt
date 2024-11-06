add_library(fruit ${CMAKE_CURRENT_SOURCE_DIR}/fruit/fruit.f90)

add_executable (xtest)
target_sources(xtest

               PRIVATE

               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_test.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test_io.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test_hpc.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test_random.f90>
)
target_link_libraries (xtest

                       PUBLIC

                       korc
                       fruit
)

set_property(TARGET xtest PROPERTY LINKER_LANGUAGE Fortran)

configure_file(${CMAKE_SOURCE_DIR}/test/egyro/korc_egyro.sh.in ${CMAKE_BINARY_DIR}/egyro_test/korc_egyro.sh)
configure_file(${CMAKE_SOURCE_DIR}/test/mars/korc_mars.sh.in ${CMAKE_BINARY_DIR}/mars_test/korc_mars.sh)
configure_file(${CMAKE_SOURCE_DIR}/test/aorsa/korc_aorsa.sh.in ${CMAKE_BINARY_DIR}/aorsa_test/korc_aorsa.sh)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/egyro_test)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/mars_test)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/aorsa_test)

if(USE_PSPLINE)
  add_test (NAME mars_test_1
            COMMAND ${CMAKE_BINARY_DIR}/mars_test/korc_mars.sh 1
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  set_tests_properties (mars_test_1 PROPERTIES PROCESSORS 1 ENVIRONMENT OMP_NUM_THREADS=1)

  #add_test (NAME mars_test_2
  #          COMMAND ${CMAKE_BINARY_DIR}/mars_test/korc_mars.sh 2
  #          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  #set_tests_properties (mars_test_2 PROPERTIES PROCESSORS 2 ENVIRONMENT OMP_NUM_THREADS=1)

  #add_test (NAME aorsa_test_1
  #          COMMAND ${CMAKE_BINARY_DIR}/aorsa_test/korc_aorsa.sh 1
  #          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  #set_tests_properties (aorsa_test_1 PROPERTIES PROCESSORS 1 ENVIRONMENT OMP_NUM_THREADS=1)

  #add_test (NAME aorsa_test_2
  #          COMMAND ${CMAKE_BINARY_DIR}/aorsa_test/korc_aorsa.sh 2
  #          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  #set_tests_properties (aorsa_test_2 PROPERTIES PROCESSORS 2 ENVIRONMENT OMP_NUM_THREADS=1)
endif()

foreach(RANK IN ITEMS 1 2 4 8)
  if(${RANK} LESS_EQUAL ${MPIEXEC_MAX_NUMPROCS})
    add_test (NAME egyro_test_${RANK}
	      COMMAND ${CMAKE_BINARY_DIR}/egyro_test/korc_egyro.sh ${RANK}
              WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    set_tests_properties (egyro_test_${RANK} PROPERTIES PROCESSORS ${RANK} ENVIRONMENT OMP_NUM_THREADS=1)

    #if(USE_PSPLINE)
    #  add_test (NAME mars_test_${RANK}
    #      COMMAND ${CMAKE_BINARY_DIR}/mars_test/korc_mars.sh ${RANK}
    #            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    #  set_tests_properties (mars_test_${RANK} PROPERTIES PROCESSORS ${RANK} ENVIRONMENT OMP_NUM_THREADS=1)

    #  add_test (NAME aorsa_test_${RANK}
    #  COMMAND ${CMAKE_BINARY_DIR}/aorsa_test/korc_aorsa.sh ${RANK}
    #        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    #  set_tests_properties (aorsa_test_${RANK} PROPERTIES PROCESSORS ${RANK} ENVIRONMENT OMP_NUM_THREADS=1)
    #endif()

    add_test (NAME unit_testing_${RANK}
              COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${RANK} ./xtest TEST_OUT
              WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
   set_tests_properties (unit_testing_${RANK} PROPERTIES PROCESSORS ${RANK} ENVIRONMENT OMP_NUM_THREADS=1)
  endif()
endforeach()

