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
)

add_test (NAME egyro_test
          COMMAND ${CMAKE_SOURCE_DIR}/test/egyro/korc_egyro.sh ${CMAKE_SOURCE_DIR}/test/egyro
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/build/bin)

add_test (NAME unit_testing
          COMMAND xtest TEST_OUT
          WORKING_DIRECTORY build/bin)
