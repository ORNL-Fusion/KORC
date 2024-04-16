add_library(korc STATIC)

target_sources(korc

               PRIVATE

               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_collisions.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_fields.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_interp.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_spatial_distribution.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_constants.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_finalize.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_ppusher.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_HDF5.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_coords.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_hammersley_generator.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_profiles.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_types.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_avalanche.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_hpc.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_rnd_numbers.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_units.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_experimental_pdf.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_initialize.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_velocity_distribution.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_c_random.cpp>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_random.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_fio_interface.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_input.f90>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/korc_pspline.f90>
)
target_link_libraries (korc

                       PUBLIC

                       korc_depends
)


add_executable (xkorc)
target_sources (xkorc

                PRIVATE

                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/main.f90>
)
target_link_libraries (xkorc

                       PUBLIC

                       korc
)


set_property(TARGET xkorc PROPERTY LINKER_LANGUAGE Fortran)

