#  Append all source files to variable korc_sources. As new files are added,
#  this must be updated.
set (korc_sources "")
list (APPEND korc_sources
    korc_collisions.f90
    korc_fields.f90
    korc_interp.f90
    korc_spatial_distribution.f90
    korc_constants.f90
    korc_finalize.f90
    korc_ppusher.f90
    korc_HDF5.f90
    korc_coords.f90
    korc_hammersley_generator.f90
    korc_profiles.f90
    korc_types.f90
    korc_avalanche.f90
    korc_hpc.f90
    korc_rnd_numbers.f90
    korc_units.f90
    korc_experimental_pdf.f90
    korc_initialize.f90
    korc_velocity_distribution.f90
    korc_c_random.cpp
    korc_random.f90
    korc_fio_interface.f90
    korc_input.f90
    main.f90
)

add_executable (xkorc ${korc_sources})
set_property(TARGET xkorc PROPERTY LINKER_LANGUAGE Fortran)

option (CORI_DIR "Use Cori directories" OFF)
if (${CORI_DIR})
    add_definitions (-DCORI)
else ()
    add_definitions (-DMAC)
    if (${KORC_TEST})
       add_library (korcsrc STATIC ${korc_sources})
    endif ()
endif ()


target_link_libraries (xkorc ${MPI_Fortran_LIBRARIES})
target_link_libraries (xkorc ${HDF5_LIBRARIES})

if (${USE_PSPLINE})
   target_link_libraries (xkorc ${EZSpline_LIBRARIES})
endif()

if(${USE_FIO})
   target_link_libraries (xkorc ${FIO_LIBRARIES} ${M3DC1_LIBRARIES})
   target_link_libraries (xkorc ${NIMFIO_LIBRARIES} ${NIMBLK_LIBRARIES} ${NIMBTY_LIBRARIES} ${NIMLOC_LIBRARIES} ${NIMITER_LIBRARIES} ${NIMMAT_LIBRARIES} ${NIMPAR_LIBRARIES} ${NIMLIB_LIBRARIES} ${NIMMPI_LIBRARIES} ${NIMLAP_LIBRARIES} ${NIMSLU_LIBRARIES} ${NIMLUD_LIBRARIES})
endif()
