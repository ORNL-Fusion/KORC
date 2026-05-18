list(APPEND fusionio_headers
    ${CMAKE_CURRENT_SOURCE_DIR}/fusion_io_c.h
    ${CMAKE_CURRENT_SOURCE_DIR}/gato_source.h
    ${CMAKE_CURRENT_SOURCE_DIR}/gpec_source.h
    ${CMAKE_CURRENT_SOURCE_DIR}/fusion_io_defs.h
    ${CMAKE_CURRENT_SOURCE_DIR}/mars_field.h
    ${CMAKE_CURRENT_SOURCE_DIR}/geqdsk_field.h
    ${CMAKE_CURRENT_SOURCE_DIR}/fusion_io_field.h
    ${CMAKE_CURRENT_SOURCE_DIR}/mars_source.h
    ${CMAKE_CURRENT_SOURCE_DIR}/compound_field.h
    ${CMAKE_CURRENT_SOURCE_DIR}/geqdsk_source.h
    ${CMAKE_CURRENT_SOURCE_DIR}/fusion_io.h
    ${CMAKE_CURRENT_SOURCE_DIR}/options.h
    ${CMAKE_CURRENT_SOURCE_DIR}/interpolate.h
    ${CMAKE_CURRENT_SOURCE_DIR}/fusion_io_source.h
    ${CMAKE_CURRENT_SOURCE_DIR}/fusion_io_series.h
    ${CMAKE_CURRENT_SOURCE_DIR}/fio_operations.h
    ${CMAKE_CURRENT_SOURCE_DIR}/fusion_io_species.h
    ${CMAKE_CURRENT_SOURCE_DIR}/m3dc1_field.h
    ${CMAKE_CURRENT_SOURCE_DIR}/gato_field.h
    ${CMAKE_CURRENT_SOURCE_DIR}/gpec_field.h
    ${CMAKE_CURRENT_SOURCE_DIR}/m3dc1_source.h
    ${CMAKE_CURRENT_SOURCE_DIR}/compound_field.h
    ${CMAKE_CURRENT_SOURCE_DIR}/interpolate.h
)
list(APPEND fusionio_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/c_interface.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gato_source.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gpec_source.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fusion_io.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mars_field.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/geqdsk_field.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fusion_io_field.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mars_source.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/compound_field.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/geqdsk_source.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/options.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/interpolate.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/push_interface.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/fio_operations.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fusion_io_species.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fusion_io_series.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/m3dc1_field.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gato_field.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gpec_field.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/m3dc1_source.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fusion_io_fortran.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/compound_field.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/interpolate.cpp
)

option(USE_NIM "Use NIMROD inputs" OFF)

configure_file(Config.h.in Config.h)

IF(USE_NIM)

  list(APPEND fusionio_headers
    ${CMAKE_CURRENT_SOURCE_DIR}/nimrod_source.h
    ${CMAKE_CURRENT_SOURCE_DIR}/nimrod_field.h
  )

  list(APPEND fusionio_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/nimrod_source.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/nimrod_field.cpp
  )
endif()

add_library(fusionio
    ${fusionio_headers}
    ${fusionio_sources}
)

target_include_directories(fusionio PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

target_compile_options(fusionio
PRIVATE
    -fPIC
    $<$<COMPILE_LANGUAGE:Fortran>:>
)
target_link_libraries(fusionio PUBLIC m3dc1)
set_target_properties(fusionio
    PROPERTIES PUBLIC_HEADER "${fusionio_headers}"
)

install(
    TARGETS fusionio
    PUBLIC_HEADER DESTINATION include
    PRIVATE_HEADER DESTINATION include
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
