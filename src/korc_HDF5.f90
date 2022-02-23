!! @note KORC module containing subroutines to read and write data in HDF5
!! files. @endnote
!! This module contains interfaces to use the HDF5 library in a more friendly
!! way. This module is intended to help developers to create new I/O
!! subroutines without having to deal with the sometimes cumbersome details
!! of the HDF5 API.
module korc_HDF5
  use korc_hpc
  use korc_types
  use korc_constants
  use HDF5

  IMPLICIT NONE

  INTEGER(HID_T), PRIVATE 	:: KORC_HDF5_REAL
  !! HDF5 real precision data type to be used in the simulation.
  INTEGER(SIZE_T), PRIVATE 	:: rp_hdf5
  !! Size of the HDF5 real precision data type used in the simulation.

  INTERFACE load_from_hdf5
     !! @note Fortran interface to subroutines loading a real or integer
     !! value from HDF5 files. @endnote
     module procedure iload_from_hdf5, rload_from_hdf5
  END INTERFACE load_from_hdf5


  INTERFACE load_array_from_hdf5
     !! @note Fortran interface to subroutines loading 2-D and 3-D arrays
     !! of real values from HDF5 files.
     module procedure rload_1d_array_from_hdf5, rload_3d_array_from_hdf5, rload_2d_array_from_hdf5
  END INTERFACE load_array_from_hdf5


  INTERFACE save_to_hdf5
     !! @note Fortran interface to subroutines saving real or integer
     !! values to HDF5 files.
     module procedure i1save_to_hdf5,i2save_to_hdf5,i4save_to_hdf5,i8save_to_hdf5,rsave_to_hdf5
  END INTERFACE save_to_hdf5

  !! @note Fortran interface to subroutines saving real and integer
  !! values to HDF5 files.
  INTERFACE save_1d_array_to_hdf5
     module procedure isave_1d_array_to_hdf5,rsave_1d_array_to_hdf5
  END INTERFACE save_1d_array_to_hdf5

  !! @note Fortran interface to subroutines saving 2-D arrays of real values to HDF5 files.
  !! @todo To code the corresponding subroutines for saving integer 2-D arrays.
  INTERFACE save_2d_array_to_hdf5
     module procedure rsave_2d_array_to_hdf5
  END INTERFACE save_2d_array_to_hdf5

  !! @note Fortran interface to subroutines saving 3-D arrays of real values to HDF5 files.
  !! @todo To include the corresponding subroutines for saving arrays of integers.
  INTERFACE save_3d_array_to_hdf5
     module procedure rsave_3d_array_to_hdf5
  END INTERFACE save_3d_array_to_hdf5

  !! @note Fortran interface to subroutines saving 1-D, 2-D or 3-D arrays of real values to HDF5 files.
  !! @todo To include the corresponding subroutines for saving arrays of integers.
  INTERFACE save_array_to_hdf5
     module procedure isave_1d_array_to_hdf5,rsave_1d_array_to_hdf5,rsave_2d_array_to_hdf5,rsave_3d_array_to_hdf5
  END INTERFACE save_array_to_hdf5

  PRIVATE :: rsave_to_hdf5,&
       isave_1d_array_to_hdf5,&
       rsave_1d_array_to_hdf5,&
       rsave_2d_array_to_hdf5,&
       iload_from_hdf5,&
       rload_from_hdf5,&
       rload_1d_array_from_hdf5,&
       rload_3d_array_from_hdf5,&
       rload_2d_array_from_hdf5,&
       i1save_to_hdf5,&
       i2save_to_hdf5,&
       i4save_to_hdf5,&
       i8save_to_hdf5

  PUBLIC :: initialize_HDF5,&
       finalize_HDF5,&
       save_simulation_parameters,&
       save_to_hdf5,&
       save_1d_array_to_hdf5,&
       save_2d_array_to_hdf5,&
       load_from_hdf5,&
       load_array_from_hdf5,&
       save_string_parameter,&
       load_time_stepping_params,&
       load_prev_time,&
       load_prev_iter,&
       save_restart_variables,&
       load_particles_ic

CONTAINS

  !! @note Initialization of HDF5 library.
  !!
  !! @param h5error HDF5 error status.
  subroutine initialize_HDF5()
    INTEGER :: h5error  ! Error flag
    call h5open_f(h5error)

#ifdef HDF5_DOUBLE_PRESICION
    call h5tcopy_f(H5T_NATIVE_DOUBLE, KORC_HDF5_REAL, h5error)
#elif HDF5_SINGLE_PRESICION
    call h5tcopy_f(H5T_NATIVE_REAL, KORC_HDF5_REAL, h5error)
#endif

    call h5tget_size_f(KORC_HDF5_REAL, rp_hdf5, h5error)
  end subroutine initialize_HDF5

  !! @note Finalization of HDF5 library.
  !!
  !! @param h5error HDF5 error status.
  subroutine finalize_HDF5()
    INTEGER :: h5error  ! Error flag
    call h5close_f(h5error)
  end subroutine finalize_HDF5

  !! @note Subroutine to load an integer datum from an HDF5 file.
  !!
  !! @todo Implement the reading of the attribute of idatum.
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the datum.
  !! @param[out] idatum Integer datum read from HDF5 file.
  !! @param[out] attr Attribute of datum read from HDF5 file.
  !! @param aname Name of idatum attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param h5error HDF5 error status.
  subroutine iload_from_hdf5(h5file_id,dset,idatum,attr)
    INTEGER(HID_T), INTENT(IN) 				:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 		:: dset
    INTEGER, INTENT(OUT) 				:: idatum
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(OUT) :: attr
    CHARACTER(4) 					:: aname = "Info"
    INTEGER(HID_T) 					:: dset_id
    INTEGER(HID_T) 					:: dspace_id
    INTEGER(HID_T) 					:: aspace_id
    INTEGER(HID_T) 					:: attr_id
    INTEGER(HID_T) 					:: atype_id
    INTEGER(HSIZE_T), DIMENSION(1) 			:: dims = (/1/)
    INTEGER(HSIZE_T), DIMENSION(1) 			:: adims = (/1/)
    INTEGER 						:: h5error

    ! * * * Read datum from file * * *

    call h5dopen_f(h5file_id, TRIM(dset), dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: iload_from_hdf5 &
            --> h5dopen_f")')
       call KORC_ABORT(14)
    end if

    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, idatum, dims, h5error)

    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: iload_from_hdf5 &
            --> h5dread_f")')
       call KORC_ABORT(14)
    end if

    call h5dclose_f(dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: iload_from_hdf5 &
            --> h5dclose_f")')
       call KORC_ABORT(14)
    end if

    if (PRESENT(attr)) then
       ! * * * Read attribute from file * * *

       ! * * * Read attribute from file * * *
    end if

    ! * * * Read datum from file * * *
  end subroutine iload_from_hdf5

  !! @note Subroutine to load a real datum from an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the datum.
  !! @param[out] rdatum Real datum read from HDF5 file and casted to
  !! KORC's real precision type.
  !! @param[out] attr Attribute of datum read from HDF5 file.
  !! @param raw_datum Datum read from HDF5 file.
  !! @param aname Name of rdatum attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param h5error HDF5 error status.
  !! @todo Implement the reading of the attribute of rdatum.
  subroutine rload_from_hdf5(h5file_id,dset,rdatum,attr)
    INTEGER(HID_T), INTENT(IN) 				:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 		:: dset
    REAL(rp), INTENT(OUT) 				:: rdatum
    REAL 						:: raw_datum
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(OUT) :: attr
    CHARACTER(4) 					:: aname = "Info"
    INTEGER(HID_T) 					:: dset_id
    INTEGER(HID_T) 					:: dspace_id
    INTEGER(HID_T) 					:: aspace_id
    INTEGER(HID_T) 					:: attr_id
    INTEGER(HID_T) 					:: atype_id
    INTEGER(HSIZE_T), DIMENSION(1) 			:: dims = (/1/)
    INTEGER(HSIZE_T), DIMENSION(1) 			:: adims = (/1/)
    INTEGER 						:: h5error

    ! * * * Read datum from file * * *

    call h5dopen_f(h5file_id, TRIM(dset), dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 &
            --> h5dopen_f")')
       call KORC_ABORT(14)
    end if

    call h5dread_f(dset_id, H5T_NATIVE_REAL, raw_datum, dims, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 &
            --> h5dread_f")')
       call KORC_ABORT(14)
    end if
    rdatum = REAL(raw_datum,rp)

    call h5dclose_f(dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 &
            --> h5dclose_f")')
       call KORC_ABORT(14)
    end if

    if (PRESENT(attr)) then
       ! * * * Read attribute from file * * *

       ! * * * Read attribute from file * * *
    end if

    ! * * * Read datum from file * * *
  end subroutine rload_from_hdf5

  !! @note Subroutine to load a 1-D array of reals from an HDF5 file.
  !! @details The dimension of the 1-D array rdata is determined by the
  !! input-output array rdata.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the data.
  !! @param[out] rdata 1-D array of real values read from HDF5 file and
  !! casted to KORC's real precision type.
  !! @param[out] attr 1-D array of attributes of rdata.
  !! @param raw_data 1-D array read from HDF5 file.
  !! @param aname Name of rdata attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param h5error HDF5 error status.
  !! @todo Implement the reading of the attributes of rdata.
  subroutine rload_1d_array_from_hdf5(h5file_id,dset,rdata,attr)
    INTEGER(HID_T), INTENT(IN) 				:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN)		:: dset
    REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) 	:: rdata
    REAL, DIMENSION(:), ALLOCATABLE 			:: raw_data
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(OUT) 	:: attr
    CHARACTER(MAX_STRING_LENGTH) 			:: aname
    INTEGER(HID_T) 					:: dset_id
    INTEGER(HID_T) 					:: dspace_id
    INTEGER(HID_T) 					:: aspace_id
    INTEGER(HID_T) 					:: attr_id
    INTEGER(HID_T) 					:: atype_id
    INTEGER(HSIZE_T), DIMENSION(1) 			:: dims
    INTEGER(HSIZE_T), DIMENSION(1) 			:: adims
    INTEGER 						:: h5error

    dims = (/ shape(rdata) /)

    ALLOCATE( raw_data(dims(1)) )

    ! * * * Read data from file * * *

    call h5dopen_f(h5file_id, TRIM(dset), dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dopen_f")')
       call KORC_ABORT(14)
    end if

    call h5dread_f(dset_id, H5T_NATIVE_REAL, raw_data, dims, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dread_f")')
       call KORC_ABORT(14)
    end if
    rdata = REAL(raw_data,rp)

    call h5dclose_f(dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dclose_f")')
       call KORC_ABORT(14)
    end if

    DEALLOCATE( raw_data )

    if (PRESENT(attr)) then
       ! * * * Read data attribute(s) from file * * *
    end if

    ! * * * Read data from file * * *
  end subroutine rload_1d_array_from_hdf5

  !! @note Subroutine to load a 2-D array of reals from an HDF5 file.
  !! @details The dimensions of the 2-D array rdata is determined by the input-output array rdata.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the data.
  !! @param[out] rdata 2-D array of real values read from HDF5 file and casted to KORC's real precision type.
  !! @param[out] attr 2-D array of attributes of rdata.
  !! @param raw_data 2-D array read from HDF5 file.
  !! @param aname Name of rdata attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param h5error HDF5 error status.
  !! @todo Implement the reading of the attributes of rdata.
  subroutine rload_2d_array_from_hdf5(h5file_id,dset,rdata,attr)
    INTEGER(HID_T), INTENT(IN) 														:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 										:: dset
    REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) 							:: rdata
    REAL, DIMENSION(:,:), ALLOCATABLE 												:: raw_data
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) 	:: attr
    CHARACTER(MAX_STRING_LENGTH) 													:: aname
    INTEGER(HID_T) 																	:: dset_id
    INTEGER(HID_T) 																	:: dspace_id
    INTEGER(HID_T) 																	:: aspace_id
    INTEGER(HID_T) 																	:: attr_id
    INTEGER(HID_T) 																	:: atype_id
    INTEGER(HSIZE_T), DIMENSION(2) 													:: dims
    INTEGER(HSIZE_T), DIMENSION(2) 													:: adims
    INTEGER 																		:: h5error

    dims = shape(rdata)

    ALLOCATE( raw_data(dims(1),dims(2)) )

    ! * * * Read data from file * * *

    call h5dopen_f(h5file_id, TRIM(dset), dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dopen_f")')
       call KORC_ABORT(14)
    end if

    call h5dread_f(dset_id, H5T_NATIVE_REAL, raw_data, dims, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dread_f")')
       call KORC_ABORT(14)
    end if
    rdata = REAL(raw_data,rp)

    call h5dclose_f(dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dclose_f")')
       call KORC_ABORT(14)
    end if

    DEALLOCATE( raw_data )

    if (PRESENT(attr)) then
       ! * * * Read data attribute(s) from file * * *
    end if

    ! * * * Read data from file * * *
  end subroutine rload_2d_array_from_hdf5

  !! @note Subroutine to load a 3-D array of reals from an HDF5 file.
  !! @details The dimensions of the 3-D array rdata is determined by the input-output array rdata.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the data.
  !! @param[out] rdata 3-D array of real values read from HDF5 file and casted to KORC's real precision type.
  !! @param[out] attr 3-D array of attributes of rdata.
  !! @param raw_data 3-D array read from HDF5 file.
  !! @param aname Name of rdata attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param h5error HDF5 error status.
  !! @todo Implement the reading of the attributes of rdata.
  subroutine rload_3d_array_from_hdf5(h5file_id,dset,rdata,attr)
    INTEGER(HID_T), INTENT(IN) 					:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 			:: dset
    REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) 	:: rdata
    REAL, DIMENSION(:,:,:), ALLOCATABLE 				:: raw_data
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) 	:: attr
    CHARACTER(MAX_STRING_LENGTH) 							:: aname
    INTEGER(HID_T) 									:: dset_id
    INTEGER(HID_T) 									:: dspace_id
    INTEGER(HID_T) 									:: aspace_id
    INTEGER(HID_T) 									:: attr_id
    INTEGER(HID_T) 									:: atype_id
    INTEGER(HSIZE_T), DIMENSION(3) 				:: dims
    INTEGER(HSIZE_T), DIMENSION(3) 				:: adims
    INTEGER 							:: h5error

    dims = shape(rdata)

    ALLOCATE( raw_data(dims(1),dims(2),dims(3)) )

    ! * * * Read data from file * * *

    call h5dopen_f(h5file_id, TRIM(dset), dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dopen_f")')
       call KORC_ABORT(14)
    end if

    call h5dread_f(dset_id, H5T_NATIVE_REAL, raw_data, dims, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dread_f")')
       call KORC_ABORT(14)
    end if
    rdata = REAL(raw_data,rp)

    call h5dclose_f(dset_id, h5error)
    if (h5error .EQ. -1) then
       write(output_unit_write,'("KORC ERROR: Something went wrong in: rload_from_hdf5 --> h5dclose_f")')
       call KORC_ABORT(14)
    end if

    DEALLOCATE( raw_data )

    if (PRESENT(attr)) then
       ! * * * Read data attribute(s) from file * * *
    end if

    ! * * * Read data from file * * *
  end subroutine rload_3d_array_from_hdf5

  !! @note Subroutine to write a 1 byte (8 bits) integer to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the datum.
  !! @param[in] idatum Integer datum read from HDF5 file.
  !! @param[in] attr Attribute of datum read from HDF5 file.
  !! @param aname Name of idatum attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param rank Number of dimensions of idatum's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param attrlen Lenght of idatum attribute's name.
  !! @param h5error HDF5 error status.
  subroutine i1save_to_hdf5(h5file_id,dset,idatum,attr)
    INTEGER(HID_T), INTENT(IN) 							:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 			:: dset
    INTEGER(KIND=1), INTENT(IN) 						:: idatum
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) 	:: attr
    CHARACTER(4) 										:: aname = "Info"
    INTEGER(HID_T) 										:: dset_id
    INTEGER(HID_T) 										:: dspace_id
    INTEGER(HID_T) 										:: aspace_id
    INTEGER(HID_T) 										:: attr_id
    INTEGER(HID_T) 										:: atype_id
    INTEGER(HSIZE_T), DIMENSION(1) 						:: dims = (/1/)
    INTEGER(HSIZE_T), DIMENSION(1) 						:: adims = (/1/)
    INTEGER 											:: rank = 1
    INTEGER 											:: arank = 1
    INTEGER(SIZE_T) 									:: attrlen
    INTEGER 											:: h5error

    ! * * * Write data to file * * *
    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_INTEGER, dspace_id, dset_id, h5error)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INT(idatum,idef), dims, h5error)

    if (PRESENT(attr)) then
       ! * * * Write attribute of data to file * * *
       attrlen = LEN_TRIM(attr)
       call h5screate_simple_f(arank,adims,aspace_id,h5error)
       call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
       call h5tset_size_f(atype_id, attrlen, h5error)
       call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
       call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

       call h5aclose_f(attr_id, h5error)
       call h5sclose_f(aspace_id, h5error)
       ! * * * Write attribute of data to file * * *
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *
  end subroutine i1save_to_hdf5

  !! @note Subroutine to write a 2 byte (16 bits) integer to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the datum.
  !! @param[in] idatum Integer datum read from HDF5 file.
  !! @param[in] attr Attribute of datum read from HDF5 file.
  !! @param aname Name of idatum attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param rank Number of dimensions of idatum's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param attrlen Lenght of idatum attribute's name.
  !! @param h5error HDF5 error status.
  subroutine i2save_to_hdf5(h5file_id,dset,idatum,attr)
    INTEGER(HID_T), INTENT(IN) 							:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 			:: dset
    INTEGER(KIND=2), INTENT(IN) 						:: idatum
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) 	:: attr
    CHARACTER(4) 										:: aname = "Info"
    INTEGER(HID_T) 										:: dset_id
    INTEGER(HID_T) 										:: dspace_id
    INTEGER(HID_T) 										:: aspace_id
    INTEGER(HID_T) 										:: attr_id
    INTEGER(HID_T) 										:: atype_id
    INTEGER(HSIZE_T), DIMENSION(1) 						:: dims = (/1/)
    INTEGER(HSIZE_T), DIMENSION(1) 						:: adims = (/1/)
    INTEGER 											:: rank = 1
    INTEGER 											:: arank = 1
    INTEGER(SIZE_T) 									:: attrlen
    INTEGER 											:: h5error

    ! * * * Write data to file * * *
    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_INTEGER, dspace_id, dset_id, h5error)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INT(idatum,idef), dims, h5error)

    if (PRESENT(attr)) then
       ! * * * Write attribute of data to file * * *
       attrlen = LEN_TRIM(attr)
       call h5screate_simple_f(arank,adims,aspace_id,h5error)
       call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
       call h5tset_size_f(atype_id, attrlen, h5error)
       call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
       call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

       call h5aclose_f(attr_id, h5error)
       call h5sclose_f(aspace_id, h5error)
       ! * * * Write attribute of data to file * * *
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *
  end subroutine i2save_to_hdf5

  !! @note Subroutine to write a 4 byte (32 bits) integer to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the datum.
  !! @param[in] idatum Integer datum read from HDF5 file.
  !! @param[in] attr Attribute of datum read from HDF5 file.
  !! @param aname Name of idatum attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param rank Number of dimensions of idatum's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param attrlen Lenght of idatum attribute's name.
  !! @param h5error HDF5 error status.
  subroutine i4save_to_hdf5(h5file_id,dset,idatum,attr)
    INTEGER(HID_T), INTENT(IN) 							:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 			:: dset
    INTEGER(KIND=4), INTENT(IN) 						:: idatum
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) 	:: attr
    CHARACTER(4) 										:: aname = "Info"
    INTEGER(HID_T) 										:: dset_id
    INTEGER(HID_T) 										:: dspace_id
    INTEGER(HID_T) 										:: aspace_id
    INTEGER(HID_T) 										:: attr_id
    INTEGER(HID_T) 										:: atype_id
    INTEGER(HSIZE_T), DIMENSION(1) 						:: dims = (/1/)
    INTEGER(HSIZE_T), DIMENSION(1) 						:: adims = (/1/)
    INTEGER 											:: rank = 1
    INTEGER 											:: arank = 1
    INTEGER(SIZE_T) 									:: attrlen
    INTEGER 											:: h5error

    ! * * * Write data to file * * *
    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_INTEGER, dspace_id, dset_id, h5error)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, INT(idatum,idef), dims, h5error)

    if (PRESENT(attr)) then
       ! * * * Write attribute of data to file * * *
       attrlen = LEN_TRIM(attr)
       call h5screate_simple_f(arank,adims,aspace_id,h5error)
       call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
       call h5tset_size_f(atype_id, attrlen, h5error)
       call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
       call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

       call h5aclose_f(attr_id, h5error)
       call h5sclose_f(aspace_id, h5error)
       ! * * * Write attribute of data to file * * *
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *
  end subroutine i4save_to_hdf5

  !! @note Subroutine to write a 8 byte (64 bits) integer to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the datum.
  !! @param[in] idatum Integer datum read from HDF5 file.
  !! @param[in] attr Attribute of datum read from HDF5 file.
  !! @param aname Name of idatum attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data read from HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param rank Number of dimensions of idatum's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param attrlen Lenght of idatum attribute's name.
  !! @param h5error HDF5 error status.
  subroutine i8save_to_hdf5(h5file_id,dset,idatum,attr)
    INTEGER(HID_T), INTENT(IN) 							:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 			:: dset
    INTEGER(KIND=8), INTENT(IN) 						:: idatum
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) 	:: attr
    CHARACTER(4) 										:: aname = "Info"
    INTEGER(HID_T) 										:: dset_id
    INTEGER(HID_T) 										:: dspace_id
    INTEGER(HID_T) 										:: aspace_id
    INTEGER(HID_T) 										:: attr_id
    INTEGER(HID_T) 										:: atype_id
    INTEGER(HSIZE_T), DIMENSION(1) 						:: dims = (/1/)
    INTEGER(HSIZE_T), DIMENSION(1) 						:: adims = (/1/)
    INTEGER 											:: rank = 1
    INTEGER 											:: arank = 1
    INTEGER(SIZE_T) 									:: attrlen
    INTEGER 											:: h5error

    ! * * * Write data to file * * *
    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_DOUBLE, dspace_id, dset_id, h5error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, REAL(idatum,8), dims, h5error)


    if (PRESENT(attr)) then
       ! * * * Write attribute of data to file * * *
       attrlen = LEN_TRIM(attr)
       call h5screate_simple_f(arank,adims,aspace_id,h5error)
       call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
       call h5tset_size_f(atype_id, attrlen, h5error)
       call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
       call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

       call h5aclose_f(attr_id, h5error)
       call h5sclose_f(aspace_id, h5error)
       ! * * * Write attribute of data to file * * *
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *
  end subroutine i8save_to_hdf5

  !! @note Subroutine to write a 1-D array of integer values to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the data.
  !! @param[in] idata Data written to HDF5 file.
  !! @param[in] attr Attributes of data written to HDF5 file.
  !! @param aname Name of idata attributes.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 data space identifier.
  !! @param aspace_id HDF5 data's attribute space identifier.
  !! @param attr_id HDF5 data's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data writen to HDF5 file.
  !! @param adims Dimensions of data's attributes written to HDF5 file.
  !! @param rank Number of dimensions of idata's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param attrlen Lenght of idata attribute's name.
  !! @param h5error HDF5 error status.
  !! @param rr Rank iterator.
  !! @param dd Dimension iterator.
  !! @bug When using a 1-D array of attributes, only the first attribute is saved.
  subroutine isave_1d_array_to_hdf5(h5file_id,dset,idata,attr)
    INTEGER(HID_T), INTENT(IN) 														:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 										:: dset
    INTEGER, DIMENSION(:), INTENT(IN) 												:: idata
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) 	:: attr
    CHARACTER(4) 																	:: aname = "Info"
    INTEGER(HID_T) 																	:: dset_id
    INTEGER(HID_T) 																	:: dspace_id
    INTEGER(HID_T) 																	:: aspace_id
    INTEGER(HID_T) 																	:: attr_id
    INTEGER(HID_T) 																	:: atype_id
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 									:: dims
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 									:: adims
    INTEGER 																		:: rank
    INTEGER 																		:: arank
    INTEGER(SIZE_T) 																:: attrlen
    INTEGER(SIZE_T) 																:: tmplen
    INTEGER 																		:: h5error
    INTEGER 																		:: rr,dd

    rank = size(shape(idata))
    ALLOCATE(dims(rank))
    dims = shape(idata)

    ! * * * Write data to file * * *
    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), H5T_NATIVE_INTEGER, dspace_id, dset_id, h5error)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, idata, dims, h5error)

    if (PRESENT(attr)) then
       arank = size(shape(attr))
       ALLOCATE(adims(arank))
       adims = shape(attr)

       ! * * * Write attribute of data to file * * *
       tmplen = 0
       attrlen = 0
       do rr=1_idef,arank
          do dd=1_idef,adims(rr)
             tmplen = LEN_TRIM(attr(dd))
             if ( tmplen .GT. attrlen) then
                attrlen = tmplen
             end if
          end do
       end do

       call h5screate_simple_f(arank,adims,aspace_id,h5error)
       call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
       call h5tset_size_f(atype_id, attrlen, h5error)
       call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
       call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

       call h5aclose_f(attr_id, h5error)
       call h5sclose_f(aspace_id, h5error)
       ! * * * Write attribute of data to file * * *

       DEALLOCATE(adims)
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *

    DEALLOCATE(dims)
  end subroutine isave_1d_array_to_hdf5

  !! @note Subroutine to write a real to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the datum.
  !! @param[in] rdatum Real datum written to HDF5 file.
  !! @param[in] attr Attribute of datum written to HDF5 file.
  !! @param aname Name of rdatum attribute.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 datum space identifier.
  !! @param aspace_id HDF5 datum's attribute space identifier.
  !! @param attr_id HDF5 datum's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data written to HDF5 file.
  !! @param adims Dimensions of data's attributes read from HDF5 file.
  !! @param rank Number of dimensions of rdatum's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param attrlen Lenght of rdatum attribute's name.
  !! @param h5error HDF5 error status.
  subroutine rsave_to_hdf5(h5file_id,dset,rdatum,attr)
    INTEGER(HID_T), INTENT(IN) 							:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 			:: dset
    REAL(rp), INTENT(IN) 								:: rdatum
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, INTENT(IN) 	:: attr
    CHARACTER(4) 										:: aname = "Info"
    INTEGER(HID_T) 										:: dset_id
    INTEGER(HID_T) 										:: dspace_id
    INTEGER(HID_T) 										:: aspace_id
    INTEGER(HID_T) 										:: attr_id
    INTEGER(HID_T) 										:: atype_id
    INTEGER(HSIZE_T), DIMENSION(1) 						:: dims = (/1/)
    INTEGER(HSIZE_T), DIMENSION(1) 						:: adims = (/1/)
    INTEGER 											:: rank = 1
    INTEGER 											:: arank = 1
    INTEGER(SIZE_T) 									:: attrlen
    INTEGER 											:: h5error

    ! * * * Write data to file * * *

    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), KORC_HDF5_REAL, dspace_id, dset_id, h5error)

    if (rp .EQ. INT(rp_hdf5)) then
       call h5dwrite_f(dset_id, KORC_HDF5_REAL, rdatum, dims, h5error)
    else
       call h5dwrite_f(dset_id, KORC_HDF5_REAL, REAL(rdatum,4), dims, h5error)
    end if

    if (PRESENT(attr)) then
       ! * * * Write attribute of data to file * * *
       attrlen = LEN_TRIM(attr)
       call h5screate_simple_f(arank,adims,aspace_id,h5error)
       call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
       call h5tset_size_f(atype_id, attrlen, h5error)
       call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
       call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

       call h5aclose_f(attr_id, h5error)
       call h5sclose_f(aspace_id, h5error)
       ! * * * Write attribute of data to file * * *
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *
  end subroutine rsave_to_hdf5

  !! @note Subroutine to write a 1-D array of real values to an HDF5 file.
  !!
  !! @bug When using a 1-D array of attributes, only the first attribute is saved.
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the data.
  !! @param[in] rdata Data written to HDF5 file.
  !! @param[in] attr Attributes of data written to HDF5 file.
  !! @param aname Name of rdata attributes.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 data space identifier.
  !! @param aspace_id HDF5 data's attribute space identifier.
  !! @param attr_id HDF5 data's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data writen to HDF5 file.
  !! @param adims Dimensions of data's attributes written to HDF5 file.
  !! @param rank Number of dimensions of rdata's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param tmplen Temporary length of rdata attribute's name.
  !! @param attrlen Lenght of rdata attribute's name.
  !! @param h5error HDF5 error status.
  !! @param rr Rank iterator.
  !! @param dd Dimension iterator.
  subroutine rsave_1d_array_to_hdf5(h5file_id,dset,rdata,attr)
    INTEGER(HID_T), INTENT(IN) 														:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 										:: dset
    REAL(rp), DIMENSION(:), INTENT(IN) 												:: rdata
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) 	:: attr
    CHARACTER(4) 																	:: aname = "Info"
    INTEGER(HID_T) 																	:: dset_id
    INTEGER(HID_T) 																	:: dspace_id
    INTEGER(HID_T) 																	:: aspace_id
    INTEGER(HID_T) 																	:: attr_id
    INTEGER(HID_T) 																	:: atype_id
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 									:: dims
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 									:: adims
    INTEGER 																		:: rank
    INTEGER 																		:: arank
    INTEGER(SIZE_T) 																:: tmplen
    INTEGER(SIZE_T) 																:: attrlen
    INTEGER 																		:: h5error
    INTEGER 																		:: rr,dd

    rank = size(shape(rdata))
    ALLOCATE(dims(rank))
    dims = shape(rdata)

    ! * * * Write data to file * * *

    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), KORC_HDF5_REAL, dspace_id, dset_id, h5error)

    if (rp .EQ. INT(rp_hdf5)) then
       call h5dwrite_f(dset_id, KORC_HDF5_REAL, rdata, dims, h5error)
    else
       call h5dwrite_f(dset_id, KORC_HDF5_REAL, REAL(rdata,4), dims, h5error)
    end if

    if (PRESENT(attr)) then
       arank = size(shape(attr))
       ALLOCATE(adims(arank))
       adims = shape(attr)

       ! * * * Write attribute of data to file * * *
       tmplen = 0
       attrlen = 0
       do rr=1_idef,arank
          do dd=1_idef,adims(rr)
             tmplen = LEN_TRIM(attr(dd))
             if ( tmplen .GT. attrlen) then
                attrlen = tmplen
             end if
          end do
       end do

       call h5screate_simple_f(arank,adims,aspace_id,h5error)
       call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, h5error)
       call h5tset_size_f(atype_id, attrlen, h5error)
       call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, h5error)
       call h5awrite_f(attr_id, atype_id, attr, adims, h5error)

       call h5aclose_f(attr_id, h5error)
       call h5sclose_f(aspace_id, h5error)
       ! * * * Write attribute of data to file * * *

       DEALLOCATE(adims)
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *

    DEALLOCATE(dims)
  end subroutine rsave_1d_array_to_hdf5

  !! @note Subroutine to write a 2-D array of real values to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the data.
  !! @param[in] rdata Data written to HDF5 file.
  !! @param[in] attr Attributes of data written to HDF5 file.
  !! @param aname Name of rdata attributes.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 data space identifier.
  !! @param aspace_id HDF5 data's attribute space identifier.
  !! @param attr_id HDF5 data's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data writen to HDF5 file.
  !! @param adims Dimensions of data's attributes written to HDF5 file.
  !! @param rank Number of dimensions of rdata's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param attrlen Lenght of rdata attribute's name.
  !! @param h5error HDF5 error status.
  !! @param rr Rank iterator.
  !! @param dd Dimension iterator.
  !! @todo Implement the writting of attributes to HDF5 file.
  subroutine rsave_2d_array_to_hdf5(h5file_id,dset,rdata,attr)
    INTEGER(HID_T), INTENT(IN) 														:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 										:: dset
    REAL(rp), DIMENSION(:,:), INTENT(IN) 											:: rdata
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) 	:: attr
    CHARACTER(4) 																	:: aname = "Info"
    INTEGER(HID_T) 																	:: dset_id
    INTEGER(HID_T) 																	:: dspace_id
    INTEGER(HID_T) 																	:: aspace_id
    INTEGER(HID_T) 																	:: attr_id
    INTEGER(HID_T) 																	:: atype_id
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 									:: dims
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 									:: adims
    INTEGER 																		:: rank
    INTEGER 																		:: arank
    INTEGER(SIZE_T) 																:: attrlen
    INTEGER 																		:: h5error
    INTEGER 																		:: rr,dd

    rank = size(shape(rdata))
    ALLOCATE(dims(rank))
    dims = shape(rdata)

    ! * * * Write data to file * * *

    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), KORC_HDF5_REAL, dspace_id, dset_id, h5error)

    if (rp .EQ. INT(rp_hdf5)) then
       call h5dwrite_f(dset_id, KORC_HDF5_REAL, rdata, dims, h5error)
    else
       call h5dwrite_f(dset_id, KORC_HDF5_REAL, REAL(rdata,4), dims, h5error)
    end if

    if (PRESENT(attr)) then
       ! * * * Write attribute of data to file * * *
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *

    DEALLOCATE(dims)
  end subroutine rsave_2d_array_to_hdf5

  !! @note Subroutine to write a 3-D array of real values to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the data.
  !! @param[in] rdata Data written to HDF5 file.
  !! @param[in] attr Attributes of data written to HDF5 file.
  !! @param aname Name of rdata attributes.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 data space identifier.
  !! @param aspace_id HDF5 data's attribute space identifier.
  !! @param attr_id HDF5 data's attribute identifier.
  !! @param atype_id Native HDF5 attribute type.
  !! @param dims Dimensions of data writen to HDF5 file.
  !! @param adims Dimensions of data's attributes written to HDF5 file.
  !! @param rank Number of dimensions of rdata's dataspace.
  !! @param arank Number of dimensions of attr's dataspace.
  !! @param attrlen Lenght of rdata attribute's name.
  !! @param h5error HDF5 error status.
  !! @param rr Rank iterator.
  !! @param dd Dimension iterator.
  !! @todo Implement the writting of attributes to HDF5 file.
  subroutine rsave_3d_array_to_hdf5(h5file_id,dset,rdata,attr)
    INTEGER(HID_T), INTENT(IN) 														:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 										:: dset
    REAL(rp), DIMENSION(:,:,:), INTENT(IN) 											:: rdata
    CHARACTER(MAX_STRING_LENGTH), OPTIONAL, DIMENSION(:), ALLOCATABLE, INTENT(IN) 	:: attr
    CHARACTER(4) 																	:: aname = "Info"
    INTEGER(HID_T) 																	:: dset_id
    INTEGER(HID_T) 																	:: dspace_id
    INTEGER(HID_T) 																	:: aspace_id
    INTEGER(HID_T) 																	:: attr_id
    INTEGER(HID_T) 																	:: atype_id
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 									:: dims
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 									:: adims
    INTEGER 																		:: rank
    INTEGER 																		:: arank
    INTEGER(SIZE_T) 																:: attrlen
    INTEGER 																		:: h5error
    INTEGER 																		:: rr,dd

    rank = size(shape(rdata))
    ALLOCATE(dims(rank))
    dims = shape(rdata)

    ! * * * Write data to file * * *

    call h5screate_simple_f(rank,dims,dspace_id,h5error)
    call h5dcreate_f(h5file_id, TRIM(dset), KORC_HDF5_REAL, dspace_id, dset_id, h5error)

    if (rp .EQ. INT(rp_hdf5)) then
       call h5dwrite_f(dset_id, KORC_HDF5_REAL, rdata, dims, h5error)
    else
       call h5dwrite_f(dset_id, KORC_HDF5_REAL, REAL(rdata,4), dims, h5error)
    end if

    if (PRESENT(attr)) then
       ! * * * Write attribute of data to file * * *
    end if

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)
    ! * * * Write data to file * * *

    DEALLOCATE(dims)
  end subroutine rsave_3d_array_to_hdf5

  !! @note Subroutine to write an array of strings to an HDF5 file.
  !!
  !! @param[in] h5file_id HDF5 file identifier.
  !! @param[in] dset String containing the name of the array of strings.
  !! @param[in] string_array Array of characters containing the strings to be written to HDF5 file.
  !! @param dset_id HDF5 data set identifier.
  !! @param dspace_id HDF5 data space identifier.
  !! @param dims Number of strings to be written to file.
  !! @param data_dims Dimensions of data written to HDF5 file. This is equal to (Maximum length of KORC string)x(Number of strings).
  !! @param str_len Size of strings to be written to file without blank spaces.
  !! @param string_type Native HDF5 string type.
  !! @param h5error HDF5 error status.
  subroutine save_string_parameter(h5file_id,dset,string_array)
    INTEGER(HID_T), INTENT(IN) 								:: h5file_id
    CHARACTER(MAX_STRING_LENGTH), INTENT(IN) 				:: dset
    CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), INTENT(IN) 	:: string_array
    INTEGER(HID_T) 											:: dset_id
    INTEGER(HID_T) 											:: dspace_id
    INTEGER(HSIZE_T), DIMENSION(1) 							:: dims
    INTEGER(HSIZE_T), DIMENSION(2) 							:: data_dims
    INTEGER(SIZE_T), DIMENSION(:), ALLOCATABLE 				:: str_len
    INTEGER(HID_T) 											:: string_type
    INTEGER 												:: h5error

    ALLOCATE(str_len(SIZE(string_array)))

    dims = (/SIZE(string_array)/)
    data_dims = (/MAX_STRING_LENGTH,SIZE(string_array)/)
    str_len = (/LEN_TRIM(string_array)/)

    call h5tcopy_f(H5T_STRING,string_type,h5error)
    call h5tset_strpad_f(string_type,H5T_STR_SPACEPAD_F,h5error)

    call h5screate_simple_f(1,dims,dspace_id,h5error)

    call h5dcreate_f(h5file_id,TRIM(dset),string_type,dspace_id,dset_id,h5error)

    call h5dwrite_vl_f(dset_id,string_type,string_array,data_dims,str_len,h5error,dspace_id)

    call h5sclose_f(dspace_id, h5error)
    call h5dclose_f(dset_id, h5error)

    DEALLOCATE(str_len)
  end subroutine save_string_parameter


  subroutine save_simulation_parameters(params,spp,F,P)
    !! @note Subroutine to save to a HDF5 file all the relevant simulation
    !! parameters. @endnote
    !! This subroutine saves to the HDF5 file "<a>simulation_parameters.h5</a>"
    !! all the relevant simulation parameters of KORC, most of them being part
    !! of the input file, but also including some derived quantities from the
    !! input parameters. This file is intended to facilitate the
    !! post-processing of KORC data using any software that supports
    !! the HDF5 software.
    TYPE(KORC_PARAMS), INTENT(IN) 				:: params
    !!Core KORC simulation parameters.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) 	:: spp
    !! An instance of KORC's derived type SPECIES containing all
    !! the information of different electron species. See [[korc_types]].
    TYPE(FIELDS), INTENT(IN) 					:: F
    !! An instance of KORC's derived type FIELDS containing all the information
    !! about the fields used in the simulation. See [[korc_types]]
    !! and [[korc_fields]].
    TYPE(PROFILES), INTENT(IN) 					:: P
    !! An instance of KORC's derived type PROFILES containing all the
    !! information about the plasma profiles used in the simulation.
    !! See [[korc_types]] and [[korc_profiles]].
    CHARACTER(MAX_STRING_LENGTH) 				:: filename
    !! String containing the name of the HDF5 file.
    CHARACTER(MAX_STRING_LENGTH) 				:: gname
    !! String containing the group name of a set of KORC parameters.
    CHARACTER(MAX_STRING_LENGTH) 				:: dset
    !! Name of data set to be saved to file.
    INTEGER(HID_T) 						:: h5file_id
    !!  HDF5 file identifier.
    INTEGER(HID_T) 						:: group_id
    !! HDF5 group identifier.
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 		:: dims
    !! Dimensions of data saved to HDF5 file.
    REAL(rp), DIMENSION(:), ALLOCATABLE 			:: rdata
    !! 1-D array of real data to be saved to HDF5 file.
    INTEGER, DIMENSION(:), ALLOCATABLE 				:: idata
    !! 1-D array of integer data to be saved to HDF5 file.
    CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
    !! An 1-D array with attributes of 1-D real or integer arrays that are
    !! passed to KORC interfaces of HDF5 I/O subroutines.
    CHARACTER(MAX_STRING_LENGTH) 				:: attr
    !!  A single attributes of real or integer data that is passed to KORC
    !! interfaces of HDF5 I/O subroutines.
    INTEGER 							:: h5error
    !! HDF5 error status.
    CHARACTER(19) 						:: tmp_str
    !! Temporary string used to manipulate various strings.
    REAL(rp) 							:: units
    !! Temporary variable used to add physical units to KORC parameters.

    ! * * * Error handling * * * !
    call h5eset_auto_f(params%HDF5_error_handling, h5error)
    ! Turn off: 0_idef. Turn on: 1_idef

    if (.NOT.(params%restart)) then

    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,'("Saving simulations parameters")')
    end if

       
       if (SIZE(params%outputs_list).GT.1_idef) then
          write(tmp_str,'(I18)') params%mpi_params%rank
          filename = TRIM(params%path_to_outputs) // "file_"  &
               // TRIM(ADJUSTL(tmp_str)) // ".h5"
          call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)
          call h5fclose_f(h5file_id, h5error)
       end if

       if (params%mpi_params%rank .EQ. 0) then
          filename = TRIM(params%path_to_outputs) // "simulation_parameters.h5"

          call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)

          ! Simulation parameters group
          gname = "simulation"
          call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

          ALLOCATE(attr_array(1))
          ALLOCATE(idata(1))

          dset = TRIM(gname) // "/field_model"
          call save_string_parameter(h5file_id,dset,(/params%field_model/))

          dset = TRIM(gname) // "/profile_model"
          call save_string_parameter(h5file_id,dset,(/params%profile_model/))

          dset = TRIM(gname) // "/simulation_time"
          attr = "Total aimed simulation time in seconds"
          call save_to_hdf5(h5file_id,dset,params%simulation_time* &
               params%cpp%time,attr)

          dset = TRIM(gname) // "/snapshot_frequency"
          attr = "Time between snapshots in seconds"
          call save_to_hdf5(h5file_id,dset,params%snapshot_frequency* &
               params%cpp%time,attr)

          dset = TRIM(gname) // "/dt"
          attr = "Time step in secs"
          call save_to_hdf5(h5file_id,dset,params%dt*params%cpp%time,attr)

          dset = TRIM(gname) // "/t_steps"
          attr_array(1) = "Number of time steps"
          idata = params%t_steps
          call save_1d_array_to_hdf5(h5file_id,dset,idata,attr_array)

          dset = TRIM(gname) // "/num_omp_threads"
          attr = "Number of omp threads"
          call save_to_hdf5(h5file_id,dset, params%num_omp_threads,attr)

          dset = TRIM(gname) // "/output_cadence"
          attr_array(1) = "Cadence of output files"
          idata = params%output_cadence
          call save_1d_array_to_hdf5(h5file_id,dset,idata,attr_array)

          dset = TRIM(gname) // "/HDF5_error_handling"
          attr_array(1) = "Error handling option: 0=OFF, 1=ON"
          idata = params%HDF5_error_handling
          call save_1d_array_to_hdf5(h5file_id,dset,idata,attr_array)

          dset = TRIM(gname) // "/restart_output_cadence"
          attr_array(1) = "Cadence of output files"
          idata = params%restart_output_cadence
          call save_1d_array_to_hdf5(h5file_id,dset,idata,attr_array)

          dset = TRIM(gname) // "/num_snapshots"
          attr_array(1) = "Number of outputs for each variable"
          idata = params%num_snapshots
          call save_1d_array_to_hdf5(h5file_id,dset,idata,attr_array)

          dset = TRIM(gname) // "/num_species"
          attr = "Number of particle species"
          call save_to_hdf5(h5file_id,dset,params%num_species,attr)

          dset = TRIM(gname) // "/nmpi"
          attr = "Number of mpi processes"
          call save_to_hdf5(h5file_id,dset,params%mpi_params%nmpi,attr)

          dset = TRIM(gname) // "/minimum_particle_energy"
          attr = "Minimum energy of simulated particles in eV"
          call save_to_hdf5(h5file_id,dset,params%minimum_particle_energy* &
               params%cpp%energy/C_E,attr)

          dset = TRIM(gname) // "/minimum_particle_g"
          attr = "Minimum relativistic factor gamma of simulated particles"
          call save_to_hdf5(h5file_id,dset,params%minimum_particle_g,attr)

          dset = TRIM(gname) // "/radiation"
          attr = "Radiation losses included in simulation"
          if(params%radiation) then
             call save_to_hdf5(h5file_id,dset,1_idef,attr)
          else
             call save_to_hdf5(h5file_id,dset,0_idef,attr)
          end if

          dset = TRIM(gname) // "/collisions"
          attr = "Collisions included in simulation"
          if(params%collisions) then
             call save_to_hdf5(h5file_id,dset,1_idef,attr)
          else
             call save_to_hdf5(h5file_id,dset,0_idef,attr)
          end if

          dset = TRIM(gname) // "/outputs_list"
          call save_string_parameter(h5file_id,dset,params%outputs_list)

          dset = TRIM(gname) // "/orbit_model"
          call save_string_parameter(h5file_id,dset,(/params%orbit_model/))

          dset = TRIM(gname) // "/field_eval"
          call save_string_parameter(h5file_id,dset,(/params%field_eval/))

          DEALLOCATE(idata)
          DEALLOCATE(attr_array)

          call h5gclose_f(group_id, h5error)


          ! Plasma species group
          gname = "species"
          call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

          ALLOCATE(attr_array(params%num_species))

          dset = TRIM(gname) // "/spatial_distribution"
          call save_string_parameter(h5file_id,dset,spp%spatial_distribution)

          dset = TRIM(gname) // "/energy_distribution"
          call save_string_parameter(h5file_id,dset,spp%energy_distribution)

          dset = TRIM(gname) // "/pitch_distribution"
          call save_string_parameter(h5file_id,dset,spp%pitch_distribution)

          dset = TRIM(gname) // "/ppp"
          attr_array(1) = "Particles per (mpi) process"
          call save_1d_array_to_hdf5(h5file_id,dset,spp%ppp,attr_array)

          dset = TRIM(gname) // "/q"
          attr_array(1) = "Electric charge"
          call save_1d_array_to_hdf5(h5file_id,dset,spp%q* &
               params%cpp%charge,attr_array)

          dset = TRIM(gname) // "/m"
          attr_array(1) = "Species mass in kg"
          call save_1d_array_to_hdf5(h5file_id,dset,spp%m* &
               params%cpp%mass,attr_array)

          dset = TRIM(gname) // "/Eo"
          attr_array(1) = "Initial (average) energy in eV"
          call save_1d_array_to_hdf5(h5file_id,dset,spp%Eo* &
               params%cpp%energy/C_E,attr_array)

          dset = TRIM(gname) // "/go"
          attr_array(1) = "Initial relativistic g factor."
          call save_1d_array_to_hdf5(h5file_id,dset,spp%go,attr_array)

          dset = TRIM(gname) // "/etao"
          attr_array(1) = "Initial pitch angle in degrees"
          call save_1d_array_to_hdf5(h5file_id,dset,spp%etao,attr_array)

          dset = TRIM(gname) // "/wc"
          attr_array(1) = "Average relativistic cyclotron frequency in Hz"
          call save_1d_array_to_hdf5(h5file_id,dset,spp%wc_r/ &
               params%cpp%time,attr_array)

          dset = TRIM(gname) // "/Ro"
          attr_array(1) = "Initial radial position of population"
          call save_1d_array_to_hdf5(h5file_id,dset,spp%Ro* &
               params%cpp%length,attr_array)

          dset = TRIM(gname) // "/PHIo"
          attr_array(1) = "Azimuthal angle in degrees."
          call save_1d_array_to_hdf5(h5file_id,dset,spp%PHIo* &
               180.0_rp/C_PI,attr_array)

          dset = TRIM(gname) // "/Zo"
          attr_array(1) = "Initial Z position of population"
          call save_1d_array_to_hdf5(h5file_id,dset,spp%Zo* &
               params%cpp%length,attr_array)

          dset = TRIM(gname) // "/ri"
          attr_array(1) = "Inner radius of initial spatial distribution"
          call save_1d_array_to_hdf5(h5file_id,dset,spp%r_inner* &
               params%cpp%length,attr_array)

          dset = TRIM(gname) // "/ro"
          attr_array(1) = "Outter radius of initial spatial distribution"
          call save_1d_array_to_hdf5(h5file_id,dset,spp%r_outter* &
               params%cpp%length,attr_array)

          dset = TRIM(gname) // "/falloff_rate"
          attr_array(1) = "Falloff of gaussian or exponential radial &
               profile in m"
          call save_1d_array_to_hdf5(h5file_id,dset,spp%falloff_rate/ &
               params%cpp%length,attr_array)

          dset = TRIM(gname) // "/shear_factor"
          attr_array(1) = "Shear factor (in case ELLIPTIC-TORUS  &
               spatial distribution is used."
          call save_1d_array_to_hdf5(h5file_id,dset,spp%shear_factor, &
               attr_array)

          dset = TRIM(gname) // "/sigmaR"
          attr_array(1) = "Variance of first dimension of 2D spatial & 
               distribution."
          call save_1d_array_to_hdf5(h5file_id,dset, &
               spp%sigmaR*params%cpp%length,attr_array)

          dset = TRIM(gname) // "/sigmaZ"
          attr_array(1) = "Variance of second dimension of 2D spatial &
               distribution."
          call save_1d_array_to_hdf5(h5file_id,dset, &
               spp%sigmaZ*params%cpp%length,attr_array)

          dset = TRIM(gname) // "/theta_gauss"
          attr_array(1) = "Angle of rotation of 2D spatial distribution."
          call save_1d_array_to_hdf5(h5file_id,dset,spp%theta_gauss,attr_array)

          dset = TRIM(gname) // "/psi_max"
          attr_array(1) = "Indicator function level of the argument of &
               the 2D gaussian exponential."
          call save_1d_array_to_hdf5(h5file_id,dset,spp%psi_max,attr_array)

          dset = TRIM(gname) // "/dth"
          attr_array(1) = "Variance of sampling normal variate for pitch angle."
          call save_1d_array_to_hdf5(h5file_id,dset,spp%dth,attr_array)

          dset = TRIM(gname) // "/dgam"
          attr_array(1) = "Variance of sampling normal variate for gamma."
          call save_1d_array_to_hdf5(h5file_id,dset,spp%dgam,attr_array)

          dset = TRIM(gname) // "/dR"
          attr_array(1) = "Variance of sampling normal variate for R."
          call save_1d_array_to_hdf5(h5file_id,dset,spp%dR*params%cpp%length, &
               attr_array)

          dset = TRIM(gname) // "/dZ"
          attr_array(1) = "Variance of sampling normal variate for Z."
          call save_1d_array_to_hdf5(h5file_id,dset,spp%dZ*params%cpp%length, &
               attr_array)
          
          call h5gclose_f(group_id, h5error)

          DEALLOCATE(attr_array)


          ! Plasma profiles group
!          if (params%collisions) then
             gname = "profiles"
             call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

             dset = TRIM(gname) // "/density_profile"
             call save_string_parameter(h5file_id,dset,(/P%ne_profile/))

             dset = TRIM(gname) // "/temperature_profile"
             call save_string_parameter(h5file_id,dset,(/P%Te_profile/))

             dset = TRIM(gname) // "/Zeff_profile"
             call save_string_parameter(h5file_id,dset,(/P%Zeff_profile/))

             dset = TRIM(gname) // "/neo"
             attr = "Density at the magnetic axis (m^-3)"
             call save_to_hdf5(h5file_id,dset,P%neo*params%cpp%density,attr)

             dset = TRIM(gname) // "/Teo"
             attr = "Temperature at the magnetic axis (eV)"
             call save_to_hdf5(h5file_id,dset,P%Teo* &
                  params%cpp%temperature/C_E,attr)

             dset = TRIM(gname) // "/Zeffo"
             attr = "Zeff at the magnetic axis"
             call save_to_hdf5(h5file_id,dset,P%Zeffo,attr)

             if (TRIM(params%profile_model) .EQ. 'ANALYTICAL') then
                dset = TRIM(gname) // "/n_ne"
                attr = "Exponent of tanh(x)^n for density profile"
                call save_to_hdf5(h5file_id,dset,P%n_ne,attr)

                dset = TRIM(gname) // "/a_ne"
                attr = "Coefficients f=ao+a1*r+a2*r^2+a3*r^3.  &
                     a_ne=[a0,a1,a2,a3]"
                call save_1d_array_to_hdf5(h5file_id,dset,P%a_ne)

                dset = TRIM(gname) // "/n_Te"
                attr = "Exponent of tanh(x)^n for density profile"
                call save_to_hdf5(h5file_id,dset,P%n_Te,attr)

                dset = TRIM(gname) // "/a_Te"
                attr = "Coefficients f=ao+a1*r+a2*r^2+a3*r^3.  &
                     a_Te=[a0,a1,a2,a3]"
                call save_1d_array_to_hdf5(h5file_id,dset,P%a_Te)

                dset = TRIM(gname) // "/n_Zeff"
                attr = "Exponent of tanh(x)^n for Zeff profile"
                call save_to_hdf5(h5file_id,dset,P%n_Zeff,attr)

                dset = TRIM(gname) // "/a_Zeff"
                attr = "Coefficients f=ao+a1*r+a2*r^2+a3*r^3.  &
                     a_Zeff=[a0,a1,a2,a3]"
                call save_1d_array_to_hdf5(h5file_id,dset,P%a_Zeff)

                if  (params%field_eval.EQ.'interp') then

                   ALLOCATE(attr_array(1))
                   dset = TRIM(gname) // "/dims"
                   attr_array(1) = "Mesh dimension of the profile (NR,NPHI,NZ)"
                   call save_1d_array_to_hdf5(h5file_id,dset,F%dims,attr_array)

                   dset = TRIM(gname) // "/R"
                   attr_array(1) = "Radial position of the magnetic field grid nodes"
                   call save_1d_array_to_hdf5(h5file_id,dset, &
                        F%X%R*params%cpp%length,attr_array)

                   if (ALLOCATED(F%X%PHI)) then
                      dset = TRIM(gname) // "/PHI"
                      attr_array(1) = "Azimuthal angle of the magnetic &
                           field grid nodes"
                      call save_1d_array_to_hdf5(h5file_id,dset,F%X%PHI,attr_array)
                   end if
                   
                   dset = TRIM(gname) // "/Z"
                   attr_array(1) = "Z position of the magnetic field grid nodes"
                   call save_1d_array_to_hdf5(h5file_id,dset,F%X%Z* &
                        params%cpp%length,attr_array)

                   dset = TRIM(gname) // "/ne"
                   units = params%cpp%density
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*P%ne_2D)

                   dset = TRIM(gname) // "/Te"
                   units = params%cpp%temperature
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*P%Te_2D)

                   dset = TRIM(gname) // "/Zeff"
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        P%Zeff_2D)

                   DEALLOCATE(attr_array)
                end if
                
             else if (params%profile_model(1:8) .EQ. 'EXTERNAL') then
                ALLOCATE(attr_array(1))
                dset = TRIM(gname) // "/dims"
                attr_array(1) = "Mesh dimension of the profiles (NR,NPHI,NZ)"
                call save_1d_array_to_hdf5(h5file_id,dset,P%dims,attr_array)

                dset = TRIM(gname) // "/R"
                attr_array(1) = "Grid nodes of profiles along the &
                     radial position"
                call save_1d_array_to_hdf5(h5file_id,dset,P%X%R* &
                     params%cpp%length,attr_array)

                if (ALLOCATED(F%X%PHI)) then
                   dset = TRIM(gname) // "/PHI"
                   attr_array(1) = "Grid nodes of profiles along the &
                        azimuthal position"
                   call save_1d_array_to_hdf5(h5file_id,dset, &
                        P%X%PHI,attr_array)
                end if

                dset = TRIM(gname) // "/Z"
                attr_array(1) = "Grid nodes of profiles along the Z position"
                call save_1d_array_to_hdf5(h5file_id,dset, &
                     P%X%Z*params%cpp%length,attr_array)

                dset = TRIM(gname) // "/ne"
                units = params%cpp%density
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*P%ne_2D)

                dset = TRIM(gname) // "/Te"
                units = params%cpp%temperature
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*P%Te_2D)

                dset = TRIM(gname) // "/Zeff"
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     P%Zeff_2D)

                if (params%profile_model(10:10).eq.'H') then

                   dset = TRIM(gname) // "/RHON"
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        P%RHON)
                   
                   dset = TRIM(gname) // "/nRE"
                   units = params%cpp%density
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*P%nRE_2D)

                   dset = TRIM(gname) // "/nAr0"
                   units = params%cpp%density
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*P%nAr0_2D)

                   dset = TRIM(gname) // "/nAr1"
                   units = params%cpp%density
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*P%nAr1_2D)

                   dset = TRIM(gname) // "/nAr2"
                   units = params%cpp%density
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*P%nAr2_2D)

                   dset = TRIM(gname) // "/nAr3"
                   units = params%cpp%density
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*P%nAr3_2D)

                   dset = TRIM(gname) // "/nD"
                   units = params%cpp%density
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*P%nD_2D)

                   dset = TRIM(gname) // "/nD1"
                   units = params%cpp%density
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*P%nD1_2D)
                   
                endif

                DEALLOCATE(attr_array)
             else if (params%profile_model .EQ. 'UNIFORM') then
                ! Something
             end if

             call h5gclose_f(group_id, h5error)
          !end if


          ! Electromagnetic fields group

          gname = "fields"
          call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

          if (TRIM(params%field_model(1:10)) .EQ. 'ANALYTICAL') then
             dset = TRIM(gname) // "/Bo"
             attr = "Toroidal field at the magnetic axis in T"
             call save_to_hdf5(h5file_id,dset,F%Bo*params%cpp%Bo,attr)

             dset = TRIM(gname) // "/current_direction"
             call save_string_parameter(h5file_id,dset, &
                  (/F%AB%current_direction/))

             dset = TRIM(gname) // "/a"
             attr = "Minor radius in m"
             call save_to_hdf5(h5file_id,dset,F%AB%a*params%cpp%length,attr)

             dset = TRIM(gname) // "/Ro"
             attr = "Magnetic axis radial position"
             call save_to_hdf5(h5file_id,dset,F%Ro*params%cpp%length,attr)

             dset = TRIM(gname) // "/Zo"
             attr = "Magnetic axis vertical position"
             call save_to_hdf5(h5file_id,dset,F%Zo*params%cpp%length,attr)
             
             dset = TRIM(gname) // "/qa"
             attr = "Safety factor at minor radius"
             call save_to_hdf5(h5file_id,dset,F%AB%qa,attr)

             dset = TRIM(gname) // "/qo"
             attr = "Safety factor at the magnetic axis"
             call save_to_hdf5(h5file_id,dset,F%AB%qo,attr)

             dset = TRIM(gname) // "/lambda"
             attr = "Parameter lamda in m"
             call save_to_hdf5(h5file_id,dset,F%AB%lambda* &
                  params%cpp%length,attr)

             dset = TRIM(gname) // "/Bpo"
             attr = "Poloidal magnetic field in T"
             call save_to_hdf5(h5file_id,dset,F%AB%Bpo*params%cpp%Bo,attr)
             
             dset = TRIM(gname) // "/Eo"
             attr = "Electric field at the magnetic axis in V/m"
             call save_to_hdf5(h5file_id,dset,F%Eo*params%cpp%Eo,attr)

             if (params%SC_E) then
                dset = TRIM(gname) // "/dt_E_SC"
                attr = "Time step for self-consistent E calculation"
                call save_to_hdf5(h5file_id,dset,F%dt_E_SC,attr)

                dset = TRIM(gname) // "/Ip_exp"
                attr = "Scaling for self-consistent current density"
                call save_to_hdf5(h5file_id,dset,F%Ip_exp,attr)

                dset = TRIM(gname) // "/Ip0"
                attr = "Total RE current normalization"
                call save_to_hdf5(h5file_id,dset,F%Ip0,attr)
                
                ALLOCATE(attr_array(1))
                dset = TRIM(gname) // "/r_1D"
                attr_array(1) = "1D minor radial mesh for &
                     self-consistent fields"
                call save_1d_array_to_hdf5(h5file_id,dset,F%r_1D,attr_array)
                DEALLOCATE(attr_array)
             end if
             
             if  (params%field_eval.EQ.'interp') then

                ALLOCATE(attr_array(1))
                dset = TRIM(gname) // "/dims"
                attr_array(1) = "Mesh dimension of the magnetic  &
                     field (NR,NPHI,NZ)"
                call save_1d_array_to_hdf5(h5file_id,dset,F%dims,attr_array)

                dset = TRIM(gname) // "/R"
                attr_array(1) = "Radial position of the magnetic field grid nodes"
                call save_1d_array_to_hdf5(h5file_id,dset, &
                     F%X%R*params%cpp%length,attr_array)

                if (ALLOCATED(F%X%PHI)) then
                   dset = TRIM(gname) // "/PHI"
                   attr_array(1) = "Azimuthal angle of the magnetic &
                        field grid nodes"
                   call save_1d_array_to_hdf5(h5file_id,dset,F%X%PHI,attr_array)
                end if

                dset = TRIM(gname) // "/Z"
                attr_array(1) = "Z position of the magnetic field grid nodes"
                call save_1d_array_to_hdf5(h5file_id,dset,F%X%Z* &
                     params%cpp%length,attr_array)

                if (ALLOCATED(F%PSIp)) then
                   dset = TRIM(gname) // "/psi_p"
                   units = params%cpp%Bo*params%cpp%length**2
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*F%PSIp)
                end if
                
                if(params%field_model(12:13).eq.'2D') then
                
                   dset = TRIM(gname) // "/BR"
                   units = params%cpp%Bo
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*F%B_2D%R)

                   dset = TRIM(gname) // "/BPHI"
                   units = params%cpp%Bo
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*F%B_2D%PHI)

                   dset = TRIM(gname) // "/BZ"
                   units = params%cpp%Bo
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*F%B_2D%Z)
                   
                   if (ALLOCATED(F%FLAG2D)) then
                      dset = TRIM(gname) // "/Flag"
                      call rsave_2d_array_to_hdf5(h5file_id, dset, &
                           F%FLAG2D)
                   end if
                
                   if  (params%orbit_model(3:5).EQ.'pre') then

                      dset = TRIM(gname) // "/gradBR"
                      units = params%cpp%Bo/params%cpp%length
                      call rsave_2d_array_to_hdf5(h5file_id, dset, &
                           units*F%gradB_2D%R)

                      dset = TRIM(gname) // "/gradBPHI"
                      units = params%cpp%Bo/params%cpp%length
                      call rsave_2d_array_to_hdf5(h5file_id, dset, &
                           units*F%gradB_2D%PHI)

                      dset = TRIM(gname) // "/gradBZ"
                      units = params%cpp%Bo/params%cpp%length
                      call rsave_2d_array_to_hdf5(h5file_id, dset, &
                           units*F%gradB_2D%Z)

                      dset = TRIM(gname) // "/curlbR"
                      units = 1./params%cpp%length
                      call rsave_2d_array_to_hdf5(h5file_id, dset, &
                           units*F%curlb_2D%R)

                      dset = TRIM(gname) // "/curlbPHI"
                      units = 1./params%cpp%length
                      call rsave_2d_array_to_hdf5(h5file_id, dset, &
                           units*F%curlb_2D%PHI)

                      dset = TRIM(gname) // "/curlbZ"
                      units =1./params%cpp%length
                      call rsave_2d_array_to_hdf5(h5file_id, dset, &
                           units*F%curlb_2D%Z)

                   end if

                else if(params%field_model(12:13).eq.'3D') then
                
                   dset = TRIM(gname) // "/BR"
                   units = params%cpp%Bo
                   call rsave_3d_array_to_hdf5(h5file_id, dset, &
                        units*F%B_3D%R)

                   dset = TRIM(gname) // "/BPHI"
                   units = params%cpp%Bo
                   call rsave_3d_array_to_hdf5(h5file_id, dset, &
                        units*F%B_3D%PHI)

                   dset = TRIM(gname) // "/BZ"
                   units = params%cpp%Bo
                   call rsave_3d_array_to_hdf5(h5file_id, dset, &
                        units*F%B_3D%Z)
                   
                   if (ALLOCATED(F%FLAG3D)) then
                      dset = TRIM(gname) // "/Flag"
                      call rsave_3d_array_to_hdf5(h5file_id, dset, &
                           F%FLAG3D)
                   end if
                
                   if  (params%orbit_model(3:5).EQ.'pre') then

                      dset = TRIM(gname) // "/gradBR"
                      units = params%cpp%Bo/params%cpp%length
                      call rsave_3d_array_to_hdf5(h5file_id, dset, &
                           units*F%gradB_3D%R)

                      dset = TRIM(gname) // "/gradBPHI"
                      units = params%cpp%Bo/params%cpp%length
                      call rsave_3d_array_to_hdf5(h5file_id, dset, &
                           units*F%gradB_3D%PHI)

                      dset = TRIM(gname) // "/gradBZ"
                      units = params%cpp%Bo/params%cpp%length
                      call rsave_3d_array_to_hdf5(h5file_id, dset, &
                           units*F%gradB_3D%Z)

                      dset = TRIM(gname) // "/curlbR"
                      units = 1./params%cpp%length
                      call rsave_3d_array_to_hdf5(h5file_id, dset, &
                           units*F%curlb_3D%R)

                      dset = TRIM(gname) // "/curlbPHI"
                      units = 1./params%cpp%length
                      call rsave_3d_array_to_hdf5(h5file_id, dset, &
                           units*F%curlb_3D%PHI)

                      dset = TRIM(gname) // "/curlbZ"
                      units =1./params%cpp%length
                      call rsave_3d_array_to_hdf5(h5file_id, dset, &
                           units*F%curlb_3D%Z)

                   end if
                end if
                   
                DEALLOCATE(attr_array)
             end if

          else if (params%field_model(1:8) .EQ. 'EXTERNAL') then
             ALLOCATE(attr_array(1))
             
             dset = TRIM(gname) // "/dims"
             attr_array(1) = "Mesh dimension of the magnetic  &
                  field (NR,NPHI,NZ)"
             call save_1d_array_to_hdf5(h5file_id,dset,F%dims,attr_array)

             dset = TRIM(gname) // "/R"
             attr_array(1) = "Radial position of the magnetic field grid nodes"
             call save_1d_array_to_hdf5(h5file_id,dset, &
                  F%X%R*params%cpp%length,attr_array)

             if (ALLOCATED(F%X%PHI)) then
                if (F%Dim2x1t) then
                   dset = TRIM(gname) // "/PHI"
                   attr_array(1) = "Azimuthal angle of the magnetic &
                        field grid nodes"
                   call save_1d_array_to_hdf5(h5file_id,dset, &
                        F%X%PHI*params%cpp%time,attr_array)
                else
                   dset = TRIM(gname) // "/PHI"
                   attr_array(1) = "Azimuthal angle of the magnetic &
                        field grid nodes"
                   call save_1d_array_to_hdf5(h5file_id,dset,F%X%PHI,attr_array)
                end if
             end if

             dset = TRIM(gname) // "/Z"
             attr_array(1) = "Z position of the magnetic field grid nodes"
             call save_1d_array_to_hdf5(h5file_id,dset,F%X%Z* &
                  params%cpp%length,attr_array)

             dset = TRIM(gname) // "/Bo"
             attr = "Toroidal field at the magnetic axis in T"
             call save_to_hdf5(h5file_id,dset,F%Bo*params%cpp%Bo,attr)

             dset = TRIM(gname) // "/Eo"
             attr = "Electric field at the magnetic axis in V/m"
             call save_to_hdf5(h5file_id,dset,F%Eo*params%cpp%Eo,attr)

             dset = TRIM(gname) // "/E_dyn"
             attr = "Magnitude of dynamic E"
             call save_to_hdf5(h5file_id,dset,F%E_dyn*params%cpp%Eo,attr)
             
             dset = TRIM(gname) // "/E_pulse"
             attr = "Magnitude of dynamic E"
             call save_to_hdf5(h5file_id,dset,F%E_pulse*params%cpp%time,attr)

             dset = TRIM(gname) // "/E_width"
             attr = "Magnitude of dynamic E"
             call save_to_hdf5(h5file_id,dset,F%E_width*params%cpp%time,attr)
             
             dset = TRIM(gname) // "/Ro"
             attr = "Radial position of magnetic axis"
             call save_to_hdf5(h5file_id,dset,F%Ro*params%cpp%length,attr)

             dset = TRIM(gname) // "/Zo"
             attr = "Radial position of magnetic axis"
             call save_to_hdf5(h5file_id,dset,F%Zo*params%cpp%length,attr)

             dset = TRIM(gname) // "/PSIP0"
             attr = "Poloidal flux at magnetic axis"
             call save_to_hdf5(h5file_id,dset, &
                  F%PSIP_min*(params%cpp%Bo*params%cpp%length**2),attr)

             dset = TRIM(gname) // "/PSIPlim"
             attr = "Poloidal flux at LCFS"
             call save_to_hdf5(h5file_id,dset, &
                  F%PSIp_lim*(params%cpp%Bo*params%cpp%length**2),attr)

             dset = TRIM(gname) // "/PSIP_conv"
             attr = "Scaling factor fof poloidal flux function"
             call save_to_hdf5(h5file_id,dset, &
                  F%psip_conv,attr)
             
             dset = TRIM(gname) // "/Axisymmetric"
             attr = "Radial position of magnetic axis"
             if(F%axisymmetric_fields) then
                call save_to_hdf5(h5file_id,dset,1_idef,attr)
             else
                call save_to_hdf5(h5file_id,dset,0_idef,attr)
             end if

             if (ALLOCATED(F%PSIp)) then
                dset = TRIM(gname) // "/psi_p"
                units = params%cpp%Bo*params%cpp%length**2
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%PSIp)
             end if

             if (ALLOCATED(F%E_3D%R)) then
                dset = TRIM(gname) // "/ER"
                units = params%cpp%Eo
                call rsave_3d_array_to_hdf5(h5file_id, dset, &
                     units*F%E_3D%R)
             end if
             
             if (ALLOCATED(F%E_3D%PHI)) then
                dset = TRIM(gname) // "/EPHI"
                units = params%cpp%Eo
                call rsave_3d_array_to_hdf5(h5file_id, dset, &
                     units*F%E_3D%PHI)
             end if

             if (ALLOCATED(F%E_3D%Z)) then
                dset = TRIM(gname) // "/EZ"
                units = params%cpp%Eo
                call rsave_3d_array_to_hdf5(h5file_id, dset, &
                     units*F%E_3D%Z)
             end if
             
             if (ALLOCATED(F%PSIp3D)) then
                dset = TRIM(gname) // "/psi_p3D"
                units = params%cpp%Bo*params%cpp%length**2
                call rsave_3d_array_to_hdf5(h5file_id, dset, &
                     units*F%PSIp3D)
             end if

             if (ALLOCATED(F%FLAG2D)) then
                dset = TRIM(gname) // "/Flag2D"
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     F%FLAG2D)
             end if

             if (ALLOCATED(F%FLAG3D)) then
                dset = TRIM(gname) // "/Flag3D"
                call rsave_3d_array_to_hdf5(h5file_id, dset, &
                     F%FLAG3D)
             end if


             if (ALLOCATED(F%B1Re_2D%R)) then
                dset = TRIM(gname) // "/BR1_Re"
                units = params%cpp%Bo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%B1Re_2D%R)
             end if

             if (ALLOCATED(F%B1Re_2D%PHI)) then
                dset = TRIM(gname) // "/BPHI1_Re"
                units = params%cpp%Bo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%B1Re_2D%PHI)
             end if

             if (ALLOCATED(F%B1Re_2D%Z)) then
                dset = TRIM(gname) // "/BZ1_Re"
                units = params%cpp%Bo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%B1Re_2D%Z)
             end if

             if (ALLOCATED(F%B1Im_2D%R)) then
                dset = TRIM(gname) // "/BR1_Im"
                units = params%cpp%Bo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%B1Im_2D%R)
             end if

             if (ALLOCATED(F%B1Im_2D%PHI)) then
                dset = TRIM(gname) // "/BPHI1_Im"
                units = params%cpp%Bo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%B1Im_2D%PHI)
             end if

             if (ALLOCATED(F%B1Im_2D%Z)) then
                dset = TRIM(gname) // "/BZ1_Im"
                units = params%cpp%Bo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%B1Im_2D%Z)
             end if

             if (ALLOCATED(F%B1Re_2DX%X)) then
                dset = TRIM(gname) // "/BX1_Re"
                units = params%cpp%Bo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%B1Re_2DX%X)
             end if

             if (ALLOCATED(F%B1Re_2DX%Y)) then
                dset = TRIM(gname) // "/BY1_Re"
                units = params%cpp%Bo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%B1Re_2DX%Y)
             end if

             if (ALLOCATED(F%B1Re_2DX%Z)) then
                dset = TRIM(gname) // "/BZ1_Re"
                units = params%cpp%Bo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%B1Re_2DX%Z)
             endif

             if (ALLOCATED(F%B1Im_2DX%X)) then
                dset = TRIM(gname) // "/BX1_Im"
                units = params%cpp%Bo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%B1Im_2DX%X)
             end if

             if (ALLOCATED(F%B1Im_2DX%Y)) then
                dset = TRIM(gname) // "/BY1_Im"
                units = params%cpp%Bo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%B1Im_2DX%Y)
             end if

             if (ALLOCATED(F%B1Im_2DX%Z)) then
                dset = TRIM(gname) // "/BZ1_Im"
                units = params%cpp%Bo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%B1Im_2DX%Z)
             endif
             
             if (ALLOCATED(F%E1Re_2DX%X)) then
                dset = TRIM(gname) // "/EX1_Re"
                units = params%cpp%Eo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%E1Re_2DX%X)
             end if

             if (ALLOCATED(F%E1Re_2DX%Y)) then
                dset = TRIM(gname) // "/EY1_Re"
                units = params%cpp%Eo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%E1Re_2DX%Y)
             end if

             if (ALLOCATED(F%E1Re_2DX%Z)) then
                dset = TRIM(gname) // "/EZ1_Re"
                units = params%cpp%Eo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%E1Re_2DX%Z)
             endif

             if (ALLOCATED(F%E1Im_2DX%X)) then
                dset = TRIM(gname) // "/EX1_Im"
                units = params%cpp%Eo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%E1Im_2DX%X)
             end if

             if (ALLOCATED(F%E1Im_2DX%Y)) then
                dset = TRIM(gname) // "/EY1_Im"
                units = params%cpp%Eo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%E1Im_2DX%Y)
             end if

             if (ALLOCATED(F%E1Im_2DX%Z)) then
                dset = TRIM(gname) // "/EZ1_Im"
                units = params%cpp%Eo
                call rsave_2d_array_to_hdf5(h5file_id, dset, &
                     units*F%E1Im_2DX%Z)
             endif
             
             if (params%SC_E) then
                dset = TRIM(gname) // "/dt_E_SC"
                attr = "Time step for self-consistent E calculation"
                call save_to_hdf5(h5file_id,dset,F%dt_E_SC,attr)

                dset = TRIM(gname) // "/Ip_exp"
                attr = "Scaling for self-consistent current density"
                call save_to_hdf5(h5file_id,dset,F%Ip_exp,attr)

                dset = TRIM(gname) // "/Ip0"
                attr = "Total RE current normalization"
                call save_to_hdf5(h5file_id,dset,F%Ip0,attr)
                
                dset = TRIM(gname) // "/PSIP_1D"
                attr_array(1) = "1D minor radial mesh for &
                     self-consistent fields"
                call save_1d_array_to_hdf5(h5file_id,dset,F%PSIP_1D,attr_array)
             end if
             
             if  (F%axisymmetric_fields.and. &
                  .not.(params%field_model(10:12).eq.'PSI')) then

                if (ALLOCATED(F%B_2D%R)) then
                   dset = TRIM(gname) // "/BR"
                   units = params%cpp%Bo
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*F%B_2D%R)
                end if

                if (ALLOCATED(F%B_2D%PHI)) then
                   dset = TRIM(gname) // "/BPHI"
                   units = params%cpp%Bo
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*F%B_2D%PHI)
                end if

                if (ALLOCATED(F%B_2D%Z)) then
                   dset = TRIM(gname) // "/BZ"
                   units = params%cpp%Bo
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*F%B_2D%Z)
                end if
                
                if (ALLOCATED(F%E_2D%R)) then
                   dset = TRIM(gname) // "/ER"
                   units = params%cpp%Eo
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*F%E_2D%R)

                   dset = TRIM(gname) // "/EPHI"
                   units = params%cpp%Eo
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*F%E_2D%PHI)

                   dset = TRIM(gname) // "/EZ"
                   units = params%cpp%Eo
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*F%E_2D%Z)
                endif

                if  (params%orbit_model(3:5).EQ.'pre') then

                   dset = TRIM(gname) // "/gradBR"
                   units = params%cpp%Bo/params%cpp%length
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*F%gradB_2D%R)

                   dset = TRIM(gname) // "/gradBPHI"
                   units = params%cpp%Bo/params%cpp%length
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*F%gradB_2D%PHI)

                   dset = TRIM(gname) // "/gradBZ"
                   units = params%cpp%Bo/params%cpp%length
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*F%gradB_2D%Z)

                   dset = TRIM(gname) // "/curlbR"
                   units = 1./params%cpp%length
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*F%curlb_2D%R)

                   dset = TRIM(gname) // "/curlbPHI"
                   units = 1./params%cpp%length
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*F%curlb_2D%PHI)

                   dset = TRIM(gname) // "/curlbZ"
                   units = 1./params%cpp%length
                   call rsave_2d_array_to_hdf5(h5file_id, dset, &
                        units*F%curlb_2D%Z)

                end if
                
             else if ((.not.F%axisymmetric_fields).and. &
                  .not.((params%field_model(10:12).eq.'PSI').or. &
                  (params%field_model(10:13).eq.'MARS').or. &
                  (params%field_model(10:14).eq.'AORSA'))) then

                dset = TRIM(gname) // "/BR"
                units = params%cpp%Bo
                call rsave_3d_array_to_hdf5(h5file_id, dset, &
                     units*F%B_3D%R)

                dset = TRIM(gname) // "/BPHI"
                units = params%cpp%Bo
                call rsave_3d_array_to_hdf5(h5file_id, dset, &
                     units*F%B_3D%PHI)

                dset = TRIM(gname) // "/BZ"
                units = params%cpp%Bo
                call rsave_3d_array_to_hdf5(h5file_id, dset, &
                     units*F%B_3D%Z)

                if  (params%orbit_model(3:5).EQ.'pre') then

                   dset = TRIM(gname) // "/gradBR"
                   units = params%cpp%Bo/params%cpp%length
                   call rsave_3d_array_to_hdf5(h5file_id, dset, &
                        units*F%gradB_3D%R)

                   dset = TRIM(gname) // "/gradBPHI"
                   units = params%cpp%Bo/params%cpp%length
                   call rsave_3d_array_to_hdf5(h5file_id, dset, &
                        units*F%gradB_3D%PHI)

                   dset = TRIM(gname) // "/gradBZ"
                   units = params%cpp%Bo/params%cpp%length
                   call rsave_3d_array_to_hdf5(h5file_id, dset, &
                        units*F%gradB_3D%Z)

                   dset = TRIM(gname) // "/curlbR"
                   units = 1./params%cpp%length
                   call rsave_3d_array_to_hdf5(h5file_id, dset, &
                        units*F%curlb_3D%R)

                   dset = TRIM(gname) // "/curlbPHI"
                   units = 1./params%cpp%length
                   call rsave_3d_array_to_hdf5(h5file_id, dset, &
                        units*F%curlb_3D%PHI)

                   dset = TRIM(gname) // "/curlbZ"
                   units = 1./params%cpp%length
                   call rsave_3d_array_to_hdf5(h5file_id, dset, &
                        units*F%curlb_3D%Z)

                end if
             end if
             
             DEALLOCATE(attr_array)
          else if (params%field_model .EQ. 'UNIFORM') then
             dset = TRIM(gname) // "/Bo"
             attr = "Magnetic field in T"
             call save_to_hdf5(h5file_id,dset,F%Bo*params%cpp%Bo,attr)

             dset = TRIM(gname) // "/Eo"
             attr = "Electric field in V/m"
             call save_to_hdf5(h5file_id,dset,F%Eo*params%cpp%Eo,attr)
          end if

          call h5gclose_f(group_id, h5error)


          ! Characteristic scales
          gname = "scales"
          call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

          dset = TRIM(gname) // "/t"
          attr = "Characteristic time in secs"
          call save_to_hdf5(h5file_id,dset,params%cpp%time,attr)

          dset = TRIM(gname) // "/m"
          attr = "Characteristic mass in kg"
          call save_to_hdf5(h5file_id,dset,params%cpp%mass,attr)

          dset = TRIM(gname) // "/q"
          attr = "Characteristic charge in Coulombs"
          call save_to_hdf5(h5file_id,dset,params%cpp%charge,attr)

          dset = TRIM(gname) // "/l"
          attr = "Characteristic length in m"
          call save_to_hdf5(h5file_id,dset,params%cpp%length,attr)

          dset = TRIM(gname) // "/v"
          attr = "Characteristic velocity in m"
          call save_to_hdf5(h5file_id,dset,params%cpp%velocity,attr)

          dset = TRIM(gname) // "/K"
          attr = "Characteristic kinetic energy in J"
          call save_to_hdf5(h5file_id,dset,params%cpp%energy,attr)

          dset = TRIM(gname) // "/n"
          attr = "Characteristic plasma density in m^-3"
          call save_to_hdf5(h5file_id,dset,params%cpp%density,attr)

          dset = TRIM(gname) // "/E"
          attr = "Characteristic electric field in V/m"
          call save_to_hdf5(h5file_id,dset,params%cpp%Eo,attr)

          dset = TRIM(gname) // "/B"
          attr = "Characteristic magnetic field in T"
          call save_to_hdf5(h5file_id,dset,params%cpp%Bo,attr)

          dset = TRIM(gname) // "/P"
          attr = "Characteristic pressure in Pa"
          call save_to_hdf5(h5file_id,dset,params%cpp%pressure,attr)

          dset = TRIM(gname) // "/T"
          attr = "Characteristic plasma temperature in J"
          call save_to_hdf5(h5file_id,dset,params%cpp%temperature,attr)

          call h5gclose_f(group_id, h5error)

          call h5fclose_f(h5file_id, h5error)
       end if

    end if
  end subroutine save_simulation_parameters


  subroutine save_simulation_outputs(params,spp,F)
    !! @note Subroutine that saves the electrons' variables specified in
    !! params::outputs_list to HDF5 files. @endnote
    TYPE(KORC_PARAMS), INTENT(IN) 				:: params
    !! Core KORC simulation parameters.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) 	:: spp
    !! An instance of KORC's derived type SPECIES containing all
    !! the information
    !! of different electron species. See [[korc_types]].
    TYPE(FIELDS), INTENT(IN)                 :: F
    CHARACTER(MAX_STRING_LENGTH) 				:: filename
    !! String containing the name of the HDF5 file.
    CHARACTER(MAX_STRING_LENGTH) 				:: gname
    !! String containing the group name of a set of KORC parameters.
    CHARACTER(MAX_STRING_LENGTH) 				:: subgname
    !! String containing the subgroup name of a set of KORC parameters.
    CHARACTER(MAX_STRING_LENGTH) 				:: dset
    !! Name of data set to be saved to file.
    INTEGER(HID_T) 						:: h5file_id
    !! HDF5 file identifier.
    INTEGER(HID_T) 						:: group_id
    !! HDF5 group identifier.
    INTEGER(HID_T) 						:: subgroup_id
    !! HDF5 subgroup identifier.
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 		:: dims
    !! Dimensions of data saved to HDF5 file.
    REAL(rp), DIMENSION(:), ALLOCATABLE 			:: rdata
    !! 1-D array of real data to be saved to HDF5 file.
    INTEGER, DIMENSION(:), ALLOCATABLE 				:: idata
    !!1-D array of integer data to be saved to HDF5 file.
    CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE       :: attr_array
    !! An 1-D array with attributes of 1-D real or integer arrays that are
    !! passed to KORC interfaces of HDF5 I/O subroutines.
    CHARACTER(MAX_STRING_LENGTH) 				:: attr
    !! A single attributes of real or integer data that is passed to KORC
    !! interfaces of HDF5 I/O subroutines.
    INTEGER 							:: h5error
    !!HDF5 error status.
    CHARACTER(19) 						:: tmp_str
    !!Temporary string used to manipulate various strings.
    REAL(rp) 						:: units
    !! Temporary variable used to add physical units to electrons' variables.
    INTEGER 						:: ss
    !! Electron species iterator.
    INTEGER 						:: jj
    !! Iterator for reading all the entried of params::outputs_list.
    LOGICAL 						:: object_exists
    !! Flag determining if a certain dataset is already present in
    !! the HDF5 output files.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE  ::YY
    !! Temporary variable get proper units on vars%Y(1,:) and vars%Y(3,:), which
    !! are lengths, while keeping vars%Y(2,:), which is an angle

    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,'("Saving snapshot: ",I15)') &
            params%it/(params%t_skip)
       !write(output_unit_write,*) 'it',params%it,'t_skip',params%t_skip,'t_SC',params%t_it_SC
       
    end if

    if (SIZE(params%outputs_list).GT.1_idef) then
       write(tmp_str,'(I18)') params%mpi_params%rank
       filename = TRIM(params%path_to_outputs) // "file_" &
            // TRIM(ADJUSTL(tmp_str)) // ".h5"
       call h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, h5file_id, h5error)

       ! Create group 'it'
       write(tmp_str,'(I18)') params%it
       gname = TRIM(ADJUSTL(tmp_str))
       call h5lexists_f(h5file_id,TRIM(gname),object_exists,h5error)

       if (.NOT.object_exists) then ! Check if group does exist.
          call h5gcreate_f(h5file_id, TRIM(gname), group_id, h5error)

          dset = TRIM(gname) // "/time"
          attr = "Simulation time in secs"
          call save_to_hdf5(h5file_id,dset,params%init_time*params%cpp%time &
               + REAL(params%it,rp)*params%dt*params%cpp%time,attr)
          
          do ss=1_idef,params%num_species

             write(tmp_str,'(I18)') ss
             subgname = "spp_" // TRIM(ADJUSTL(tmp_str))
             call h5gcreate_f(group_id, TRIM(subgname), subgroup_id, h5error)

             do jj=1_idef,SIZE(params%outputs_list)
                SELECT CASE (TRIM(params%outputs_list(jj)))
                CASE ('X')
                   dset = "X"
                   units = params%cpp%length
                   call rsave_2d_array_to_hdf5(subgroup_id, dset, &
                        units*spp(ss)%vars%X)
                CASE ('Y')
                   dset = "Y"
                   units = params%cpp%length

                   YY=spp(ss)%vars%Y
                   YY(:,1)=units*YY(:,1)
                   YY(:,3)=units*YY(:,3)

                   call rsave_2d_array_to_hdf5(subgroup_id, dset, &
                        YY)

                   DEALLOCATE(YY)

                CASE('V')
                   dset = "V"
                   if (params%orbit_model(1:2).eq.'FO') then
                      units = params%cpp%velocity
                      call rsave_2d_array_to_hdf5(subgroup_id, dset, &
                           units*spp(ss)%vars%V)
                   else if (params%orbit_model(1:2).eq.'GC') then
                      YY=spp(ss)%vars%V

                      YY(:,1)=YY(:,1)*params%cpp%mass*params%cpp%velocity
                      YY(:,2)=YY(:,2)*params%cpp%mass* &
                           (params%cpp%velocity)**2/params%cpp%Bo

                      call rsave_2d_array_to_hdf5(subgroup_id, dset, &
                           YY)
                      DEALLOCATE(YY)

                   end if
                CASE('RHS')
                   dset = "RHS"
                   YY=spp(ss)%vars%RHS

                   units = params%cpp%length/params%cpp%time
                   YY(:,1)=YY(:,1)*units
                   YY(:,2)=YY(:,2)*units
                   YY(:,3)=YY(:,3)*units
                   units = params%cpp%mass*params%cpp%velocity/params%cpp%time
                   YY(:,4)=YY(:,4)*units
                   YY(:,5)=YY(:,5)*params%cpp%mass* &
                           (params%cpp%velocity)**2/params%cpp%Bo/params%cpp%time

                   call rsave_2d_array_to_hdf5(subgroup_id, dset, &
                        YY)
                   DEALLOCATE(YY)
                   
                CASE('Rgc')
                   dset = "Rgc"
                   units = params%cpp%length
                   call rsave_2d_array_to_hdf5(subgroup_id, dset, &
                        units*spp(ss)%vars%Rgc)
                CASE('g')
                   dset = "g"
                   call save_1d_array_to_hdf5(subgroup_id, dset, &
                        spp(ss)%vars%g)
                CASE('eta')
                   dset = "eta"
                   call save_1d_array_to_hdf5(subgroup_id, dset, &
                        spp(ss)%vars%eta)
                CASE('mu')
                   dset = "mu"
                   units = params%cpp%mass*params%cpp%velocity**2/params%cpp%Bo
                   call save_1d_array_to_hdf5(subgroup_id, dset, &
                        units*spp(ss)%vars%mu)
                CASE('Prad')
                   dset = "Prad"
                   units = params%cpp%mass*(params%cpp%velocity**3)/ &
                        params%cpp%length
                   call save_1d_array_to_hdf5(subgroup_id, dset, &
                        units*spp(ss)%vars%Prad)
                CASE('Pin')
                   dset = "Pin"
                   units = params%cpp%mass*(params%cpp%velocity**3)/ &
                        params%cpp%length
                   call save_1d_array_to_hdf5(subgroup_id, dset, &
                        units*spp(ss)%vars%Pin)
                CASE('flagCon')
                   dset = "flagCon"
                   call save_1d_array_to_hdf5(subgroup_id,dset, &
                        INT(spp(ss)%vars%flagCon,idef))
                CASE('flagCol')
                   dset = "flagCol"
                   call save_1d_array_to_hdf5(subgroup_id,dset, &
                        INT(spp(ss)%vars%flagCol,idef))
                CASE('flagRE')
                   dset = "flagRE"
                   call save_1d_array_to_hdf5(subgroup_id,dset, &
                        INT(spp(ss)%vars%flagRE,idef))
                CASE('B')
                   dset = "B"
                   units = params%cpp%Bo
                   call rsave_2d_array_to_hdf5(subgroup_id, dset, &
                        units*spp(ss)%vars%B)
                CASE('gradB')
                   if (params%orbit_model(3:5).eq.'pre') then
                      dset = "gradB"
                      units = params%cpp%Bo/params%cpp%length
                      call rsave_2d_array_to_hdf5(subgroup_id, dset, &
                           units*spp(ss)%vars%gradB)
                   end if
                CASE('curlb')
                   if (params%orbit_model(3:5).eq.'pre') then
                      dset = "curlb"
                      units = 1./params%cpp%length
                      call rsave_2d_array_to_hdf5(subgroup_id, dset, &
                           units*spp(ss)%vars%curlb)
                   end if
                CASE('E')
                   dset = "E"
                   units = params%cpp%Eo
                   call rsave_2d_array_to_hdf5(subgroup_id, dset, &
                        units*spp(ss)%vars%E)
                   
                CASE('PSIp')
                   
                   dset = "PSIp"
                   if (.not.params%field_model.eq.'M3D_C1') then
                      units = params%cpp%Bo*params%cpp%length**2
                   else
                      units = 1._rp
                   end if
                   call save_1d_array_to_hdf5(subgroup_id, dset, &
                        units*spp(ss)%vars%PSI_P)


                   
                CASE('AUX')
                   dset = "AUX"
                   call save_1d_array_to_hdf5(subgroup_id, dset, &
                        spp(ss)%vars%AUX)
                CASE ('ne')
                   dset = "ne"
                   units = params%cpp%density
                   call save_1d_array_to_hdf5(subgroup_id, dset, &
                        units*spp(ss)%vars%ne)
                CASE ('nimp')
                   dset = "nimp"
                   units = params%cpp%density
                   call save_2d_array_to_hdf5(subgroup_id, dset, &
                        units*spp(ss)%vars%nimp)
                CASE ('Te')
                   dset = "Te"
                   units = params%cpp%temperature
                   call save_1d_array_to_hdf5(subgroup_id, dset, &
                        units*spp(ss)%vars%Te/C_E)
                CASE ('Zeff')
                   dset = "Zeff"
                   call save_1d_array_to_hdf5(subgroup_id, dset, &
                        spp(ss)%vars%Zeff)                   
                   
                CASE ('J_SC')

                   dset = "J_SC"
                   if (params%SC_E) then
                      call save_1d_array_to_hdf5(subgroup_id, dset, &
                           F%J1_SC_1D%PHI)
                   end if


                CASE ('A_SC')

                   dset = "A_SC"
                   if (params%SC_E) then
                      call save_1d_array_to_hdf5(subgroup_id, dset, &
                           F%A1_SC_1D%PHI)
                   end if

                CASE ('E_SC')

                   dset = "E_SC"
                   units = params%cpp%Eo
                   if (params%SC_E) then
                      call save_1d_array_to_hdf5(subgroup_id, dset, &
                           units*F%E_SC_1D%PHI)
                   end if

                CASE DEFAULT

                   
                END SELECT
             end do

             call h5gclose_f(subgroup_id, h5error)
          end do

          call h5gclose_f(group_id, h5error)
       end if ! Check if group does exist.

       call h5fclose_f(h5file_id, h5error)
    end if
  end subroutine save_simulation_outputs



  subroutine save_restart_variables(params,spp,F)
    !! @note Subroutine that saves all the variables that KORC needs for
    !! restarting a simulation. These variables are saved to "restart_file.h5".
    TYPE(KORC_PARAMS), INTENT(IN) 				:: params
    !! params Core KORC simulation parameters.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(IN) 	:: spp
    !! An instance of KORC's derived type SPECIES containing
    !! all the information of different electron species. See [[korc_types]].
    TYPE(FIELDS), INTENT(IN)      :: F
    REAL(rp), DIMENSION(:), ALLOCATABLE :: send_buffer_rp, receive_buffer_rp
    !! Temporary buffer to be used by MPI to gather different electrons'
    !! variables.
    !! Temporary buffer to be used by MPI to gather different electrons'
    !! variables.
    INTEGER(is), DIMENSION(:), ALLOCATABLE :: send_buffer_is, receive_buffer_is
    !! Temporary buffer to be used by MPI to gather different electrons'
    !! variables.
    !! Temporary buffer to be used by MPI to gather different electrons'
    !! variables.
    REAL(rp), DIMENSION(:,:), ALLOCATABLE 			:: X
    REAL(rp), DIMENSION(:,:), ALLOCATABLE 			:: V
    REAL(rp), DIMENSION(:), ALLOCATABLE 			:: g
    REAL(rp), DIMENSION(:), ALLOCATABLE 			:: J0_SC
    REAL(rp), DIMENSION(:), ALLOCATABLE 			:: J1_SC
    REAL(rp), DIMENSION(:), ALLOCATABLE 			:: J2_SC
    REAL(rp), DIMENSION(:), ALLOCATABLE 			:: J3_SC
    REAL(rp), DIMENSION(:), ALLOCATABLE 			:: E_SC
    INTEGER(is), DIMENSION(:), ALLOCATABLE :: flagCon,flagCol,flagRE
    CHARACTER(MAX_STRING_LENGTH) 				:: filename
    !! String containing the name of the HDF5 file.
    CHARACTER(MAX_STRING_LENGTH) 				:: gname
    !! String containing the group name of a set of KORC parameters.
    CHARACTER(MAX_STRING_LENGTH) 				:: subgname
    CHARACTER(MAX_STRING_LENGTH) 				:: dset
    !! Name of data set to be saved to file.
    INTEGER(HID_T) 						:: h5file_id
    !! HDF5 file identifier.
    INTEGER(HID_T) 						:: group_id
    !! HDF5 group identifier.
    INTEGER(HID_T) 						:: subgroup_id
    !! HDF5 subgroup identifier.
    INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE 		:: dims
    !!  Dimensions of data saved to HDF5 file.
    REAL(rp), DIMENSION(:), ALLOCATABLE 			:: rdata
    !! 1-D array of real data to be saved to HDF5 file.
    INTEGER, DIMENSION(:), ALLOCATABLE 				:: idata
    !! 1-D array of integer data to be saved to HDF5 file.
    CHARACTER(MAX_STRING_LENGTH), DIMENSION(:), ALLOCATABLE :: attr_array
    !! An 1-D array with attributes of 1-D real or integer arrays that
    !! are passed to KORC interfaces of HDF5 I/O subroutines.
    CHARACTER(MAX_STRING_LENGTH) 				:: attr
    !! A single attributes of real or integer data that is passed to KORC
    !! interfaces of HDF5 I/O subroutines.
    INTEGER 							:: h5error
    !! HDF5 error status.
    CHARACTER(19) 						:: tmp_str
    !! Temporary string used to manipulate various strings.
    REAL(rp) 							:: units
    !! Temporary variable used to add physical units to restart variables.
    INTEGER 							:: ss,jj
    !! Electron species iterator.
    !! Iterator for reading all the entried of params::outputs_list.
    INTEGER 							:: mpierr
    !! MPI error status.
    INTEGER 					:: numel_send, numel_receive
    !! Variable used by MPI to count the amount of data sent by each MPI
    !! procces.
    !! Variable used by MPI to count the amount of data received by the main
    !! MPI procces.


!    if ( MODULO(params%it,params%restart_output_cadence) .EQ. 0_ip ) then 
    if (params%mpi_params%rank.EQ.0_idef) then

       write(output_unit_write,'("Saving restart: ",I15)') &
            params%it/(params%t_skip*params%t_it_SC)

       filename = TRIM(params%path_to_outputs) // "restart_file.h5"
       call h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, h5file_id, h5error)

       dset = "it"
       attr = "Iteration"
       call save_to_hdf5(h5file_id,dset,params%it,attr)

       dset = "time"
       attr = "Current simulation time in secs"
       call save_to_hdf5(h5file_id,dset,params%init_time*params%cpp%time &
            + REAL(params%it,rp)*params%dt*params%cpp%time,attr)

       dset = "simulation_time"
       attr = "Total simulation time in secs"
       call save_to_hdf5(h5file_id,dset,params%simulation_time* &
            params%cpp%time,attr)

       dset = "snapshot_frequency"
       attr = "Snapshot frequency in secs"
       call save_to_hdf5(h5file_id,dset,params%snapshot_frequency* &
            params%cpp%time,attr)

       dset = "dt"
       attr = "Time step in secs"
       call save_to_hdf5(h5file_id,dset,params%dt*params%cpp%time,attr)

       dset = "t_steps"
       attr = "Time steps in simulation"
       call save_to_hdf5(h5file_id,dset,params%t_steps,attr)

       dset = "output_cadence"
       attr = "Output cadence"
       call save_to_hdf5(h5file_id,dset,params%output_cadence,attr)

       dset = "restart_output_cadence"
       attr = "Restart output cadence"
       call save_to_hdf5(h5file_id,dset,params%restart_output_cadence,attr)

       dset = "num_snapshots"
       attr = "Number of snapshots in time for saving simulation variables"
       call save_to_hdf5(h5file_id,dset,params%num_snapshots,attr)

       if (F%ReInterp_2x1t) then
          dset = "ind_2x1t"
          attr = "ReInterp_2x1t iteration"
          call save_to_hdf5(h5file_id,dset,F%ind_2x1t,attr)
       end if

       dset = "num_mpi"
       attr = "Number of mpi nodes used in simulation"
       call save_to_hdf5(h5file_id,dset,params%mpi_params%nmpi,attr)
       
    end if

    do ss=1_idef,params%num_species
       numel_send = 3_idef*spp(ss)%ppp
       numel_receive = 3_idef*spp(ss)%ppp*params%mpi_params%nmpi

       if (params%mpi_params%rank.EQ.0_idef) then
          ALLOCATE(X(spp(ss)%ppp*params%mpi_params%nmpi,3))
          ALLOCATE(V(spp(ss)%ppp*params%mpi_params%nmpi,3))
          ALLOCATE(g(spp(ss)%ppp*params%mpi_params%nmpi))
          ALLOCATE(flagCon(spp(ss)%ppp*params%mpi_params%nmpi))
          ALLOCATE(flagCol(spp(ss)%ppp*params%mpi_params%nmpi))
          ALLOCATE(flagRE(spp(ss)%ppp*params%mpi_params%nmpi))
       end if

       ALLOCATE(send_buffer_rp(numel_send))
       ALLOCATE(receive_buffer_rp(numel_receive))

       if (params%orbit_model(1:2).EQ.'FO') then             
          send_buffer_rp = RESHAPE(spp(ss)%vars%X,(/numel_send/))
       else if (params%orbit_model(1:2).EQ.'GC') then
          send_buffer_rp = RESHAPE(spp(ss)%vars%Y,(/numel_send/))
       end if
       receive_buffer_rp = 0.0_rp
       CALL MPI_GATHER(send_buffer_rp,numel_send,MPI_REAL8, &
            receive_buffer_rp,numel_send,MPI_REAL8,0,MPI_COMM_WORLD, &
            mpierr)
       if (params%mpi_params%rank.EQ.0_idef) then
          X = RESHAPE(receive_buffer_rp,(/spp(ss)%ppp* &
               params%mpi_params%nmpi,3/))
       end if

       send_buffer_rp = RESHAPE(spp(ss)%vars%V,(/numel_send/))
       receive_buffer_rp = 0.0_rp
       CALL MPI_GATHER(send_buffer_rp,numel_send,MPI_REAL8, &
            receive_buffer_rp,numel_send,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
       if (params%mpi_params%rank.EQ.0_idef) then
          V = RESHAPE(receive_buffer_rp,(/spp(ss)%ppp* &
               params%mpi_params%nmpi,3/))
       end if

       DEALLOCATE(send_buffer_rp)
       DEALLOCATE(receive_buffer_rp)

       numel_send = spp(ss)%ppp
       numel_receive = spp(ss)%ppp*params%mpi_params%nmpi

       ALLOCATE(send_buffer_is(numel_send))
       ALLOCATE(receive_buffer_is(numel_receive))

       send_buffer_is = spp(ss)%vars%flagRE
       receive_buffer_is = 0_is
       CALL MPI_GATHER(send_buffer_is,numel_send,MPI_INTEGER1, &
            receive_buffer_is,numel_send,&
            MPI_INTEGER1,0,MPI_COMM_WORLD,mpierr)
       if (params%mpi_params%rank.EQ.0_idef) then
          flagRE = receive_buffer_is
       end if
       
       send_buffer_is = spp(ss)%vars%flagCon
       receive_buffer_is = 0_is
       CALL MPI_GATHER(send_buffer_is,numel_send,MPI_INTEGER1, &
            receive_buffer_is,numel_send,&
            MPI_INTEGER1,0,MPI_COMM_WORLD,mpierr)
       if (params%mpi_params%rank.EQ.0_idef) then
          flagCon = receive_buffer_is
       end if

       send_buffer_is = spp(ss)%vars%flagCol
       receive_buffer_is = 0_is
       CALL MPI_GATHER(send_buffer_is,numel_send,MPI_INTEGER1, &
            receive_buffer_is,numel_send,&
            MPI_INTEGER1,0,MPI_COMM_WORLD,mpierr)
       if (params%mpi_params%rank.EQ.0_idef) then
          flagCol = receive_buffer_is
       end if
       
       DEALLOCATE(send_buffer_is)
       DEALLOCATE(receive_buffer_is)

       ALLOCATE(send_buffer_rp(numel_send))
       ALLOCATE(receive_buffer_rp(numel_receive))

       send_buffer_rp = spp(ss)%vars%g
       receive_buffer_rp = 0_rp
       CALL MPI_GATHER(send_buffer_rp,numel_send,MPI_REAL8, &
            receive_buffer_rp,numel_send,&
            MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
       if (params%mpi_params%rank.EQ.0_idef) then
          g = receive_buffer_rp
       end if

       DEALLOCATE(send_buffer_rp)
       DEALLOCATE(receive_buffer_rp)

       if (params%mpi_params%rank.EQ.0_idef) then
          write(tmp_str,'(I18)') ss
          subgname = "spp_" // TRIM(ADJUSTL(tmp_str))
          call h5gcreate_f(h5file_id, TRIM(subgname), group_id, h5error)

          dset = "X"
          call rsave_2d_array_to_hdf5(group_id, dset, X)

          dset = "V"
          call rsave_2d_array_to_hdf5(group_id, dset, V)

          dset = "flagCon"
          call save_1d_array_to_hdf5(group_id,dset, INT(flagCon,idef))

          dset = "flagCol"
          call save_1d_array_to_hdf5(group_id,dset, INT(flagCol,idef))

          dset = "flagRE"
          call save_1d_array_to_hdf5(group_id,dset, INT(flagRE,idef))


          dset = "g"
          call save_1d_array_to_hdf5(group_id,dset, g)

          if (params%SC_E) then

             ALLOCATE(J0_SC(F%dim_1D))
             ALLOCATE(J1_SC(F%dim_1D))
             ALLOCATE(J2_SC(F%dim_1D))
             ALLOCATE(J3_SC(F%dim_1D))
             ALLOCATE(E_SC(F%dim_1D))
             J0_SC=F%J1_SC_1D%PHI/F%Ip0
             J1_SC=F%J1_SC_1D%PHI
             J2_SC=F%J2_SC_1D%PHI
             J3_SC=F%J3_SC_1D%PHI
             E_SC=F%E_SC_1D%PHI

             dset = "J0_SC"
             call save_1d_array_to_hdf5(group_id,dset,J0_SC)
             dset = "J1_SC"
             call save_1d_array_to_hdf5(group_id,dset,J1_SC)
             dset = "J2_SC"
             call save_1d_array_to_hdf5(group_id,dset,J2_SC)
             dset = "J3_SC"
             call save_1d_array_to_hdf5(group_id,dset,J3_SC)
             dset = "E_SC"
             call save_1d_array_to_hdf5(group_id,dset,E_SC)
             
             DEALLOCATE(J0_SC)
             DEALLOCATE(J1_SC)
             DEALLOCATE(J2_SC)
             DEALLOCATE(J3_SC)
             DEALLOCATE(E_SC)
          end if
          
          call h5gclose_f(group_id, h5error)
       end if

       if (params%mpi_params%rank.EQ.0_idef) then
          DEALLOCATE(X)
          DEALLOCATE(V)
          DEALLOCATE(g)
          DEALLOCATE(flagCon)
          DEALLOCATE(flagCol)
          DEALLOCATE(flagRE)
       end if
    end do

    
    if (params%mpi_params%rank.EQ.0_idef) then
       call h5fclose_f(h5file_id, h5error)
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    
!    end if
  end subroutine save_restart_variables

  ! * * * * * * * * * * * * * * * * * * * * * * * * * !
  ! * * * SUBROUTINES FOR RESTARTING SIMULATION * * * !
  ! * * * * * * * * * * * * * * * * * * * * * * * * * !

  subroutine load_time_stepping_params(params)
    !! @note Subroutine that loads KORC parameters that control the time
    !! stepping in [[main]].    
    TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
    !! Core KORC simulation parameters.
    CHARACTER(MAX_STRING_LENGTH) 		:: filename
    !! String containing the name of the HDF5 file.
    CHARACTER(MAX_STRING_LENGTH) 		:: dset
    !! Name of data set to be read from file.
    INTEGER(HID_T) 						:: h5file_id
    !! HDF5 file identifier.
    REAL(KIND=8) 						:: real_number
    !! A temporary real number.
    CHARACTER(19) 						:: tmp_str
    !! Temporary string used to manipulate various strings.
    INTEGER 							:: h5error
    !! HDF5 error status.
    INTEGER 							:: mpierr
    !!  MPI error status.
    INTEGER 							:: ss
    !! Electron species iterator.

    if (params%mpi_params%rank.EQ.0_idef) then
       filename = TRIM(params%path_to_outputs) // "restart_file.h5"
       call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
       if (h5error .EQ. -1) then
          write(output_unit_write,'("KORC ERROR: Something went wrong in: &
               &load_particles_ic --> h5fopen_f")')
          call KORC_ABORT(14)
       end if

       dset = "/it"
       call load_from_hdf5(h5file_id,dset,real_number)
       params%ito = INT(real_number,ip) + 1_ip

       dset = "/dt"
       call load_from_hdf5(h5file_id,dset,params%dt)

       dset = "/t_steps"
       call load_from_hdf5(h5file_id,dset,real_number)
       params%t_steps = INT(real_number,ip)

       dset = "/simulation_time"
       call load_from_hdf5(h5file_id,dset,params%simulation_time)

       dset = "/snapshot_frequency"
       call load_from_hdf5(h5file_id,dset,params%snapshot_frequency)

       dset = "/output_cadence"
       call load_from_hdf5(h5file_id,dset,real_number)
       params%output_cadence = INT(real_number,ip)

       dset = "/restart_output_cadence"
       call load_from_hdf5(h5file_id,dset,real_number)
       params%restart_output_cadence = INT(real_number,ip)

       dset = "/num_snapshots"
       call load_from_hdf5(h5file_id,dset,real_number)
       params%num_snapshots = INT(real_number,ip)       

       call h5fclose_f(h5file_id, h5error)
    end if

    CALL MPI_BCAST(params%ito,1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)

    CALL MPI_BCAST(params%dt,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

    CALL MPI_BCAST(params%t_steps,1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)

    CALL MPI_BCAST(params%simulation_time,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

    CALL MPI_BCAST(params%snapshot_frequency,1,MPI_REAL8,0,MPI_COMM_WORLD, &
         mpierr)

    CALL MPI_BCAST(params%output_cadence,1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)

    CALL MPI_BCAST(params%restart_output_cadence,1,MPI_INTEGER8,0, &
         MPI_COMM_WORLD,mpierr)

    CALL MPI_BCAST(params%num_snapshots,1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
  end subroutine load_time_stepping_params

  subroutine load_prev_time(params)
    !! @note Subroutine that loads KORC parameters that control the time
    !! stepping in [[main]].    
    TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
    !! Core KORC simulation parameters.
    CHARACTER(MAX_STRING_LENGTH) 		:: filename
    !! String containing the name of the HDF5 file.
    CHARACTER(MAX_STRING_LENGTH) 		:: dset
    !! Name of data set to be read from file.
    INTEGER(HID_T) 						:: h5file_id
    !! HDF5 file identifier.
    REAL(KIND=8) 						:: real_number
    !! A temporary real number.
    CHARACTER(19) 						:: tmp_str
    !! Temporary string used to manipulate various strings.
    INTEGER 							:: h5error
    !! HDF5 error status.
    INTEGER 							:: mpierr
    !!  MPI error status.
    INTEGER 							:: ss
    !! Electron species iterator.

    if (params%mpi_params%rank.EQ.0_idef) then
       filename = TRIM(params%path_to_outputs) // "restart_file.h5"
       call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
       if (h5error .EQ. -1) then
          write(output_unit_write,'("KORC ERROR: Something went wrong in: &
               &load_prev_time --> h5fopen_f")')
          call KORC_ABORT(14)
       end if

       dset = "/time"
       call load_from_hdf5(h5file_id,dset,params%init_time)      

       dset = "/num_mpi"
       call load_from_hdf5(h5file_id,dset,real_number)
       params%mpi_params%nmpi_prev = INT(real_number,ip)
       
       call h5fclose_f(h5file_id, h5error)
    end if
    
    CALL MPI_BCAST(params%init_time,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    
    CALL MPI_BCAST(params%mpi_params%nmpi_prev,1,MPI_INTEGER8,0, &
         MPI_COMM_WORLD,mpierr)

    ! Not sure why, but params%mpi_params%rank is reset to zero in the above
    ! call to MPI_BCAST (but not other mpi_params values). Added the following
    ! MPI_COMM_RANK to reintialize
    call MPI_COMM_RANK(MPI_COMM_WORLD, params%mpi_params%rank, mpierr)
    

  end subroutine load_prev_time

    subroutine load_prev_iter(params)
    !! @note Subroutine that loads KORC parameters that control the time
    !! stepping in [[main]].    
    TYPE(KORC_PARAMS), INTENT(INOUT) 	:: params
    !! Core KORC simulation parameters.
    CHARACTER(MAX_STRING_LENGTH) 		:: filename
    !! String containing the name of the HDF5 file.
    CHARACTER(MAX_STRING_LENGTH) 		:: dset
    !! Name of data set to be read from file.
    INTEGER(HID_T) 						:: h5file_id
    !! HDF5 file identifier.
    REAL(KIND=8) 						:: real_number
    !! A temporary real number.
    CHARACTER(19) 						:: tmp_str
    !! Temporary string used to manipulate various strings.
    INTEGER 							:: h5error
    !! HDF5 error status.
    INTEGER 							:: mpierr
    !!  MPI error status.
    INTEGER 							:: ss
    !! Electron species iterator.

    if (params%mpi_params%rank.EQ.0_idef) then
       filename = TRIM(params%path_to_outputs) // "restart_file.h5"
       call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
       if (h5error .EQ. -1) then
          write(output_unit_write,'("KORC ERROR: Something went wrong in: &
               &load_prev_iter --> h5fopen_f")')
          call KORC_ABORT(14)
       end if

       dset = "/ind_2x1t"
       call load_from_hdf5(h5file_id,dset,real_number)
       params%prev_iter_2x1t = INT(real_number,ip)

       call h5fclose_f(h5file_id, h5error)
    end if

    CALL MPI_BCAST(params%prev_iter_2x1t,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

  end subroutine load_prev_iter


  subroutine load_particles_ic(params,spp,F)
    !! @note Subroutine that loads all the electrons' data from
    !! "restart_file.h5" to restart a simulation.
    TYPE(KORC_PARAMS), INTENT(INOUT) 			:: params
    !! Core KORC simulation parameters.
    TYPE(SPECIES), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: spp
    !! An instance of KORC's derived type SPECIES containing all the
    !! information of different electron species. See korc_types.f90.
    TYPE(FIELDS), INTENT(INOUT) 			:: F
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: X_send_buffer,X_send_buffer_tmp,X_send_buffer_tmp1
    !! Temporary buffer used by MPI for scattering the electrons' position
    !! to different MPI processes.
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: X_receive_buffer
    !! Temporary buffer used by MPI for scattering the electrons' position
    !! among MPI processes.
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: V_send_buffer,V_send_buffer_tmp,V_send_buffer_tmp1
    !! Temporary buffer used by MPI for scattering the electrons' velocity
    !! among MPI processes.
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: V1_send_buffer_tmp,V2_send_buffer_tmp,V3_send_buffer_tmp
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: V1_send_buffer_tmp1,V2_send_buffer_tmp1,V3_send_buffer_tmp1
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: V_receive_buffer
    !! Temporary buffer used by MPI for scattering the electrons' velocity
    !! among MPI processes.
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: AUX_send_buffer
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: AUX_send_buffer_tmp
    !!  Temporary buffer used by MPI to scatter various electrons' variables
    !! among MPI processes.
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: AUX_receive_buffer
    !! Temporary buffer used by MPI to scatter various electrons' variables
    !! among MPI processes.
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: JSC0_buffer
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: JSC1_buffer
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: JSC2_buffer
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: JSC3_buffer
    REAL(rp), DIMENSION(:), ALLOCATABLE 		:: ESC_buffer


    CHARACTER(MAX_STRING_LENGTH) 			:: filename
    !! String containing the name of the HDF5 file.
    CHARACTER(MAX_STRING_LENGTH) 			:: dset
    !! Name of data set to be saved to file.
    INTEGER(HID_T) 					:: h5file_id
    !! HDF5 file identifier.
    CHARACTER(19) 					:: tmp_str
    !! Temporary string used to manipulate various strings.
    INTEGER 						:: h5error
    !! HDF5 error status.
    INTEGER 						:: mpierr
    !! Electron species iterator.
    INTEGER 						:: ss,ii,jj
    !! MPI error status.
    INTEGER  :: recieve_num,send_num,nmpi_ratio
    
    do ss=1_idef,params%num_species

       !write(6,*) params%mpi_params%nmpi,params%mpi_params%nmpi_prev
       
       if (mod(spp(ss)%ppp,params%mpi_params%nmpi/ &
            params%mpi_params%nmpi_prev).ne.0) then
          write(output_unit_write,'("ppp must be divisible by factor increase in MPI nodes")')
          call KORC_ABORT(14)
       endif

       nmpi_ratio=params%mpi_params%nmpi/params%mpi_params%nmpi_prev
       
       send_num=spp(ss)%ppp*params%mpi_params%nmpi_prev
       recieve_num=spp(ss)%ppp/nmpi_ratio
       
       ALLOCATE(X_send_buffer(3*send_num))
       ALLOCATE(X_send_buffer_tmp(3*send_num))
       ALLOCATE(X_send_buffer_tmp1(3*send_num))
       ALLOCATE(X_receive_buffer(3*recieve_num))

       ALLOCATE(V_send_buffer(3*send_num))
       ALLOCATE(V_send_buffer_tmp(3*send_num))
       ALLOCATE(V_send_buffer_tmp1(3*send_num))
       ALLOCATE(V_receive_buffer(3*recieve_num))

       ALLOCATE(V1_send_buffer_tmp(send_num))
       ALLOCATE(V2_send_buffer_tmp(send_num))
       ALLOCATE(V3_send_buffer_tmp(send_num))
       ALLOCATE(V1_send_buffer_tmp1(send_num))
       ALLOCATE(V2_send_buffer_tmp1(send_num))
       ALLOCATE(V3_send_buffer_tmp1(send_num))
       
       ALLOCATE(AUX_send_buffer(send_num))
       ALLOCATE(AUX_send_buffer_tmp(send_num))
       ALLOCATE(AUX_receive_buffer(recieve_num))
       
       if (params%mpi_params%rank.EQ.0_idef) then
          
          filename = TRIM(params%path_to_outputs) // "restart_file.h5"
          call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
          if (h5error .EQ. -1) then
             write(output_unit_write,'("KORC ERROR: Something went wrong in: &
                  &load_particles_ic --> h5fopen_f")')
             call KORC_ABORT(14)
          end if

          write(tmp_str,'(I18)') ss

          dset = "/spp_" // TRIM(ADJUSTL(tmp_str)) // "/X"
          call load_array_from_hdf5(h5file_id,dset,X_send_buffer)

          call h5fclose_f(h5file_id, h5error)

          !write(6,*) 'shape of X_send_buffer',shape(X_send_buffer)
          !write(6,*) 'X_send_buffer',X_send_buffer
          !write(6,*) 'recieve_num',recieve_num

#if DBG_CHECK
          write(6,*) 'X_send_buffer',X_send_buffer*params%cpp%length
#endif
          
          if (params%load_balance) then
             do ii=0,params%mpi_params%nmpi_prev-1

                V1_send_buffer_tmp(ii*spp(ss)%ppp+1:(ii+1)*spp(ss)%ppp)= &
                     X_send_buffer(3*ii*spp(ss)%ppp+1:(3*ii+1)*spp(ss)%ppp)
                V2_send_buffer_tmp(ii*spp(ss)%ppp+1:(ii+1)*spp(ss)%ppp)= &
                     X_send_buffer((3*ii+1)*spp(ss)%ppp+1:(3*ii+2)*spp(ss)%ppp)
                V3_send_buffer_tmp(ii*spp(ss)%ppp+1:(ii+1)*spp(ss)%ppp)= &
                     X_send_buffer((3*ii+2)*spp(ss)%ppp+1:(3*ii+3)*spp(ss)%ppp)
                
             end do
             do ii=1,params%mpi_params%nmpi_prev*spp(ss)%ppp

                V1_send_buffer_tmp1((mod(ii-1,params%mpi_params%nmpi))* &
                     spp(ss)%ppp/nmpi_ratio+(ii-1)/params%mpi_params%nmpi+1) = &
                     V1_send_buffer_tmp(ii)
                V2_send_buffer_tmp1((mod(ii-1,params%mpi_params%nmpi))* &
                     spp(ss)%ppp/nmpi_ratio+(ii-1)/params%mpi_params%nmpi+1) = &
                     V2_send_buffer_tmp(ii)
                V3_send_buffer_tmp1((mod(ii-1,params%mpi_params%nmpi))* &
                     spp(ss)%ppp/nmpi_ratio+(ii-1)/params%mpi_params%nmpi+1) = &
                     V3_send_buffer_tmp(ii)
                                                
             end do
             do ii=0,params%mpi_params%nmpi_prev-1

                X_send_buffer_tmp1(3*ii*spp(ss)%ppp+1:(3*ii+1)*spp(ss)%ppp)= &
                     V1_send_buffer_tmp1(ii*spp(ss)%ppp+1:(ii+1)*spp(ss)%ppp)
                X_send_buffer_tmp1((3*ii+1)*spp(ss)%ppp+1:(3*ii+2)*spp(ss)%ppp)= &
                     V2_send_buffer_tmp1(ii*spp(ss)%ppp+1:(ii+1)*spp(ss)%ppp)
                X_send_buffer_tmp1((3*ii+2)*spp(ss)%ppp+1:(3*ii+3)*spp(ss)%ppp)= &
                     V3_send_buffer_tmp1(ii*spp(ss)%ppp+1:(ii+1)*spp(ss)%ppp)


             end do
          else
             X_send_buffer_tmp1=X_send_buffer
          end if
          
          do jj=0_idef,params%mpi_params%nmpi_prev-1_idef
             do ii=0_idef,nmpi_ratio-1_idef

                X_send_buffer_tmp(jj*spp(ss)%ppp*3+recieve_num*(3*ii)+1: &
                     jj*spp(ss)%ppp*3+recieve_num*(3*ii+1))= &
                     X_send_buffer_tmp1(jj*spp(ss)%ppp*3+recieve_num*ii+1: &
                     jj*spp(ss)%ppp*3+recieve_num*(ii+1))

                X_send_buffer_tmp(jj*spp(ss)%ppp*3+recieve_num*(3*ii+1)+1: &
                     jj*spp(ss)%ppp*3+recieve_num*(3*ii+2))=&
                     X_send_buffer_tmp1(spp(ss)%ppp*(3*jj+1)+recieve_num*ii+1: &
                     spp(ss)%ppp*(3*jj+1)+recieve_num*(ii+1))

                X_send_buffer_tmp(jj*spp(ss)%ppp*3+recieve_num*(3*ii+2)+1: &
                     jj*spp(ss)%ppp*3+recieve_num*(3*ii+3))=&
                     X_send_buffer_tmp1(spp(ss)%ppp*(3*jj+2)+recieve_num*ii+1: &
                     spp(ss)%ppp*(3*jj+2)+recieve_num*(ii+1))
             end do
          end do

#if DBG_CHECK
          write(6,*) 'X_send_buffer_tmp1',X_send_buffer_tmp1*params%cpp%length
          write(6,*) 'X_send_buffer_tmp',X_send_buffer_tmp*params%cpp%length
#endif

          
       end if

       X_receive_buffer = 0.0_rp
       CALL MPI_SCATTER(X_send_buffer_tmp,3*recieve_num,MPI_REAL8, &
            X_receive_buffer,3*recieve_num,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
       if (params%orbit_model(1:2).EQ.'FO') then  
          spp(ss)%vars%X(1:recieve_num,:) = &
               RESHAPE(X_receive_buffer,(/recieve_num,3/))
          !set dummy values for new particles, so EZspline doesn't throw errors

          do ii=1_idef,spp(ss)%ppp-recieve_num
             spp(ss)%vars%X(recieve_num+ii,:)=spp(ss)%vars%X(1,:)
          end do
          
       else if (params%orbit_model(1:2).EQ.'GC') then
          spp(ss)%vars%Y(1:recieve_num,:) = &
               RESHAPE(X_receive_buffer,(/recieve_num,3/))

          do ii=1_idef,spp(ss)%ppp-recieve_num
             spp(ss)%vars%Y(recieve_num+ii,:)=spp(ss)%vars%Y(1,:)
          end do
          
          spp(ss)%vars%Y(:,2)=modulo(spp(ss)%vars%Y(:,2),2*C_PI)
       end if

       if (params%mpi_params%rank.EQ.0_idef) then
          filename = TRIM(params%path_to_outputs) // "restart_file.h5"
          call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
          if (h5error .EQ. -1) then
             write(output_unit_write,'("KORC ERROR: Something went wrong in: &
                  &load_particles_ic --> h5fopen_f")')
             call KORC_ABORT(14)
          end if

          write(tmp_str,'(I18)') ss

          dset = "/spp_" // TRIM(ADJUSTL(tmp_str)) // "/V"
          call load_array_from_hdf5(h5file_id,dset,V_send_buffer)

          call h5fclose_f(h5file_id, h5error)

#if DBG_CHECK
          write(6,*) 'V_send_buffer',V_send_buffer
#endif
          
          if (params%load_balance) then
             do ii=0,params%mpi_params%nmpi_prev-1

                V1_send_buffer_tmp(ii*spp(ss)%ppp+1:(ii+1)*spp(ss)%ppp)= &
                     V_send_buffer(3*ii*spp(ss)%ppp+1:(3*ii+1)*spp(ss)%ppp)
                V2_send_buffer_tmp(ii*spp(ss)%ppp+1:(ii+1)*spp(ss)%ppp)= &
                     V_send_buffer((3*ii+1)*spp(ss)%ppp+1:(3*ii+2)*spp(ss)%ppp)
                V3_send_buffer_tmp(ii*spp(ss)%ppp+1:(ii+1)*spp(ss)%ppp)= &
                     V_send_buffer((3*ii+2)*spp(ss)%ppp+1:(3*ii+3)*spp(ss)%ppp)
                
             end do
             do ii=1,params%mpi_params%nmpi_prev*spp(ss)%ppp

                V1_send_buffer_tmp1((mod(ii-1,params%mpi_params%nmpi))* &
                     spp(ss)%ppp/nmpi_ratio+(ii-1)/params%mpi_params%nmpi+1) = &
                     V1_send_buffer_tmp(ii)
                V2_send_buffer_tmp1((mod(ii-1,params%mpi_params%nmpi))* &
                     spp(ss)%ppp/nmpi_ratio+(ii-1)/params%mpi_params%nmpi+1) = &
                     V2_send_buffer_tmp(ii)
                V3_send_buffer_tmp1((mod(ii-1,params%mpi_params%nmpi))* &
                     spp(ss)%ppp/nmpi_ratio+(ii-1)/params%mpi_params%nmpi+1) = &
                     V3_send_buffer_tmp(ii)
                                                
             end do
             do ii=0,params%mpi_params%nmpi_prev-1

                V_send_buffer_tmp1(3*ii*spp(ss)%ppp+1:(3*ii+1)*spp(ss)%ppp)= &
                     V1_send_buffer_tmp1(ii*spp(ss)%ppp+1:(ii+1)*spp(ss)%ppp)
                V_send_buffer_tmp1((3*ii+1)*spp(ss)%ppp+1:(3*ii+2)*spp(ss)%ppp)= &
                     V2_send_buffer_tmp1(ii*spp(ss)%ppp+1:(ii+1)*spp(ss)%ppp)
                V_send_buffer_tmp1((3*ii+2)*spp(ss)%ppp+1:(3*ii+3)*spp(ss)%ppp)= &
                     V3_send_buffer_tmp1(ii*spp(ss)%ppp+1:(ii+1)*spp(ss)%ppp)


             end do
          else
             V_send_buffer_tmp1=V_send_buffer
          end if
          
          
          do jj=0_idef,params%mpi_params%nmpi_prev-1_idef
             do ii=0_idef,nmpi_ratio-1_idef

                V_send_buffer_tmp(jj*spp(ss)%ppp*3+recieve_num*(3*ii)+1: &
                     jj*spp(ss)%ppp*3+recieve_num*(3*ii+1))= &
                     V_send_buffer_tmp1(jj*spp(ss)%ppp*3+recieve_num*ii+1: &
                     jj*spp(ss)%ppp*3+recieve_num*(ii+1))

                V_send_buffer_tmp(jj*spp(ss)%ppp*3+recieve_num*(3*ii+1)+1: &
                     jj*spp(ss)%ppp*3+recieve_num*(3*ii+2))=&
                     V_send_buffer_tmp1(spp(ss)%ppp*(3*jj+1)+recieve_num*ii+1: &
                     spp(ss)%ppp*(3*jj+1)+recieve_num*(ii+1))

                V_send_buffer_tmp(jj*spp(ss)%ppp*3+recieve_num*(3*ii+2)+1: &
                     jj*spp(ss)%ppp*3+recieve_num*(3*ii+3))=&
                     V_send_buffer_tmp1(spp(ss)%ppp*(3*jj+2)+recieve_num*ii+1: &
                     spp(ss)%ppp*(3*jj+2)+recieve_num*(ii+1))
             end do
          end do


#if DBG_CHECK
          write(6,*) 'V_send_buffer_tmp1',V_send_buffer_tmp1
          write(6,*) 'V_send_buffer_tmp',V_send_buffer_tmp
#endif
          
       end if

       

       V_receive_buffer = 0.0_rp
       CALL MPI_SCATTER(V_send_buffer_tmp,3*recieve_num,MPI_REAL8, &
            V_receive_buffer,3*recieve_num,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
       spp(ss)%vars%V(1:recieve_num,:) = &
            RESHAPE(V_receive_buffer,(/recieve_num,3/))

       do ii=1_idef,spp(ss)%ppp-recieve_num
          spp(ss)%vars%V(recieve_num+ii,:)=spp(ss)%vars%V(1,:)
       end do
       
       if (params%mpi_params%rank.EQ.0_idef) then
          filename = TRIM(params%path_to_outputs) // "restart_file.h5"
          call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
          if (h5error .EQ. -1) then
             write(output_unit_write,'("KORC ERROR: Something went wrong in: &
                  &load_particles_ic --> h5fopen_f")')
             call KORC_ABORT(14)
          end if

          write(tmp_str,'(I18)') ss

          dset = "/spp_" // TRIM(ADJUSTL(tmp_str)) // "/flagCon"
          call load_array_from_hdf5(h5file_id,dset,AUX_send_buffer)

          call h5fclose_f(h5file_id, h5error)

          if (params%load_balance) then
             do ii=1,params%mpi_params%nmpi_prev*spp(ss)%ppp

                !write(6,*) 'send index',ii
                !write(6,*) 'tmp index',(mod(ii-1,params%mpi_params%nmpi))* &
                !     spp(ss)%ppp/nmpi_ratio+(ii-1)/params%mpi_params%nmpi+1

                AUX_send_buffer_tmp((mod(ii-1,params%mpi_params%nmpi))* &
                     spp(ss)%ppp/nmpi_ratio+(ii-1)/params%mpi_params%nmpi+1) = &
                     AUX_send_buffer(ii)
             end do
          else
             AUX_send_buffer_tmp=AUX_send_buffer
          end if
          
       end if
       
       AUX_receive_buffer = 0.0_rp
       CALL MPI_SCATTER(AUX_send_buffer_tmp,recieve_num,MPI_REAL8, &
            AUX_receive_buffer,recieve_num,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
       spp(ss)%vars%flagCon(1:recieve_num) = INT(AUX_receive_buffer,is)

       if (params%mpi_params%rank.EQ.0_idef) then
          filename = TRIM(params%path_to_outputs) // "restart_file.h5"
          call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
          if (h5error .EQ. -1) then
             write(output_unit_write,'("KORC ERROR: Something went wrong in: &
                  &load_particles_ic --> h5fopen_f")')
             call KORC_ABORT(14)
          end if

          write(tmp_str,'(I18)') ss

          dset = "/spp_" // TRIM(ADJUSTL(tmp_str)) // "/flagCol"
          call load_array_from_hdf5(h5file_id,dset,AUX_send_buffer)

          call h5fclose_f(h5file_id, h5error)

          if (params%load_balance) then
             do ii=1,params%mpi_params%nmpi_prev*spp(ss)%ppp

                !write(6,*) 'send index',ii
                !write(6,*) 'tmp index',(mod(ii-1,params%mpi_params%nmpi))* &
                !     spp(ss)%ppp/nmpi_ratio+(ii-1)/params%mpi_params%nmpi+1

                AUX_send_buffer_tmp((mod(ii-1,params%mpi_params%nmpi))* &
                     spp(ss)%ppp/nmpi_ratio+(ii-1)/params%mpi_params%nmpi+1) = &
                     AUX_send_buffer(ii)
             end do
          else
             AUX_send_buffer_tmp=AUX_send_buffer
          end if
          
       end if
       
       AUX_receive_buffer = 0.0_rp
       CALL MPI_SCATTER(AUX_send_buffer_tmp,recieve_num,MPI_REAL8, &
            AUX_receive_buffer,recieve_num,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
       spp(ss)%vars%flagCol(1:recieve_num) = INT(AUX_receive_buffer,is)

       if (params%mpi_params%rank.EQ.0_idef) then
          filename = TRIM(params%path_to_outputs) // "restart_file.h5"
          call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
          if (h5error .EQ. -1) then
             write(output_unit_write,'("KORC ERROR: Something went wrong in: &
                  &load_particles_ic --> h5fopen_f")')
             call KORC_ABORT(14)
          end if

          write(tmp_str,'(I18)') ss

          dset = "/spp_" // TRIM(ADJUSTL(tmp_str)) // "/flagRE"
          call load_array_from_hdf5(h5file_id,dset,AUX_send_buffer)

          call h5fclose_f(h5file_id, h5error)

          !write(6,*) 'ppp',spp(ss)%ppp
          !write(6,*) 'nmpi_prev',params%mpi_params%nmpi_prev
          !write(6,*) 'nmpi',params%mpi_params%nmpi
          !write(6,*) 'nmpi_ratio',nmpi_ratio
          
          if (params%load_balance) then
             do ii=1,params%mpi_params%nmpi_prev*spp(ss)%ppp

                !write(6,*) 'send index',ii
                !write(6,*) 'tmp index',(mod(ii-1,params%mpi_params%nmpi))* &
                !     spp(ss)%ppp/nmpi_ratio+(ii-1)/params%mpi_params%nmpi+1

                AUX_send_buffer_tmp((mod(ii-1,params%mpi_params%nmpi))* &
                     spp(ss)%ppp/nmpi_ratio+(ii-1)/params%mpi_params%nmpi+1) = &
                     AUX_send_buffer(ii)
             end do
          else
             AUX_send_buffer_tmp=AUX_send_buffer
          end if

#if DBG_CHECK
          write(6,*) 'flagRE_send_buffer',AUX_send_buffer
          write(6,*) 'flagRE_send_buffer_tmp',AUX_send_buffer_tmp
#endif
          
       end if             
       
       AUX_receive_buffer = 0.0_rp
       CALL MPI_SCATTER(AUX_send_buffer_tmp,recieve_num,MPI_REAL8, &
            AUX_receive_buffer,recieve_num,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
       spp(ss)%vars%flagRE(1:recieve_num) = INT(AUX_receive_buffer,is)
       
       if (params%mpi_params%rank.EQ.0_idef) then
          filename = TRIM(params%path_to_outputs) // "restart_file.h5"
          call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
          if (h5error .EQ. -1) then
             write(output_unit_write,'("KORC ERROR: Something went wrong in: &
                  &load_particles_ic --> h5fopen_f")')
             call KORC_ABORT(14)
          end if

          write(tmp_str,'(I18)') ss

          dset = "/spp_" // TRIM(ADJUSTL(tmp_str)) // "/g"
          call load_array_from_hdf5(h5file_id,dset,AUX_send_buffer)

          call h5fclose_f(h5file_id, h5error)

          if (params%load_balance) then
             do ii=1,params%mpi_params%nmpi_prev*spp(ss)%ppp

                !write(6,*) 'send index',ii
                !write(6,*) 'tmp index',(mod(ii-1,params%mpi_params%nmpi))* &
                !     spp(ss)%ppp/nmpi_ratio+(ii-1)/params%mpi_params%nmpi+1

                AUX_send_buffer_tmp((mod(ii-1,params%mpi_params%nmpi))* &
                     spp(ss)%ppp/nmpi_ratio+(ii-1)/params%mpi_params%nmpi+1) = &
                     AUX_send_buffer(ii)
             end do
          else
             AUX_send_buffer_tmp=AUX_send_buffer
          end if
          
       end if

       AUX_receive_buffer = 0.0_rp
       CALL MPI_SCATTER(AUX_send_buffer_tmp,recieve_num,MPI_REAL8, &
            AUX_receive_buffer,recieve_num,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
       spp(ss)%vars%g(1:recieve_num) = AUX_receive_buffer

       if (params%SC_E) then
          
          ALLOCATE(JSC0_buffer(F%dim_1D))
          ALLOCATE(JSC1_buffer(F%dim_1D))
          ALLOCATE(JSC2_buffer(F%dim_1D))
          ALLOCATE(JSC3_buffer(F%dim_1D))
          ALLOCATE(ESC_buffer(F%dim_1D))

          filename = TRIM(params%path_to_outputs) // "restart_file.h5"
          call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file_id, h5error)
          if (h5error .EQ. -1) then
             write(output_unit_write,'("KORC ERROR: Something went wrong in: &
                  &load_particles_ic --> h5fopen_f")')
             call KORC_ABORT(14)
          end if

          write(tmp_str,'(I18)') ss

          dset = "/spp_" // TRIM(ADJUSTL(tmp_str)) // "/J0_SC"
          call load_array_from_hdf5(h5file_id,dset,JSC0_buffer)
          
          dset = "/spp_" // TRIM(ADJUSTL(tmp_str)) // "/J1_SC"
          call load_array_from_hdf5(h5file_id,dset,JSC1_buffer)

          dset = "/spp_" // TRIM(ADJUSTL(tmp_str)) // "/J2_SC"
          call load_array_from_hdf5(h5file_id,dset,JSC2_buffer)

          dset = "/spp_" // TRIM(ADJUSTL(tmp_str)) // "/J3_SC"
          call load_array_from_hdf5(h5file_id,dset,JSC3_buffer)

          dset = "/spp_" // TRIM(ADJUSTL(tmp_str)) // "/E_SC"
          call load_array_from_hdf5(h5file_id,dset,ESC_buffer)
          
          call h5fclose_f(h5file_id, h5error)

          F%J0_SC_1D%PHI=JSC0_buffer
          F%J1_SC_1D%PHI=JSC1_buffer
          F%J2_SC_1D%PHI=JSC2_buffer
          F%J3_SC_1D%PHI=JSC3_buffer
          F%E_SC_1D%PHI=ESC_buffer/params%cpp%Eo
          
          DEALLOCATE(JSC0_buffer)
          DEALLOCATE(JSC1_buffer)
          DEALLOCATE(JSC2_buffer)
          DEALLOCATE(JSC3_buffer)
          DEALLOCATE(ESC_buffer)
          
       end if
       
       DEALLOCATE(X_send_buffer)
       DEALLOCATE(X_receive_buffer)

       DEALLOCATE(V_send_buffer)
       DEALLOCATE(V_receive_buffer)

       DEALLOCATE(AUX_send_buffer)
       DEALLOCATE(AUX_receive_buffer)
    end do

    if (params%orbit_model(1:2).EQ.'GC') then
       params%GC_coords=.TRUE.
    end if

#if DBG_CHECK    
    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    
    if (params%mpi_params%rank.eq.0_idef) then
       write(6,*) 'mpi',params%mpi_params%rank
       write(6,*) 'Y_R',spp(1)%vars%Y(:,1)*params%cpp%length
       write(6,*) 'Y_PHI',spp(1)%vars%Y(:,2)*params%cpp%length
       write(6,*) 'Y_Z',spp(1)%vars%Y(:,3)*params%cpp%length
       write(6,*) 'V_PLL',spp(1)%vars%V(:,1)*params%cpp%velocity*params%cpp%mass
       write(6,*) 'V_MU',spp(1)%vars%V(:,2)*params%cpp%energy/params%cpp%Bo
       write(6,*) 'flagCon',spp(1)%vars%flagCon
       write(6,*) 'flagCol',spp(1)%vars%flagCol
       write(6,*) 'flagRE',spp(1)%vars%flagRE
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    
    if (params%mpi_params%rank.eq.1_idef) then
       write(6,*) 'mpi',params%mpi_params%rank
       write(6,*) 'Y_R',spp(1)%vars%Y(:,1)*params%cpp%length
       write(6,*) 'Y_PHI',spp(1)%vars%Y(:,2)*params%cpp%length
       write(6,*) 'Y_Z',spp(1)%vars%Y(:,3)*params%cpp%length
       write(6,*) 'V_PLL',spp(1)%vars%V(:,1)*params%cpp%velocity*params%cpp%mass
       write(6,*) 'V_MU',spp(1)%vars%V(:,2)*params%cpp%energy/params%cpp%Bo
       write(6,*) 'flagCon',spp(1)%vars%flagCon
       write(6,*) 'flagCol',spp(1)%vars%flagCol
       write(6,*) 'flagRE',spp(1)%vars%flagRE
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    
    if (params%mpi_params%rank.eq.2_idef) then
       write(6,*) 'mpi',params%mpi_params%rank
       write(6,*) 'Y_R',spp(1)%vars%Y(:,1)*params%cpp%length
       write(6,*) 'Y_PHI',spp(1)%vars%Y(:,2)*params%cpp%length
       write(6,*) 'Y_Z',spp(1)%vars%Y(:,3)*params%cpp%length
       write(6,*) 'V_PLL',spp(1)%vars%V(:,1)*params%cpp%velocity*params%cpp%mass
       write(6,*) 'V_MU',spp(1)%vars%V(:,2)*params%cpp%energy/params%cpp%Bo
       write(6,*) 'flagCon',spp(1)%vars%flagCon
       write(6,*) 'flagCol',spp(1)%vars%flagCol
       write(6,*) 'flagRE',spp(1)%vars%flagRE
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    
    if (params%mpi_params%rank.eq.3_idef) then
       write(6,*) 'mpi',params%mpi_params%rank
       write(6,*) 'Y_R',spp(1)%vars%Y(:,1)*params%cpp%length
       write(6,*) 'Y_PHI',spp(1)%vars%Y(:,2)*params%cpp%length
       write(6,*) 'Y_Z',spp(1)%vars%Y(:,3)*params%cpp%length
       write(6,*) 'V_PLL',spp(1)%vars%V(:,1)*params%cpp%velocity*params%cpp%mass
       write(6,*) 'V_MU',spp(1)%vars%V(:,2)*params%cpp%energy/params%cpp%Bo
       write(6,*) 'flagCon',spp(1)%vars%flagCon
       write(6,*) 'flagCol',spp(1)%vars%flagCol
       write(6,*) 'flagRE',spp(1)%vars%flagRE
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

    
  end subroutine load_particles_ic

end module korc_HDF5
