#ifdef FIO
!*******************************************************************************
!  @file korc_fio_interface.f90
!  @brief Interface for the fio interpolation library.
!*******************************************************************************

MODULE korc_fio
  USE, INTRINSIC :: iso_c_binding
  USE korc_types
  USE korc_input
  USE korc_HDF5
  USE mpi

  IMPLICIT NONE

  INTEGER (C_INT), PARAMETER :: FIO_SUCCESS        = 0
  INTEGER (C_INT), PARAMETER :: FIO_OUT_OF_BOUNDS  = 10002
  INTEGER (C_INT), PARAMETER :: FIO_NO_DATA        = 10006

  INTEGER (C_INT), PARAMETER :: FIO_M3DC1_SOURCE   = 3
  INTEGER (C_INT), PARAMETER :: FIO_NIMROD_SOURCE   = 5

  INTEGER (C_INT), PARAMETER :: FIO_TIMESLICE      = 1

  INTEGER (C_INT), PARAMETER :: FIO_SPECIES        = 3

  INTEGER (C_INT), PARAMETER :: FIO_ELECTRON       = 1
  INTEGER (C_INT), PARAMETER :: FIO_MAIN_ION       = -1

  INTEGER (C_INT), PARAMETER :: FIO_DENSITY        = 102
  INTEGER (C_INT), PARAMETER :: FIO_TEMPERATURE    = 103

  INTEGER (C_INT), PARAMETER :: FIO_ELECTRIC_FIELD = 1001
  INTEGER (C_INT), PARAMETER :: FIO_MAGNETIC_FIELD = 1003
  INTEGER (C_INT), PARAMETER :: FIO_VECTOR_POTENTIAL = 1002

  INTEGER (C_INT), PARAMETER :: FIO_TIME = 7001

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_add_field(icfield, ifield, op, fac)      &
          BIND(C, NAME='fio_add_field')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: icfield
       INTEGER (C_INT), VALUE, INTENT(IN) :: ifield
       INTEGER (C_INT), VALUE, INTENT(IN) :: op
       REAL (C_DOUBLE), VALUE, INTENT(IN) :: fac
     END FUNCTION fio_add_field
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_allocate_search_hint(isrc, hint)         &
          BIND(C, NAME='fio_allocate_search_hint')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: isrc
       TYPE (C_PTR), INTENT(OUT)          :: hint
     END FUNCTION fio_allocate_search_hint
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_close_field(ifield)                      &
          BIND(C, NAME='fio_close_field')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: ifield
     END FUNCTION fio_close_field
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_close_series(iseries)                    &
          BIND(C, NAME='fio_close_series')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: iseries
     END FUNCTION fio_close_series
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_close_source(isrc)                       &
          BIND(C, NAME='fio_close_source')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: isrc
     END FUNCTION fio_close_source
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_create_compound_field(ifield)            &
          BIND(C, NAME='fio_create_compound_field')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), INTENT(IN) :: ifield
     END FUNCTION fio_create_compound_field
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_deallocate_search_hint(isrc, hint)       &
          BIND(C, NAME='fio_deallocate_search_hint')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: isrc
       TYPE (C_PTR), INTENT(INOUT)        :: hint
     END FUNCTION fio_deallocate_search_hint
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_eval_field(ifield, x, v, hint)         &
          BIND(C, NAME='fio_eval_field')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: ifield
       REAL (C_DOUBLE), INTENT(IN)        :: x
       REAL (C_DOUBLE), INTENT(OUT)       :: v
       TYPE (C_PTR), VALUE, INTENT(IN)    :: hint
     END FUNCTION fio_eval_field
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_eval_field_deriv(ifield, x, v, hint)   &
          BIND(C, NAME='fio_eval_field_deriv')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: ifield
       REAL (C_DOUBLE), INTENT(IN)        :: x
       REAL (C_DOUBLE), INTENT(OUT)       :: v
       TYPE (C_PTR), VALUE, INTENT(IN)    :: hint
     END FUNCTION fio_eval_field_deriv
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_eval_series(iseries, x, v)             &
          BIND(C, NAME='fio_eval_series')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: iseries
       REAL (C_DOUBLE), INTENT(IN)        :: x
       REAL (C_DOUBLE), INTENT(OUT)       :: v
     END FUNCTION fio_eval_series
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_get_options(isrc)                        &
          BIND(C, NAME='fio_get_options')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN) :: isrc
     END FUNCTION fio_get_options
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_get_available_fields(isrc, n, f)         &
          BIND(C, NAME='fio_get_available_fields')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN)         :: isrc
       INTEGER (C_INT), INTENT(OUT)               :: n
       INTEGER (C_INT), DIMENSION(:), INTENT(OUT) :: f
     END FUNCTION fio_get_available_fields
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_open_source(itype, filename, handle)     &
          BIND(C, NAME='fio_open_source')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(IN)        :: itype
       CHARACTER (kind=C_CHAR,len=1), INTENT(IN) :: filename
       INTEGER (C_INT), INTENT(OUT)              :: handle
     END FUNCTION fio_open_source
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_get_field(isrc, type, handle)            &
          BIND(C, NAME='fio_get_field')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(in) :: isrc
       INTEGER (C_INT), VALUE, INTENT(IN) :: type
       INTEGER (C_INT), INTENT(INOUT)     :: handle
     END FUNCTION fio_get_field
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_set_int_option(iopt, v)                  &
          BIND(C, NAME='fio_set_int_option')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE, INTENT(in) :: iopt
       INTEGER (C_INT), VALUE, INTENT(in) :: v
     END FUNCTION fio_set_int_option
  END INTERFACE

  INTERFACE
     INTEGER (C_INT) FUNCTION fio_get_real_field_parameter(ifield, t, p) &
          BIND(C, NAME='fio_get_real_field_parameter')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       integer(c_int), intent(in), value :: ifield, t
       real(c_double) :: p
     END FUNCTION fio_get_real_field_parameter
  END INTERFACE

CONTAINS

  SUBROUTINE initialize_m3d_c1(params, F, P, spp,init)

    TYPE(KORC_PARAMS), INTENT(INOUT)           :: params
    TYPE(FIELDS), INTENT(INOUT)                :: F
    TYPE(PROFILES), INTENT(INOUT)              :: P
    TYPE(SPECIES), DIMENSION(:), INTENT(INOUT) :: spp
    LOGICAL, INTENT(IN)  :: init
    
    INTEGER                                    :: ii
    INTEGER                                    :: pp
    INTEGER                                    :: status
    INTEGER                                    :: isrc
    real(c_double)  ::  time0,time1
    INTEGER (C_INT)                         :: FIO_tmp
    TYPE(C_PTR) :: hint_tmp
    real(rp), DIMENSION(3) :: x
    REAL(rp), DIMENSION(3)         :: Btmp
    
    if (init) then

       status = fio_open_source(FIO_M3DC1_SOURCE,           &
            TRIM(params%magnetic_field_filename)            &
            // C_NULL_CHAR, F%isrc)
       
       F%Efield = Efield
       F%PSIp_lim=PSIp_lim
       F%PSIp_0=PSIp_0
       F%ReInterp_2x1t=ReInterp_2x1t

       if (params%restart.OR.params%proceed) then

          isrc=F%isrc
          status = fio_get_options(isrc)
          status = fio_set_int_option(FIO_TIMESLICE, ind0_2x1t)
          status = fio_get_field(isrc, FIO_MAGNETIC_FIELD, F%FIO_B)
          
          hint_tmp=c_null_ptr
          x(1)=spp(1)%Ro
          x(2)=spp(1)%PHIo
          x(3)=spp(1)%Zo

          status = fio_eval_field(F%FIO_B, x(1),                      &
               Btmp(1),hint_tmp)

          F%Bo = Btmp(2)
          F%Eo = 1.0
          F%Ro = 1.0
          F%Zo = 1.0

          status = fio_close_field(F%FIO_B)
          
          call load_prev_iter(params)
          F%ind0_2x1t=params%prev_iter_2x1t+1
          
       else
          F%ind0_2x1t=ind0_2x1t
       end if

       F%ind_2x1t=F%ind0_2x1t



    else
       status = fio_close_field(F%FIO_B)
       status = fio_close_field(F%FIO_B+1)
       status = fio_close_field(F%FIO_A)

       status = fio_close_field(P%FIO_ne)
       status = fio_close_field(P%FIO_te)
       status = fio_close_field(P%FIO_ni)
       
    end if

    isrc=F%isrc
    
    status = fio_get_options(isrc)
       

    if (.not.F%ReInterp_2x1t) then
       status = fio_set_int_option(FIO_TIMESLICE, params%time_slice)
    else
       status = fio_set_int_option(FIO_TIMESLICE, F%ind_2x1t)
    end if

    
    status = fio_get_field(isrc, FIO_MAGNETIC_FIELD, F%FIO_B)
    status = fio_get_field(isrc, FIO_ELECTRIC_FIELD, F%FIO_E)
    status = fio_get_field(isrc, FIO_VECTOR_POTENTIAL, F%FIO_A)

    status = fio_get_real_field_parameter(F%FIO_B, FIO_TIME, time0)

    write(output_unit_write,*) 'FIO present time index',F%ind_2x1t
    write(output_unit_write,*) 'FIO present time',time0

    if (F%ReInterp_2x1t) then
       status = fio_set_int_option(FIO_TIMESLICE, F%ind_2x1t+1)

       status = fio_get_field(isrc, FIO_MAGNETIC_FIELD, FIO_tmp)
       
       status = fio_get_real_field_parameter(FIO_tmp, FIO_TIME, time1)
       write(output_unit_write,*) 'FIO next time index',F%ind_2x1t+1
       write(output_unit_write,*) 'FIO next time',time1

       status = fio_set_int_option(FIO_TIMESLICE, F%ind_2x1t)

       params%snapshot_frequency=time1-time0

       !write(6,*) 'snapshot_frequency',params%snapshot_frequency
       !write(6,*) 'dt',params%dt
       
       if (.not.init) then
          
          params%t_skip = FLOOR(params%snapshot_frequency/params%cpp%time/ &
               params%dt,ip)
          
       end if
              
    end if    
    
    if (.not.F%Efield) F%FIO_E=-1

    status = fio_set_int_option(FIO_SPECIES, FIO_ELECTRON);
    status = fio_get_field(isrc, FIO_DENSITY, P%FIO_ne);
    status = fio_get_field(isrc, FIO_TEMPERATURE, P%FIO_te);

    status = fio_set_int_option(FIO_SPECIES, FIO_MAIN_ION);
    status = fio_get_field(isrc, FIO_DENSITY, P%FIO_ni);


    if (init) then
       if (.NOT.(params%restart.OR.params%proceed)) then
          hint_tmp=c_null_ptr
          x(1)=spp(1)%Ro
          x(2)=spp(1)%PHIo
          x(3)=spp(1)%Zo

          status = fio_eval_field(F%FIO_B, x(1),                      &
               Btmp(1),hint_tmp)

          F%Bo = Btmp(2)
          F%Eo = 1.0
          F%Ro = 1.0
          F%Zo = 1.0
       end if

       do ii = 1, params%num_species

          do pp = 1, spp(ii)%ppp
             status = fio_allocate_search_hint(isrc, spp(ii)%vars%hint(pp))
             !spp(ii)%vars%hint(pp)=c_null_ptr
          end do

          spp(ii)%vars%cart = .false.
       end do
    end if

    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,*) 'Calculate B',F%FIO_B
       write(output_unit_write,'("Magnetic field: ",E17.10)') F%Bo
       write(output_unit_write,*) 'Calculate E',F%FIO_E
       write(output_unit_write,'("Electric field: ",E17.10)') F%Eo
       write(output_unit_write,*) 'Calculate A',F%FIO_A
       write(output_unit_write,*) 'Calculate ne',P%FIO_ne
       write(output_unit_write,*) 'Calculate Te',P%FIO_te
       write(output_unit_write,*) 'Calculate ni',P%FIO_ni
    end if
       
  END SUBROUTINE initialize_m3d_c1
  
  SUBROUTINE initialize_nimrod(params, F, P, spp,init)

    TYPE(KORC_PARAMS), INTENT(INOUT)           :: params
    TYPE(FIELDS), INTENT(INOUT)                :: F
    TYPE(PROFILES), INTENT(INOUT)              :: P
    TYPE(SPECIES), DIMENSION(:), INTENT(INOUT) :: spp
    LOGICAL, INTENT(IN)  :: init
    
    INTEGER                                    :: ii
    INTEGER                                    :: pp
    INTEGER                                    :: status
    INTEGER                                    :: isrc
    real(c_double)  ::  time0,time1
    INTEGER (C_INT)                         :: FIO_tmp
    TYPE(C_PTR) :: hint_tmp
    real(rp), DIMENSION(3) :: x
    REAL(rp), DIMENSION(3)         :: Btmp
    CHARACTER(150) :: filename
    INTEGER, PARAMETER :: list_file = 102
    INTEGER :: num_files,ios


    if (init) then
       F%Efield = Efield
       F%PSIp_lim=PSIp_lim
       F%PSIp_0=PSIp_0
       F%ReInterp_2x1t=ReInterp_2x1t
       
       if (F%ReInterp_2x1t) then
          params%magnetic_field_directory=magnetic_field_directory
          params%magnetic_field_list=magnetic_field_list

          filename=TRIM(params%magnetic_field_directory)//TRIM(params%magnetic_field_list)
          
          OPEN(UNIT=list_file,FILE=filename,STATUS='OLD')

          READ(UNIT=list_file,FMT='(I4)',IOSTAT=ios) num_files
          if (ios/=0) then
             write(6,*) 'Error reading in magnetic_file_list for initialize_nimrod'
             call KORC_abort(22)
          endif

          !write(6,*) 'num_files',num_files
          
          ALLOCATE(params%magnetic_field_filenames(num_files))
          ALLOCATE(params%time_of_filenames(num_files))

          do ii=1,num_files
             READ(UNIT=list_file,FMT='(A15,F12.10)',IOSTAT=ios) params%magnetic_field_filenames(ii),params%time_of_filenames(ii)

             params%magnetic_field_filenames(ii)=ADJUSTL(params%magnetic_field_filenames(ii))
             
             !write(6,*) 'ii',ii,'filename ',TRIM(params%magnetic_field_filenames(ii)), &
             !     ' time',params%time_of_filenames(ii)
             
          enddo

          CLOSE(UNIT=list_file)
          

          
       endif       
       
       if (params%restart.or.params%proceed) then

          if (.not.F%ReInterp_2x1t) then
             filename=TRIM(params%magnetic_field_filename)  
          else
             filename=TRIM(params%magnetic_field_directory)//TRIM(params%magnetic_field_filenames(ind0_2x1t))  
          end if           
       
          status = fio_open_source(FIO_NIMROD_SOURCE, &
               filename // C_NULL_CHAR, F%isrc)
          
          isrc=F%isrc
          status = fio_get_options(isrc)
          status = fio_get_field(isrc, FIO_MAGNETIC_FIELD, F%FIO_B)
          
          hint_tmp=c_null_ptr
          x(1)=spp(1)%Ro
          x(2)=spp(1)%PHIo
          x(3)=spp(1)%Zo

          status = fio_eval_field(F%FIO_B, x(1),                      &
               Btmp(1),hint_tmp)

          F%Bo = Btmp(2)
          F%Eo = 1.0
          F%Ro = 1.0
          F%Zo = 1.0

          status = fio_close_field(F%FIO_B)
          status=fio_close_source(F%isrc)
          
          call load_prev_iter(params)
          F%ind0_2x1t=params%prev_iter_2x1t+1
       else
          F%ind0_2x1t=ind0_2x1t
       end if

       F%ind_2x1t=F%ind0_2x1t

       if (.not.F%ReInterp_2x1t) then
          filename=TRIM(params%magnetic_field_filename)  
       else
          filename=TRIM(params%magnetic_field_directory)//TRIM(params%magnetic_field_filenames(F%ind_2x1t))  
       end if
       
       status = fio_open_source(FIO_NIMROD_SOURCE, &
            filename // C_NULL_CHAR, F%isrc)

    else
       status = fio_close_field(F%FIO_B)
       status = fio_close_field(F%FIO_B+1) ! actually for E, to account for F%Efield=-1
!       status = fio_close_field(F%FIO_A)

       status = fio_close_field(P%FIO_ne)
!       status = fio_close_field(P%FIO_te)
!       status = fio_close_field(P%FIO_ni)

       status=fio_close_source(F%isrc)

       if (.not.F%ReInterp_2x1t) then
          filename=TRIM(params%magnetic_field_filename)  
       else
          filename=TRIM(params%magnetic_field_directory)//TRIM(params%magnetic_field_filenames(F%ind_2x1t))  
       end if
       
       status = fio_open_source(FIO_NIMROD_SOURCE, &
            filename // C_NULL_CHAR, F%isrc)       
             
    end if

    isrc=F%isrc
    
    status = fio_get_options(isrc)


    !if (.not.F%ReInterp_2x1t) then
    !   status = fio_set_int_option(FIO_TIMESLICE, params%time_slice)
    !else
    !   status = fio_set_int_option(FIO_TIMESLICE, F%ind_2x1t)
    !end if
    ! For NIMROD dump*, status=fio_get_options(isrc)=0, can't set FIO_TIMESLICE
    ! as was the case for C1.h5 having information on multiple timeslices 

    
    status = fio_get_field(isrc, FIO_MAGNETIC_FIELD, F%FIO_B)
    status = fio_get_field(isrc, FIO_ELECTRIC_FIELD, F%FIO_E)
    !status = fio_get_field(isrc, FIO_VECTOR_POTENTIAL, F%FIO_A)
    F%FIO_A=-1

    !write(6,*) 'time_of_filenames',time_of_filenames
    !write(6,*) 'ind_2x1t',F%ind_2x1t
    


    if (F%ReInterp_2x1t) then

       !status = fio_get_real_field_parameter(F%FIO_B, FIO_TIME, time0)
       time0=params%time_of_filenames(F%ind_2x1t)

       write(output_unit_write,*) 'FIO present time index',F%ind_2x1t
       write(output_unit_write,*) 'FIO present time',time0
       
       !status = fio_set_int_option(FIO_TIMESLICE, F%ind_2x1t+1)

       !status = fio_get_field(isrc, FIO_MAGNETIC_FIELD, FIO_tmp)
       
       !status = fio_get_real_field_parameter(FIO_tmp, FIO_TIME, time1)
       time1=params%time_of_filenames(F%ind_2x1t+1)
       write(output_unit_write,*) 'FIO next time index',F%ind_2x1t+1
       write(output_unit_write,*) 'FIO next time',time1

       !status = fio_set_int_option(FIO_TIMESLICE, F%ind_2x1t)

       params%snapshot_frequency=time1-time0

       !write(6,*) 'snapshot_frequency',params%snapshot_frequency
       !write(6,*) 'dt',params%dt
       
       if (.not.init) then
          
          params%t_skip = FLOOR(params%snapshot_frequency/params%cpp%time/ &
               params%dt,ip)
          
       end if
              
    end if    
    
    if (.not.F%Efield) F%FIO_E=-1

    
    status = fio_set_int_option(FIO_SPECIES, FIO_ELECTRON);
       
    
    status = fio_get_field(isrc, FIO_DENSITY, P%FIO_ne);
    !status = fio_get_field(isrc, FIO_TEMPERATURE, P%FIO_te);

    !status = fio_set_int_option(FIO_SPECIES, FIO_MAIN_ION);
    !status = fio_get_field(isrc, FIO_DENSITY, P%FIO_ni);


    if (init) then
       if (.NOT.(params%restart.OR.params%proceed)) then
          hint_tmp=c_null_ptr
          x(1)=spp(1)%Ro
          x(2)=spp(1)%PHIo
          x(3)=spp(1)%Zo

          status = fio_eval_field(F%FIO_B, x(1),                      &
               Btmp(1),hint_tmp)

          F%Bo = Btmp(2)
          F%Eo = 1.0
          F%Ro = 1.0
          F%Zo = 1.0
       end if    

       do ii = 1, params%num_species

          do pp = 1, spp(ii)%ppp
             status = fio_allocate_search_hint(isrc, spp(ii)%vars%hint(pp))
             !spp(ii)%vars%hint(pp)=c_null_ptr
          end do

          spp(ii)%vars%cart = .false.
       end do
    end if

    if (params%mpi_params%rank .EQ. 0) then
       write(output_unit_write,*) 'Calculate B',F%FIO_B
       write(output_unit_write,*) 'Calculate E',F%FIO_E
       write(output_unit_write,*) 'Calculate ne',P%FIO_ne
    end if
       
  END SUBROUTINE initialize_nimrod
  
  SUBROUTINE initialize_m3d_c1_imp(params,F,P,num_imp,init)

    TYPE(KORC_PARAMS), INTENT(IN)           :: params
    TYPE(FIELDS), INTENT(IN)                :: F
    TYPE(PROFILES), INTENT(INOUT)              :: P
    INTEGER, INTENT(IN)				:: num_imp
    LOGICAL, INTENT(IN)  :: init

    INTEGER                                    :: ii
    INTEGER                                    :: status
    INTEGER                                    :: isrc
    INTEGER,ALLOCATABLE,DIMENSION(:)          :: Zo
    INTEGER  :: A,Zo1

    !status = fio_open_source(FIO_M3DC1_SOURCE,           &
    !     TRIM(params%magnetic_field_filename)            &
    !     // C_NULL_CHAR, isrc)

    isrc=F%isrc
    
    status = fio_get_options(isrc)
       
    
    if (.not.F%ReInterp_2x1t) then
       status = fio_set_int_option(FIO_TIMESLICE, params%time_slice)
    else
       status = fio_set_int_option(FIO_TIMESLICE, F%ind_2x1t)
    end if

    if (init) then
       ALLOCATE(P%FIO_nimp(num_imp))
    else
       do ii=1,num_imp
          status = fio_close_field(P%FIO_nimp(ii))
       end do          
    endif
       
    !write(6,*) size(params%Zj)
    !write(6,*) size(params%Zj(ubound(params%Zj)))
    Zo=int(params%Zj(ubound(params%Zj)))
    Zo1=Zo(1)

    if (Zo1.eq.18) then
       A=40
    else if (Zo1.eq.10) then
       A=20
    else if (Zo1.eq.6) then
       A=12
    end if
    
    do ii=1,num_imp
       status = fio_set_int_option(FIO_SPECIES, &
            FIO_MAKE_SPECIES(A, Zo1, Zo1+1-ii));
       status = fio_get_field(isrc, FIO_DENSITY, P%FIO_nimp(ii));
    end do

    if (params%mpi_params%rank .EQ. 0) then
       do ii=1,num_imp
          write(output_unit_write,*) 'Calculate nimp_',ii,P%FIO_nimp(ii)
       end do
    end if
       
  END SUBROUTINE initialize_m3d_c1_imp

  SUBROUTINE finalize_fio(params, F, P)
    TYPE(KORC_PARAMS), INTENT(IN)           :: params
    TYPE(FIELDS), INTENT(IN)                :: F
    TYPE(PROFILES), INTENT(INOUT)              :: P
    INTEGER                                    :: status
    INTEGER                                    :: ii

    status = fio_close_field(F%FIO_B)
    status = fio_close_field(F%FIO_B+1)
    if (F%FIO_A.gt.0) then
       status = fio_close_field(F%FIO_A)
    endif
       
    status = fio_close_field(P%FIO_ne)
    if (P%FIO_te.gt.0) then
       status = fio_close_field(P%FIO_te)
    endif
    if (P%FIO_ni.gt.0) then
       status = fio_close_field(P%FIO_ni)
    endif
       
    if (params%collisions) then
       do ii=1,params%num_impurity_species
          status = fio_close_field(P%FIO_nimp(ii))
       end do
    end if

    status=fio_close_source(F%isrc)
    
  end SUBROUTINE FINALIZE_FIO

  
  FUNCTION FIO_MAKE_SPECIES(m, p, e)
    INTEGER, INTENT(IN) :: m
    INTEGER, INTENT(IN) :: p
    INTEGER, INTENT(IN) :: e
    INTEGER(C_INT) :: FIO_MAKE_SPECIES
    
    FIO_MAKE_SPECIES = e + p*256 + (m-p)*65536
  END FUNCTION FIO_MAKE_SPECIES


END MODULE korc_fio
#endif
