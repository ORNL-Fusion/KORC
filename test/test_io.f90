module test_io
  implicit none

  ! * * * * * * * * * * * * * * * * * * * * !
  ! * * * Real and integer precisions * * * !
  ! * * * * * * * * * * * * * * * * * * * * !

  INTEGER, PUBLIC, PARAMETER 	:: is = KIND(INT(1,1)) 
  !! Definition of 1 Byte (8 bits) Fortran KORC integer type.
  INTEGER, PUBLIC, PARAMETER 	:: ip = KIND(INT(1,8)) 
  !! Definition of 8 Bytes (64 bits) Fortran KORC integer type.
  INTEGER, PUBLIC, PARAMETER 	:: idef = KIND(1) 
  !! Definition of the default KORC integer type on the system where
  !! KORC is compiled.
  INTEGER, PUBLIC, PARAMETER 	:: rdef = KIND(1.0) 
  !! Definition of the default KORC real type on the system where
  !! KORC is compiled.
#ifdef DOUBLE_PRECISION
  INTEGER, PUBLIC, PARAMETER 	:: rp = KIND(0.d0) 
  !! Definition of the KORC double precision real type.
#elif SINGLE_PRECISION
  INTEGER, PUBLIC, PARAMETER 	:: rp = KIND(1.0) 
  !! Definition of the KORC single precision real type.
#endif
  INTEGER, PUBLIC, PARAMETER 	:: MX_STRING_LENGTH = 1000 
  !! Default length of a KORC_STRING variable.
  INTEGER, PUBLIC, PARAMETER 	:: default_unit_open = 101 
  INTEGER, PUBLIC, PARAMETER 	:: test_unit_write = 203 


contains

  subroutine set_paths(path_to_outputs)
    !! @note Subroutine that sets the input/output paths.@endnote
    INTEGER 				:: argn,read_stat,exei
    !! Number of command line inputs. The default value is 
    !! two: the input files path and the outputs path.
    CHARACTER(MX_STRING_LENGTH) :: ctmp
    CHARACTER(MX_STRING_LENGTH),intent(out) :: path_to_outputs

    argn = command_argument_count()

    if (argn .EQ. 1_idef) then
       call get_command_argument(1,path_to_outputs)
    else
       write(6,*) 'Error with command line arguments'
       call exit(2)
    end if

    OPEN(UNIT=test_unit_write, &
         FILE=TRIM(path_to_outputs)//"output.test", &
         STATUS='UNKNOWN',FORM='FORMATTED',POSITION='REWIND')

    write(test_unit_write,'(/,"* * * * * PATHS * * * * *")')    
    write(test_unit_write,*) 'The output folder is: ',&
         TRIM(path_to_outputs)
    write(test_unit_write,'("* * * * * * * * * * * * *",/)')      
  end subroutine set_paths

end module test_io
