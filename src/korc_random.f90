MODULE korc_random

  USE, INTRINSIC :: iso_c_binding
  USE korc_types
  USE korc_hpc, ONLY : get_thread_number, get_max_threads

  IMPLICIT NONE

!*******************************************************************************
!  Interface binding for the c++ random functions.
!*******************************************************************************

  INTERFACE
     TYPE (C_PTR) FUNCTION random_construct_U(seed) BIND(C, NAME='random_construct_U')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE :: seed

     END FUNCTION random_construct_U
  END INTERFACE

  INTERFACE
     TYPE (C_PTR) FUNCTION random_construct_N(seed) BIND(C, NAME='random_construct_N')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       INTEGER (C_INT), VALUE :: seed

     END FUNCTION random_construct_N
  END INTERFACE

  INTERFACE
     REAL (rp) FUNCTION random_get_number_U(r) BIND(C, NAME='random_get_number_U')
       USE korc_types, ONLY : rp
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       TYPE (C_PTR), VALUE :: r

     END FUNCTION random_get_number_U
  END INTERFACE

  INTERFACE
     REAL (rp) FUNCTION random_get_number_N(r) BIND(C, NAME='random_get_number_N')
       USE korc_types, ONLY : rp
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       TYPE (C_PTR), VALUE :: r

     END FUNCTION random_get_number_N
  END INTERFACE

  INTERFACE
     SUBROUTINE random_destroy_U(r) BIND(C, NAME='random_destroy_U')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       TYPE (C_PTR), VALUE :: r

     END SUBROUTINE random_destroy_U
  END INTERFACE

  INTERFACE
     SUBROUTINE random_destroy_N(r) BIND(C, NAME='random_destroy_N')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       TYPE (C_PTR), VALUE :: r

     END SUBROUTINE random_destroy_N
  END INTERFACE
  
  INTERFACE
     SUBROUTINE random_set_dist_U(r, low, high) BIND(C, NAME='random_set_dist_U')
       USE korc_types, ONLY : rp
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       TYPE (C_PTR), VALUE :: r
       REAL (rp), VALUE    :: low
       REAL (rp), VALUE    :: high
     END SUBROUTINE random_set_dist_U
  END INTERFACE

  INTERFACE
     SUBROUTINE random_set_dist_N(r, low, high) BIND(C, NAME='random_set_dist_N')
       USE korc_types, ONLY : rp
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       TYPE (C_PTR), VALUE :: r
       REAL (rp), VALUE    :: low
       REAL (rp), VALUE    :: high
     END SUBROUTINE random_set_dist_N
  END INTERFACE

  INTERFACE
     SUBROUTINE random_set_seed_U(r, seed) BIND(C, NAME='random_set_seed_U')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       TYPE (C_PTR), VALUE    :: r
       INTEGER (C_INT), VALUE :: seed
     END SUBROUTINE random_set_seed_U
  END INTERFACE

  INTERFACE
     SUBROUTINE random_set_seed_N(r, seed) BIND(C, NAME='random_set_seed_N')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       TYPE (C_PTR), VALUE    :: r
       INTEGER (C_INT), VALUE :: seed
     END SUBROUTINE random_set_seed_N
  END INTERFACE

!*******************************************************************************
!  Class Defintions
!*******************************************************************************

  TYPE :: random_base
     TYPE(C_PTR), DIMENSION(:), POINTER :: states => null()
  END TYPE

  TYPE, EXTENDS(random_base) :: random_U_context
  CONTAINS
     PROCEDURE :: get => random_U_get_random
     PROCEDURE :: get_array => random_U_get_randoms
     PROCEDURE :: set => random_U_set_dist
     PROCEDURE :: seed => random_U_set_seed
     FINAL     :: random_U_context_destruct
  END TYPE

  TYPE, EXTENDS(random_base) :: random_N_context
  CONTAINS
     PROCEDURE :: get => random_N_get_random
     PROCEDURE :: get_array => random_N_get_randoms
     PROCEDURE :: set => random_N_set_dist
     PROCEDURE :: seed => random_N_set_seed
     FINAL     :: random_N_context_destruct
  END TYPE

  TYPE :: random_context
     CLASS (random_U_context), POINTER :: uniform => null()
     CLASS (random_N_context), POINTER :: normal => null()
  CONTAINS
     FINAL :: random_context_destruct
  END TYPE

CONTAINS

!*******************************************************************************
!  Constructors
!*******************************************************************************

  FUNCTION random_U_context_construct(seed, mpi_rank)
  IMPLICIT NONE

!  Arguments
  CLASS (random_U_context), POINTER :: random_U_context_construct
  INTEGER, INTENT(IN)               :: seed
  INTEGER, INTENT(IN)               :: mpi_rank

!  Local Variables
  INTEGER                           :: num_threads
  INTEGER                           :: thread_num

!  Start of executable code.
  ALLOCATE(random_U_context_construct)

  num_threads = get_max_threads()
  ALLOCATE(random_U_context_construct%states(0:num_threads))

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(thread_num)
  thread_num = get_thread_number()
  random_U_context_construct%states(thread_num) =                              &
     random_construct_U(seed + mpi_rank*num_threads + thread_num)
!$OMP END PARALLEL
  END FUNCTION

  FUNCTION random_N_context_construct(seed, mpi_rank)
  IMPLICIT NONE

!  Arguments
  CLASS (random_N_context), POINTER :: random_N_context_construct
  INTEGER, INTENT(IN)               :: seed
  INTEGER, INTENT(IN)               :: mpi_rank

!  Local Variables
  INTEGER                           :: num_threads
  INTEGER                           :: thread_num

!  Start of executable code.
  ALLOCATE(random_N_context_construct)

  num_threads = get_max_threads()
  ALLOCATE(random_N_context_construct%states(0:num_threads))

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(thread_num)
  thread_num = get_thread_number()
  random_N_context_construct%states(thread_num) =                              &
     random_construct_N(seed + mpi_rank*num_threads + thread_num)
!$OMP END PARALLEL
  END FUNCTION

  FUNCTION random_context_construct(seed, mpi_rank)

  IMPLICIT NONE

!  Arguments
  CLASS (random_context), POINTER :: random_context_construct
  INTEGER, INTENT(IN)             :: seed
  INTEGER, INTENT(IN)             :: mpi_rank

!  Start of executable code.
  ALLOCATE(random_context_construct)

  random_context_construct%uniform => random_U_context_construct(seed, mpi_rank)
  random_context_construct%normal => random_N_context_construct(seed, mpi_rank)

  END FUNCTION

!*******************************************************************************
!  Destructors
!*******************************************************************************

  SUBROUTINE random_U_context_destruct(this)
  IMPLICIT NONE

!  Arguments
  TYPE (random_U_context), INTENT(inout) :: this

!  Start of executable code
  IF (ASSOCIATED(this%states)) THEN
!$OMP PARALLEL DEFAULT(SHARED)
     CALL random_destroy_U(this%states(get_thread_number()))
!$OMP END PARALLEL
     DEALLOCATE(this%states)
     this%states => null()
  END IF

  END SUBROUTINE

  SUBROUTINE random_N_context_destruct(this)
  IMPLICIT NONE

!  Arguments
  TYPE (random_N_context), INTENT(inout) :: this

!  Start of executable code
  IF (ASSOCIATED(this%states)) THEN
!$OMP PARALLEL DEFAULT(SHARED)
     CALL random_destroy_N(this%states(get_thread_number()))
!$OMP END PARALLEL
     DEALLOCATE(this%states)
     this%states => null()
  END IF

  END SUBROUTINE

  SUBROUTINE random_context_destruct(this)
  IMPLICIT NONE

!  Arguments
  TYPE (random_context), INTENT(inout) :: this

!  Start of executable code
  IF (ASSOCIATED(this%uniform)) THEN
     DEALLOCATE(this%uniform)
     this%uniform => null()
  END IF

  IF (ASSOCIATED(this%normal)) THEN
     DEALLOCATE(this%normal)
     this%normal => null()
  END IF

  END SUBROUTINE

!*******************************************************************************
!  Getters
!*******************************************************************************
  FUNCTION random_U_get_random(this)
  IMPLICIT NONE

!  Arguments
  REAL(rp)                             :: random_U_get_random
  CLASS (random_U_context), INTENT(in) :: this

!  Start of executable code
  random_U_get_random = random_get_number_U(this%states(get_thread_number()))

  END FUNCTION

  FUNCTION random_N_get_random(this)
  IMPLICIT NONE

!  Arguments
  REAL(rp)                             :: random_N_get_random
  CLASS (random_N_context), INTENT(in) :: this

!  Start of executable code
  random_N_get_random = random_get_number_N(this%states(get_thread_number()))

  END FUNCTION

  SUBROUTINE random_U_get_randoms(this, nums)
  IMPLICIT NONE

!  Arguments
  CLASS (random_U_context), INTENT(in) :: this
  REAL(rp), DIMENSION(:), INTENT(out)  :: nums

!  Local Variables
  INTEGER                                    :: i

!  Start of executable code
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
  DO i = 1, SIZE(nums)
     nums(i) = this%get()
  END DO
!$OMP END PARALLEL DO
  END SUBROUTINE

  SUBROUTINE random_N_get_randoms(this, nums)
  IMPLICIT NONE

!  Arguments
  CLASS (random_N_context), INTENT(in) :: this
  REAL(rp), DIMENSION(:), INTENT(out)  :: nums

!  Local Variables
  INTEGER                                    :: i

!  Start of executable code
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
  DO i = 1, SIZE(nums)
     nums(i) = this%get()
  END DO
!$OMP END PARALLEL DO
  END SUBROUTINE

!*******************************************************************************
!  Setters
!*******************************************************************************
  SUBROUTINE random_U_set_dist(this, low, high)
  IMPLICIT NONE

!  Arguments
  CLASS (random_U_context), INTENT(in) :: this
  REAL(rp), INTENT(IN)                 :: low
  REAL(rp), INTENT(IN)                 :: high

!  Start of executable code
!$OMP PARALLEL DEFAULT(SHARED)
  CALL random_set_dist_U(this%states(get_thread_number()), low, high)
!$OMP END PARALLEL

  END SUBROUTINE

  SUBROUTINE random_N_set_dist(this, low, high)
  IMPLICIT NONE

!  Arguments
  CLASS (random_N_context), INTENT(in) :: this
  REAL(rp), INTENT(IN)                 :: low
  REAL(rp), INTENT(IN)                 :: high

!  Start of executable code
!$OMP PARALLEL DEFAULT(SHARED)
  CALL random_set_dist_N(this%states(get_thread_number()), low, high)
!$OMP END PARALLEL

  END SUBROUTINE

  SUBROUTINE random_U_set_seed(this, seed, mpi_rank)
  IMPLICIT NONE

!  Arguments
  CLASS (random_U_context), INTENT(in) :: this
  INTEGER, INTENT(in)                  :: seed
  INTEGER, INTENT(IN)                  :: mpi_rank

!  Local Variables
  INTEGER                              :: thread_num

!  Start of executable code
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(thread_num)
  thread_num = get_thread_number()
  CALL random_set_seed_U(this%states(thread_num),                             &
                         seed + mpi_rank*get_max_threads() + thread_num)
!$OMP END PARALLEL

  END SUBROUTINE

  SUBROUTINE random_N_set_seed(this, seed, mpi_rank)
  IMPLICIT NONE

!  Arguments
  CLASS (random_N_context), INTENT(in) :: this
  INTEGER, INTENT(in)                  :: seed
  INTEGER, INTENT(IN)                  :: mpi_rank

!  Local Variables
  INTEGER                              :: thread_num

!  Start of executable code
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(thread_num)
   thread_num = get_thread_number()
   CALL random_set_seed_N(this%states(thread_num),                             &
                          seed + mpi_rank*get_max_threads() + thread_num)
!$OMP END PARALLEL

   END SUBROUTINE

END MODULE korc_random
