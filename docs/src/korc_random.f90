!include 'mkl_vsl.f90'

MODULE korc_random

  USE, INTRINSIC :: iso_c_binding
  USE korc_types
  !    use mkl_vsl_type
  !    use mkl_vsl

  IMPLICIT NONE

  TYPE(C_PTR), DIMENSION(:), ALLOCATABLE , PRIVATE :: states
  TYPE(C_PTR),  PRIVATE :: state
  !    TYPE(VSL_STREAM_STATE), PRIVATE :: stream

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
     REAL (C_DOUBLE) FUNCTION random_get_number_U(r) BIND(C, NAME='random_get_number_U')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       TYPE (C_PTR), VALUE :: r

     END FUNCTION random_get_number_U
  END INTERFACE

  INTERFACE
     REAL (C_DOUBLE) FUNCTION random_get_number_N(r) BIND(C, NAME='random_get_number_N')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       TYPE (C_PTR), VALUE :: r

     END FUNCTION random_get_number_N
  END INTERFACE

  INTERFACE
     REAL (C_DOUBLE) FUNCTION random_get_number(r) BIND(C, NAME='random_get_number')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       TYPE (C_PTR), VALUE :: r

     END FUNCTION random_get_number
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
     SUBROUTINE random_set_dist(r, low, high) BIND(C, NAME='random_set_dist')
       USE, INTRINSIC :: iso_c_binding

       IMPLICIT NONE

       TYPE (C_PTR), VALUE    :: r
       REAL (C_DOUBLE), VALUE :: low
       REAL (C_DOUBLE), VALUE :: high
     END SUBROUTINE random_set_dist
  END INTERFACE
  
  PUBLIC :: initialize_random

CONTAINS

  SUBROUTINE initialize_random(seed)
    USE omp_lib
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: seed
    INTEGER             :: num_threads
    INTEGER             :: thread_num

    num_threads = OMP_GET_MAX_THREADS()
    IF (.NOT. ALLOCATED(states)) THEN
       ALLOCATE(states(0:num_threads - 1))
    END IF

    !$OMP PARALLEL PRIVATE(thread_num)
    thread_num = OMP_GET_THREAD_NUM()
    states(thread_num) = random_construct_U(seed + thread_num)
    !$OMP END PARALLEL
  END SUBROUTINE initialize_random

  SUBROUTINE initialize_random_U(seed)
    USE omp_lib
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: seed

    state = random_construct_U(seed)

  END SUBROUTINE initialize_random_U

  SUBROUTINE initialize_random_N(seed)
    USE omp_lib
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: seed

    state = random_construct_N(seed)

  END SUBROUTINE initialize_random_N
  
  FUNCTION get_random()
    USE omp_lib
    IMPLICIT NONE

    REAL(rp)            :: get_random

    get_random = random_get_number(states(OMP_GET_THREAD_NUM()))
  END FUNCTION get_random

  FUNCTION get_random_U()
    USE omp_lib
    IMPLICIT NONE

    REAL(rp)            :: get_random_U

    get_random_U = random_get_number_U(state)
  END FUNCTION get_random_U

  FUNCTION get_random_N()
    USE omp_lib
    IMPLICIT NONE

    REAL(rp)            :: get_random_N

    get_random_N = random_get_number_N(state)
  END FUNCTION get_random_N

  SUBROUTINE get_randoms(nums)
    USE omp_lib
    IMPLICIT NONE

    REAL(rp), DIMENSION(:), INTENT(OUT) :: nums

    INTEGER                             :: i

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
    DO i = 1, SIZE(nums)
       nums(i) = get_random()
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE get_randoms

  SUBROUTINE set_random_dist(low, high)
    USE omp_lib
    IMPLICIT NONE

    REAL(rp), INTENT(IN) :: low
    REAL(rp), INTENT(IN) :: high

    !$OMP PARALLEL DEFAULT(SHARED)
    CALL random_set_dist(states(OMP_GET_THREAD_NUM()), low, high)
    !$OMP END PARALLEL

  END SUBROUTINE set_random_dist
  
  !SUBROUTINE initialize_random_mkl(seed)
  !    USE omp_lib
  !    IMPLICIT NONE

  !    INTEGER, INTENT(IN) :: seed
  !    INTEGER             :: num_threads
  !    INTEGER             :: thread_num
  !    INTEGER         :: errcode
  !    integer         :: brng

  !    brng=VSL_BRNG_MT19937
  !   brng=VSL_BRNG_MT2203

  !    num_threads = OMP_GET_MAX_THREADS()
  !    IF (.NOT. ALLOCATED(streams)) THEN
  !        ALLOCATE(states(0:num_threads - 1))
  !    END IF

!!$OMP PARALLEL PRIVATE(thread_num)
  !    thread_num = OMP_GET_THREAD_NUM()
  !    errcode=vslnewstream(streams(thread_num),brng,seed)
  !    errcode=vslnewstream(stream,brng,seed)
!!$OMP END PARALLEL
  !END SUBROUTINE

  !FUNCTION get_random_mkl()
  !    USE omp_lib
  !    IMPLICIT NONE
  !
  !    INTEGER, PARAMETER         :: n=8_idef
  !    REAL(rp),DIMENSION(n)            :: get_random_mkl
  !    INTEGER         :: errcode
  !    INTEGER         :: method

  !    real(rp)        :: mu,sigma

  !    method=VSL_RNG_METHOD_GAUSSIAN_ICDF

  !    mu=0._rp
  !    sigma=1._rp

  !    errcode = vdrnggaussian(method,streams(OMP_GET_THREAD_NUM()),n, &
  !         get_random_mkl,mu,sigma)
  !END FUNCTION

  !FUNCTION get_random_mkl_N(mu,sigma)
  !    USE omp_lib
  !    IMPLICIT NONE
  !
  !    INTEGER, PARAMETER         :: n=8_idef
  !    REAL(rp)            :: get_random_mkl_N
  !    REAL(rp),dimension(2)            :: buffer
  !    INTEGER         :: errcode
  !    INTEGER         :: method

  !    real(rp),intent(in)        :: mu,sigma

  !    method=VSL_RNG_METHOD_GAUSSIAN_BOXMULLER

  !    mu=0._rp
  !    sigma=1._rp

  !    errcode = vdrnggaussian(method,stream,2_idef, &
  !         buffer,mu,sigma)
  !    get_random_mkl_N=buffer(1)
  !END FUNCTION

  !FUNCTION get_random_mkl_U()
  !    USE omp_libvdrnggaussian
  !    IMPLICIT NONE
  !
  !    INTEGER, PARAMETER         :: n=8_idef
  !    REAL(rp)            :: get_random_mkl_U
  !    REAL(rp),dimension(2)           :: buffer
  !    INTEGER         :: errcode
  !    INTEGER         :: method


  !    method=VSL_RNG_METHOD_UNIFORM_STD_ACCURATE

  !    mu=0._rp
  !    sigma=1._rp

  !    errcode = vdrnguniform(method,stream,2_idef, &
  !         buffer,0._rp,1._rp)
  !    get_random_mkl_U=buffer(2)
  !END FUNCTION

END MODULE korc_random
