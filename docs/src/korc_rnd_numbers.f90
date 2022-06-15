module korc_rnd_numbers
  !! @note Module with subrotuines for generating integer 
  !! and real random numbers.@endnote
  !! This subroutines were taken from Numerical Recipes in Fortran 90, 
  !! and provide a way for generating random numbers of 'better quality' 
  !! in a faster way than build-in Fortran random generators (at least  
  !! better than those of Fortran 77). For more details we refer the user 
  !! to Numerical Recipes in Fortran 90.

#ifdef __INTEL_COMPILER
  use ifport
#endif
  use korc_types

  IMPLICIT NONE

  ! Parameters and variables used in generator of uniform random numbers
  INTEGER(8), PARAMETER  :: iv = 4101842887655102017_8
  INTEGER(8), PARAMETER  :: iw = 1_8
  INTEGER(8), PARAMETER  :: a = 4294957665_8
  INTEGER(8), PARAMETER  :: b = 4294967295_8
  INTEGER(8), PARAMETER  :: d = 2862933555777941757_8
  INTEGER(8), PARAMETER  :: e = 7046029254386353087_8
  REAL(rp), PARAMETER    :: rcoeff = 5.42101086242752217E-20_rp


  TYPE, PRIVATE :: URAND
     INTEGER(8) :: u
     INTEGER(8) :: v
     INTEGER(8) :: w
  END TYPE URAND


  TYPE(URAND), PRIVATE :: urand_vars
  ! Parameters and variables used in generator of uniform random numbers


  INTERFACE u_random
     module procedure rand_int64,rand_int32,rand_real,rand_real_array
  END INTERFACE u_random


  PUBLIC :: init_u_random, u_random
  PRIVATE :: rand_int64,rand_int32,rand_real

CONTAINS

subroutine init_u_random(seed)
	INTEGER(8), INTENT(IN) :: seed
	INTEGER(8)             :: dummy_int64

	urand_vars%u = seed**urand_vars%v
	call rand_int64()
	urand_vars%v = urand_vars%u
	call rand_int64()
	urand_vars%w = urand_vars%v
	call rand_int64()
end subroutine init_u_random


subroutine rand_int64(irand)
	INTEGER(8), OPTIONAL, INTENT(OUT)  :: irand
	INTEGER(8)                         :: x


	urand_vars%u = urand_vars%u*d + e !u=u* 2862933555777941757LL + 7046029254386353087LL;

	urand_vars%v = IEOR(urand_vars%v,ISHFT(urand_vars%v,-17)) ! v^=v>>17;
	urand_vars%v = IEOR(urand_vars%v,ISHFT(urand_vars%v,31)) ! v^=v<<31;
	urand_vars%v = IEOR(urand_vars%v,ISHFT(urand_vars%v,-8)) ! v^=v>>8;

	urand_vars%w = a*IAND(urand_vars%w,b) + ISHFT(urand_vars%w,-32) ! w = 4294957665U*(w & 0xffffffff) + (w >> 32);

	x = IEOR(urand_vars%u,ISHFT(urand_vars%u,21)) ! Ullong x=u^(u<< 21);
	x = IEOR(x,ISHFT(x,-35)) ! x ^= x >> 35;
	x = IEOR(x,ISHFT(x,4)) ! x ^= x << 4;

	if (PRESENT(irand)) then
		irand = IEOR(x + urand_vars%v,urand_vars%w)
	end if
end subroutine rand_int64


subroutine rand_int32(irand32)
	INTEGER(4), INTENT(OUT)    :: irand32
	INTEGER(8)                 :: irand64

	call rand_int64(irand64)

	irand32 = INT(irand64,4)
end subroutine rand_int32


subroutine rand_real_array(rrand)
	REAL(rp), DIMENSION(:), INTENT(INOUT)  :: rrand
	INTEGER(8)                             :: irand64
    INTEGER                                :: ii ! Iterator

    do ii=1_idef,SIZE(rrand)
        	call rand_int64(irand64)
        rrand(ii) = rcoeff*REAL(irand64,rp) + 0.5_rp
    end do
end subroutine rand_real_array


subroutine rand_real(rrand)
	REAL(rp), INTENT(OUT)  :: rrand
	INTEGER(8)             :: irand64

 	call rand_int64(irand64)
    rrand = rcoeff*REAL(irand64,rp) + 0.5_rp
end subroutine rand_real


subroutine init_random_seed()
#ifdef PARALLEL_RANDOM
  use korc_random
#endif
  INTEGER, allocatable       :: seed(:)
  INTEGER(8), DIMENSION(8)   :: dt
  INTEGER(8)                 :: i
  INTEGER(8)                 :: istat
  INTEGER(8)                 :: pid
  INTEGER(4)                 :: n
  INTEGER(8)                 :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(default_unit_open, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(default_unit_open) seed
     close(default_unit_open)
     
  else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970_8) * 365_8 * 24_8 * 60_8 * 60_8 * 1000_8 &
             + dt(2) * 31_8 * 24_8 * 60_8 * 60_8 * 1000_8 &
             + dt(3) * 24_8 * 60_8 * 60_8 * 1000_8 &
             + dt(5) * 60_8 * 60_8 * 1000_8 &
             + dt(6) * 60_8 * 1000_8 &
             + dt(7) * 1000_8 &
             + dt(8)
     end if
     pid = getpid()
     write(output_unit_write,'("PID: ",I15)') pid
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
#ifdef PARALLEL_RANDOM
  call initialize_random(seed(1))
  call initialize_random_U(seed(1))
  call initialize_random_N(seed(1))

!  call initialize_random_mkl(seed(1))
#else
  call random_seed(put=seed)
#endif
contains

  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    INTEGER :: lcg
    INTEGER(8) :: s
    if (s == 0) then
       s = 104729_8
    else
       s = mod(s, 4294967296_8)
    end if
    s = mod(s * 279470273_8, 4294967291_8)
    lcg = int(mod(s, int(huge(0), 8)), kind(0))
  end function lcg
end subroutine init_random_seed

end module korc_rnd_numbers
