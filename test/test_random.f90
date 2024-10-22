module test_random
  use korc_random
  use mpi
  use fruit
  implicit none

contains

  SUBROUTINE test_random_auto
    IMPLICIT NONE

!  Local Variables
    CLASS (random_U_context), POINTER   :: uniform => null()
    INTEGER                             :: mpierr
    INTEGER                             :: i
    INTEGER                             :: rank
    INTEGER                             :: size
    INTEGER                             :: localsize
    REAL(rp), DIMENSION(:), ALLOCATABLE :: buffer
    INTEGER, DIMENSION(:), ALLOCATABLE  :: counts
    INTEGER, DIMENSION(:), ALLOCATABLE  :: offsets
    REAL(rp)                            :: base
    REAL(rp)                            :: test

!  Local parameters
    INTEGER, PARAMETER                  :: totalsize = 10000
    INTEGER, PARAMETER                  :: window = totalsize/2

!  Start of executable code.
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, mpierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpierr)

    IF (rank .eq. 0) THEN
      WRITE (*,*)
      WRITE (*,*) "Starting Random Test"
    END IF

    uniform => random_U_context_construct(0, rank)
    CALL uniform%set(-1.0d0, 1.0d0)

    localsize = totalsize/size
    IF (rank .lt. MOD(totalsize, size)) THEN
      localsize = localsize + 1
    ENDIF

    IF (rank .eq. 0) THEN
      ALLOCATE(buffer(totalsize))
      ALLOCATE(counts(size))
      ALLOCATE(offsets(size))
      CALL uniform%get_array(buffer(1:localsize))
    ELSE
      ALLOCATE(buffer(localsize))
      CALL uniform%get_array(buffer)
    END IF

    CALL MPI_GATHER(localsize, 1, MPI_INTEGER,                                 &
                    counts, 1, MPI_INTEGER,                                    &
                    0, MPI_COMM_WORLD, mpierr)
    IF (rank .eq. 0) THEN
      offsets(1) = 0
      DO i = 2, size
        offsets(i) = offsets(i - 1)+counts(i - 1)
      END DO
    END IF
    CALL MPI_GATHERV(buffer, localsize, MPI_DOUBLE,                            &
                     buffer, counts, offsets, MPI_DOUBLE,                      &
                     0, MPI_COMM_WORLD, mpierr)

    IF (rank .eq. 0) THEN
      DEALLOCATE(counts)
      DEALLOCATE(offsets)

!  Check against 10% of the first peak. We could go lower but I don't want to
!  trigger a test failure on noise.
      base = autocorrelation(buffer, 0)*0.1
      DO i = 1, window
        test = autocorrelation(buffer, i)
        CALL assert_equals(test .gt. base, .false.)
      END DO
    END IF

    DEALLOCATE(buffer)

  END SUBROUTINE

  FUNCTION autocorrelation(sequence, offset)

    IMPLICIT NONE

!  Declare arguments.
    REAL(rp) :: autocorrelation
    REAL(rp), DIMENSION(:) :: sequence
    INTEGER                :: offset

!  Start of executable code.
    autocorrelation = DOT_PRODUCT(sequence(:SIZE(sequence) - offset),          &
                                  sequence(offset + 1:))                       &
                    / (SIZE(sequence) - offset)
  END FUNCTION

end module
