module test_hpc
  use mpi
  use fruit
  use korc_types
  implicit none

contains

  subroutine test_mpi_initialization(params)
    TYPE(KORC_PARAMS), INTENT(INOUT) :: params
    integer :: ierror
    integer :: size, rank
    integer :: size_k, rank_k

    size_k=params%mpi_params%nmpi
    rank_k=params%mpi_params%rank

    call MPI_COMM_SIZE (MPI_COMM_WORLD, size, ierror)
    call MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierror)

    call assert_equals(size_k,size)
    call assert_equals(rank_k,rank)
    
  end subroutine test_mpi_initialization

end module test_hpc



