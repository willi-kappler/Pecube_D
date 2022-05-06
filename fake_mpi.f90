! Fake module if you want to compile without MPI support
! TODO: Add more fake functions

module mpi

  private

  integer(4), parameter :: MPI_SUCCESS = 1
  integer(4), parameter :: MPI_STATUS_SIZE = 1
  integer(4), parameter :: MPI_ANY_SOURCE = 1
  integer(4), parameter :: MPI_ANY_TAG = 1
  integer(4), parameter :: MPI_COMM_WORLD = 1
  integer(4), parameter :: MPI_LOGICAL = 1
  integer(4), parameter :: MPI_DOUBLE_PRECISION = 1
  integer(4), parameter :: MPI_INTEGER = 1
  integer(4), parameter :: MPI_INTEGER8 = 1
  integer(4), parameter :: MPI_SOURCE = 1




  public MPI_SUCCESS
  public MPI_STATUS_SIZE
  public MPI_ANY_SOURCE
  public MPI_ANY_TAG
  public MPI_COMM_WORLD
  public MPI_LOGICAL
  public MPI_DOUBLE_PRECISION
  public MPI_INTEGER
  public MPI_INTEGER8
  public MPI_SOURCE

  public MPI_INIT



contains

subroutine MPI_INIT(error_result)
  implicit none

  integer(4), intent(out) :: error_result

  error_result = MPI_SUCCESS

end subroutine MPI_INIT

end module mpi
