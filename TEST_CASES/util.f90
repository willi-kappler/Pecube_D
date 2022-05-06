module m_util
  implicit none

  type :: test_case
    character(256) :: name = "unknown"
    procedure(), nopass, pointer :: run_test => null()
  end type test_case

  interface assert
      module procedure assert_string
      module procedure assert_integer
      module procedure assert_logical
      module procedure assert_real8
  end interface assert

  private

  real(8), parameter :: precision = 0.00000001

  public assert
  public test_case

contains

  subroutine assert_string(actual, expected, result)
      implicit none

      character(*), intent(in) :: actual, expected

      logical, intent(inout) :: result
      logical :: original_result

      original_result = result
      result = (actual == expected)

      if (.not. result) then
          print "(A)", "assertion (string) failed: actual: '" // trim(adjustl(actual)) // &
          "', expected: '" // trim(adjustl(expected)) // "'"
      endif

      result = result .and. original_result

  end subroutine assert_string

  subroutine assert_integer(actual, expected, result)
      implicit none

      integer(4), intent(in) :: actual, expected

      logical, intent(inout) :: result
      logical :: original_result

      original_result = result
      result = (actual == expected)

      if (.not. result) then
          print "(A, I5, A, I5, A)", "assertion (integer) failed: actual: '", actual, "', expected: '", expected, "'"
      endif

      result = result .and. original_result

  end subroutine assert_integer

  subroutine assert_logical(actual, result)
      implicit none

      logical, intent(in) :: actual

      logical, intent(inout) :: result
      logical :: original_result

      original_result = result
      result = actual

      if (.not. result) then
          print "(A, L1, A, L1, A)", "assertion (logical) failed: actual: '", actual, "'"
      endif

      result = result .and. original_result

  end subroutine assert_logical

  subroutine assert_real8(actual, expected, result)
      implicit none

      real(8), intent(in) :: actual, expected

      logical, intent(inout) :: result
      logical :: original_result

      original_result = result
      result = (abs(actual - expected) < precision)

      if (.not. result) then
          print "(A, F30.10, A, F30.10, A)", "assertion (real8) failed: actual: '", actual, "', expected: '", expected, "'"
      endif

      result = result .and. original_result

  end subroutine assert_real8

end module m_util
