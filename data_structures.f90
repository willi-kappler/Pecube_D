module m_data_structures
  implicit none

  type vector3D_t
    real(8) :: x, y, z
  end type vector3D_t

  ! WK: data structure for system clock information
  ! Note: we use the more accurate system_clock suboutine, since
  ! cpu_time does not count the real time. If we use multiple
  ! CPUs / cores a call to cpu_time returns a different result from
  ! what we want
  type :: system_clock_info_t
      integer(4) :: sys_count, sys_count_rate, sys_count_max
  end type system_clock_info_t

  type key_value_t
    character(10) :: key
    real(8) :: value
  end type key_value_t

  type velocity_info_t
    real(8) :: zl, x1f, y1f, x2f, y2f, def, dif, theta, phi, mft_ratein
    real(8) :: mbt_ratein, mct_ratein, stf_ratein, Peclet, Peclet2, dt
    real(8), dimension(:), allocatable :: vz
    integer(4) :: geoflag, ntime
  end type velocity_info_t

  type age_info_t
    real(8), dimension(:, :), allocatable :: all_ages, ftld
    real(8), dimension(:), allocatable :: ftldmean, ftldsd
  end type age_info_t

contains

subroutine allocate_age_info(age_info_array, num_of_steps, num_of_points, num_of_flags)
  implicit none

  integer(4), intent(in) :: num_of_steps, num_of_points, num_of_flags
  integer(4) :: i

  type(age_info_t), dimension(:), allocatable :: age_info_array

  if (allocated(age_info_array)) then
    ! Nothing to do, already allocated
    return
  endif

  allocate(age_info_array(num_of_steps))

  do i = 1, num_of_steps
    allocate(age_info_array(i)%all_ages(num_of_points, num_of_flags))
    age_info_array(i)%all_ages = 0.0

    allocate(age_info_array(i)%ftld(num_of_points, 17))
    age_info_array(i)%ftld = 0.0

    allocate(age_info_array(i)%ftldmean(num_of_points))
    age_info_array(i)%ftldmean = 0.0

    allocate(age_info_array(i)%ftldsd(num_of_points))
    age_info_array(i)%ftldsd = 0.0
  enddo
end subroutine allocate_age_info

subroutine deallocate_age_info(age_info_array, num_of_steps)
  implicit none

  integer(4), intent(in) :: num_of_steps
  integer(4) :: i

  type(age_info_t), dimension(:), allocatable :: age_info_array

  if (.not. allocated(age_info_array)) then
    ! Nothing to do, not allocated
    return
  endif

  do i = 1, num_of_steps
    deallocate(age_info_array(i)%all_ages)
    deallocate(age_info_array(i)%ftld)
    deallocate(age_info_array(i)%ftldmean)
    deallocate(age_info_array(i)%ftldsd)
  enddo

  deallocate(age_info_array)
end subroutine deallocate_age_info

end module m_data_structures
