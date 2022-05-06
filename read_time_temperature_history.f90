module m_read_time_temperature_history
contains

subroutine read_time_temperature_history(num_recs, temperature_history_file, &
      outer_step1, num_of_points, ztime_local, ztemp_local, vi_pos, vi_velo)
  use m_logger
  use m_data_structures

  implicit none

  integer(4), intent(in) :: num_recs, outer_step1, num_of_points

  character(*), intent(in) :: temperature_history_file

  integer(4), parameter :: file_unit_temperature_history = 85
  integer(4) :: start_index, current_step, i, j, ntime
  integer(4) :: current_time_record, current_temp_record, sub_step
  integer(4) :: num_of_points_file, outer_step2

  real(8), dimension(:), allocatable, intent(out) :: ztime_local
  real(8), dimension(:,:), allocatable, intent(out) :: ztemp_local

  type(vector3D_t), dimension(:,:), allocatable, intent(out) :: vi_pos, vi_velo


  allocate(ztime_local(num_recs), ztemp_local(num_of_points, num_recs), &
    vi_pos(num_of_points, num_recs), vi_velo(num_of_points, num_recs))

  open(file_unit_temperature_history, file=trim(temperature_history_file)//".bin", &
      status="old", form="unformatted", access="stream")

  read(file_unit_temperature_history) num_of_points_file, outer_step2

  if (outer_step1 /= outer_step2) then
    call log_message("Error reading in temperature history: outer_step1 != outer_step2")
    call log_message("outer_step1: " + outer_step1 + ", outer_step2: " + outer_step2)
    call log_message("error code: ee4ef42d65a196ed026a1ae73850cfc9")
    error stop 1
  endif

  current_time_record = 0
  current_temp_record = 0
  ntime = 0
  start_index = num_recs

  do current_step = outer_step1, 1, -1
    read(file_unit_temperature_history) i, ntime

    ! calculate correct index. The order of the temperature history on disk is not what is
    ! needed for the age calculation, so we need to reorder it here:
    start_index         = start_index - ntime
    current_time_record = start_index
    current_temp_record = start_index

    if (i /= current_step) then
      call log_message("Error reading in temperature history: current_step != i")
      call log_message("current_step: " + current_step + ", i: " + i)
      call log_message("error code: 327df24acae0486adad5f770edde4a1e")
      error stop 1
    endif

    do sub_step = 1, ntime
      current_time_record = current_time_record + 1

      read(file_unit_temperature_history) j, ztime_local(current_time_record)

      if (j /= sub_step) then
        call log_message("Error reading in temperature history: j != sub_step")
        call log_message("j: " + j + ", sub_step: " + sub_step)
        call log_message("error code: 570984767d5d0d36426e5a087a42add3")
        error stop 1
      endif
    enddo

    do sub_step = 1, ntime
      current_temp_record = current_temp_record + 1

      do i = 1, num_of_points
        read(file_unit_temperature_history) j, ztemp_local(i, current_temp_record), &
          vi_pos(i, current_temp_record)%x, vi_pos(i, current_temp_record)%y, &
          vi_pos(i, current_temp_record)%z, vi_velo(i, current_temp_record)%x, &
          vi_velo(i, current_temp_record)%y, vi_velo(i, current_temp_record)%z

        if (j /= i) then
          call log_message("Error reading in temperature history: j != i")
          call log_message("j: " + j + ", i: " + i)
          call log_message("error code: 9bb8ec51607f8c5e1f709fb5975c217c")
          error stop 1
        endif
      enddo
    enddo
  enddo ! current_step = outer_step1, 1, -1

  if (current_time_record /= ntime) then
    call log_message("Error reading in temperature history: current_time_record != ntime")
    call log_message("current_time_record: " + current_time_record + ", ntime: " + ntime)
    call log_message("error code: 713bb63f46d666a0105e3e598c9b929b")
    error stop 1
  endif

  if (current_temp_record /= ntime) then
    call log_message("Error reading in temperature history: current_temp_record != num_recs")
    call log_message("current_temp_record: " + current_temp_record + ", ntime: " + ntime)
    call log_message("error code: 17efe049ef7dc0f06c5224a31ebb026f")
    error stop 1
  endif

  close(file_unit_temperature_history)

  ! Smooth temperature curve a bit to handle numeric instabilities
  do i = 1, num_of_points
    do j = 3, num_recs - 2
       ztemp_local(i, j) = (ztemp_local(i, j - 2) + ztemp_local(i, j - 1) + ztemp_local(i, j) + ztemp_local(i, j + 1) + ztemp_local(i, j + 2)) / 5.0
    enddo
  enddo




end subroutine read_time_temperature_history

end module m_read_time_temperature_history
