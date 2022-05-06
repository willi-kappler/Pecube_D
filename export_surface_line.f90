module m_export_surface_line
contains

subroutine export_surface_line(outer_step, num_of_points, num_recs, &
    model_depth, temperature_history_file, line_points)
  use m_logger
  use m_data_structures
  use m_read_time_temperature_history

  implicit none

  integer(4), intent(in) :: outer_step, num_of_points, num_recs

  character(*), intent(in) :: temperature_history_file

  real(8), intent(in) :: model_depth
  real(8), dimension(4), intent(in) :: line_points

  integer(4), parameter :: file_unit_temperature_history = 85
  integer(4) :: i, j

  real(8), dimension(:), allocatable :: ztime_local
  real(8), dimension(:,:), allocatable :: ztemp_local

  type(vector3D_t), dimension(:,:), allocatable :: vi_pos, vi_velo


  call read_time_temperature_history(num_recs, temperature_history_file, &
        outer_step, num_of_points, ztime_local, ztemp_local, vi_pos, vi_velo)

  open(file_unit_temperature_history, file="export_surface_line.txt", status="unknown")
  write(file_unit_temperature_history, *) "# point id, sub time step, time [myr], temperature [deg C], px, py, pz, vx, vy, vz"
  close(file_unit_temperature_history)

  call log_message("points of line: " + line_points)

  do i = 1, num_of_points
    do j = 1, num_recs
        ! TODO: filter points on the line
        write(file_unit_temperature_history, *) i, j, ztime_local(j), ztemp_local(i, j), &
          vi_pos(i, j)%x, vi_pos(i, j)%y, vi_pos(i, j)%z - model_depth, &
          vi_velo(i, j)%x, vi_velo(i, j)%y, vi_velo(i, j)%z
    enddo
  enddo

  deallocate(ztemp_local, ztime_local)
  deallocate(vi_pos, vi_velo)
end subroutine export_surface_line
end module m_export_surface_line
