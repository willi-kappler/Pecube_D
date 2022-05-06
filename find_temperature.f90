module m_find_temperature

type point_info
  real(8) :: distance
  real(8) :: factor
  real(8) :: temperature
  integer(4) :: index
end type point_info

contains

subroutine find_temperature(px, py, pz, tx, ty, tz, t_in, number_of_elements, radius, t_out, debug)
  use m_logger

  implicit none

  ! WK: TODO: interpolate onto a regular 3d grid and do a tri-linear interpolation for the given point
  ! this will be much faster.

  integer(4), intent(in) :: number_of_elements

  real(8), intent(in) :: px, py, pz, radius
  real(8), intent(in) :: tx(number_of_elements), ty(number_of_elements), tz(number_of_elements), t_in(number_of_elements)

  real(8), intent(out) :: t_out

  logical, intent(in) :: debug

  integer(4), parameter :: max_points = 10

  real(8) :: d, factor_sum
  integer(4) :: i, j, k

  type(point_info) :: p_info(max_points)

  do i = 1, max_points
    p_info%distance = huge(d)
    p_info%factor = 0.0
    p_info%temperature = 0.0
    p_info%index = 1
  enddo

  factor_sum = 0.0
  t_out = 0.0

  do i = 1, number_of_elements
      d = sqrt(((px - tx(i))**2) + ((py - ty(i))**2) + ((pz - tz(i))**2))

      do j = 1, max_points
        if (d < p_info(j)%distance) then
          do k = max_points, j + 1, -1
            p_info(k) = p_info(k - 1)
          enddo

          p_info(j)%distance = d
          p_info(j)%temperature = t_in(i)
          p_info(j)%index = i

          exit
        endif
      enddo
  enddo

  do i = 1, max_points
    p_info(i)%factor = exp(-p_info(i)%distance / radius)
    !if (p_info(i)%distance > radius) then
    !  p_info(i)%factor = 0.0
    !else
    !  p_info(i)%factor = 1 - (p_info(i)%distance / radius)
    !endif
    factor_sum = factor_sum + p_info(i)%factor
    t_out = t_out + (p_info(i)%temperature * p_info(i)%factor)
  enddo

  t_out = t_out / factor_sum

  if (debug .and. px == 0.0 .and. py == 0.0) then
    call log_message("find_temperature.f90:")
    call log_message("px: " + px + ", py: " + py + ", pz: " + pz)

    do i = 1, max_points
      call log_message("i: " + i + ", tx: " + tx(p_info(i)%index) + ", ty: " + ty(p_info(i)%index) + ", tz: " + tz(p_info(i)%index))
      call log_message("index: " + p_info(i)%index + ", distance: " + p_info(i)%distance + ", factor: " + p_info(i)%factor + ", temperature: " + p_info(i)%temperature)
    enddo

    call log_message("t_out: " + t_out + ", factor_sum: " + factor_sum)
  endif

end subroutine find_temperature

end module m_find_temperature
