module m_bivar
    contains

  subroutine idbvip ( num_input, xin, yin, zin, num_output, xout, yout, zout )
! Better and easier interpolation: more stable, better to understand and extend

  implicit none

  integer(4) :: num_input, num_output

  real(8), dimension(num_input), intent(in) :: xin, yin, zin
  real(8), dimension(num_output), intent(in) :: xout, yout
  real(8), dimension(num_output), intent(out) :: zout

  ! local variables
  integer(4) :: i, j, k
  integer(4), dimension(4) :: current_index

  real(8) :: dist, dx2, dy2
  real(8), dimension(4) :: current_dist

  do i=1,num_output
     current_index = 1
     current_dist = huge(dist)

     do j=1,num_input
        dx2 = (xout(i) - xin(j))**2
        dy2 = (yout(i) - yin(j))**2
        dist = sqrt(dx2 + dy2)

        if (dist < current_dist(1)) then
           do k=4,2,-1
              current_dist(k) = current_dist(k-1)
              current_index(k) = current_index(k-1)
           enddo

           current_dist(1) = dist
           current_index(1) = j
        else if (dist < current_dist(2)) then
           do k=4,3,-1
              current_dist(k) = current_dist(k-1)
              current_index(k) = current_index(k-1)
           enddo

           current_dist(2) = dist
           current_index(2) = j

        else if (dist < current_dist(3)) then
           current_dist(4) = current_dist(3)
           current_index(4) = current_index(3)
           current_dist(3) = dist
           current_index(3) = j
        else if (dist < current_dist(4)) then
           current_dist(4) = dist
           current_index(4) = j
        endif
     enddo ! num_input

     zout(i) = 0.0

     do k=1,4
        zout(i) = zout(i) + zin(current_index(k))
     enddo

     zout(i) = zout(i) / 4.0

  enddo ! num_output
end subroutine idbvip

end module m_bivar
