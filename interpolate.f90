  module m_interpolate
    use m_logger
    implicit none

    contains
  subroutine interpolate (t,tp,z,zp,n)

      implicit none

! interpolation routine

      integer(4) :: n, ip, i
      real(8) :: t(n), tp(n), z(n), zp(n)
      real(8) :: xp, x0, x1, r

      do ip=2,n-1
        xp=zp(ip)
        do i=1,n-1
          x0=z(i)
          x1=z(i+1)
            if ((xp-x0) * (xp-x1) <= 0.0) then
              r = (xp-x0) / (x1-x0)
              tp(ip) = t(i) + (t(i+1) - t(i)) * r
              goto 1
            endif
        enddo
        call log_message("problem in interpolation")
        call log_message("ip: " + ip + ", i: " + i + ", x0: " + x0 + ", x1: " + x1)
        stop
1       continue
      enddo

      return
  end subroutine interpolate
  end module m_interpolate

