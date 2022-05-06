module m_find_dt
    use m_logger
contains

subroutine find_dt (zl, diffusivity, nsurf, zsurf, zsurfp, &
                 Pecletz, timesurf, timesurfp, istep, eps, &
                 dt, ntime, istatic)

      implicit none

      real(8), intent(in) :: zsurf(nsurf), zsurfp(nsurf)
      real(8), intent(in) :: zl, diffusivity, Pecletz, timesurf, timesurfp, eps
      real(8), intent(out) :: dt
      integer(4), intent(in) :: istep, nsurf
      integer(4), intent(out) :: ntime, istatic

      real(8) :: dt1, dt2, dt3, dzmax, Pesurf
      integer(4) :: i

      ! Check if time values are valid

      if (istep /= 0) then
        if (timesurf == timesurfp) then
            call log_message("Error in find_dt: timesurf and timesurfp are identical: " + timesurf + ", " + timesurfp)
            call log_message("Please check your input file (time steps in input 12a)")
            call log_message("itep: " + istep + ", zl: " + zl + ", ntime: " + ntime)
            call logger_flush()
            stop
        endif
      endif

! first constraint on time step from conduction

      dt1 = zl**2 / diffusivity / 100.0

      !call log_message("find_dt.f90: dt1: " + dt1)

! second constrain on time step from advection

      dt2 = dt1
      if (abs(Pecletz) > eps) then
          dt2 = zl /abs(Pecletz) / 100.0
      endif

      !call log_message("find_dt.f90: dt2: " + dt2)

! third constrain on time step from surface lowering

      dt3 = dt1
      if (istep /= 0) then
          dzmax = 0.0
          do i = 1, nsurf
              dzmax = max(dzmax, zsurfp(i) - zsurf(i))
          enddo
          Pesurf = dzmax / (timesurf - timesurfp)
         if (abs(Pesurf)> eps) then
            dt3 =zl / Pesurf / 5.0
         endif
      endif

      !call log_message("find_dt.f90: dzmax: " + dzmax)
      !call log_message("find_dt.f90: timesurf: " + timesurf)
      !call log_message("find_dt.f90: timesurfp: " + timesurfp)
      !call log_message("find_dt.f90: Pesurf: " + Pesurf)
      !call log_message("find_dt.f90: dt3: " + dt3)


! find optimum time step and number of steps

      dt = min(dt1, dt2)
      dt = min(dt, dt3)
      if (istep /= 0) then
          dt = min(dt, timesurf - timesurfp)
      endif

      !call log_message("find_dt.f90: dt: " + dt)

      ntime = int((timesurf - timesurfp) / dt)
      ntime = ntime + 1
      dt = (timesurf - timesurfp) / ntime
      istatic = 0

      if (istep == 0) then
        ntime = 1
        dt = 0.0
        istatic = 1
      endif

      !call log_message("ntime/dt/dt1/dt2/dt3= " + ntime + ", " + dt + ", " + dt1 + ", " + dt2 + ", " + dt3)
      !call log_message("timesurf: " + timesurf + ", timesurfp: " + timesurfp)


      return
      end subroutine find_dt
end module m_find_dt
