module m_find_velo
    contains


! Subroutine that calls function to calculate velocities in the x, y, and z directions

    subroutine find_velo (x, y, z, vx, vy, vz, velo_info, id)

      use m_geometry
      use m_logger
      use m_data_structures

      implicit none

      real(8), intent(in) :: x, y, z
      real(8), intent(out) :: vx, vy, vz

      real(8) :: xx0, yy0, zz0, xx, yy, zz, xnorm, xxn, yyn, zzn

      integer(4) :: v_niter
      integer(4), intent(in) :: id

      type(velocity_info_t), intent(in) :: velo_info

! this extra layer between Pecube and the velocity function geometry ensures
! that the advection of the "lagrangian particles" for which we track the thermal
! history (to calculate the age) is second order accurate; it uses a mid-point
! algorithm.
!      implicit real*8 (a-h,o-z)

    if (x /= x) then
        call log_message("find_velo.f90: x is NaN")
        call log_message("x:" + x + ", y:" + y + ", z:" + z)
        stop
    endif

    if (y /= y) then
        call log_message("find_velo.f90: y is NaN")
        call log_message("x:" + x + ", y:" + y + ", z:" + z)
        stop
    endif

    if (z /= z) then
        call log_message("find_velo.f90: z is NaN")
        call log_message("x:" + x + ", y:" + y + ", z:" + z)
        stop
    endif

      xx0 = x
      yy0 = y
      zz0 = z
      xx = xx0
      yy = yy0
      zz = zz0
      v_niter = 0


1     call geometry (xx, yy, zz, vx, vy, vz, velo_info, id)
      xxn = xx0 + velo_info%dt * vx / 2.0
      yyn = yy0 + velo_info%dt * vy / 2.0
      zzn = zz0 + velo_info%dt * vz / 2.0


    ! if (vx /= vx) then
    !     print *, "vx is NaN (find_velo, line 69)"
    !     stop
    ! endif

    xnorm = sqrt((xxn-xx)**2 + (yyn-yy)**2 + (zzn-zz)**2)
    xx = xxn
    yy = yyn
    zz = zzn
    v_niter = v_niter + 1

    ! WK: debug
    ! if (vx < 0.0) then
    !   call log_message("find_velo.f90")
    !   call log_message("x: " + x + ", y: " + y + ", z:" + z)
    !   call log_message("vx: " + vx + ", vy: " + vy + ", vz:" + vz)
    ! endif


! In some cases (low resolution mesh) this mid-point algorithm cannot converge
! we assume a simple explicit estimate for the velocity
    if (v_niter > 10) then
      call geometry (xx0, yy0, zz0, vx, vy, vz, velo_info, id)
      ! if (vx /= vx) then
      !   print *, "vx is NaN (find_velo, line 86)"
      !   stop
      ! endif

      return
    endif

    if (xnorm > velo_info%zl * 1.e-6) goto 1
    return
  end subroutine find_velo
end module m_find_velo
