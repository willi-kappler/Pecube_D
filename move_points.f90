module m_move_points
contains
  subroutine move_points(px, py, pz, velo_info, file_unit, id, xmax, ymax)
    use m_find_velo
    use m_logger
    use m_data_structures

    implicit none

    real(8), intent(inout) :: px, py, pz
    real(8), intent(in) :: xmax, ymax
    real(8) :: vx, vy, vz

    integer(4), intent(in) :: id, file_unit

    type(velocity_info_t), intent(inout) :: velo_info

    ! if (px /= px) then
    !   call log_message("move_points.f90, px is NaN")
    !   call log_message("px: " + px + ", py: " + py + ", pz: " + pz)
    ! endif
    !
    ! if (py /= py) then
    !   call log_message("move_points.f90, py is NaN")
    !   call log_message("px: " + px + ", py: " + py + ", pz: " + pz)
    ! endif
    !
    ! if (pz /= pz) then
    !   call log_message("move_points.f90, pz is NaN")
    !   call log_message("px: " + px + ", py: " + py + ", pz: " + pz)
    ! endif

    call find_velo(px, py, pz, vx, vy, vz, velo_info, id)

     if (id == 1) then
       call log_message("move_points.f90: id: " + id + ", velo_info%dt: " + velo_info%dt)
       call log_message("velo_info%peclet: " + velo_info%peclet + ", velo_info%ntime: " + velo_info%ntime)
       call log_message("move_points.f90: px: " + px + ", py: " + py + ", pz: " + pz)
       call log_message("move_points.f90: vx: " + vx + ", vy: " + vy + ", vz: " + vz)
     endif

    px = px - (vx * velo_info%dt * velo_info%ntime)
    py = py - (vy * velo_info%dt * velo_info%ntime)
    pz = pz - (vz * velo_info%dt * velo_info%ntime)

    ! Important:
    ! Write new position after the point has been moved!
    ! Back tracking starts from this new moved position and moves the points up.
    write(file_unit) id, px, py, pz, vx, vy, vz

    ! Update velo_info

    velo_info%vz(id) = vz

     if (id == 1) then
       call log_message("move_points.f90: px2: " + px + ", py2: " + py + ", pz2: " + pz)
     endif

    ! Check bounds
    if (px < 0.0) then
      px = 0.0
    elseif (px > xmax) then
      px = xmax
    endif

    if (py < 0.0) then
      py = 0.0
    elseif (py > ymax) then
      py = ymax
    endif

    if (pz < 0.0) then
      pz = 0.0
    endif
  end subroutine move_points
end module m_move_points
