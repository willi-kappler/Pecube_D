module m_tec_mat_output
contains


! Subroutine to output x, y, and z postions along with node number,
! temperature value, and x, y, and z velocities

      subroutine tec_mat_output (x,y,z,t,nnode,icon,nelem,mpe,run,&
                     file_id1,xlonmin,xlatmin,nrun,&
                     dynamic_thermal_conductivity,&
                     has_dynamic_thermal_conductivity,velo_info)

! This subroutine formats the output from Pecube so that it can easily
! be read into Tecplot or Matlab.

    use m_find_velo
    use m_logger
    use m_dynamic_thermal_conductivity
    use m_data_structures

! Declare variables
      implicit none

      integer(4), intent(in) :: nrun, nnode, nelem, mpe
      integer(4), intent(in), dimension(mpe, nelem) :: icon

      real(8), intent(in) :: xlonmin, xlatmin
      real(8), intent(in), dimension(nnode) :: x, y, z, t

      type(thermal_conductivity_t), intent(in) :: dynamic_thermal_conductivity

      logical, intent(in) :: has_dynamic_thermal_conductivity

      character, intent(in) ::  run*100, file_id1*4

      character :: mpe_char*10

      real(8) :: vx, vy, vz, z_mod, current_diffusivity
      integer(4) :: i, ie, k

      type(velocity_info_t), intent(in) :: velo_info

! Open files to create
      open(77,file=run(1:nrun)//'/Temps_tec'//file_id1//'.dat',status='unknown')
!      open(78,file=run//'/Temps_mat.dat',status='unknown')


! Write headers to files
      write (77,*) 'TITLE = "Pecube temperature output"'
      if (has_dynamic_thermal_conductivity) then
        write (77,'(a)') 'VARIABLES = "x (km)" "y (km)" "z (km)" "node" "temperature (C)" "U (mm/yr)" "V (mm/yr)" "W (mm/yr)" "diffusivity"'
      else
        write (77,'(a)') 'VARIABLES = "x (km)" "y (km)" "z (km)" "node" "temperature (C)" "U (mm/yr)" "V (mm/yr)" "W (mm/yr)"'
      endif
      write (77,*) 'ZONE T = "Pecube"'
      write (77,*) 'n=',nnode,', e=',nelem,'et=brick,f=fepoint'

! Write temperatures, velocities, and positions to files
      call log_message("tec_mat_output.f90")
      call log_message("xlonmin: " + xlonmin + ", xlatmin: " + xlatmin + ", zl: " + velo_info%zl)
      do i = 1, nnode

        if (x(i) /= x(i)) then
             print *, "x(i) is NaN (tec_mat_output)"
            stop
        endif

        call find_velo(x(i), y(i), z(i), vx, vy, vz, velo_info, 0)

! Calculates elevation from model thickness to node elevation
        z_mod = z(i) - velo_info%zl
! Write thermal field data to file
        if (has_dynamic_thermal_conductivity) then
          call get_conductivity(dynamic_thermal_conductivity, &
            velo_info%zl - z(i), current_diffusivity)
          write (77,'(9f16.3)') x(i) + xlonmin, y(i) + xlatmin, z_mod, &
            dble(i), t(i), vx, vy, vz, current_diffusivity
        else
          write (77,'(8f16.3)') x(i) + xlonmin, y(i) + xlatmin, z_mod, &
            dble(i), t(i), vx, vy, vz
        endif

!        stop

      enddo

      write (mpe_char,'(i10)') mpe

! Write out node-element connectivities for tecplot version
      do ie = 1, nelem
        write (77,'('//mpe_char//'i10)') (icon(k,ie),k=1,mpe)
      enddo

! Close files
      close(77)
      end subroutine
end module m_tec_mat_output
