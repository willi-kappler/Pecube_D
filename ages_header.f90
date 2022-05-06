module m_ages_header
  private

  integer(4), parameter :: file_unit_ages = 92

  public ages_header_surface
  public ages_header_borehole

  contains
  subroutine ages_header_surface (num_of_points, xlonmin, xsurf, &
                xlatmin, ysurf, zl, xdepth, ydepth, zdepth, &
                age_info, file_prefix, nstep, &
                header_info, thermflag, age_flags, &
                nx, num_of_elements, num_of_flags)

    ! Subroutine for Pecube that writes ages from every time step to
    ! Tecplot formatted files called Ages_tec000.dat, Ages_tec001.dat, etc

    use m_logger
    use m_data_structures

    implicit none

    integer(4) :: num_of_flags, num_of_points, nstep, i, k

    real(8) :: xlonmin,xlatmin,zl
    real(8) :: xsurf(num_of_points),ysurf(num_of_points)
    real(8), intent(in) :: header_info(6)
    real(8), dimension(nstep, num_of_points) :: xdepth, ydepth, zdepth
    integer(4), intent(in) :: age_flags(num_of_flags), nx, num_of_elements
    integer(4), intent(in) :: thermflag(nstep + 1)
    integer(4) :: counter

    type(age_info_t) :: age_info(nstep)
    character(*) :: file_prefix
    character(4) :: file_id

    do k = 1, nstep
      call log_message("ages_header.f90, thermflag(" + k + ") = " + &
          thermflag(k))
      if (thermflag(k) == 1) then
        write(file_id,"(i4.4)") k

        open(file_unit_ages, file=trim(file_prefix)//file_id//".dat", status="unknown")

        call log_message("write ages: "//trim(file_prefix)//file_id//".dat")

        call write_header(header_info, age_flags, num_of_flags, num_of_points, num_of_elements, &
          'VARIABLES = "surface point id" "surface coord x (km)" "surface coord y (km)" "real world x (km)" "real world y (km)" "real topography z (km)"')

        do i = 1, num_of_points
          write (file_unit_ages,'(I10, 5F12.4)',advance="no") i, xsurf(i), ysurf(i), xdepth(k, i) + xlonmin, ydepth(k, i) + xlatmin, zdepth(k, i) - zl

          call write_content(num_of_flags, age_flags, age_info(k)%all_ages(i, :))
        enddo

        counter=0
        do i=1,num_of_elements
          counter=counter+1
          if (mod(counter,nx).eq.0) counter=counter+1
          write (file_unit_ages,*) counter+nx,counter+nx+1,counter+1,counter
        enddo

        close(file_unit_ages)

      endif
    enddo
  end subroutine ages_header_surface

  subroutine ages_header_borehole(num_of_points, offset_x, offset_y, offset_z, &
        borehole_ages_points, age_info, file_prefix, num_of_steps, header_info, &
        therm_his_val, age_flags, num_of_flags)
    use m_data_structures

    implicit none

    integer(4), intent(in) :: num_of_points, num_of_steps, num_of_flags
    integer(4), intent(in) :: therm_his_val(num_of_steps + 1), age_flags(num_of_flags)

    real(8), intent(in) :: offset_x, offset_y, offset_z, header_info(6)

    character(*), intent(in) :: file_prefix

    type(age_info_t), intent(in) :: age_info(num_of_steps)
    type(vector3D_t), intent(in) :: borehole_ages_points(num_of_points)

    integer(4) :: i, j
    character(4) :: file_id


    do i = 1, num_of_steps
      if (therm_his_val(i) == 1) then
        write(file_id, "(i4.4)") i

        open(file_unit_ages, file=trim(file_prefix)//file_id//".dat", status="unknown")

        call write_header(header_info, age_flags, num_of_flags, num_of_points, 0, &
          'VARIABLES = "borehole point id" "borehole x (km)" "borhole y (km)" "borehole z (km)"')

        do j = 1, num_of_points
          write (file_unit_ages,'(I10, 3F12.4)',advance="no") j, borehole_ages_points(j)%x + offset_x, &
              borehole_ages_points(j)%y + offset_y, borehole_ages_points(j)%z - offset_z

          call write_content(num_of_flags, age_flags, age_info(i)%all_ages(j, :))
        enddo

        close(file_unit_ages)
      endif
    enddo ! i = 1, num_of_steps
  end subroutine ages_header_borehole

  subroutine write_header(header_info, age_flags, num_of_flags, num_of_points, num_of_elements, header_text)
    use m_age_algorithms

    implicit none

    integer(4), intent(in) :: num_of_flags, age_flags(num_of_flags)
    integer(4), intent(in) :: num_of_points, num_of_elements

    real(8), intent(in) :: header_info(6)

    character(*), intent(in) :: header_text

    character(5) :: as, am, al, eUl, eUm, eUh

    ! Writes the various grain sizes and radiation levels to
    ! character arrays printed to ages output Tecplot file
    write(as, '(f5.1)') header_info(1)
    write(am, '(f5.1)') header_info(2)
    write(al, '(f5.1)') header_info(3)
    write(eUl, '(f5.1)') header_info(4)
    write(eUm, '(f5.1)') header_info(5)
    write(eUh, '(f5.1)') header_info(6)

    ! Writes Tecplot header in Ages_tec*.dat file
    write(file_unit_ages, *) 'TITLE = "Pecube Ages"'
    write(file_unit_ages, '(A)', advance="no") header_text

    if (has_ahe_1(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "AHe Age - Farley, 2000 (Ma)"'
    if (has_ahe_2(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "AHe Age - a='//as//' um (Ma)"'
    if (has_ahe_3(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "AHe Age - a='//am//' um (Ma)"'
    if (has_ahe_4(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "AHe Age - a='//al//' um (Ma)"'
    if (has_ahe_5(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "AHe Age - eU='//eUl//' ppm (Ma)"'
    if (has_ahe_6(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "AHe Age - eU='//eUm//' ppm (Ma)"'
    if (has_ahe_7(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "AHe Age - eU='//eUh//' ppm (Ma)"'
    if (has_aft(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "AFT Age (Ma)"'
    if (has_kfeldar(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "KFeldAr Age (Ma)"'
    if (has_zhe(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "ZHe Age (Ma)"'
    if (has_zhe_low_damage(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "ZHe low Age (Ma)"'
    if (has_zhe_med_damage(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "ZHe med Age (Ma)"'
    if (has_zhe_high_damage(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "ZHe high Age (Ma)"'
    if (has_zft(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "ZFT Age (Ma)"'
    if (has_muscar(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "MuscAr Age (Ma)"'
    if (has_bioar(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "BioAr Age (Ma)"'
    if (has_biotite(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "BioAr2 Age (Ma)"'
    if (has_mar(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "MAr Age (Ma)"'
    if (has_apb(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "APb Age (Ma)"'
    if (has_hornar(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "HornAr Age (Ma)"'
    if (has_rupb(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "RUPb Age (Ma)"'
    if (has_tupb(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "TUPb Age (Ma)"'
    if (has_zupb(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "ZUPb Age (Ma)"'
    if (has_tuth_he(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "TUTh Age (Ma)"'
    if (has_RDAAM_ApUThHe(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "RDAAM ApUThHe (Ma)"'
    if (has_RDAAM_ZirUThHe(age_flags)) write(file_unit_ages,'(A)',advance="no") ' "RDAAM ZirUThHe (Ma)"'

    write(file_unit_ages, *)
    write(file_unit_ages, *) 'ZONE T="Ages"'
    write(file_unit_ages, '(A2,i10)', advance="no") 'n=', num_of_points
    write(file_unit_ages, '(A4,i10)', advance="no") ', e=', num_of_elements
    write(file_unit_ages, *) ', et=quadrilateral, f=fepoint'
  end subroutine write_header

  subroutine write_content(num_of_flags, age_flags, all_ages)
    use m_age_algorithms

    implicit none

    integer(4), intent(in) :: num_of_flags, age_flags(num_of_flags)
    real(8), intent(in) :: all_ages(num_of_flags)

    integer(4) :: i

    do i = 1, 7
      if (age_flags(i) == 1) then
        write (file_unit_ages,'(F12.4)',advance="no") all_ages(i)
      endif
    enddo

    if (has_aft(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_aft)
    if (has_kfeldar(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_kfeldar)
    if (has_zhe(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_zhe)
    if (has_zhe_low_damage(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_zhe_low_damage)
    if (has_zhe_med_damage(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_zhe_med_damage)
    if (has_zhe_high_damage(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_zhe_high_damage)
    if (has_zft(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_zft)
    if (has_muscar(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_muscar)
    if (has_bioar(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_bioar)
    if (has_biotite(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_biotite)
    if (has_mar(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_mar)
    if (has_apb(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_apb)
    if (has_hornar(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_hornar)
    if (has_rupb(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_rupb)
    if (has_tupb(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_tupb)
    if (has_zupb(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_zupb)
    if (has_tuth_he(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_tuth_he)
    if (has_RDAAM_ApUThHe(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_RDAAM_ApUThHe)
    if (has_RDAAM_ZirUThHe(age_flags)) write(file_unit_ages,'(f12.4)',advance="no") all_ages(index_RDAAM_ZirUThHe)

    write(file_unit_ages, *)
  end subroutine write_content
end module m_ages_header
