module m_calculate_ages
  use m_logger
  use m_mad_he
  use m_mad_he2
  use m_mad_trax_zirkon
  use m_muscovite
  use m_zft
  use m_mad_trax
  use m_age_algorithms
  use m_data_structures
  use m_read_time_temperature_history
  use m_pecube_config

contains

subroutine calculate_ages(config, num_recs, temperature_history_file, &
      outer_step, num_of_points, age_flags, age_info, header_info, model_depth, export_history)
  implicit none

  type(config_t), intent(in) :: config

  integer(4), intent(in) :: num_recs, outer_step, num_of_points
  integer(4), intent(in) :: age_flags(config%num_of_age_flags)

  character(*), intent(in) :: temperature_history_file

  real(8), intent(in) :: header_info(6), model_depth

  integer(4), parameter :: file_unit_temperature_history = 85
  integer(4) :: i, j, sum_ages

  logical, intent(in) :: export_history

  real(8), parameter :: AGE_ZERO_THRESHOLD = 0.01
  real(8), parameter :: AGE_COMPARE_THRESHOLD = 1.03
  real(8), dimension(:), allocatable :: all_ages_local, ztime_local
  real(8), dimension(:,:), allocatable :: ztemp_local
  real(4), dimension(:), allocatable :: ketch_time, ketch_temp
  real(8) :: ketch_alo, ketch_final_age, ketch_oldest_age, ketch_fmean, ketch_fdist(41)
  real(8) :: time_now

  type(age_info_t), intent(inout) :: age_info
  type(vector3D_t), dimension(:,:), allocatable :: vi_pos, vi_velo

  allocate(all_ages_local(config%num_of_age_flags))
  allocate(ketch_time(num_recs), ketch_temp(num_recs))

  all_ages_local = 0.0
  ketch_alo = 16.0
  ketch_final_age = 0.0
  ketch_oldest_age = 0.0
  ketch_fmean = 0.0
  ketch_fdist = 0.0
  ketch_time = 0.0
  ketch_temp = 0.0
  time_now = 0.0


  call read_time_temperature_history(num_recs, temperature_history_file, &
        outer_step, num_of_points, ztime_local, ztemp_local, vi_pos, vi_velo)

  if (export_history) then
    open(file_unit_temperature_history, file=trim(temperature_history_file)//".txt", status="unknown")

    write(file_unit_temperature_history, *) "# point id, sub time step, time [myr], temperature [deg C], px, py, pz, vx, vy, vz"

    do i = 1, num_of_points
      do j = 1, num_recs
          write(file_unit_temperature_history, *) i, j, ztime_local(j), ztemp_local(i, j), &
            vi_pos(i, j)%x, vi_pos(i, j)%y, vi_pos(i, j)%z - model_depth, &
            vi_velo(i, j)%x, vi_velo(i, j)%y, vi_velo(i, j)%z
      enddo
    enddo

    close(file_unit_temperature_history)

  endif ! export_to_text_file

  ! Only needed for exporting to text file, so no longer needed
  deallocate(vi_pos, vi_velo)

  do i = 1, num_of_points
      ! values for age flags:
      ! 1: AHe Age - Farley, 2000 (Ma) (a)
      ! 2: AHe (b)
      ! 3: AHe (c)
      ! 4: AHe (d)
      ! 5: AHe (e)
      ! 6: AHe (f)
      ! 7: AHe (g)
      ! 8: ZHe (h)
      ! 9: AFT (i)
      ! 10: ZFT (j)
      ! 11: MAr (k)
      ! 12: K-feldspar Ar (l)
      ! 13: Biotite Ar (m)
      ! 14: Muscovite Ar (n)
      ! 15: Hornblende Ar (o)
      ! 16: Apatite U-Th / Pb (APb) (p)
      ! 17: Biotite (q)
      ! 18: Rutile U-Pb (RUPb) (r)
      ! 19: Titanite U-Pb (TUPb) (s)
      ! 20: Zircon U-Pb (ZUPb) (t)
      ! 21: Titanite (U-Th) / He (u)
      ! 22: ZHe low damage (v)
      ! 23: ZHe med damage (w)
      ! 24: ZHe high damage (x)
      ! 25: RDAAM Apatite U-Th/He (y)
      ! 26: RDAAM Zircon U-Th/He (z)

      sum_ages = 0
      do j = 1, 8
          sum_ages = sum_ages + age_flags(j)
      enddo

      do j = 12, 15
          sum_ages = sum_ages + age_flags(j)
      enddo

      do j = 21, 24
          sum_ages = sum_ages + age_flags(j)
      enddo

      if (sum_ages > 0) then
         all_ages_local = age_info%all_ages(i,:)
         call Mad_He(ztime_local, ztemp_local(i,:), num_recs, all_ages_local, i, &
              header_info, age_flags, config%num_of_age_flags)
         age_info%all_ages(i,:) = all_ages_local
         ! 2013.02.12, WK: call code from Jean Braun here
         !call Mad_He2 (ztime_local, ztemp_local, num_recs, all_ages(i,1), 1)
         !call Mad_He2 (ztime_local, ztemp_local, num_recs, all_ages(i,8), 2)
      endif

      if (has_aft(age_flags)) then
          if (config%use_aft_ketcham) then
            ! Use Ketcham model
            time_now = ztime_local(num_recs)
            ketch_alo = 16.0
            ketch_final_age = 0.0
            ketch_oldest_age = 0.0
            ketch_fmean = 0.0
            ketch_fdist = 0.0
            ketch_time = sngl(ztime_local(num_recs:1:-1))
            ketch_temp = sngl(ztemp_local(i, num_recs:1:-1))

            call ketch_main(num_recs, ketch_time, ketch_temp, ketch_alo, ketch_final_age, ketch_oldest_age, ketch_fmean, ketch_fdist)

            age_info%all_ages(i, index_aft) = ketch_final_age
          else
            call Mad_Trax(ztime_local, ztemp_local(i,:), num_recs, 1, 2 ,&
                 age_info%all_ages(i, index_aft), age_info%ftldmean(i), &
                 age_info%ftldsd(i), age_info%ftld(i, :))
            ! 2013.06.13, WK: call aft code from Terra.
            ! Note that the terra code expects the time in reverse order:
            ! for ex. 1.0My, 2.0My, 3.0My, 4.0My, ...
            ! The order of the temperature must not be changed!
            !call get_aft_age(ztime_local(num_recs : 0 : -1), ztemp_local, num_recs, all_ages_local)
          endif
      endif

      if (has_zft(age_flags)) then
          !call ZFT_old (ztemp_local,ztime_local,num_recs,age_info%all_ages(i,10))
          !call ZFT(ztemp_local(i,:), ztime_local, num_recs, age_info%all_ages(i, index_zft), outer_step, i)
          call Mad_Zirc(ztime_local, ztemp_local(i,:), num_recs, 0, 1, age_info%all_ages(i, index_zft))

          !  call log_message("num_recs: " + num_recs)
          !  call log_message("age: " + all_ages(i, 10))
          !  call log_message("temperature: " + ztemp_local)
          !  call log_message("time: " + ztime_local)
          !
          !  stop
      endif

      if (has_mar(age_flags)) then
         ! MAr age
          call Muscovite(ztemp_local(i,:), ztime_local, num_recs, age_info%all_ages(i, index_mar), outer_step, i)
      endif

      if (has_apb(age_flags)) then
          ! Apatite U-Th / Pb
          call Apatite_U_Th_Pb(ztemp_local(i,:), ztime_local, num_recs, age_info%all_ages(i, index_apb), outer_step, i)
      end if

      if (has_biotite(age_flags)) then
          ! Biotite
          call Biotite(ztemp_local(i,:), ztime_local, num_recs, age_info%all_ages(i, index_biotite), outer_step, i)
      end if

      if (has_rupb(age_flags)) then
          ! Rutile
          call Rutile_U_Pb(ztemp_local(i,:), ztime_local, num_recs, age_info%all_ages(i, index_rupb), outer_step, i)
      end if

      if (has_tupb(age_flags)) then
          ! Titanite U-Pb
          call Titanite_U_Pb(ztemp_local(i,:), ztime_local, num_recs, age_info%all_ages(i, index_tupb), outer_step, i)
      end if

      if (has_zupb(age_flags)) then
          ! Zircon U-Pb
          call Zircon_U_Pb(ztemp_local(i,:), ztime_local, num_recs, age_info%all_ages(i, index_zupb), outer_step, i)
      end if

      if (has_zhe_low_damage(age_flags)) then
        ! TODO
      end if

      if (has_zhe_med_damage(age_flags)) then
        ! TODO
      end if

      if (has_zhe_high_damage(age_flags)) then
        ! TODO
      end if

      ! You can either have RDAAM or ZRDAAM but not both in the same simulation
      if (has_RDAAM_ApUThHe(age_flags)) then
        !call log_message("RDAAM1: " + age_info%all_ages(i, index_RDAAM_ApUThHe))
        !call rdaam_fortran(num_recs, ztime_local, ztemp_local(i, :), age_info%all_ages(i, index_RDAAM_ApUThHe))
        call rdaam_fortran(num_recs, ztime_local(num_recs:1:-1), ztemp_local(i, num_recs:1:-1), age_info%all_ages(i, index_RDAAM_ApUThHe))
        !call log_message("RDAAM2: " + age_info%all_ages(i, index_RDAAM_ApUThHe))
      else if (has_RDAAM_ZirUThHe(age_flags)) then
        !call log_message("RDAAM3: " + age_info%all_ages(i, index_RDAAM_ZirUThHe))
        call rdaam_fortran(num_recs, ztime_local(num_recs:1:-1), ztemp_local(i, num_recs:1:-1), age_info%all_ages(i, index_RDAAM_ZirUThHe))
        !call log_message("RDAAM4: " + age_info%all_ages(i, index_RDAAM_ZirUThHe))
      end if

      ! if (i == 1) then
      !    call log_message("m: " + m + ", i: " + i + ", all_ages(i): " + all_ages(i, :))
      ! end if

      ! Consistency checks for age values:

      ! Check for zero age:
      ! do j = 1, config%num_of_age_flags
      !   if ((age_flags(j) == 1) .and. (age_info%all_ages(i, j) < AGE_ZERO_THRESHOLD)) then
      !     call log_message("----------")
      !     call log_message("age is zero: " + j + " (type of age) " + age_name_from_index(j))
      !     call log_message("age: " + age_info%all_ages(i, j))
      !     call log_message("surface point id: " + i)
      !     call log_message("time step: " + outer_step)
      !     call log_message("----------")
      !   endif
      ! enddo ! j = 1, num_of_age_flags

      ! Check if age values overlap:
      ! Age(MAr) < Age(ZFT) < Age(ZHE)
      ! if (has_zft(age_flags) .and. has_mar(age_flags)) then
      !   if ((get_zft(age_info, i) * AGE_COMPARE_THRESHOLD) < get_mar(age_info, i)) then
      !     call log_message("age(zft) < age(mar): " + get_zft(age_info, i) + ", " + get_mar(age_info, i) + &
      !       ", i: " + i)
      !   endif
      ! endif

      ! if (has_zft(age_flags) .and. has_zhe(age_flags)) then
      !   if (get_zft(age_info, i) > (get_zhe(age_info, i) * AGE_COMPARE_THRESHOLD)) then
      !     call log_message("age(zft) > age(zhe): " + get_zft(age_info, i) + ", " + get_zhe(age_info, i) + &
      !       ", i: " + i)
      !   endif
      ! endif
      !stop
  enddo ! i = 1, num_of_points

  deallocate(all_ages_local)
  deallocate(ztemp_local, ztime_local)
  deallocate(ketch_time, ketch_temp)
end subroutine calculate_ages

function age_name_from_index(i)
  implicit none

  character(30) :: age_name_from_index
  integer, intent(in) :: i

  select case(i)
  case(index_ahe_1:index_ahe_7)
      age_name_from_index = "AHe"
    case(index_zhe)
      age_name_from_index = "ZHe"
    case(index_aft)
      age_name_from_index = "AFT"
    case(index_zft)
      age_name_from_index = "ZFT"
    case(index_mar)
      age_name_from_index = "MAr"
    case(index_kfeldar)
      age_name_from_index = "K-feldspar Ar"
    case(index_bioar)
      age_name_from_index = "Biotite Ar"
    case(index_muscar)
      age_name_from_index = "Muscovite Ar"
    case(index_hornar)
      age_name_from_index = "Hornblende Ar"
    case(index_apb)
      age_name_from_index = "Apatite U-Th / Pb (APb)"
    case(index_biotite)
      age_name_from_index = "Biotite"
    case(index_rupb)
      age_name_from_index = "Rutile U-Pb (RUPb)"
    case(index_tupb)
      age_name_from_index = "Titanite U-Pb (TUPb)"
    case(index_zupb)
      age_name_from_index = "Zircon U-Pb (ZUPb)"
    case(index_tuth_he)
      age_name_from_index = "Titanite (U-Th)"
    case(22)
      age_name_from_index = "ZHe low damage"
    case(23)
      age_name_from_index = "ZHe med damage"
    case(24)
      age_name_from_index = "ZHe high damage"
    case default
      age_name_from_index = "unknown"
  end select
end function age_name_from_index
end module m_calculate_ages
