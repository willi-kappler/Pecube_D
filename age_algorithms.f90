module m_age_algorithms
  use m_logger
  use m_data_structures

  implicit none
  save

! Written by Willi Kappler
! willi.kappler@uni-tuebingen.de

  !> Make everything private
  private

  !> Universal constants
  real(8), parameter :: PI = 3.1415926535897932384626433
  real(8), parameter :: univ_gas_const = 8.314462175 !> [J / (Mol * Kelvin)]
  real(8), parameter :: days_per_year = 365.25
  real(8), parameter :: hours_per_day = 24.0
  real(8), parameter :: sec_per_hour = 3600.0
  real(8), parameter :: mil_year = 1.0e6
  real(8), parameter :: sec_per_mil_year = sec_per_hour * hours_per_day * days_per_year * mil_year
  real(8), parameter :: kelvin_to_deg_cel = -273.15
  real(8), parameter :: deg_cel_to_kelvin = +273.15
  real(8), parameter :: um_to_m = 1.0e-6
  real(8), parameter :: kJ_to_J = 1.0e3

  !> Constants for age indices
  integer(4), parameter :: index_ahe_1 = 1
  integer(4), parameter :: index_ahe_2 = 2
  integer(4), parameter :: index_ahe_3 = 3
  integer(4), parameter :: index_ahe_4 = 4
  integer(4), parameter :: index_ahe_5 = 5
  integer(4), parameter :: index_ahe_6 = 6
  integer(4), parameter :: index_ahe_7 = 7
  integer(4), parameter :: index_zhe = 8
  integer(4), parameter :: index_aft = 9
  integer(4), parameter :: index_zft = 10
  integer(4), parameter :: index_mar = 11
  integer(4), parameter :: index_kfeldar = 12
  integer(4), parameter :: index_bioar = 13
  integer(4), parameter :: index_muscar = 14
  integer(4), parameter :: index_hornar = 15
  integer(4), parameter :: index_apb = 16
  integer(4), parameter :: index_biotite = 17
  integer(4), parameter :: index_rupb = 18
  integer(4), parameter :: index_tupb = 19
  integer(4), parameter :: index_zupb = 20
  integer(4), parameter :: index_tuth_he = 21
  integer(4), parameter :: index_zhe_low_damage = 22
  integer(4), parameter :: index_zhe_med_damage = 23
  integer(4), parameter :: index_zhe_high_damage = 24
  integer(4), parameter :: index_RDAAM_ApUThHe = 25
  integer(4), parameter :: index_RDAAM_ZirUThHe = 26

  character(30) :: method

  !> Exported subroutines
  !public age_with_closure
  public Muscovite
  public ZFT
  public Apatite_U_Th_Pb
  public Biotite
  public Rutile_U_Pb
  public Titanite_U_Pb
  public Zircon_U_Pb

  public index_ahe_1
  public index_ahe_2
  public index_ahe_3
  public index_ahe_4
  public index_ahe_5
  public index_ahe_6
  public index_ahe_7
  public index_zhe
  public index_aft
  public index_zft
  public index_mar
  public index_kfeldar
  public index_bioar
  public index_muscar
  public index_hornar
  public index_apb
  public index_biotite
  public index_rupb
  public index_tupb
  public index_zupb
  public index_tuth_he
  public index_zhe_low_damage
  public index_zhe_med_damage
  public index_zhe_high_damage
  public index_RDAAM_ApUThHe
  public index_RDAAM_ZirUThHe

  public has_ahe_1
  public has_ahe_2
  public has_ahe_3
  public has_ahe_4
  public has_ahe_5
  public has_ahe_6
  public has_ahe_7
  public has_zhe
  public has_aft
  public has_zft
  public has_mar
  public has_kfeldar
  public has_bioar
  public has_muscar
  public has_hornar
  public has_apb
  public has_biotite
  public has_rupb
  public has_tupb
  public has_zupb
  public has_tuth_he
  public has_zhe_low_damage
  public has_zhe_med_damage
  public has_zhe_high_damage
  public has_RDAAM_ApUThHe
  public has_RDAAM_ZirUThHe

  public get_zhe
  public get_zft
  public get_mar

  public sec_per_mil_year

  contains

    !> Calculates the muscovite argon age given the time / temperature history path
    !! @param temperature Array of temperature values temperature path [degC]
    !! @param time Array of time values / time path [My]
    !! @param n_elems Number of elements for the above arrays
    !! @param age Output parameter for the calculated age value
    subroutine Muscovite(temperature, time, n_elems, age, m, i)
      implicit none

      !> Subroutine parameters
      integer(4), intent(in) :: n_elems, m, i

      real(8), intent(in), dimension(n_elems) :: temperature, time

      real(8), intent(out) :: age

      !> local constants
      real(8), parameter :: D0 = 4.0e-8_8
      real(8), parameter :: energy = 180.0_8
      real(8), parameter :: grain_size = 750.0_8
      real(8), parameter :: geometry_factor = 27.0_8 ! cylinder geometry, changes from Byron Adams

      method = "Muscovite Ar"
      call age_with_closure(dble(temperature), dble(time), n_elems, D0, energy, grain_size, geometry_factor, age, m, i)
    end subroutine Muscovite

    !> Calculates the zircon fission track age given the time / temperature history path
    !! @param temperature Array of temperature values temperature path [degC]
    !! @param time Array of time values / time path [My]
    !! @param n_elems Number of elements for the above arrays
    !! @param age Output parameter for the calculated age value
    subroutine ZFT(temperature, time, n_elems, age, m, i)
      implicit none

      !> Subroutine parameters
      integer(4), intent(in) :: n_elems, m, i

      real(8), intent(in), dimension(n_elems) :: temperature, time

      real(8), intent(out) :: age

      !> local constants
      real(8), parameter :: D0 = 0.001_8
      real(8), parameter :: energy = 208.32768_8
      real(8), parameter :: grain_size = 3.158_8
      real(8), parameter :: geometry_factor = 55.0_8

      method = "ZFT"
      call age_with_closure(dble(temperature), dble(time), n_elems, D0, energy, grain_size, geometry_factor, age, m, i)
    end subroutine ZFT

    !> Calculates the apatite U-Th/Pb age given the time / temperature history path
    !! @param temperature Array of temperature values temperature path [degC]
    !! @param time Array of time values / time path [My]
    !! @param n_elems Number of elements for the above arrays
    !! @param age Output parameter for the calculated age value
    subroutine Apatite_U_Th_Pb(temperature, time, n_elems, age, m, i)
      implicit none

      !> Subroutine parameters
      integer(4), intent(in) :: n_elems, m, i

      real(8), intent(in), dimension(n_elems) :: temperature, time

      real(8), intent(out) :: age

      !> local constants
      real(8), parameter :: D0 = 2.0e-8_8
      real(8), parameter :: energy = 230.0_8
      real(8), parameter :: grain_size = 50.0_8
      real(8), parameter :: geometry_factor = 55.0_8

      method = "Apatite U-Th / Pb (APb)"
      call age_with_closure(dble(temperature), dble(time), n_elems, D0, energy, grain_size, geometry_factor, age, m, i)
    end subroutine Apatite_U_Th_Pb

    !> Calculates the Rutile U-Pb (RUPb) age given the time / temperature history path
    !! @param temperature Array of temperature values temperature path [degC]
    !! @param time Array of time values / time path [My]
    !! @param n_elems Number of elements for the above arrays
    !! @param age Output parameter for the calculated age value
    subroutine Rutile_U_Pb(temperature, time, n_elems, age, m, i)
      implicit none

      !> Subroutine parameters
      integer(4), intent(in) :: n_elems, m, i

      real(8), intent(in), dimension(n_elems) :: temperature, time

      real(8), intent(out) :: age

      !> local constants
      real(8), parameter :: D0 = 1.6e-10_8
      real(8), parameter :: energy = 243.0_8
      real(8), parameter :: grain_size = 250.0_8
      real(8), parameter :: geometry_factor = 55.0_8 ! spherical

      method = "Rutile U-Pb (RUPb)"
      call age_with_closure(dble(temperature), dble(time), n_elems, D0, energy, grain_size, geometry_factor, age, m, i)
    end subroutine Rutile_U_Pb

    !> Calculates the Titanite U-Pb (TUPb) age given the time / temperature history path
    !! @param temperature Array of temperature values temperature path [degC]
    !! @param time Array of time values / time path [My]
    !! @param n_elems Number of elements for the above arrays
    !! @param age Output parameter for the calculated age value
    subroutine Titanite_U_Pb(temperature, time, n_elems, age, m, i)
      implicit none

      !> Subroutine parameters
      integer(4), intent(in) :: n_elems, m, i

      real(8), intent(in), dimension(n_elems) :: temperature, time

      real(8), intent(out) :: age

      !> local constants
      real(8), parameter :: D0 = 1.1e-4_8
      real(8), parameter :: energy = 331.0_8
      real(8), parameter :: grain_size = 500.0_8
      real(8), parameter :: geometry_factor = 55.0_8 ! spherical

      method = "Titanite U-Pb (TUPb)"
      call age_with_closure(dble(temperature), dble(time), n_elems, D0, energy, grain_size, geometry_factor, age, m, i)
    end subroutine Titanite_U_Pb

    !> Calculates the Zircon U-Pb (TUPb) age given the time / temperature history path
    !! @param temperature Array of temperature values temperature path [degC]
    !! @param time Array of time values / time path [My]
    !! @param n_elems Number of elements for the above arrays
    !! @param age Output parameter for the calculated age value
    subroutine Zircon_U_Pb(temperature, time, n_elems, age, m, i)
      implicit none

      !> Subroutine parameters
      integer(4), intent(in) :: n_elems, m, i

      real(8), intent(in), dimension(n_elems) :: temperature, time

      real(8), intent(out) :: age

      !> local constants
      real(8), parameter :: D0 = 7.8e-3_8
      real(8), parameter :: energy = 544.0_8
      real(8), parameter :: grain_size = 50.0_8
      real(8), parameter :: geometry_factor = 55.0_8 ! spherical

      method = "Zircon U-Pb (ZUPb)"
      call age_with_closure(dble(temperature), dble(time), n_elems, D0, energy, grain_size, geometry_factor, age, m, i)
    end subroutine Zircon_U_Pb

    !> Calculates the biotite age given the time / temperature history path
    !! @param temperature Array of temperature values temperature path [degC]
    !! @param time Array of time values / time path [My]
    !! @param n_elems Number of elements for the above arrays
    !! @param age Output parameter for the calculated age value
    subroutine Biotite(temperature, time, n_elems, age, m, i)
      implicit none

      !> Subroutine parameters
      integer(4), intent(in) :: n_elems, m, i

      real(8), intent(in), dimension(n_elems) :: temperature, time

      real(8), intent(out) :: age

      !> local constants
      real(8), parameter :: D0 = 2.0e-13_8
      real(8), parameter :: energy = 105.0_8
      real(8), parameter :: grain_size = 500.0_8
      real(8), parameter :: geometry_factor = 27.0_8 ! cylinder geometry, changes from Byron Adams

      method = "Biotite"
      call age_with_closure(dble(temperature), dble(time), n_elems, D0, energy, grain_size, geometry_factor, age, m, i)
    end subroutine Biotite

    !> Calculates the age using the closure temperature
    !! @param temperature Array of temperature values temperature path [degC]
    !! @param time Array of time values / time path [My]
    !! @param n_elems Number of elements for the above arrays
    !! @param D0 Diffusivity at infinite temperature [m² / s]
    !! @param energy Activation energy [kJ / mol]
    !! @param grain_size Grain size [µm]
    !! @param geometry_factor Geometry factor for plane sheets
    !! @param age Output parameter for the calculated age value
    subroutine age_with_closure(temperature, time, n_elems, D0, energy, grain_size, geometry_factor, age, m, i)
      implicit none

      integer(4), parameter :: file_unit_closure = 81

      !> Subroutine parameters
      integer(4), intent(in) :: n_elems, m, i

      real(8), intent(in) :: D0, energy, grain_size, geometry_factor

      real(8), intent(in), dimension(n_elems) :: temperature, time

      real(8), intent(out) :: age

      !> local variables
      real(8), dimension(n_elems) :: local_time
      real(8) :: cooling_rate !> [Kelvin / s]
      real(8) :: closure_temp !> [degC]
      real(8) :: prev_closure_temp !> [degC]
      real(8) :: prev_temp !> [degC]
      real(8) :: tau
      real(8) :: diff
      real(8) :: ratio
      real(8) :: energy_in_J !> Activation energy in [J]

      integer(4) :: j

      !> What is diff ?
      !> diff [1/s]
      diff = sec_per_mil_year * (D0 / ((grain_size * um_to_m)**2.0))

      !> Adjust time to local time
      local_time = time - time(n_elems)

      !> Initial age value
      age = local_time(1)

      !> Convert activation energy from kJ to J
      !> [J / Mol]
      energy_in_J = energy * kJ_to_J

      !> Initial value for variables
      cooling_rate = 0.0 !> [Kelvin / s]
      closure_temp = 0.0 !> [degC]
      prev_closure_temp = 0.0 !> [degC]
      prev_temp = 0.0 !> [degC]
      tau = 0.0 !> [s]
      ratio = 0.0

      !> Iterative age calculation
      !> Taken from Muscovite.f
      do j = n_elems, 2, -1
          !> For the temperature difference it doesn't matter if it's in degC or Kelvin:
          !> temp1 - temp2 == (temp1 + deg_cel_to_kelvin) - (temp2 + deg_cel_to_kelvin)
          if (j == 1) then
              cooling_rate = (temperature(j + 1) - temperature(j)) / (local_time(j + 1) - local_time(j))
          else if (j == n_elems) then
              cooling_rate = (temperature(j) - temperature(j - 1)) / (local_time(j) - local_time(j - 1))
          else
              cooling_rate = (temperature(j + 1) - temperature(j - 1)) / (local_time(j + 1) - local_time(j - 1))
          end if

          !> Cooling rate must be a valid value
          cooling_rate = max(cooling_rate, 0.1_8)

          !> What is tau ?
          !> tau [s]
          tau = univ_gas_const * ((temperature(j) + deg_cel_to_kelvin)**2.0) / (energy_in_J * cooling_rate)
          !> closure_temp [degC]
          closure_temp = kelvin_to_deg_cel + (energy_in_J / (log(geometry_factor * tau * diff) * univ_gas_const))

          if (temperature(j) > closure_temp) then
              ratio = (prev_closure_temp - prev_temp) / (prev_closure_temp - prev_temp + temperature(j) - closure_temp)
              age = local_time(j) + (ratio * (local_time(j-1) - local_time(j)))

              return
          end if

          prev_closure_temp = closure_temp
          prev_temp = temperature(j)

          write (file_unit_closure, *) m, i, method, prev_temp, prev_closure_temp

      end do
    end subroutine age_with_closure

    logical function has_ahe_1(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_ahe_1 = flags(index_ahe_1) == 1
    end function has_ahe_1

    logical function has_ahe_2(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_ahe_2 = flags(index_ahe_2) == 1
    end function has_ahe_2

    logical function has_ahe_3(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_ahe_3 = flags(index_ahe_3) == 1
    end function has_ahe_3

    logical function has_ahe_4(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_ahe_4 = flags(index_ahe_4) == 1
    end function has_ahe_4

    logical function has_ahe_5(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_ahe_5 = flags(index_ahe_5) == 1
    end function has_ahe_5

    logical function has_ahe_6(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_ahe_6 = flags(index_ahe_6) == 1
    end function has_ahe_6

    logical function has_ahe_7(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_ahe_7 = flags(index_ahe_7) == 1
    end function has_ahe_7

    logical function has_zhe(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_zhe = flags(index_zhe) == 1
    end function has_zhe

    logical function has_aft(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_aft = flags(index_aft) == 1
    end function has_aft

    logical function has_zft(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_zft = flags(index_zft) == 1
    end function has_zft

    logical function has_mar(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_mar = flags(index_mar) == 1
    end function has_mar

    logical function has_kfeldar(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_kfeldar = flags(index_kfeldar) == 1
    end function has_kfeldar

    logical function has_bioar(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_bioar = flags(index_bioar) == 1
    end function has_bioar

    logical function has_muscar(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_muscar = flags(index_muscar) == 1
    end function has_muscar

    logical function has_hornar(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_hornar = flags(index_hornar) == 1
    end function has_hornar

    logical function has_apb(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_apb = flags(index_apb) == 1
    end function has_apb

    logical function has_biotite(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_biotite = flags(index_biotite) == 1
    end function has_biotite

    logical function has_rupb(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_rupb = flags(index_rupb) == 1
    end function has_rupb

    logical function has_tupb(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_tupb = flags(index_tupb) == 1
    end function has_tupb

    logical function has_zupb(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_zupb = flags(index_zupb) == 1
    end function has_zupb

    logical function has_tuth_he(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_tuth_he = flags(index_tuth_he) == 1
    end function has_tuth_he

    logical function has_RDAAM_ApUThHe(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      ! Apatite U-Th/He
      has_RDAAM_ApUThHe = flags(index_RDAAM_ApUThHe) == 1
    end function has_RDAAM_ApUThHe

    logical function has_RDAAM_ZirUThHe(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      ! Zircon U-Th/He
      has_RDAAM_ZirUThHe = flags(index_RDAAM_ZirUThHe) == 1
    end function has_RDAAM_ZirUThHe

    logical function has_zhe_low_damage(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_zhe_low_damage = flags(index_zhe_low_damage) == 1
    end function has_zhe_low_damage

    logical function has_zhe_med_damage(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_zhe_med_damage = flags(index_zhe_med_damage) == 1
    end function has_zhe_med_damage

    logical function has_zhe_high_damage(flags)
      implicit none

      integer(4), dimension(:), intent(in) :: flags

      has_zhe_high_damage = flags(index_zhe_high_damage) == 1
    end function has_zhe_high_damage

    real(8) function get_zhe(age_info, i)
      implicit none

      type(age_info_t), intent(in) :: age_info
      integer(4), intent(in) :: i

      get_zhe = age_info%all_ages(i, index_zhe)
    end function get_zhe

    real(8) function get_zft(age_info, i)
      implicit none

      type(age_info_t), intent(in) :: age_info
      integer(4), intent(in) :: i

      get_zft = age_info%all_ages(i, index_zft)
    end function get_zft

    real(8) function get_mar(age_info, i)
      implicit none

      type(age_info_t), intent(in) :: age_info
      integer(4), intent(in) :: i

      get_mar = age_info%all_ages(i, index_mar)
    end function get_mar

end module m_age_algorithms
