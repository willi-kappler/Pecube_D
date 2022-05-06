module m_test_ages

    use m_age_algorithms
    use m_util

    implicit none

    contains
        subroutine test_ages()
          print *, "Begin test_ages"
          call system("date")

          call test_01()
          call test_02()

          print *, "Finished test_ages"
          call system("date")
        end subroutine test_ages

        subroutine test_01
          implicit none

          integer(4), parameter :: n_elems = 54

          real(8), dimension(n_elems) :: temperature, time
          real(8) :: age

          logical :: result

          temperature = (/ 116.9231, 115.8599, 114.7967, 113.7335, 112.6703, 111.6072, 110.5440, 109.4808, 108.4176, 107.3544, 106.2912, &
                           105.2280, 104.1648, 103.1017, 102.0385, 100.9753, 99.9121, 98.8489, 97.7857, 96.7225, 95.6594, 94.5962, 93.5330, &
                           92.4698, 91.4066, 90.3434, 89.2802, 88.2170, 87.1539, 86.0907, 85.0275, 83.9643, 82.9011, 81.8379, 80.7747, 79.7115, &
                           77.6442, 75.5769, 73.5096, 71.4423, 69.3750, 67.3077, 65.2404, 63.1731, 61.1058, 59.0385, 56.9712, 54.9039, 52.8365, &
                           50.7692, 48.7019, 46.6346, 44.5673, 42.5000 /)

          time = (/ 30.0000, 29.4286, 28.8571, 28.2857, 27.7143, 27.1429, 26.5714, 26.0000, 25.4286, 24.8571, 24.2857, 23.7143, 23.1429, 22.5714, &
                    22.0000, 21.4286, 20.8571, 20.2857, 19.7143, 19.1429, 18.5714, 18.0000, 17.4286, 16.8571, 16.2857, 15.7143, 15.1429, 14.5714, &
                    14.0000, 13.4286, 12.8571, 12.2857, 11.7143, 11.1429, 10.5714, 10.0000, 9.4444, 8.8889, 8.3333, 7.7778, 7.2222, 6.6667, 6.1111, &
                    5.5556, 5.0000, 4.4444, 3.8889, 3.3333, 2.7778, 2.2222, 1.6667, 1.1111, 0.5556, 0.0000 /)

          result = .true.

          call Muscovite(temperature, time, n_elems, age, 0, 0)
          call assert(age, 30.0_8, result)

          call ZFT(temperature, time, n_elems, age, 0, 0)
          call assert(age, 30.0_8, result)

          call Apatite_U_Th_Pb(temperature, time, n_elems, age, 0, 0)
          call assert(age, 30.0_8, result)

          call Rutile_U_Pb(temperature, time, n_elems, age, 0, 0)
          call assert(age, 30.0_8, result)

          call Titanite_U_Pb(temperature, time, n_elems, age, 0, 0)
          call assert(age, 30.0_8, result)

          call Zircon_U_Pb(temperature, time, n_elems, age, 0, 0)
          call assert(age, 30.0_8, result)

          call Biotite(temperature, time, n_elems, age, 0, 0)
          call assert(age, 30.0_8, result)

          if (result) then
            print *, "test_01 successfull!"
          else
            print *, "test_01 had errors"
          endif
        end subroutine test_01

        subroutine test_02
          implicit none

          integer(4), parameter :: n_elems = 5

          real(8), dimension(n_elems) :: temperature, time
          real(8) :: age

          logical :: result

          temperature = (/ 120.0, 119.0, 118.0, 117.0, 116.0 /)

          time = (/ 30.0000, 29.0, 28.0, 27.0, 26.0 /)

          result = .true.

          call Muscovite(temperature, time, n_elems, age, 0, 0)
          call assert(age, 4.0_8, result)

          call ZFT(temperature, time, n_elems, age, 0, 0)
          call assert(age, 4.0_8, result)

          call Apatite_U_Th_Pb(temperature, time, n_elems, age, 0, 0)
          call assert(age, 4.0_8, result)

          call Rutile_U_Pb(temperature, time, n_elems, age, 0, 0)
          call assert(age, 4.0_8, result)

          call Titanite_U_Pb(temperature, time, n_elems, age, 0, 0)
          call assert(age, 4.0_8, result)

          call Zircon_U_Pb(temperature, time, n_elems, age, 0, 0)
          call assert(age, 4.0_8, result)

          call Biotite(temperature, time, n_elems, age, 0, 0)
          call assert(age, 4.0_8, result)

          if (result) then
            print *, "test_02 successfull!"
          else
            print *, "test_02 had errors"
          endif
        end subroutine test_02

end module m_test_ages
