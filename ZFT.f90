module m_zft
contains

subroutine ZFT_old (temp,time,ntime,age)

      implicit none

      integer(4) :: ntime, i, kmod

      real(8) :: temp(ntime),time(ntime), age, a, B, closure, cooling, diff
      real(8) :: time_local(ntime), energy, energy_cal, geom, r, ratio, tau
      real(8) :: tempp, closurep, D0

      ! Kinetic parameters to use:
      ! 1=Batt et al., 2001
      ! 2=
      ! 3=
      kmod=1

      if (kmod.eq.1) then
        energy_cal=49.77    ! Activation energy [kcal/mol]
        B=3.16e-22      ! 1/(D0/a^2) [My]
        ! D0 = 100
        ! a =  9.9865
        !
        ! D0 = 10
        ! a =  3.158
      elseif (kmod.eq.2) then
        energy=0.       ! Activation energy [kJ/mol]
        D0=0.           ! D0 [cm^2/s]
        a=0.            ! a [um]
      elseif (kmod.eq.3) then
        energy=0.       ! Activation energy [kJ/mol]
        D0=0.           ! D0 [cm^2/s]
        a=0.            ! a [um]
      endif

      ! Geometry factor for plane sheet
      geom=55.

      ! Conversion factors
      if (kmod.eq.1) then
        energy=energy_cal*4.184*1.e3
        diff=1./B
      else
        energy=energy*1.e3
        diff=((D0*0.01*0.01)/(a*1e-6)**2)*365.25*24.*3600.*1.e6
      endif

      cooling = 0.0
      closure = 0.0
      closurep = 0.0
      tempp = 0.0

      r=8.314

      ! 2013.02.01, WK: adjust time to model time
      time_local = time - time(ntime)
      age=real(time_local(1))
        do i=ntime,2,-1
          if (i.eq.1) then
              cooling=real((temp(i+1)-temp(i))/(time_local(i+1)-time_local(i)))
          elseif (i.eq.ntime) then
              cooling=real((temp(i)-temp(i-1))/(time_local(i)-time_local(i-1)))
          else
              cooling=real((temp(i+1)-temp(i-1))/(time_local(i+1)-time_local(i-1)))
          endif
! note that cooling cannot be nil and is therefore forced to
! be at least 1deg/10My
        cooling=max(cooling,1./10.)
        tau=r*(temp(i)+273.)**2/energy/cooling
        closure=energy/r/log(geom*tau*diff)-273.

          if (temp(i).gt.closure) then
              ratio=real((closurep-tempp)/(closurep-tempp+temp(i)-closure))
              age=real(time_local(i)+(time_local(i-1)-time_local(i))*ratio)
              return
          endif
        closurep=closure
        tempp=temp(i)
        enddo

      return
      end subroutine

end module m_zft
