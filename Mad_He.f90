module m_mad_he

contains

! Subroutine to calculate Helium ages of grains with variable
! size and radtion levels
! Calulcates 7 Helium ages (1 to 7) along with a Zircon Helium (8) age
!
!
! WK, 2013.11.13: Added code from Jean Braun for more age parameters:
!
! age_flag(1): He in Apatite
! age_flag(2 to 7): He in Apatite various sizes
! age_flag(8): He in Zircon
! age_flag(12): Ar in K-feldspar
! age_flag(13): Ar in Biotite
! age_flag(14): Ar in Muscovite
! age_flag(15): Ar in Hornblende
! age_flag(21): Titanite (U-Th) / He
! age_flag(22): ZHe low damage
! age_flag(23): ZHe med damage
! age_flag(24): ZHe high damage

! WK: output in array output_age, work array agei

! WK: TODO: rewrite this and move it into age_algorithms

      subroutine Mad_He (time,temperature,ntime,output_age,isurf,header_info,global_age_flags,num_of_age_flags)
        use m_tridag

      implicit none

      ! Declare space for passed/new arrays
      real(8), intent(in) :: time(ntime),temperature(ntime)
      real(8),dimension(:,:),allocatable::age,diag,sup,inf,f
      real(8),dimension(:),allocatable::agei,f1,f2,Ft
      real(8), dimension(:), allocatable :: inf_temp, diag_temp, sup_temp
      real(8), dimension(:), allocatable :: f_temp, age_temp
      integer(4), intent(in) :: num_of_age_flags, ntime
      real(8),dimension(1:num_of_age_flags):: output_age
      integer(4), intent(in) :: isurf, global_age_flags(num_of_age_flags)
      integer(4) :: age_flags(num_of_age_flags)
      real(8) :: header_info(6), a, al, alpha, am, ard, as, beta, D0, D0a2
      real(8) :: D0a2rd, D0al2, D0am2, D0as2, D0z, D0za2, Da2now, Da2now_h
      real(8) :: Da2now_l, Da2now_m, Da2then, Da2then_h, Da2then_l, Da2then_m
      real(8) :: Dal2now, Dal2then, Dam2now, Dam2then, Das2now, Das2then
      real(8) :: dr, dt, dt0, Dza2now, Dza2then, Ea, Eard, Eaz, Et, eUh
      real(8) :: eUl, eUm, fact, fstep, He_h, He_l, He_m, P_h, P_l, P_m, pi
      real(8) :: psi, R, S, temp, tempf, temps, Thh, Thhc, Thl, Thlc, Thm
      real(8) :: Thmc, Uh, Uhc, Ul, Ulc, Um, Umc

      integer(4) :: i, istep, itime, k, kmod, n, nstep

      real(8), dimension(:), allocatable :: Da2now_ar, Da2then_ar, &
            D0a2_ar, Ea_ar

      real(8), parameter :: SEC_PER_MIL_YEARS = 3600.0 * 24.0 * 365.25e6

      age_flags = 0

      if (global_age_flags(1) == 1) age_flags(1) = 1 ! He Apatite
      if (global_age_flags(2) == 1) age_flags(2) = 1 ! He Apatite 2
      if (global_age_flags(3) == 1) age_flags(3) = 1 ! He Apatite 3
      if (global_age_flags(4) == 1) age_flags(4) = 1 ! He Apatite 4
      if (global_age_flags(5) == 1) age_flags(5) = 1 ! He Apatite 5
      if (global_age_flags(6) == 1) age_flags(6) = 1 ! He Apatite 6
      if (global_age_flags(7) == 1) age_flags(7) = 1 ! He Apatite 7
      if (global_age_flags(8) == 1) age_flags(8) = 1 ! He Zircon
      if (global_age_flags(12) == 1) age_flags(12) = 1 ! Ar K-feldspar
      if (global_age_flags(13) == 1) age_flags(13) = 1 ! Ar Biotite
      if (global_age_flags(14) == 1) age_flags(14) = 1 ! Ar Muscovite
      if (global_age_flags(15) == 1) age_flags(15) = 1 ! Ar Hornblende
      if (global_age_flags(21) == 1) age_flags(21) = 1 ! Titanite (U-Th) / He

      ! He Zircon
      if (global_age_flags(22) == 1) age_flags(22) = 1 ! He Zirkon LOW Radiation Damage: D0z_low = 4.6 [cm^2 / sec]
      if (global_age_flags(23) == 1) age_flags(23) = 1 ! He Zirkon MEDIUM Radiation Damage: D0z_med = 0.3 [cm^2 / sec]
      if (global_age_flags(24) == 1) age_flags(24) = 1 ! He Zirkon HIGH Radiation Damage: D0z_high = 4.6E+05 [cm^2 / sec]




      ! Define constants
      pi=3.141592       ! pi
      R=8.314           ! Universal gas constant [J/(K*mol)]
      S=20.             ! Alpha stopping distance [um]

      ! Allocates alpha ejection correction array
      allocate (Ft(num_of_age_flags))

      Ft = 0.0

      ! Default diffusion kinetic model
      ! 1 = Wolf et al., 1996
      ! 2 = Farley, 2000
      kmod=2

      if (kmod.eq.1) then
        ! Kinetics from Wolf et al., 1996
        D0a2=10.**7.7*SEC_PER_MIL_YEARS              ! [1/My]
        Ea=36.2e3*4.184                               ! Activation energy in undamaged apatite grain [J/mol]
      elseif (kmod.eq.2) then
        ! Kinetic parameters from Farley, 2000 (added by dwhipp 03/07)
        ! 2014.06.10, WK: Changed according to table from Byron Adams
        D0=90.0                                       ! Diffusivity at infinite temp. [cm^2/s]
        a=50.0                                        ! Diffusion domain [um]
        Ea=142.0e3                                    ! Activation energy in undamaged apatite grain [J/mol]
        Ft(1)=1.-(3.*S)/(4.*a)+(S**3.)/(16.*a**3.)    ! Alpha ejection correction factor []
        ! Convert diffusion kinetics to same units as Wolf et al., 1996 (above)
        D0=D0*0.01*0.01                               ! [m^2/s]
        a=a*1e-6                                      ! [m]
        D0a2=(D0/(a*a))*SEC_PER_MIL_YEARS            ! D0 over a^2 [1/My]
      endif

      ! Zircon He diffusion kinetics (Reiners et al., 2004)
      D0z = 0.46                  ! Diffusivity of zircon grain [cm^2/s]
      Eaz = 169.*1000.            ! Activation energy in undamaged zircon grain [J/mol]

      ! Radiation damaged apatite diffusion kinetics (Schuster et al., 2006)
      Eard=120.e3               ! Activation energy of undamaged apatite [J/mol]
      Et=29.e3                  ! Activation energy in damaged zones of crystal [J/mol]
      psi=1.26e-4               ! Radiation damage constant [gm/nmol]
      D0a2rd=1.58e4*SEC_PER_MIL_YEARS   ! D0/a^2 [1/My]
      ard=60.                   ! Radiation damaged grain size *****ONLY FOR FT CORRECTION***** [um]
      !D0a2rd=0.        ! D0/a^2 [1/My]

      ! User-input parameters (this and conversions below added by dwhipp 03/07)
      as=20.            ! Diffusion domain (small grain) [um]
      am=40.            ! Diffusion domain (medium grain)[um]
      al=70.            ! Diffusion domain (large grain) [um]
      Ul=10.            ! Uranium concentration (low damage) [ppm]
      Um=50.            ! Uranium concentration (moderate damage) [ppm]
      Uh=200.           ! Uranium concentration (high damage) [ppm]
      Thl=10.           ! Thorium concentration (low damage) [ppm]
      Thm=50.           ! Thorium concentration (moderate damage) [ppm]
      Thh=200.          ! Thorium concentration (high damage) [ppm]

      ! Calculate effective Uranium for Tecplot output header (Schuster et al., 2006)
      ! Old version below when all thermal histories were used
      ! Commented out by dwhipp (10/07)
      !if(step.eq.1.and.isurf.eq.1) then
      if(isurf.eq.1) then
        eUl=Ul+0.235*Thl
        eUm=Um+0.235*Thm
        eUh=Uh+0.235*Thh
        ! Write input parameters to Tecplot header array
        header_info(1)=as
        header_info(2)=am
        header_info(3)=al
        header_info(4)=eUl
        header_info(5)=eUm
        header_info(6)=eUh
      endif

      ! Other parameter conversions
      D0z = D0z*0.01*0.01                                 ! [m^2/s]

      Ft(2)=1.-(3.*S)/(4.*as)+(S**3.)/(16.*as**3.)      ! Alpha ejection correction factor, small grain []
      Ft(3)=1.-(3.*S)/(4.*am)+(S**3.)/(16.*am**3.)      ! Alpha ejection correction factor, medium grain []
      Ft(4)=1.-(3.*S)/(4.*al)+(S**3.)/(16.*al**3.)      ! Alpha ejection correction factor, large grain []
      Ft(5)=1.-(3.*S)/(4.*ard)+(S**3.)/(16.*ard**3.)    ! Alpha ejection correction factor, low r.d. []
      Ft(6)=Ft(5)                                       ! Alpha ejection correction factor, medium r.d. []
      Ft(7)=Ft(5)                                       ! Alpha ejection correction factor, high r.d. []
      Ft(8)=1.                                          ! No alpha ejection correction for ZHe
      as=as*1e-6                                ! [m]
      am=am*1e-6                                ! [m]
      al=al*1e-6                                ! [m]
      D0as2=(D0/(as*as))*SEC_PER_MIL_YEARS     ! D0 over a^2 [1/My]
      D0am2=(D0/(am*am))*SEC_PER_MIL_YEARS     ! D0 over a^2 [1/My]
      D0al2=(D0/(al*al))*SEC_PER_MIL_YEARS     ! D0 over a^2 [1/My]

      D0za2 = (D0z/(a*a))*SEC_PER_MIL_YEARS            ! D0 over a^2 [1/My]

      ! Convert U/Th concentrations [ppm]-->[nmol/g]
      Ulc=Ul*(1/238.03)*1e3
      Umc=Um*(1/238.03)*1e3
      Uhc=Uh*(1/238.03)*1e3
      Thlc=Thl*(1/232.038)*1e3
      Thmc=Thm*(1/232.038)*1e3
      Thhc=Thh*(1/232.038)*1e3

      ! He production rate (low, med, high damage) [nmol/(gm*My)]
      ! *** OLD VERSION - WRONG RATE BECAUSE OF GEOMETRIC FACTOR (4/3)*pi
      !P_l=(4./3.)*pi*(8.*4.95e-18*Ulc+6.*1.56e-18*Thlc)*SEC_PER_MIL_YEARS
      !P_m=(4./3.)*pi*(8.*4.95e-18*Umc+6.*1.56e-18*Thmc)*SEC_PER_MIL_YEARS
      !P_h=(4./3.)*pi*(8.*4.95e-18*Uhc+6.*1.56e-18*Thhc)*SEC_PER_MIL_YEARS

      ! He production rate (low, med, high damage) [nmol/(gm*My)]
      P_l=(8.*4.95e-18*Ulc+6.*1.56e-18*Thlc)*SEC_PER_MIL_YEARS
      P_m=(8.*4.95e-18*Umc+6.*1.56e-18*Thmc)*SEC_PER_MIL_YEARS
      P_h=(8.*4.95e-18*Uhc+6.*1.56e-18*Thhc)*SEC_PER_MIL_YEARS

      ! Define number of iterations and time step
      dt0=0.1
      n=100

      ! Allocate arrays
      allocate (agei(num_of_age_flags),age(num_of_age_flags,n),&
                diag(num_of_age_flags,n),&
                sup(num_of_age_flags,n),inf(num_of_age_flags,n),f(num_of_age_flags,n))
      allocate (f1(num_of_age_flags),f2(num_of_age_flags))

      allocate (Da2now_ar(num_of_age_flags))
      allocate (Da2then_ar(num_of_age_flags))
      allocate (D0a2_ar(num_of_age_flags))
      allocate (Ea_ar(num_of_age_flags))

      allocate(inf_temp(n))
      allocate(diag_temp(n))
      allocate(sup_temp(n))
      allocate(f_temp(n))
      allocate(age_temp(n))

      ! Initialize age and He concentrations (added by dwhipp 03/07)
      age=0.0            ! []
      agei=0.0           ! []
      He_l=0.0           ! [nmol/gm]
      He_m=0.0           ! [nmol/gm]
      He_h=0.0           ! [nmol/gm]
      diag = 0.0
      sup = 0.0
      inf = 0.0
      f = 0.0
      f1 = 0.0
      f2 = 0.0

      Da2now_ar = 0.0
      Da2then_ar = 0.0
      D0a2_ar = 0.0
      Ea_ar = 0.0

      ! KFeldAr: 12
      D0a2_ar(12) = 5.6 * SEC_PER_MIL_YEARS
      Ea_ar(12) = 120.0e3

      ! BioAr: 13
      !D0a2_ar(13) = 160.0 * SEC_PER_MIL_YEARS
      !Ea_ar(13) = 211.0e3
      ! WK: 2017.04.18, New values from Mary Schultz
      ! These numbers are for iron rich biotite, reference for this is Hodges 2014 Tr. on Geochem,
      ! where the numbers are based on Grove and Harrison American Mineralogist 1996.
      ! The closure T for this is 370 C @ 10 C/Myr cooling rate.  For crustal rocks, the iron rich
      ! biotite is most appropriate.  For Magnesium rich biotite the numbers are different, and these
      ! values should be used for mafic rocks.
      ! D0a2_ar(13) is the same as above: D0a2_ar(13) = 160.0 * SEC_PER_MIL_YEARS
      D0a2_ar(13) = (4.0e-5 / ((500.0 * 1.0e-6)**2.0)) * SEC_PER_MIL_YEARS
      Ea_ar(13) = 213.7e3

      ! MuscAr: 14
      D0a2_ar(14) = 13.2 * SEC_PER_MIL_YEARS
      Ea_ar(14) = 183.0e3

      ! HornAr: 15
      ! 2014.06.10, WK: changes from Byron Adams
      D0a2_ar(15) = 9.6 * SEC_PER_MIL_YEARS
      Ea_ar(15) = 268.0e3

      ! Titanite (U-Th) / He: 21
      ! D0 = 5.9e-3 [m²/s]
      ! E = 188 [kJ/mol]
      ! a = 250 [µm]
      D0a2_ar(21) = (5.9e-3 / ((250.0 * 1.0e-6)**2.0)) * SEC_PER_MIL_YEARS
      Ea_ar(21) = 188.0e3

      ! New zircon helium for radiation damage:
      ! low
      D0a2_ar(22) = ((4.6 * 0.01 * 0.01) / (a*a)) * SEC_PER_MIL_YEARS
      Ea_ar(22) = Eaz

      ! mid
      D0a2_ar(23) = ((0.3 * 0.01 * 0.01) / (a*a)) * SEC_PER_MIL_YEARS
      Ea_ar(23) = Eaz

      ! high
      D0a2_ar(24) = ((4.6e+05 * 0.01 * 0.01) / (a*a)) * SEC_PER_MIL_YEARS
      Ea_ar(24) = Eaz



      ! Loop through iteratively to calculate ages
      do itime=1,ntime-1

      ! Update time step
      dt=dt0
      nstep=max(1,int((time(itime)-time(itime+1)+tiny(dt))/dt))
      dt=(time(itime)-time(itime+1))/nstep
      alpha=0.5
      dr=1./(n-1)

      ! beta determines the geometry: 2 = spherical
      !                               1 = cylindrical
      beta=2.

      ! Grab current temperature info
      temps=temperature(itime)
      tempf=temperature(itime+1)

      ! Calculate new D0/a^2 (added by cspath 3/12/07)
      Da2now=D0a2*exp(-Ea/R/(temps+273.))                                               ! [1/My]
      Das2now=D0as2*exp(-Ea/R/(temps+273.))                                             ! [1/My]
      Dam2now=D0am2*exp(-Ea/R/(temps+273.))                                             ! [1/My]
      Dal2now=D0al2*exp(-Ea/R/(temps+273.))                                             ! [1/My]
      Da2now_l=(D0a2rd*exp(-Eard/R/(temps+273.)))/((psi*He_l*exp(Et/R/(temps+273.)))+1) ! [1/My]
      Da2now_m=(D0a2rd*exp(-Eard/R/(temps+273.)))/((psi*He_m*exp(Et/R/(temps+273.)))+1) ! [1/My]
      Da2now_h=(D0a2rd*exp(-Eard/R/(temps+273.)))/((psi*He_h*exp(Et/R/(temps+273.)))+1) ! [1/My]

      Dza2now = D0za2*exp(-Eaz/R/(temps+273.))                                          ! [1/My]

      do i=12,num_of_age_flags
          Da2now_ar(i) = D0a2_ar(i) * exp(-Ea_ar(i) / R / (temps + 273.0) )
      end do

        ! Loop through all timesteps in cooling history to calculate cooling ages
        do istep=1,nstep

        ! Determine normalized time and temperature gradient
        fstep=float(istep)/(nstep)
        temp=temps+(tempf-temps)*fstep

        ! Calculate various diffusion parameters (added by cspath 03/07)
        ! Calculate new diffusion parameters (standard grain)
        if (age_flags(1) == 1) then
          Da2then=Da2now
          Da2now=D0a2*exp(-Ea/R/(temp+273.))
          f1(1)=alpha*dt*Da2now/dr**2
          f2(1)=(1.-alpha)*dt*Da2then/dr**2
        endif

        ! Calculate new diffusion parameters (small grain)
        if (age_flags(2) == 1) then
          Das2then=Das2now
          Das2now=D0as2*exp(-Ea/R/(temp+273.))
          f1(2)=alpha*dt*Das2now/dr**2
          f2(2)=(1.-alpha)*dt*Das2then/dr**2
        endif

        ! Calculate new diffusion parameters (medium grain)
        if (age_flags(3) == 1) then
          Dam2then=Dam2now
          Dam2now=D0am2*exp(-Ea/R/(temp+273.))
          f1(3)=alpha*dt*Dam2now/dr**2
          f2(3)=(1.-alpha)*dt*Dam2then/dr**2
        endif

        ! Calculate new diffusion parameters (large grain)
        if (age_flags(4) == 1) then
          Dal2then=Dal2now
          Dal2now=D0al2*exp(-Ea/R/(temp+273.))
          f1(4)=alpha*dt*Dal2now/dr**2
          f2(4)=(1.-alpha)*dt*Dal2then/dr**2
        endif

        ! Calculate new diffusion parameters (low damage)
        if (age_flags(5) == 1) then
          Da2then_l=Da2now_l
          Da2now_l=(D0a2rd*exp(-Eard/R/(temp+273.)))/((psi*He_l*exp(Et/R/(temp+273.)))+1)
          f1(5)=alpha*dt*Da2now_l/dr**2
          f2(5)=(1.-alpha)*dt*Da2then_l/dr**2
        endif

        ! Calculate new diffusion parameters (moderate damage)
        if (age_flags(6) == 1) then
          Da2then_m=Da2now_m
          Da2now_m=(D0a2rd*exp(-Eard/R/(temp+273.)))/((psi*He_m*exp(Et/R/(temp+273.)))+1)
          f1(6)=alpha*dt*Da2now_m/dr**2
          f2(6)=(1.-alpha)*dt*Da2then_m/dr**2
        endif

        ! Calculate new diffusion parameters (high damage)
        if (age_flags(7) == 1) then
          Da2then_h=Da2now_h
          Da2now_h=(D0a2rd*exp(-Eard/R/(temp+273.)))/((psi*He_h*exp(Et/R/(temp+273.)))+1)
          f1(7)=alpha*dt*Da2now_h/dr**2
          f2(7)=(1.-alpha)*dt*Da2then_h/dr**2
        endif

        ! Calculate new diffusion parameters (zircon grain)
        if (age_flags(8) == 1) then
          Dza2then = Dza2now
          Dza2now = D0za2*exp(-Eaz/R/(temp+273.0))
          f1(8) = alpha*dt*Dza2now/dr**2
          f2(8) = (1.0-alpha)*dt*Dza2then/dr**2
        endif


        do i=12,num_of_age_flags
            if (age_flags(i) == 1) then
                Da2then_ar(i) = Da2now_ar(i)
                Da2now_ar(i) = D0a2_ar(i) * exp(-Ea_ar(i) / R / (temp + 273.0))
                f1(i) = alpha * dt * Da2now_ar(i) / dr**2
                f2(i) = (1.0 - alpha) * dt * Da2then_ar(i) / dr**2
            endif
        end do

       ! Fill 2-D arrays to hold the seven age calculations (added by cspath 03/07)
       do k=1,num_of_age_flags
         if (age_flags(k) == 1) then
           do i=2,n-1
             diag(k,i)=1.+2.*f1(k)
             sup(k,i)=-f1(k)*(1.+beta/(i-1)/2.)
             inf(k,i)=-f1(k)*(1.-beta/(i-1)/2.)
             f(k,i)=age(k,i)+f2(k)* &
                    ((age(k,i+1)-2.*age(k,i)+age(k,i-1)) + &
                    beta*(age(k,i+1)-age(k,i-1))/(i-1)/2.)+dt
           enddo
         endif
       enddo

       ! Calculate the seven different ages (added by cspath 03/07)
       do k=1,num_of_age_flags
         if (age_flags(k) == 1) then
           diag(k,1)=1.
           sup(k,1)=-1.
           f(k,1)=0.
           diag(k,n)=1.
           inf(k,n)=0.
           f(k,n)=0.

           inf_temp = inf(k,:)
           diag_temp = diag(k,:)
           sup_temp = sup(k,:)
           f_temp = f(k,:)
           age_temp = age(k,:)

           ! Invert tri-diagonal matrices
           call tridag (inf_temp,diag_temp,sup_temp,f_temp,age_temp,n)
           !call tridag (inf,diag,sup,f,age,n,k,8)

           ! Calculate apparent age for current time step
           agei(k)=0.
           do i=1,n
             age(k,i) = age_temp(i)
             fact=1.
             if (i.eq.1.or.i.eq.n) fact=0.5
             agei(k)=agei(k)+age(k,i)*fact*dr**3*(i-1)*(i-1)
           enddo
           agei(k)=3.*agei(k)!*Ft(k)

         endif
       enddo

       ! Update grain He concentrations (added by cspath 03/07)
       !He_l=agei(5)*P_l        ! He concentration (low damage) [nmol/g]
       !He_m=agei(6)*P_m        ! He concentration (moderate damage) [nmol/g]
       !He_h=agei(7)*P_h        ! He concentration (high damage) [nmol/g]

       ! Calculates radiation damage factor as function of total alpha decay
        He_l=He_l+dt*P_l
        He_m=He_m+dt*P_m
        He_h=He_h+dt*P_h

        !He_m=maxval(time)*P_m

       ! Debugging (added by dwhipp 03/07)

       enddo ! istep

      enddo ! itime

      ! Store ages array in new array
      output_age=agei

      deallocate (age,diag,sup,inf,f)
      deallocate (agei,f1,f2)
      deallocate (Da2now_ar)
      deallocate (Da2then_ar)
      deallocate (D0a2_ar)
      deallocate (Ea_ar)

      deallocate(inf_temp)
      deallocate(diag_temp)
      deallocate(sup_temp)
      deallocate(f_temp)
      deallocate(age_temp)

      return
      end subroutine Mad_He
end module m_mad_he
