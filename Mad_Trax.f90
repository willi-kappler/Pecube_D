module m_mad_trax

contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                           c
!  MMM        MMM                TTTTTTTTTT                 c
!  MMMM      MMMM                TTTTTTTTTT                 c
!  MMMMM    MMMMM                    TT                     c
!  MM  MM  MM  MM                    TT                     c
!  MM   MMMM   MM             d      TT                     c
!  MM    MM    MM             d      TT                     c
!  MM          MM  aaa a  ddd d      TT r rrr  aaa a x   x  c
!  MM          MM a   aa d   dd      TT rrr   a   aa  x x   c
!  MM          MM a   aa d   dd      TT r     a   aa   x    c
!  MM          MM a   aa d   dd      TT r     a   aa  x x   c
!  MM          MM  aaa a  ddd d      TT r      aaa a x   x  c
!                                                           c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c

      subroutine Mad_Trax (time_temp,temp_i,n,out_flag,param_flag, &
                          fta, ftldmean, ftldsd, ftld)


! subroutine Mad_Trax to calculate fission track age, track length
! distribution and statistics from a given thermal history
!c
! in input:
! real*4  time_i(n)  :   the time values (in Myr) in descending order
!                        at which the thermal history is given
!                        (ex: 100,50,20,10,0); the last value
!                        should always be 0; the first value
!                        should be smaller than 1000.
! real*4  temp_i(n)  :   the thermal history in degree Celsius
! integer n          :   the number of time-temperature pairs used
!                        to describe the temperature history
! integer out_flag   :   =0 only calculate fission track age
!                        =1 also calculate track length distribution
!                           and statistics
! integer param_flag :   =1 uses Laslett et al, 1987 parameters
!                        =2 uses Crowley et al., 1991 Durango parameters
!                        =3 uses Crowley et al., 1991 F-apatite parameters
!c
! in output:
! real*4  fta        :   fission track age in Myr
! real*4  ftld(17)   :   normalised track length distribution where
!                        ftld(k) is the percentage of track with
!                        length between k-0.5 and k+0.5 microns
! real*4  ftldmean   :   mean track length in microns
! real*4  ftldsd     :   track length standard deviation in microns
!c
! This subroutine is based on the subroutine "ftmod.pas" provided by
! Peter vanderBeek in December 1995. The algorithm is explained in
! Peter's PhD thesis and is based on the work by Lutz and Omar (1991)
!c
! References:
!c
! VanderBeek, P., 1995.Tectonic evolution of continental rifts, PhD Thesis,
!       Faculty of Earth Sicences, Free University, Amsterdam.
!c
! Lutz, T.M. and Omar, G.I., 1991. An inverse method of modeling thermal
!       histories from apatite fission-track data. EPSL, 104, 181-195.
!c
! Laslett, G.M., Green, P.F., Duddy, I.R. and Gleadow, A.J.W., 1987. Thermal
!       annealing of fission tracks in apatite 2. A quantitative analysis.
!       Chem. Geol. (Isot. Geosci. Sect.) 65, 1-13.
!c
! Crowley, K.D., Cameron, M. and Schaefer, R.L., 1991. Experimental studies
!       of annealing of etched fission tracks in fluorapatite. Geochim.
!       Cosmochim. Acta, 55, 1449-1465.
!c

      implicit none

      real(8) :: temp_i(n),time_i(n),time_temp(n)
      real(8) :: r(1000),prob(101)
      real(8) :: rp, dj, x, sumprob, sumftld, devftld, sumdj
      real(8) :: tempp, temp, a, b, time, tempm, teq
      real(8) :: gr, dfr, dt, fr, c0, c1, c2, c3, deltat, h
      real(8) :: rr, time_interval, xind, xfct

      real(8), intent(out) :: ftldmean, ftldsd, ftld(17), fta

      integer(4) :: out_flag,param_flag, n, i, imax, imin, j, l
      integer(4) :: nstep

! Adjusts time values to relfect current model run time
    time_i=time_temp
    time_i=time_i-time_i(n)

      if (param_flag.eq.1) then

! from Laslett et al., 1987

      a=0.35
      b=2.7
      c0=-4.87
      c1=0.000168
      c2=0.00472416
      c3=0.

      elseif (param_flag.eq.2) then

! from Crowley et al., 1991 (Durango)

      a=0.49
      b=3.0
      c0=-3.202
      c1=0.0000937
      c2=0.001839
      c3=0.00042

      elseif (param_flag.eq.3) then

! from Crowley et al., 1991 (F-apatite)
! note: modified to use correct beta value for fluorapatite from
! Crowley et al., 1991 - dwhipp 02/07

      a=0.76
      b=4.30
      c0=-1.508
      c1=0.00002076
      c2=0.0002143
      c3=0.0009967

      else

! lc,mod fanning Arrhenius model from Ketcham et al., 1999 table 5e
! *** DO NOT USE: Not yet functional ***
! added by dwhipp (08/07)

      a=-0.05771
      b=-13.218
      c0=-9.0722
      c1=0.00029896
! Original c2 of Ketcham et al., 1999 modified for version of track length reduction
! equation used below
      c2=-15.846
      c3=0.00076370

      endif

! unannealed fission track length

      xind=16.

! mean length of spontaneous tracks in standards

      xfct=14.5

! calculate the number of time steps assuming 1My time step length
! if run time > 100 My, else take 100 time steps

      nstep=int(time_i(1))
!     if (nstep.gt.8000) stop 'Your final time is greater than '//
!    &                        'the age of the Universe...'
!     if (nstep.gt.4500) stop 'Your final time is greater than '//
!    &                        'the age of the Earth...'
!     if (nstep.gt.1000) stop 'Fission track does not work very well '//
!    &                        'for time spans greater than 1Byr...'
      if (nstep.gt.1000) stop
      time_interval=1.
      if (nstep.lt.100) then
        nstep=100
        time_interval=time_i(1)/100.
      endif
      deltat=time_interval*1.e6*365.24*24.*3600.

! calculate final temperature

      tempp=temperature(temp_i,time_i,0.0_8,n)+273.0_8
      rp=0.5_8

! begining of time stepping

        do i=1,nstep

        time=float(i)*time_interval

! calculate temperature by linear interpolation

        temp=temperature(temp_i,time_i,time,n)+273.0_8

! calculate mean temperature over the time step

        tempm=(temp+tempp)/2.0_8

! calculate the "equivalent time", teq

        teq=exp((-c2/c1)+((g(rp,a,b)-c0)/c1)*(1./tempm-c3))
        if (i.eq.1) teq=0.0_8

! check if we are not getting too close to r=0
! in which case r remains 0 for all following time steps

          if (dlog(teq+deltat).ge. &
          (expos(1.0_8/b,a)-a*c0-1.0_8)/a/c1*(1.0_8/tempm-c3)-c2/c1) then

            do j=i,nstep
            r(j)=0.
            enddo
          nstep=i
          goto 90

          endif

! otherwise calculate reduction in length, r, over the time step, dt

        dt=teq+deltat
        gr=c0+((c1*dlog(dt)+c2)/((1./tempm)-c3))
! equation for decrease in normalized track length from Ketcham et al., 1999
! note: this is slightly different than above, and may require other code mods to implement
! also note: this is currently not functional
! dwhipp - (08/07)
!        gr=c0+c1*((alog(dt)+c2)/((1./tempm)-c3))
        r(i)=xinv(gr,a,b)

! update variables for next time step

        tempp=temp
        rp=r(i)

!       print*,i,time,temp,r(i)

        enddo

90    continue

! all reduction factors for all time steps have been calculated
! now estimate the fission track age by simple summation
! (here it helps to use 1Myr time steps)

      sumdj=0.

        do i=1,nstep

          if (r(i).le.0.35) then
          dj=0.
          elseif (r(i).le.0.66) then
          dj=2.15*r(i)-0.76
          else
          dj=r(i)
          endif

        sumdj=sumdj+dj

        enddo

      fta=real((xind/xfct)*sumdj*time_interval)

! now (if out_flag.ne.0) let's do some statistics

      if (out_flag.eq.0) return

! first, calculate probability density function using Luts and Omar (1991)
! method and assuming a Gaussian distribution

      sumprob=0.

        do j=1,101

        rr=float(j-1)/100.

          if (rr.le..43) then
          h=2.53
          elseif (rr.le..67) then
          h=5.08-5.93*rr
          else
          h=1.39-.61*rr
          endif

        fr=0.

          do i=1,nstep

          x=(rr-r(i))*xind/h
          dfr=xk(x)/h
          fr=fr+dfr

          enddo

        prob(j)=fr/nstep
        sumprob=sumprob+prob(j)

        enddo

! now let's rescale the track length distribution, its mean and standard
! deviation

      ftld(17)=100.
      imin=1

        do l=1,16

        imax=int(l*100./xind)
        ftld(l)=0.

          do i=imin,imax
          ftld(l)=ftld(l)+prob(i)
          enddo

        ftld(l)=(ftld(l)*100./sumprob)
        ftld(17)=ftld(17)-ftld(l)

        imin=imax+1

        enddo

      sumftld=0.

        do l=1,17
        sumftld=sumftld+ftld(l)*(float(l)-0.5)
        enddo

      ftldmean=sumftld/100.

      devftld=0.

        do l=1,17
        devftld=devftld+ftld(l)*(float(l)-0.5-ftldmean)**2
        enddo

      ftldsd=sqrt(devftld/100.)

      return
      end subroutine

!c---
      function g(r,a,b)

      implicit none

      real(8) r, a, b, g
      g=(expos((1.0_8-expos(r,b))/b,a)-1.0_8)/a

      return
      end function

!c---
      function xinv (gr,a,b)

      implicit none

      real(8) gr, a, b, xinv

      xinv=expos(1.0_8-b*expos(a*gr+1.0_8,1.0_8/a),1.0_8/b)

      return
      end function

!c---
      function temperature (temp,time,t,n)

      implicit none

      integer(4) :: n, i

      real(8)   temp(n),time(n)
      real(8)   temperature, t, rat

      temperature=temp(1)

        do i=1,n-1
          if ((t-time(i))*(t-time(i+1)).le.0.) then
          rat=(t-time(i))/(time(i+1)-time(i))
          temperature=temp(i)+rat*(temp(i+1)-temp(i))
          return
          endif
        enddo

      temperature=temp(n)

      return
      end function

!c---
      function expos (x,a)

      implicit none

      real(8) x, a, expos
      expos=exp(a*dlog(x))

      return
      end function

!c---
      function xk (x)

      implicit none

      real(8) x, xk
      xk=0.39894228*exp(-(x**2)/2.)

      return
      end function
end module m_mad_trax
