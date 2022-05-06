module m_mad_trax_zirkon
contains

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!

      subroutine Mad_Zirc (time_i,temp_i,n,out_flag,param_flag,&
                          fta)

! subroutine Mad_Trax to calculate zircon fission track age from a given thermal history
!
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
! NB integer out_flag=1 not implemented yet for Zircon !

! integer param_flag :   =1 uses parameters for alpha-damaged zircon
!                        =2 uses parameters for zero-damage zircon
!                        (parameters from Rahn et al. 2004)
!
! in output:
! real*4  fta        :   fission track age in Myr
! real*4  ftld(17)   :   normalised track length distribution where
!                        ftld(k) is the percentage of track with
!                        length between k-0.5 and k+0.5 microns
! real*4  ftldmean   :   mean track length in microns
! real*4  ftldsd     :   track length standard deviation in microns
!
! This subroutine is based on the subroutine "ftmod.pas" provided by
! Peter van der Beek in December 1995. The algorithm is explained in
! Peter's PhD thesis and is based on the work by Lutz and Omar (1991)
! This adaptation of the program for zircon fission-track annealing was
! written by Peter in August/September 2006 and is based on algorithms
! given by Galbraith & Laslett (1997), Tagami et al. (1998) and
! Rahn et al. (2004)
!
! References:
!
! van der Beek, P., 1995.Tectonic evolution of continental rifts, PhD Thesis,
!       Faculty of Earth Sicences, Free University, Amsterdam.
!
! Lutz, T.M. and Omar, G.I., 1991. An inverse method of modeling thermal
!       histories from apatite fission-track data. EPSL, 104, 181-195.
!
! Galbraith, R. F., and G. M. Laslett (1997), Statistical modelling of thermal
!       annealing of fission tracks in zircon, Chemical Geology, 140, 123-135.
!
! Tagami, T., et al. (1998), Revised annealing kinetics of fission tracks in
!        zircon and geological implications, in Advances in Fission-Track Geochronology,
!        edited by P. Van den haute and F. De Corte, pp. 99-112, Kluwer Academic
!        Publishers, Dordrecht, Netherlands.
!
! Rahn, M. K., et al. (2004), A zero-damage model for fission track annealing in zircon,
!        American Mineralogist, 89, 473-484.


      real(8) :: temp_i(n), time_i(n), a, b, c, fta
      real(8) :: local_time(n), time_interval, deltat, tempp, time, temp
      real(8) :: r(1000), prob(101), tempm, teq, dt, gr, rp, dj, x
      real(8) :: dfr, sumprob, sumdj, sumftld, devftld, fr
      real(8) :: ftldmean, ftldsd
      real(8) :: ftld(17)

      integer  out_flag,param_flag

      local_time = time_i - time_i(n)
!      local_time = time_i

      if (param_flag.eq.1) then

! alpha-damaged zircon

      a=-10.77
      b=2.599E-4
      c=1.026E-2

      else

! zero-damage zircon

      a=-11.57
      b=2.755E-4
      c=1.075E-2

      endif

! unannealed fission track length

      xind=10.8

! mean length of spontaneous tracks in standards

      xfct=10.8

! calculate the number of time steps assuming 1My time step length
! if run time > 100 My, else take 100 time steps

      nstep=int(local_time(1))
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
        time_interval=local_time(1)/100.
      endif
      deltat=time_interval*1.e6*365.24*24.*3600.

! calculate final temperature

      tempp = temperaturezr(temp_i,local_time,0.0_8,n) + 273.0
      rp=.5

! beginning of time stepping

        do i=1,nstep

        time=float(i)*time_interval

! calculate temperature by linear interpolation

        temp=temperaturezr(temp_i,local_time,time,n)+273.

! calculate mean temperature over the time step

        tempm=(temp+tempp)/2.

! calculate the "equivalent time", teq

        if (i.eq.1)then
          teq=0.
        else
          teq=exp((gzr(rp)-a-(c*tempm))/(b*tempm))
        endif


! check if we are not getting too close to r=0
! in which case r remains 0 for all following time steps

!        if (param_flag.lt.4) then
!           rcrit=(expos(1./b,a)-a*c0-1.)/a/c1*(1./tempm-c3)-c2/c1
!        else
!           rcrit=(expos(1./b,a)-a*c0-1.)/a/c1*(alog(1./tempm)-c3)-c2/c1
!        endif

!          if (alog(teq+deltat).ge.rcrit) then
!            do j=i,nstep
!            r(j)=0.
!            enddo
!          nstep=i
!          goto 90

!          endif

! otherwise calculate reduction in length, r, over the time step, dt

        dt=teq+deltat

        gr=a+((b*tempm)*log(dt))+(c*tempm)

        r(i)=xinvzr(gr)
! stop calculation if r<0.4
        if (r(i).le.0.4) then
          nstep=i
          goto 90
        endif

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

          if (r(i).le.0.4) then
          dj=0.
          elseif (r(i).le.0.66) then
          dj=2.15*r(i)-0.76
          else
          dj=r(i)
          endif

        sumdj=sumdj+dj

        enddo

      fta=(xind/xfct)*sumdj*time_interval

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
          dfr=xkzr(x)/h
          fr=fr+dfr

          enddo

        prob(j)=fr/nstep
        sumprob=sumprob+prob(j)

        enddo

! now let's rescale the track length distribution, its mean and standard
! deviation

      ftld(11)=100.
      imin=1

        do l=1,10

        imax=int(l*100./xind)
        ftld(l)=0.

          do i=imin,imax
          ftld(l)=ftld(l)+prob(i)
          enddo

        ftld(l)=(ftld(l)*100./sumprob)
        ftld(11)=ftld(11)-ftld(l)

        imin=imax+1

        enddo

      sumftld=0.

        do l=1,11
        sumftld=sumftld+ftld(l)*(float(l)-0.5)
        enddo

      ftldmean=sumftld/100.

      devftld=0.

        do l=1,11
        devftld=devftld+ftld(l)*(float(l)-0.5-ftldmean)**2
        enddo

      ftldsd=sqrt(devftld/100.)

      return
      end subroutine

! ---
      function gzr(r)
        real(8) :: gzr, r

        gzr = log(1.0 - r)

      return
      end function

! ---
      function xinvzr (gr)
        real(8) :: gr, xinvzr

        xinvzr = 1 - exp(gr)

      return
      end function

! ---
      function temperaturezr (temp,time,t,n)

      real(8) :: temp(n), time(n), t
      real(8) :: temperaturezr, rat

      temperaturezr = temp(1)

        do i=1,n-1
          if ((t-time(i))*(t-time(i+1)).le.0.) then
          rat=(t-time(i))/(time(i+1)-time(i))
          temperaturezr=temp(i)+rat*(temp(i+1)-temp(i))
          return
          endif
        enddo

      temperaturezr = temp(n)

      return
      end function


! ---
      function xkzr (x)
        real(8) :: xkzr, x

        xkzr = 0.39894228 * exp(-(x**2) / 2.0)

      return
      end function
end module m_mad_trax_zirkon
