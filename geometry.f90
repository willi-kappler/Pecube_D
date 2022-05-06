  module m_geometry
      use m_logger
      use m_data_structures

      implicit none

      contains

      subroutine geometry (xin, yin, zin, vxp, vyp, vzp, velo_info, id)

! Subroutine that calculates velocities in the x, y, and z directions
! based on the velocity field geometry specified in Pecube.in
! 1 = Vertical motion
! 2 = Diagonal motion
! 3 = Listric fault
! 4 = Nepalese Himalaya geometry for rotated model (dwhipp - 08/07)
! 5 = Parabolic uplift field
! 6 = Ellipsoid Uplift Field (mschmiddunser - 07/10)
! 7 = Ellipsoid Uplift Field with smoothed edges (mschmiddunser - 07/10)
! 8 = free velocity definition for each time step
! For given values of xp,yp,zp this routine should return a
! velocity field in the form of a triplet vxp,vyp,vzp

! The size of the box is given in terms of zl
! note that zp is given from the bottom of the box up
! (it is not the depth but the height of the box)

      use m_global_velocities

      real*8 x0,depth,list,rad_theta,rad_phi,kb15s,kb15i,kb45s,kb45i
      real*8 mft_dipin,mbt_dipin,mct_dipin,vin,stf_dipin,xx
      real*8 yy,x,y,z,mht1_dipin,mht2_dipin,mht3_dipin
      real*8 mht4_dipin,mfts,mfti,mbts,mbti,mcts,mcti,mht1s
      real*8 mht1i,mht2s,mht2i,mht3s,mht3i,mht4s,mht4i,mht4dip,kb0s
      real*8 kb0i,kb1s,kb1i,kb2s,kb2i,kb3s,kb3i,kb4s,kb4i,kb5s,kb5i
      real*8 kb6s,kb6i,kb7as,kb7ai,kb7bs,kb7bi,kb55as,kb55ai,kb55bs
      real*8 kb55bi,SHAa,SHAb,SHBA,LHAa,LHAb,LHAc,LHBA,GHSWl,GHSEl,GHSWw
      real*8 GHSEw,GHSBw,GHSBe,vup_sh,vup_lh,vup_gh,mftrate,mbtrate
      real*8 mctrate,stfrate,utrate,mbtdip,stfdip,mht1dip,mht2dip
      real*8 mht3dip,vergin,ut_ratein
      real*8 pie,halfpi,mftdip,mctdip,verge
      real*8 axis,azimuth,lambda,lambda_corr,k,dx,dy,phi1,phi2,h,d
      real*8 a1,a2,dist,vz,veloin,veloout,amed,adif

      real(8) :: alpha, d_corr, df, dip, f, stf_seg, stfai, stfas, stfbi
      real(8) :: stfbs, vnorm, vx, vy, x1, x1p, x2, x2p, xb, xp, y1, y1p
      real(8) :: y2, y2p, yp, zb, zp

      real(8), intent(in) :: xin, yin, zin
      real(8), intent(out) :: vxp, vyp, vzp

      integer(4) :: i
      integer(4), intent(in) :: id

      type(velocity_info_t), intent(in) :: velo_info

! reset velocities

      vxp = 0.0
      vyp = 0.0
      vzp = 0.0

! x and y velocities set to 0
! z velocity set to 1.0 * Peclet
      if (velo_info%geoflag == 1) then

         vxp = 0.d0
         vyp = 0.d0
         vzp = velo_info%Peclet

         return

! Structure is shifting on angle theta from horizontal
! Angle phi is the angle relative to the xy plane
      else if (velo_info%geoflag == 2) then

        rad_theta = velo_info%theta * (3.141592654/180)
        rad_phi = velo_info%phi * (3.141592654/180)

        vxp = cos(rad_theta) * sin(rad_phi)
        vyp = cos(rad_theta) * cos(rad_phi)
        vzp = sin(rad_theta) * velo_info%Peclet

        return

! Assume a listric fault and a velocity field with velocities parallel
! to the fault
      else if (velo_info%geoflag == 3) then
        x1 = velo_info%x1f
        y1 = velo_info%y1f
        x2 = velo_info%x2f
        y2 = velo_info%y2f

        alpha = atan2(x1 - x2, y2 - y1)

        xp = xin*cos(alpha) + yin*sin(alpha)
        yp =-xin*sin(alpha) + yin*cos(alpha)
        zp = zin
        x1p = x1*cos(alpha) + y1*sin(alpha)
        y1p =-x1*sin(alpha) + y1*cos(alpha)
        x2p = x2*cos(alpha) + y2*sin(alpha)
        y2p =-x2*sin(alpha) + y2*cos(alpha)

        ! x0 is the position of the fault tip in the rotated system
        x0 = x1p

        ! depth is the fault soling depth
        depth = velo_info%def

        ! list is the listricity of the fault (surface dip=depth/list)
        dip = velo_info%dif
        list = depth / tan(dip * 3.141592654d0/180.d0)

        if (list /= list) then
            !print *, "list is NaN (geometry line 117)"
            stop
        endif

        xx = x0 - xp
        yy = yp
        z = velo_info%zl - zp

        ! calculates velocities according to algorithm described in
        ! notes/manuals

        if (xx.lt.0.) return
        if (yy.lt.y1p .or. yy.gt.y2p) return
        if (z.gt.depth*(1.-exp(-xx/list))) return

        xb = xx
        if (xb /= xb) then
            !print *, "xb is NaN (geometry line 136)"
            !print *, xx, x0, xp
            stop
        endif

        do i = 1, 3
          f=xx-xb+depth/list*(z-depth)*exp(-xb/list)+depth**2/list*exp(-2.*xb/list)
          df=-1.-depth/list**2*(z-depth)*exp(-xb/list)-2.*depth**2/list**2*exp(-2.*xb/list)
          xb = xb - f / df


          ! if dip > 68.5 xb becomes NaN
          ! Thus we limit xb to values >= -300
          ! 2016.05.11, WK, Karl
          if (xb < -300.0) then
              xb = -300.0
          endif


          !print *, "i, f, df, xb, x, y, z:", i, f, df, xb, xin, yin, zin
          if (xb /= xb) then
              ! print *, "xb is NaN (geometry line 146)"
              ! print *, i, z, f, df, xx, depth, list
              stop
          endif
        enddo
      zb=depth*(1.-exp(-xb/list))
      vnorm=sqrt(1.+depth**2/list**2*exp(-2.*xb/list))
      vx=1./vnorm
        if (vx /= vx) then
            !print *, "vx is NaN (geometry line 144)"
            !print *, vnorm, xb, list, depth
            stop
        endif

      vy=0.d0
      vz=depth/list*exp(-xb/list)/vnorm

! 'un'change of coordinate system
! Note: Possible coordinate system unrotation syntax problem
! Changed vyp equation to replace vy with vx
      vxp=vx*cos(alpha)-vy*sin(alpha)

        if (vxp /= vxp) then
            !print *, "vxp is NaN (geometry line 152)"
            !print *, vx, vy, alpha
            stop
        endif



      vyp=vx*sin(alpha)+vy*cos(alpha)

      vzp = vz * velo_info%Peclet

      return

      ! New Nepal model geometry using rotated grid - dwhipp (08/07)
    else if (velo_info%geoflag == 4) then

        !***
        ! Geometry variables (don't mess with these unless you know what you're doing)
        !***
        vergin=180          ! Convergence direction from GPS (Jouanne et al., 2004)
        mft_dipin=-40           ! From Avouac cross-section (degrees)
        mbt_dipin=-40           ! From Avouac cross-section (degrees)
        ! Original MCT dip commented out, see below
        !mct_dipin=-28          ! From Searle & Godin, 2003 - (22-28 deg is range they list)
        mct_dipin=-40                   ! From argument that high relief, physiographic transition zone should be south of MCT kink band plane
                                        !   Also, DeCelles et al., 2001 lists a range of 30-48 deg for MCT dip
        stf_dipin=-28           ! Assume same dip as MCT - Searle & Godin, 2003 say the foliation dips less than 30 deg NNE
        mht1_dipin=-4           ! From Avouac cross-section (degrees)
        mht2_dipin=-7           ! From Avouac cross-section (degrees)
        mht3_dipin=-19          ! From Avouac cross-section (degrees)
        mht4_dipin=-10          ! From Avouac cross-section (degrees)
        vin=20              ! From Lave & Avouac, 2000 - Consistent with convergence rates along MFT during Holocene
        x = xin             ! Rename xin
        y = yin             ! Rename yin
        z = zin-velo_info%zl          ! Shift z-points down into elevations ASL (km)
        stf_seg=51.3            ! Distance from left edge of model to tear in STF (km)

    ! Fault plane slopes/intercepts
    mfts=-0.70          ! MFT slope - From Avouac cross-section (unitless)
    mfti=21.25          ! MFT z-intercept - From Avouac cross-section (km)
    mbts=-0.70          ! MBT slope - From Avouac cross-section (unitless)
    mbti=40.10          ! MBT z-intercept - From Avouac cross-section (km)
        ! Original MCT slope/intercept commented out, see below
    !mcts=-0.49         ! MCT slope - From Searle & Godin, 2003 (unitless)
    !mcti=59.34         ! MCT z-intercept - Calculated from Searle & Godin, 2003 slope with
                    ! Lave & Avouac, 2000 cross-section (km)
        mcts=-0.70                      ! MCT slope - From argument that high relief, physiographic transition zone should be south of MCT kink band plane (unitless)
                                        !   Also, DeCelles et al., 2001 lists a range of 30-48 deg for MCT dip
        mcti=82.09                      ! MCT z-intercept - Calculated from above with Avouac cross-section (km)
    stfas=-0.49         ! STF-W slope - From Searle & Godin, 2003 (unitless)
    stfai=70.42         ! STF-W z-intercept - From Avouac cross-section (km)
        stfbs=-0.49         ! STF-E slope - From Searle & Godin, 2003 (unitless)
        stfbi=88.99         ! STF-E z-intercept - From Avouac cross-section (km)
    mht1s=-0.07         ! MHT section 1 slope - From Avouac cross-section (unitless)
    mht1i=-0.67         ! MHT section 1 z-intercept - From Avouac cross-section (km)
    mht2s=-0.12         ! MHT section 2 slope - From Avouac cross-section (unitless)
    mht2i=2.72          ! MHT section 2 z-intercept - From Avouac cross-section (km)
    mht3s=-0.33         ! MHT section 3 slope - From Avouac cross-section (unitless)
    mht3i=24.46         ! MHT section 3 z-intercept - From Avouac cross-section (km)
    mht4s=-0.17         ! MHT section 4 slope - From Avouac cross-section (unitless)
    mht4i=0.94          ! MHT section 4 z-intercept - From Avouac cross-section (km)

    ! Kink band plane slopes/intercepts
    kb0s=2.75           ! Kink band plane 0 slope - From Avouac cross-section (unitless)
    kb0i=-84.62         ! Kink band plane 0 z-intercept - From Avouac cross-section (km)
    kb1s=2.48           ! Kink band plane 1 slope - From Avouac cross-section (unitless)
    kb1i=-89.46         ! Kink band plane 1 z-intercept - From Avouac cross-section (km)
    kb15s=2.75          ! Kink band plane 1.5 slope - From Avouac cross-section (unitless)
    kb15i=-155.95           ! Kink band place 1.5 z-intercept - From Avouac cross-section (km)
    kb2s=10.39          ! Kink band plane 2 slope - From Avouac cross-section (unitless)
    kb2i=-679.11            ! Kink band plane 2 z-intercept - From Avouac cross-section (km)
    kb3s=2.30           ! Kink band plane 3 slope - From Avouac cross-section (unitless)
    kb3i=-154.44            ! Kink band plane 3 z-intercept - From Avouac cross-section (km)
    kb4s=4.33           ! Kink band plane 4 slope - From Avouac cross-section (unitless)
    kb4i=-459.99            ! Kink band plane 4 z-intercept - From Avouac cross-section (km)
        ! Original kb4.5 slope and intercept commented out, see above for reason why (line ~158)
    !kb45s=4.01         ! Kink band plane 4.5 slope - From Avouac cross-section (unitless)
    !kb45i=-478.35          ! Kink band place 4.5 z-intercept - From Avouac cross-section (km)
        kb45s=2.75                      ! Kink band plane 4.5 slope - From Avouac cross-section (unitless)
        kb45i=-326.90                   ! Kink band place 4.5 z-intercept - From Avouac cross-section (km)
    kb5s=3.87           ! Kink band plane 5 slope - From Avouac cross-section (unitless)
    kb5i=-588.57            ! Kink band plane 5 z-intercept - From Avouac cross-section (km)
        kb55as=4.01         ! Kink band plane 5.5a slope - From Avouac cross-section (unitless)
        kb55ai=-536.51          ! Kink band plane 5.5a z-intercept - From Avouac cross-section (km)
        kb55bs=4.01         ! Kink band plane 5.5b slope - From Avouac cross-section (unitless)
        kb55bi=-688.92          ! Kink band plane 5.5b z-intercept - From Avouac cross-section (km)
        ! Original kb6 slope and intercept commented out, see above for reason why (line ~158)
    !kb6s=2.90          ! Kink band plane 6 slope - From Avouac cross-section (unitless)
    !kb6i=-571.37           ! Kink band plane 6 z-intercept - From Avouac cross-section (km)
        kb6s=2.14                       ! Kink band plane 6 slope - From Avouac cross-section (unitless)
        kb6i=-358.20                    ! Kink band plane 6 z-intercept - From Avouac cross-section (km)
    kb7as=2.90          ! Kink band plane 7a slope - From Avouac cross-section (unitless)
    kb7ai=-660.71           ! Kink band plane 7a z-intercept - From Avouac cross-section (km)
    kb7bs=2.90          ! Kink band plane 7b slope - From Avouac cross-section (unitless)
    kb7bi=-848.62           ! Kink band plane 7b z-intercept - From Avouac cross-section (km)

        !***
        !  Convert units of input variables / Additional math
        !***
    vin = vin + velo_info%stf_ratein                      ! Account for extension on stf in keeping constant convergence rate
    ut_ratein = (vin - velo_info%mft_ratein - velo_info%mbt_ratein - velo_info%mct_ratein)    ! Underthrusting rate is equal to that not accounted for by the MFT and MCT
    mftrate=-velo_info%mft_ratein                     ! Rename fault convergence rates (for simplicity)
    mbtrate=-velo_info%mbt_ratein                     ! Rename fault convergence rates (for simplicity)
    mctrate=-velo_info%mct_ratein                     ! Rename fault convergence rates (for simplicity)
    stfrate=velo_info%stf_ratein                      ! Rename fault convergence rates (for simplicity)
    utrate=ut_ratein                        ! Rename fault convergence rates (for simplicity)
    pie=3.14159265                          ! define pi
    halfpi=pie/2                            ! define pi over two
    mftdip=(pie/180.0)*mft_dipin                    ! convert from deg to rad
    mbtdip=(pie/180.0)*mbt_dipin                    ! convert from deg to rad
    mctdip=(pie/180.0)*mct_dipin                    ! convert from deg to rad
    stfdip=(pie/180.0)*stf_dipin                    ! convert from deg to rad
    mht1dip=(pie/180.0)*mht1_dipin                  ! convert from deg to rad
    mht2dip=(pie/180.0)*mht2_dipin                  ! convert from deg to rad
    mht3dip=(pie/180.0)*mht3_dipin                  ! convert from deg to rad
    mht4dip=(pie/180.0)*mht4_dipin                  ! convert from deg to rad
    verge = (pie/180.0)*vergin                  ! convert from deg to rad

        ! Underplating calculations (assuming 100% underplating, where needed)
        SHAa=5.79       ! Surface distance between KB0 and KB1 in Sub-Himalaya [km]
        SHAb=20.76      ! Surface distance between KB1 and KB1.5 in Sub-Himalaya [km]
        SHBA=34.        ! Basal distance for underplating into Sub-Himalaya [km]
        LHAa=10.62      ! Surface distance between KB1.5 and KB3 in Lesser Himalaya [km]
        LHAb=39.01      ! Surface distance between KB3 and KB4 in Lesser Himalaya [km]
        LHAc=12.22      ! Surface distance between KB4 and KB4.5 in Lesser Himalaya [km]
        ! Original LHBA commented out, see above for reason why (line ~158)
        !LHBA=121.      ! Basal distance for underplating into Lesser Himalaya [km]
        LHBA=90.                ! Basal distance for underplating into Lesser Himalaya [km]
        GHSWl=19.76     ! Surface distance of western side of Greater Himalaya [km]
        GHSEl=52.59     ! Surface distance of eastern side of Greater Himalaya [km]
        GHSWw=150.      ! Width of Greater Himalaya south of tear in STF [km]
        GHSEw=51.3      ! Width of Greater Himalaya north of tear in STF [km]
        ! Original GHSB commented out, see above for reason why (line ~158)
        !GHSBw=29.      ! Basal distance for underplating into western Greater Himalaya [km]
        !GHSBe=90.      ! Basal distance for underplating into eastern Greater Himalaya [km]
        GHSBw=60.       ! Basal distance for underplating into western Greater Himalaya [km]
        GHSBe=121.      ! Basal distance for underplating into eastern Greater Himalaya [km]

        ! Additional velocities required for underplating
        ! Sub-Himalaya
        vup_sh=(SHAa*mftrate*tan(mftdip)+SHAb*mftrate*tan(mht1dip))
        !vup_sh=(SHAa*mftrate*tan(mftdip)+SHAb*mftrate*tan(mht1dip))/SHBA
        !vup_sh=0.
        ! Lesser Himalaya
        vup_lh=(LHAa*(mftrate+mbtrate)*tan(mbtdip)+LHAb*(mftrate+mbtrate)*tan(mht2dip)+LHAc*(mftrate+mbtrate)*tan(mht3dip))
        !vup_lh=(LHAa*(mftrate+mbtrate)*tan(mbtdip)+LHAb*(mftrate+mbtrate)*tan(mht2dip)+LHAc*(mftrate+mbtrate)*tan(mht3dip))/LHBA
        !vup_lh=0.
        ! Greater Himalaya
        !vup_gh=(((GHSWl*(mftrate+mbtrate+mctrate)*tan(mctdip))/GHSBw)+((GHSEw/GHSWw)*((GHSEl*(mftrate+mbtrate+mctrate)*tan(mctdip))/GHSBe)))
        !vup_gh=(((GHSWl*(mftrate+mbtrate+mctrate)*tan(mctdip))/GHSBw)+((GHSEw/GHSWw)*((GHSEl*(mftrate+mbtrate+mctrate)*tan(mctdip))/GHSBe)))
        vup_gh=0.

        !***
        ! Define velocity field
        !***
        vxp=0.                                  ! Always set x-component of velocity field to zero

        ! STF is only active fault
        if (velo_info%mft_ratein.eq.0. .and. velo_info%mbt_ratein.eq.0. .and. velo_info%mct_ratein.eq.0.) then
          if (x.le.stf_seg) then                                                ! Point is west of STF tear
            if (y.le.(z-kb55ai)/kb55as) then                                    ! If true, point is south of KB5.5a
              vyp=utrate
              vzp=0.
            else if (y.le.(z-kb7ai)/kb7as) then                                 ! Point is south of KB7a
              if (z.gt.stfas*y+stfai) then                                      ! Point is in STFa hanging wall
                vyp=stfrate
                vzp=stfrate*tan(stfdip)
              else                                                              ! Point is in STFa footwall
                vyp=utrate
                vzp=utrate*tan(stfdip)
              endif
            else                                                                ! Point is north of KB7a
              if (z.gt.mht4s*y+mht4i) then                                      ! Point is in MHT4 hanging wall
                vyp=(stfrate)
                vzp=(stfrate)*tan(mht4dip)
              else                                                              ! Point is in MHT4 footwall
                vyp=utrate
                vzp=utrate*tan(mht4dip)
              endif
            endif
          else                                                                  ! Point is east of STF tear
            if (y.le.(z-kb55bi)/kb55bs) then                                    ! If true, point is south of KB5.5b
              vyp=utrate
              vzp=0.
            else if (y.le.(z-kb7bi)/kb7bs) then                                 ! Point is south of KB7b
              if (z.gt.stfbs*y+stfbi) then                                      ! Point is in STFb hanging wall
                vyp=stfrate
                vzp=stfrate*tan(stfdip)
              else                                                              ! Point is in STFb footwall
                vyp=utrate
                vzp=utrate*tan(stfdip)
              endif
            else                                                                ! Point is north of KB7b
              if (z.gt.mht4s*y+mht4i) then                                      ! Point is in MHT4 hanging wall
                vyp=(stfrate)
                vzp=(stfrate)*tan(mht4dip)
              else                                                              ! Point is in MHT4 footwall
                vyp=utrate
                vzp=utrate*tan(mht4dip)
              endif
            endif
          endif

        ! MCT active (possibly w/ STF)
        else if (velo_info%mft_ratein.eq.0. .and. velo_info%mbt_ratein.eq.0.) then
          if (x.le.stf_seg) then                                                ! Point is west of STF tear
            if (y.le.(z-kb45i)/kb45s) then                                      ! If true, point is south of KB4.5
              vyp=utrate
              vzp=0.
            else if (y.le.(z-kb6i)/kb6s) then                                   ! Point is south of KB6
              if (z.gt.stfas*y+stfai) then                                      ! Point is in STFa hanging wall
                vyp=(mctrate+stfrate)
                if (velo_info%stf_ratein.eq.0.) then                                      ! Ignore STF if inactive
                  vzp=mctrate*tan(mctdip)
                else                                                            ! Move parallel to STF if active
                  vzp=(mctrate+stfrate)*tan(stfdip)
                endif
              else if (z.gt.mcts*y+mcti) then                                   ! Point is in MCT hanging wall
                vyp=mctrate
                vzp=mctrate*tan(mctdip)
              else                                                              ! Point is in MCT footwall
                vyp=utrate
                vzp=utrate*tan(mctdip)
              endif
            else if (y.le.(z-kb7ai)/kb7as) then                                 ! Point is south of KB7a
              if (z.gt.stfas*y+stfai) then                                      ! Point is in STFa hanging wall
                vyp=(mctrate+stfrate)
                if (velo_info%stf_ratein.eq.0.) then                                      ! Ignore STF if inactive
                  vzp=mctrate*tan(mht4dip)
                else                                                            ! Move parallel to STF if active
                  vzp=(mctrate+stfrate)*tan(stfdip)
                endif
              else if (z.gt.mht4s*y+mht4i) then                                 ! Point is in MHT4 hanging wall
                vyp=mctrate
                vzp=mctrate*tan(mht4dip)
              else                                                              ! Point is in MHT4 footwall
                vyp=utrate
                vzp=utrate*tan(mht4dip)+vup_gh
              endif
            else                                                                ! Point is north of KB7a
              if (z.gt.mht4s*y+mht4i) then                                      ! Point is in MHT4 hanging wall
                vyp=(mctrate+stfrate)
                vzp=(mctrate+stfrate)*tan(mht4dip)
              else                                                              ! Point is in MHT4 footwall
                vyp=utrate
                vzp=utrate*tan(mht4dip)
              endif
            endif
          else                                                                  ! Point is east of STF tear
            if (y.le.(z-kb45i)/kb45s) then                                      ! If true, point is south of KB4.5
              vyp=utrate
              vzp=0.
            else if (y.le.(z-kb6i)/kb6s) then                                   ! Point is south of KB6
              if (z.gt.stfbs*y+stfbi) then                                      ! Point is in STFb hanging wall
                vyp=(mctrate+stfrate)
                if (velo_info%stf_ratein.eq.0.) then                                      ! Ignore STF if inactive
                  vzp=mctrate*tan(mctdip)
                else                                                            ! Move parallel to STF if active
                  vzp=(mctrate+stfrate)*tan(stfdip)
                endif
              else if (z.gt.mcts*y+mcti) then                                   ! Point is in MCT hanging wall
                vyp=mctrate
                vzp=mctrate*tan(mctdip)
              else                                                              ! Point is in MCT footwall
                vyp=utrate
                vzp=utrate*tan(mctdip)
              endif
            else if (y.le.(z-kb7bi)/kb7bs) then                                 ! Point is south of KB7b
              if (z.gt.stfbs*y+stfbi) then                                      ! Point is in STFb hanging wall
                vyp=(mctrate+stfrate)
                if (velo_info%stf_ratein.eq.0.) then                                      ! Ignore STF if inactive
                  vzp=mctrate*tan(mht4dip)
                else                                                            ! Move parallel to STF if active
                  vzp=(mctrate+stfrate)*tan(stfdip)
                endif
              else if (z.gt.mht4s*y+mht4i) then                                 ! Point is in MCT hanging wall
                vyp=mctrate
                vzp=mctrate*tan(mht4dip)
              else                                                              ! Point is in MCT footwall
                vyp=utrate
                vzp=utrate*tan(mht4dip)
                if (velo_info%stf_ratein.ne.0.) vzp=vzp+vup_gh                            ! Underplate if needed (to balance erosional removal)
              endif
            else                                                                ! Point is north of KB7b
              if (z.gt.mht4s*y+mht4i) then                                      ! Point is in MHT4 hanging wall
                vyp=(mctrate+stfrate)
                vzp=(mctrate+stfrate)*tan(mht4dip)
              else                                                              ! Point is in MHT4 footwall
                vyp=utrate
                vzp=utrate*tan(mht4dip)
              endif
            endif
          endif

        ! MBT active (possibly w/ MCT)
        else if (velo_info%mft_ratein.eq.0.) then
          if (y.le.(z-kb15i)/kb15s) then                                        ! Point is south of KB1.5
            vyp=utrate
            vzp=0.
          else if (y.le.(z-kb3i)/kb3s) then                                     ! Point is south of KB3
            if (z.gt.mbts*y+mbti) then                                          ! Point is in MBT hanging wall
              vyp=mbtrate
              vzp=mbtrate*tan(mbtdip)
            else                                                                ! Point is in MBT footwall
              vyp=utrate
              vzp=utrate*tan(mbtdip)
            endif
          else if (y.le.(z-kb4i)/kb4s)  then                                    ! Point is south of KB4
            if (z.gt.mht2s*y+mht2i) then                                        ! Point is in MHT2 (MBT) hanging wall
              vyp=(mbtrate)
              vzp=(mbtrate)*tan(mht2dip)
            else                                                                ! Point is in MHT2 footwall
              vyp=utrate
              vzp=utrate*tan(mht2dip)
              !**************************************
              ! DAVE - ADD IN UNDERPLATING HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !**************************************
              ! +vup_lh
            endif
          else if (x.le.stf_seg) then                                           ! Point is west of STF tear
            if (y.le.(z-kb5i)/kb5s) then                                        ! Point is south of KB5
              if (z.gt.stfas*y+stfai) then                                      ! Point is in STFa hanging wall
                vyp=(mbtrate+mctrate+stfrate)
                if (velo_info%stf_ratein.eq.0. .and. velo_info%mct_ratein.eq.0.) then               ! Ignore STF and MCT if inactive
                  vzp=(mbtrate+mctrate+stfrate)*tan(mht3dip)
                else if (velo_info%stf_ratein.eq.0.) then                                 ! Ignore STF if inactive
                  vzp=(mbtrate+mctrate+stfrate)*tan(mctdip)
                else                                                            ! Move parallel to STF
                  vzp=(mbtrate+mctrate+stfrate)*tan(stfdip)
                endif
              else if (z.gt.mcts*y+mcti) then                                   ! Point is in MCT hanging wall
                vyp=(mbtrate+mctrate)
                if (velo_info%mct_ratein.eq.0.) then                                      ! Ignore MCT if inactive
                  vzp=(mbtrate+mctrate)*tan(mht3dip)
                else                                                            ! Move parallel to MCT
                  vzp=(mbtrate+mctrate)*tan(mctdip)
                endif
              else if (z.gt.mht3s*y+mht3i) then                                 ! Point is in MHT3 hanging wall
                vyp=mbtrate
                vzp=mbtrate*tan(mht3dip)
              else                                                              ! Point is in MHT3 footwall
                vyp=utrate
                vzp=utrate*tan(mht3dip)+vup_lh
              endif
            else if (y.le.(z-kb6i)/kb6s) then                                   ! Point is south of KB6
              if (z.gt.stfas*y+stfai) then                                      ! Point is in STFa hanging wall
                vyp=(mbtrate+mctrate+stfrate)
                if (velo_info%stf_ratein.eq.0. .and. velo_info%mct_ratein.eq.0.) then               ! Ignore MCT and STF if inactive
                  vzp=(mbtrate+mctrate+stfrate)*tan(mht4dip)
                else if (velo_info%stf_ratein.eq.0.) then                                 ! Ignore STF if inactive
                  vzp=(mbtrate+mctrate+stfrate)*tan(mctdip)
                else                                                            ! Move parallel to STF
                  vzp=(mbtrate+mctrate+stfrate)*tan(stfdip)
                endif
              else if (z.gt.mcts*y+mcti) then                                   ! Point is in MCT hanging wall
                vyp=(mbtrate+mctrate)
                if (velo_info%mct_ratein.eq.0.) then                                      ! Ignore MCT if inactive
                  vzp=(mbtrate+mctrate)*tan(mht4dip)
                else                                                            ! Move parallel to MCT
                  vzp=(mbtrate+mctrate)*tan(mctdip)
                endif
              else if (z.gt.mht4s*y+mht4i) then                                 ! Point is in MHT4 hanging wall
                vyp=mbtrate
                vzp=mbtrate*tan(mht4dip)+vup_lh
              else                                                              ! Point is in MCT footwall
                vyp=utrate
                vzp=utrate*tan(mht4dip)
              endif
            else if (y.le.(z-kb7ai)/kb7as) then                                 ! Point is south of KB7a
              if (z.gt.stfas*y+stfai) then                                      ! Point is in STFa hanging wall
                vyp=(mbtrate+mctrate+stfrate)
                if (velo_info%stf_ratein.eq.0.) then                                      ! Ignore STF if inactive
                  vzp=(mbtrate+mctrate+stfrate)*tan(mht4dip)
                else                                                            ! Move parallel to STF
                  vzp=(mbtrate+mctrate+stfrate)*tan(stfdip)
                endif
              else if (z.gt.mht4s*y+mht4i) then                                 ! Point is in MHT4 hanging wall
                vyp=(mbtrate+mctrate)
                vzp=(mbtrate+mctrate)*tan(mht4dip)
              else                                                              ! Point is in MHT4 footwall
                vyp=utrate
                vzp=utrate*tan(mht4dip)+vup_gh
              endif
            else                                                                ! Point is north of KB7a
              if (z.gt.mht4s*y+mht4i) then                                      ! Point is in MHT4 hanging wall
                vyp=(mbtrate+mctrate+stfrate)
                vzp=(mbtrate+mctrate+stfrate)*tan(mht4dip)
              else                                                              ! Point is in MHT4 footwall
                vyp=utrate
                vzp=utrate*tan(mht4dip)
              endif
            endif
          else                                                                  ! Point is east of STF tear
            if (y.le.(z-kb5i)/kb5s) then                                        ! Point is south of KB5
              if (z.gt.mcts*y+mcti) then                                        ! Point is in MCT hanging wall
                vyp=(mbtrate+mctrate)
                if (velo_info%mct_ratein.eq.0.) then                                      ! Ignore MCT if inactive
                  vzp=(mbtrate+mctrate)*tan(mht3dip)
                else                                                            ! Move parallel to MCT
                  vzp=(mbtrate+mctrate)*tan(mctdip)
                endif
              else if (z.gt.mht3s*y+mht3i) then                                 ! Point is in MHT3 hanging wall
                vyp=mbtrate
                vzp=mbtrate*tan(mht3dip)
              else                                                              ! Point is in MHT3 footwall
                vyp=utrate
                vzp=utrate*tan(mht3dip)+vup_lh
              endif
            else if (y.le.(z-kb6i)/kb6s) then                                   ! Point is south of KB6
              if (z.gt.stfbs*y+stfbi) then                                      ! Point is in STFb hanging wall
                vyp=(mbtrate+mctrate+stfrate)
                if (velo_info%stf_ratein.eq.0. .and. velo_info%mct_ratein.eq.0.) then               ! Ignore MCT and STF if inactive
                  vzp=(mbtrate+mctrate+stfrate)*tan(mht4dip)
                elseif (stfrate.eq.0.) then                                     ! Ignore STF if inactive
                  vzp=(mbtrate+mctrate+stfrate)*tan(mctdip)
                else                                                            ! Move parallel to STF
                  vzp=(mbtrate+mctrate+stfrate)*tan(stfdip)
                endif
              else if (z.gt.mcts*y+mcti) then                                   ! Point is in MCT hanging wall
                vyp=(mbtrate+mctrate)
                if (velo_info%mct_ratein.eq.0.) then                                      ! Ignore MCT if inactive
                  vzp=(mbtrate+mctrate)*tan(mht4dip)
                else                                                            ! Move parallel to MCT
                  vzp=(mbtrate+mctrate)*tan(mctdip)
                endif
              else if (z.gt.mht4s*y+mht4i) then                                 ! Point is in MHT4 hanging wall
                vyp=mbtrate
                vzp=mbtrate*tan(mht4dip)
              else                                                              ! Point is in MHT4 footwall
                vyp=utrate
                vzp=utrate*tan(mht4dip)+vup_lh
              endif
            else if (y.le.(z-kb7bi)/kb7bs) then                                 ! Point is south of KB7b
              if (z.gt.stfbs*y+stfbi) then                                      ! Point is in STFb hanging wall
                vyp=(mbtrate+mctrate+stfrate)
                if (velo_info%stf_ratein.eq.0.) then                                      ! Ignore STF if inactive
                  vzp=(mbtrate+mctrate+stfrate)*tan(mht4dip)
                else                                                            ! Move parallel to STF
                  vzp=(mbtrate+mctrate+stfrate)*tan(stfdip)
                endif
              else if (z.gt.mht4s*y+mht4i) then                                 ! Point is in MHT4 hanging wall
                vyp=(mbtrate+mctrate)
                vzp=(mbtrate+mctrate)*tan(mht4dip)
              else                                                              ! Point is in MHT4 footwall
                vyp=utrate
                vzp=utrate*tan(mht4dip)+vup_gh
              endif
            else                                                                ! Point is north of KB7b
              if (z.gt.mht4s*y+mht4i) then                                      ! Point is in MHT4 hanging wall
                vyp=(mbtrate+mctrate+stfrate)
                vzp=(mbtrate+mctrate+stfrate)*tan(mht4dip)
              else                                                              ! Point is in MHT4 footwall
                vyp=utrate
                vzp=utrate*tan(mht4dip)
              endif
            endif
          endif

        ! MFT active (perhaps w/ MBT and MCT)
        else
          if (y.le.(z-kb0i)/kb0s) then                                          ! If true, point is in front of thrust belt
            vyp=utrate
            vzp=0.
          else if (y.le.(z-kb1i)/kb1s) then                                     ! Point is south of KB1
            if (z.gt.mfts*y+mfti) then                                          ! Point is in MFT hanging wall
              vyp=mftrate
              vzp=mftrate*tan(mftdip)
            else                                                                ! Point is in MFT footwall
              vyp=utrate
              if (velo_info%mct_ratein.ne.0.) then
                vzp=utrate*tan(mftdip)+(vup_sh+vup_lh)/(SHBA+LHBA)
              else
                vzp=utrate*tan(mftdip)
              endif
            endif
          else if (y.le.(z-kb2i)/kb2s) then                                     ! Point is south of KB2
            if (z.gt.mbts*y+mbti) then                                          ! Point in above MBT
              vyp=(mftrate+mbtrate)
              if (velo_info%mbt_ratein.eq.0.) then                                        ! Ignore MBT if inactive
                vzp=(mftrate+mbtrate)*tan(mht1dip)
              else                                                              ! Move parallel to MBT
                vzp=(mftrate+mbtrate)*tan(mbtdip)
              endif
            else if (z.gt.mht1s*y+mht1i) then                                   ! Point is above MHT1
              vyp=(mftrate)
              vzp=(mftrate)*tan(mht1dip)
            else                                                                ! Point is below MHT1
              vyp=(utrate)
              if (velo_info%mct_ratein.ne.0.) then
                vzp=(utrate)*tan(mht1dip)+(vup_sh+vup_lh)/(SHBA+LHBA)
              else
                vzp=(utrate)*tan(mht1dip)
              endif
            endif
          else if (y.le.(z-kb3i)/kb3s) then                                     ! Point is south of KB3
            vyp=(mftrate+mbtrate)
            if (velo_info%mbt_ratein.eq.0.) then                                          ! Ignore MBT if inactive (Move parallel to MHT2)
              vzp=(mftrate)*tan(mht2dip)
            else                                                                ! Move parallel to MBT
              vzp=(mftrate+mbtrate)*tan(mbtdip)
            endif
          else if (y.le.(z-kb4i)/kb4s) then                                     ! Point is south of KB4
            if (z.gt.mht2s*y+mht2i) then                                        ! Point is in MHT2 (MBT) hanging wall
              vyp=(mftrate+mbtrate)
              vzp=(mftrate+mbtrate)*tan(mht2dip)
            else                                                                ! Point is in MHT2 footwall
              vyp=utrate
              if (velo_info%mct_ratein.ne.0.) then
                vzp=utrate*tan(mht2dip)+(vup_sh+vup_lh)/(SHBA+LHBA)
              else
                vzp=utrate*tan(mht2dip)
              endif
            endif
          else if (x.le.stf_seg) then                                           ! Point is west of STF tear
            if (y.le.(z-kb5i)/kb5s) then                                        ! Point is south of KB5
              if (z.gt.stfas*y+stfai) then                                      ! Point is in STFa hanging wall
                vyp=(mftrate+mbtrate+mctrate+stfrate)
                if (velo_info%stf_ratein.eq.0. .and. velo_info%mct_ratein.eq.0.) then               ! Ignore STF and MCT if inactive
                  vzp=(mftrate+mbtrate+mctrate+stfrate)*tan(mht3dip)
                else if (velo_info%stf_ratein.eq.0.) then                                 ! Ignore STF if inactive
                  vzp=(mftrate+mbtrate+mctrate+stfrate)*tan(mctdip)
                else                                                            ! Move parallel to STF
                  vzp=(mftrate+mbtrate+mctrate+stfrate)*tan(stfdip)
                endif
              else if (z.gt.mcts*y+mcti) then                                   ! Point is in MCT hanging wall
                vyp=(mftrate+mbtrate+mctrate)
                if (velo_info%mct_ratein.eq.0.) then                                      ! Ignore MCT if inactive
                  vzp=(mftrate+mbtrate+mctrate)*tan(mht3dip)
                else                                                            ! Move parallel to MCT
                  vzp=(mftrate+mbtrate+mctrate)*tan(mctdip)
                endif
              else if (z.gt.mht3s*y+mht3i) then                                 ! Point is in MHT3 hanging wall
                vyp=(mftrate+mbtrate)
                vzp=(mftrate+mbtrate)*tan(mht3dip)
              else                                                              ! Point is in MHT3 footwall
                vyp=utrate
                if (velo_info%mct_ratein.ne.0.) then
                  vzp=utrate*tan(mht3dip)+(vup_sh+vup_lh)/(SHBA+LHBA)
                else
                  vzp=utrate*tan(mht3dip)
                endif
              endif
            else if (y.le.(z-kb6i)/kb6s) then                                   ! Point is south of KB6
              if (z.gt.stfas*y+stfai) then                                      ! Point is in STFa hanging wall
                vyp=(mftrate+mbtrate+mctrate+stfrate)
                if (velo_info%stf_ratein.eq.0. .and. velo_info%mct_ratein.eq.0.) then               ! Ignore MCT and STF if inactive
                  vzp=(mftrate+mbtrate+mctrate+stfrate)*tan(mht4dip)
                else if (velo_info%stf_ratein.eq.0.) then                                 ! Ignore STF if inactive
                  vzp=(mftrate+mbtrate+mctrate+stfrate)*tan(mctdip)
                else                                                            ! Move parallel to STF
                  vzp=(mftrate+mbtrate+mctrate+stfrate)*tan(stfdip)
                endif
              else if (z.gt.mcts*y+mcti) then                                   ! Point is in MCT hanging wall
                vyp=(mftrate+mbtrate+mctrate)
                if (velo_info%mct_ratein.eq.0.) then                                      ! Ignore MCT if inactive
                  vzp=(mftrate+mbtrate+mctrate)*tan(mht4dip)
                else                                                            ! Move parallel to MCT
                  vzp=(mftrate+mbtrate+mctrate)*tan(mctdip)
                endif
              else if (z.gt.mht4s*y+mht4i) then                                 ! Point is in MHT4 hanging wall
                vyp=(mftrate+mbtrate)
                vzp=(mftrate+mbtrate)*tan(mht4dip)
              else                                                              ! Point is in MHT4 footwall
                vyp=utrate
                if (velo_info%mct_ratein.ne.0.) then
                  vzp=utrate*tan(mht4dip)+(vup_sh+vup_lh)/(SHBA+LHBA)
                else
                  vzp=utrate*tan(mht4dip)
                endif
              endif
            else if (y.le.(z-kb7ai)/kb7as) then                                 ! Point is south of KB7a
              if (z.gt.stfas*y+stfai) then                                      ! Point is in STFa hanging wall
                vyp=(mftrate+mbtrate+mctrate+stfrate)
                if (velo_info%stf_ratein.eq.0.) then                                      ! Ignore STF if inactive
                  vzp=(mftrate+mbtrate+mctrate+stfrate)*tan(mht4dip)
                else                                                            ! Move parallel to STF
                  vzp=(mftrate+mbtrate+mctrate+stfrate)*tan(stfdip)
                endif
              else if (z.gt.mht4s*y+mht4i) then                                 ! Point is in MHT4 hanging wall
                vyp=(mftrate+mbtrate+mctrate)
                vzp=(mftrate+mbtrate+mctrate)*tan(mht4dip)
              else                                                              ! Point is in MHT4 footwall
                vyp=utrate
                vzp=utrate*tan(mht4dip)+vup_gh
              endif
            else                                                                ! Point is north of KB7a
              if (z.gt.mht4s*y+mht4i) then                                      ! Point is in MHT4 hanging wall
                vyp=(mftrate+mbtrate+mctrate+stfrate)
                vzp=(mftrate+mbtrate+mctrate+stfrate)*tan(mht4dip)
              else                                                              ! Point is in MHT4 footwall
                vyp=utrate
                vzp=utrate*tan(mht4dip)
              endif
            endif
          else                                                                  ! Point is east of STF tear
            if (y.le.(z-kb5i)/kb5s) then                                        ! Point is south of KB5
              if (z.gt.mcts*y+mcti) then                                        ! Point is in MCT hanging wall
                vyp=(mftrate+mbtrate+mctrate)
                if (velo_info%mct_ratein.eq.0.) then                                      ! Ignore MCT if inactive
                  vzp=(mftrate+mbtrate+mctrate)*tan(mht3dip)
                else                                                            ! Move parallel to MCT
                  vzp=(mftrate+mbtrate+mctrate)*tan(mctdip)
                endif
              else if (z.gt.mht3s*y+mht3i) then                                 ! Point is in MHT3 hanging wall
                vyp=(mftrate+mbtrate)
                vzp=(mftrate+mbtrate)*tan(mht3dip)
              else                                                              ! Point is in MHT3 footwall
                vyp=utrate
                if (velo_info%mct_ratein.ne.0.) then
                  vzp=utrate*tan(mht3dip)+(vup_sh+vup_lh)/(SHBA+LHBA)
                else
                  vzp=utrate*tan(mht3dip)
                endif
              endif
            else if (y.le.(z-kb6i)/kb6s) then                                   ! Point is south of KB6
              if (z.gt.stfbs*y+stfbi) then                                      ! Point is in STFb hanging wall
                vyp=(mftrate+mbtrate+mctrate+stfrate)
                if (velo_info%stf_ratein.eq.0. .and. velo_info%mct_ratein.eq.0.) then               ! Ignore MCT and STF if inactive
                  vzp=(mftrate+mbtrate+mctrate+stfrate)*tan(mht4dip)
                else if (velo_info%stf_ratein.eq.0.) then                                 ! Ignore STF if inactive
                  vzp=(mftrate+mbtrate+mctrate+stfrate)*tan(mctdip)
                else                                                            ! Move parallel to STF
                  vzp=(mftrate+mbtrate+mctrate+stfrate)*tan(stfdip)
                endif
              else if (z.gt.mcts*y+mcti) then                                   ! Point is in MCT hanging wall
                vyp=(mftrate+mbtrate+mctrate)
                if (velo_info%mct_ratein.eq.0.) then                                      ! Ignore MCT if inactive
                  vzp=(mftrate+mbtrate+mctrate)*tan(mht4dip)
                else                                                            ! Move parallel to MCT
                  vzp=(mftrate+mbtrate+mctrate)*tan(mctdip)
                endif
              else if (z.gt.mht4s*y+mht4i) then                                 ! Point is in MHT4 hanging wall
                vyp=(mftrate+mbtrate)
                vzp=(mftrate+mbtrate)*tan(mht4dip)
              else                                                              ! Point is in MCT footwall
                vyp=utrate
                if (velo_info%mct_ratein.ne.0.) then
                  vzp=utrate*tan(mht4dip)+(vup_sh+vup_lh)/(SHBA+LHBA)
                else
                  vzp=utrate*tan(mht4dip)
                endif
              endif
            else if (y.le.(z-kb7bi)/kb7bs) then                                 ! Point is south of KB7b
              if (z.gt.stfbs*y+stfbi) then                                      ! Point is in STFb hanging wall
                vyp=(mftrate+mbtrate+mctrate+stfrate)
                if (velo_info%stf_ratein.eq.0.) then                                      ! Ignore STF if inactive
                  vzp=(mftrate+mbtrate+mctrate+stfrate)*tan(mht4dip)
                else                                                            ! Move parallel to STF
                  vzp=(mftrate+mbtrate+mctrate+stfrate)*tan(stfdip)
                endif
              else if (z.gt.mht4s*y+mht4i) then                                 ! Point is in MHT4 hanging wall
                vyp=(mftrate+mbtrate+mctrate)
                vzp=(mftrate+mbtrate+mctrate)*tan(mht4dip)
              else                                                              ! Point is in MHT4 footwall
                vyp=utrate
                vzp=utrate*tan(mht4dip)+vup_gh
              endif
            else                                                                ! Point is north of KB7b
              if (z.gt.mht4s*y+mht4i) then                                      ! Point is in MHT4 hanging wall
                vyp=(mftrate+mbtrate+mctrate+stfrate)
                vzp=(mftrate+mbtrate+mctrate+stfrate)*tan(mht4dip)
              else                                                              ! Point is in MHT4 footwall
                vyp=utrate
                vzp=utrate*tan(mht4dip)
              endif
            endif
          endif
        endif

        vzp = vzp * velo_info%Peclet

        return

! Assume a paraboic uplift field
! Note the following variable substitutions:
! def = maximum rock uplift rate
! dif = minimum rock uplift rate
! The velocity field is dimensionalized in this subroutine.  The wavelength of the parabola
! is hardcoded below (lambda) and will require that the Pecube be recompiled to adjust to other
! wavelengths.
! Edited by S. Olen, 6/18/10.
      else if (velo_info%geoflag == 5) then
    !  Set lateral velocities to zero
        vxp = 0
        vyp = 0

        !  Define the axis of maximum rock uplift; velocity field will be calculated
        !  based on distance from the axis.

        axis = ((velo_info%y2f - velo_info%y1f) / (velo_info%x2f - velo_info%x1f)) * (xin - velo_info%x1f) + velo_info%y1f

        !  Determine the distance of the current node from the axis line

        azimuth = 50 * (3.14159 / 180)

        dx = abs(velo_info%x1f - xin)
        dy = abs(velo_info%y1f - yin)
        h=SQRT(dx**2+dy**2)
        phi1=atan(dx/dy)
        phi2=azimuth-phi1
        d=h*sin(phi2)
        d_corr=(h*sin(phi2))/(2*3.14159)


        !  Calculate the velocity field using a cosine function
        !  Def = Minimum velocity (mmyr-1)
        !  Dif = Maximum velocity (mmyr-1)

        lambda=60
        lambda_corr=lambda/(2*3.14159)
        k=(2*3.14159)/lambda_corr

        vzp=(velo_info%dif-velo_info%def)*(cos(k*d_corr)+1)/2+velo_info%def


        if (d.lt.(-lambda/2)) then
            vzp=velo_info%def
        else if (d.gt.(lambda/2)) then
            vzp=velo_info%def
        else if (vzp.lt.velo_info%def) then
            vzp=velo_info%def
        endif

        vzp = vzp * velo_info%Peclet

        return

!Ellipsoid Uplift Field - M.Schmiddunser 07/11
!
    else if (velo_info%geoflag == 6) then
        x1 = velo_info%x1f   !Coordinates of the ellipse's two foci
        y1 = velo_info%y1f
        x2 = velo_info%x2f
        y2 = velo_info%y2f

        a1 = velo_info%def  !Inner and outer semi major axis
        a2 = velo_info%dif

        !Calculate velocities within the inner and outside the outer ellipse
        !Normalize to 1
        veloin = 1.d0

        if(velo_info%Peclet .eq. 0) then            !Calculate background velocity (veloout) as a fraction of inner velocity (veloin)
            veloout = 0.d0        !Outside velocity cannot be specified if inner velocity is 0
        else                  !(usage of Peclet in other subroutines to calculate actual velocity)
            veloout = velo_info%Peclet2 / velo_info%Peclet
        end if

        !Calculate distances to the two foci
        dist = sqrt((xin-x1)**2+(yin-y1)**2) + &
            sqrt((xin-x2)**2+(yin-y2)**2)

        !Assign vertical velocity, interpolate between inner and outer ellipse
        if (dist .le. 2*a1) then
            vz = veloin
        else if (dist .ge. 2*a2) then
            vz = veloout
        else
            vz = veloin + (dist - 2*a1)/(2*(a2-a1))*(veloout-veloin)
        end if

        !Set lateral velocities to 0
        vxp = 0.d0
        vyp = 0.d0
        vzp = vz * velo_info%Peclet

        return




!Ellipsoid Uplift Field with smoothed edges - M.Schmiddunser 07/11
!
    else if (velo_info%geoflag == 7) then
        x1 = velo_info%x1f   !Coordinates of the ellipse's two foci
        y1 = velo_info%y1f
        x2 = velo_info%x2f
        y2 = velo_info%y2f

        amed = velo_info%dif + velo_info%def  !Sum of inner and outer semi major axis
        adif = (velo_info%dif - velo_info%def) / 3.d0 !Difference of inner and outer semi major axis, division factor determines slope

        !Calculate velocities within the inner and outside the outer ellipse
        !Normalize to 1
        veloin = 1.d0

        if (velo_info%Peclet == 0) then            !Calculate background velocity (veloout) as a fraction of inner velocity (veloin)
            veloout = 0.d0        !Outside velocity cannot be specified if inner velocity is 0
        else                  !(usage of Peclet in other subroutines to calculate actual velocity)
            veloout = velo_info%Peclet2 / velo_info%Peclet
        end if

        !Calculate distances to the two foci
        dist = sqrt((xin-x1)**2+(yin-y1)**2) + &
            sqrt((xin-x2)**2+(yin-y2)**2)

        !Calculate vertical velocity using a Fermi-Distribution
        vz = veloout + (veloin - veloout)*1.05d0 / ( 1 + exp( (dist-amed)/adif ) )


        !Set lateral velocities to 0
        vxp = 0.d0
        vyp = 0.d0
        vzp = vz * velo_info%Peclet

        return

    else if (velo_info%geoflag == 8) then
      ! arbitrary velocity vector field, written by Willi Kappler
      ! input coordiantes: xin, yin, zin
      ! that means we have to provide the velocity of that point
      ! output velocity: vxp, vyp, vzp
      call get_global_velocity(xin, yin, zin, vxp, vyp, vzp)

    ! Final statement if the geometry value specified is not between 1 and 8
    else
      call log_message("Geometry flag value: " + velo_info%geoflag + ", id: " + id)
      call log_message("Geometry flag value is not between 1-8, please correct problem and start again")
      error stop 1
    endif

    return

  end subroutine geometry
end module m_geometry
