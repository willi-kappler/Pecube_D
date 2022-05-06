module m_detrital_mc
    use m_compiler
    use m_logger

    contains

! detrital_mc.f90
!
! This fortran code uses a monte carlo simulation method to randomly grab n
! samples from a distribution of model predicted cooling ages, assign a
! designated uncertainty to the samples, generate a PDF of the model predicted
! cooling ages, compare that PDF to an observed age PDF using a Kuiper test,
! and finally record the results of that test.  This process is repeated a
! large number of times (~10000).
!
! dwhipp - 04/08

!       program detrital_mc
      subroutine detrital_mc(obasin,olc,oage,oageu,&
                             page,perate,plc,run,x_basin,y_basin,nrun)

        use m_init_random_seed
        use m_kptwo

! Variable declaration
      real*8,dimension(:),allocatable :: on,opsum,op,opdf,opdfv
      real*8,dimension(:),allocatable :: pagesc,pageusc
      real*8,dimension(:),allocatable :: pagemc,pageumc,pnmc,ppsummc,ppmc,ppdfmc
      real*8,dimension(:),allocatable :: ppdfvmc
      integer,dimension(:),allocatable :: kuiper_res
      real*8 :: d,prob,pagemu,pagemed,pagesd,pdfvsc,mc_iterf,jf,kpct
      real*8 :: dx,osum,randflt,psummc,pi,peratemin,peratescl
      integer :: olc,onum,oamin,oamax,h,i,j,k,l,mc_iter
      integer :: plc,plcsc,cnt,rint,paminmc,pamaxmc,pnummc,cnt2,cnt3,hm,hm2,cnt4
      character :: obasin*300,jc*3,hc*5,run*100
      real*8 :: oage(olc),oageu(olc),age_temp,oageu_percent(olc)
      real*8 :: page(plc),perate(plc),pageu(plc),peratesc(plc)
      real*8 :: x_basin,y_basin
      character :: x_basin_char*20,y_basin_char*20
      integer :: num_xchars,num_ychars,nrun

      ! may be uninitialized
      num_xchars = 0
      num_ychars = 0

      call log_message('Performing Monte Carlo runs')
      call log_message('This may take awhile')
      call log_message('Every 500th iteration will be announced')

! Initialize random number generator
      call init_random_seed()
      !call random_seed()

! Variable initialization/declaration
      pi=atan(1.)*4.                                                            ! Define pi
!       basnum=21                                                                 ! Number of basins to analyze
      mc_iter=10000                                                             ! Number of iterations in the Monte Carlo simulation
      pdfvsc=50.                                                                 ! Approximate number of values in scaled PDF vectors

!       open (82,file='inyobasin_agedist_for_Pecube.dat',status='unknown')
!       do i=1,plc
!         read (82,*) page(i)
!       enddo
!       close (82)

      do j=1,20
        x_basin_char(j:j)=' '
        y_basin_char(j:j)=' '
      enddo
      write (x_basin_char,'(f20.2)') x_basin
      write (y_basin_char,'(f20.2)') y_basin
      do j=1,20
        if (x_basin_char(j:j).ne.'0' .and. x_basin_char(j:j).ne.'1' .and. x_basin_char(j:j).ne.'2' &
            .and. x_basin_char(j:j).ne.'3' .and. x_basin_char(j:j).ne.'4' .and. x_basin_char(j:j).ne.'5' &
            .and. x_basin_char(j:j).ne.'6' .and. x_basin_char(j:j).ne.'7' .and. x_basin_char(j:j).ne.'8' &
            .and. x_basin_char(j:j).ne.'9' .and. x_basin_char(j:j).ne.'.') num_xchars=j+1
        if (y_basin_char(j:j).ne.'0' .and. y_basin_char(j:j).ne.'1' .and. y_basin_char(j:j).ne.'2' &
            .and. y_basin_char(j:j).ne.'3' .and. y_basin_char(j:j).ne.'4' .and. y_basin_char(j:j).ne.'5' &
            .and. y_basin_char(j:j).ne.'6' .and. y_basin_char(j:j).ne.'7' .and. y_basin_char(j:j).ne.'8' &
            .and. y_basin_char(j:j).ne.'9' .and. y_basin_char(j:j).ne.'.') num_ychars=j+1
      enddo


! Read in basin info
!       open(10,file='basin_summary_info.txt',status='old')
      call sys_command('if [ ! -d '//run(1:nrun)//'/mc_pdf_output ]; then mkdir '//run(1:nrun)//'/mc_pdf_output; fi')
      call sys_command('if [ ! -d '//run(1:nrun)//'/mc_pdf_output/Basin_X_'//x_basin_char(num_xchars:20)//'_Y_'//y_basin_char(num_ychars:20)//' ]; then mkdir '//run(1:nrun)//'/mc_pdf_output/Basin_X'//x_basin_char(num_xchars:20)//'_Y'//y_basin_char(num_ychars:20)//'; fi')
      call sys_command('if [ ! -d '//run(1:nrun)//'/mc_pdf_output/Basin_X_'//x_basin_char(num_xchars:20)//'_Y_'//y_basin_char(num_ychars:20)//'/pass_mc_age_PDF ]; then mkdir '//run(1:nrun)//'/mc_pdf_output/Basin_X'//x_basin_char(num_xchars:20)//'_Y'//y_basin_char(num_ychars:20)//'/pass_mc_age_PDF; fi')
      call sys_command('if [ ! -d '//run(1:nrun)//'/mc_pdf_output/Basin_X_'//x_basin_char(num_xchars:20)//'_Y_'//y_basin_char(num_ychars:20)//'/mc_age_PDF ]; then mkdir '//run(1:nrun)//'/mc_pdf_output/Basin_X'//x_basin_char(num_xchars:20)//'_Y'//y_basin_char(num_ychars:20)//'/mc_age_PDF; fi')

      open(21,file=run(1:nrun)//'/mc_pdf_output/Basin_X'//x_basin_char(num_xchars:20)//'_Y'//y_basin_char(num_ychars:20)//'/kmc_percent_pass_summary.dat',status='unknown')             ! Open basin summary results file
      write(21,'(a37)') '                      Percent passing'                           ! Write header
      write(21,'(a35)') 'Basin                   Kuiper test'


      oageu_percent=oageu/oage

      pagemu=sum(oageu_percent)/olc

      do i=1,olc
        sumsd = real(sumsd + (oageu_percent(i)-pagemu)**2.0)
      enddo
      pagesd = sqrt(sumsd/(olc-1))

      do i=1,olc-1
        do j=i+1,olc
          if (oageu_percent(j).lt.oageu_percent(i)) then
            age_temp=oageu_percent(i)
            oageu_percent(i)=oageu_percent(j)
            oageu_percent(j)=age_temp
          endif
        enddo
      enddo

      if (mod(olc,2).eq.0) then
        pagemed=(oageu_percent(olc/2)+oageu_percent(olc/2+1))/2.
      else
        pagemed=page(olc/2+1)
      endif


! Read in observed cooling ages and generate PDF
          do j=1,plc
            !pageu(j)=page(j)*(pagemu/100.)                                       ! Assign uncertainty (using mean measured value for dataset)
!             pageu(j)=page(j)*.11
            pageu(j)=page(j)*(pagemed/100.)                                       ! Assign uncertainty (using median measured value for dataset)
          enddo

!           do j=1,plc
            !peratesc(j)=nint(perate(j)*100.)                                      ! Generate scaling factors by multiplying erosion rates by 100 and converting to integers
!           enddo

          peratemin=minval(perate)                                              ! Determine minimum erosion rate in model domain

          if (peratemin.ne.0) then
            peratescl=10./peratemin                                               ! Set scaling value to ensure at least 10 occurances of min rate ages
          else
            peratescl=100.
          endif

          peratesc=nint(perate*peratescl)                                       ! Generate scaling factors by multiplying erosion rates by 100 and converting to integers
! Generate data PDF
          dx=0.001                                                                ! Specify x spacing for data PDF generation
          oamin=nint(minval(oage)-10*maxval(oageu))                                ! Find range of ages + uncertainties
          oamax=nint(maxval(oage)+10*maxval(oageu))
          onum = nint((oamax-oamin)/dx)                                                   ! Find number of values in PDF arrays
          !do j=1,olc                                                              ! Scale uncertainties using optimal scaling factor of 0.6
          !  oageu(j)=oageu(j)*0.6                                                 ! See Brandon, M., Probability Density Plot for Fission-Track Grain-Age Samples,
          !enddo                                                                   ! Radiation Measurements, Vol. 26, No. 5, pp. 663-676, 1996
          allocate(on(onum+1),opsum(onum+1),opdf(onum+1),op(onum+1))              ! Allocate data PDF arrays
          do j=1,onum+1
            on(j)=oamin+(j-1)*dx                                                  ! Fill age range array
          enddo
          opsum=0.
          do j=1,olc
            do k=1,onum+1
              op(k)=(1./(oageu(j)*sqrt(2.*pi)))*exp(-0.5*((on(k)-oage(j))/&       ! Fill probability array
                    (oageu(j)))**2.)
              opsum(k)=opsum(k)+op(k)                                             ! Fill sum array to check area under array curve
            enddo
          enddo

          osum=0.
          do j=1,onum+1
            !opdf(j)=(opsum(j)/olc)*dx                                             ! Scale PDF array to normalize area under PDF curve
            opdf(j)=(opsum(j)/olc)                                                ! Scale PDF array to normalize area under PDF curve
            osum=osum+opdf(j)                                                     ! Calculate area under curve
          enddo

! Generate data PDF vector
          cnt2=0
          do j=1,onum+1
            hm=nint(pdfvsc*opdf(j))                                                ! Set number of occurances of given age at current probability
            do k=1,hm
              cnt2=cnt2+1                                                         ! Count total number of ages in PDF vector array for allocation below
            enddo
          enddo
          allocate(opdfv(cnt2))                                                    ! Allocate PDF vector array
          cnt2=0
          opdfv=0.
          do j=1,onum+1
            hm=nint(pdfvsc*opdf(j))
            do k=1,hm
              cnt2=cnt2+1
              opdfv(cnt2)=on(j)                                                    ! Fill PDF vector array scaling number of ages by the probability they occur at given age
            enddo
          enddo

! Generate erosion-rate-scaled model age distribution
          plcsc = nint(sum(peratesc))                                                     ! Determine number of ages in scaled model age distribution
          allocate(pagesc(plcsc),pageusc(plcsc))                                  ! Allocate scaled age/uncertainty arrays
!           allocate(pagesc(plc),pageusc(plc))
          cnt=0
          do j=1,plc
            do k=1,int(peratesc(j))
              cnt=cnt+1
              pagesc(cnt)=page(j)                                                 ! Fill scaled age array
              pageusc(cnt)=pageu(j)                                               ! Fill scaled age uncertainty array
            enddo
          enddo

!!!
! Start monte carlo runs for testing model/data fits with Kuiper test
!!!!
          allocate(kuiper_res(mc_iter),pagemc(olc),pageumc(olc))                  ! Allocate kuiper test array and arrays of size olc for random subsampling of the true
          !allocate(pmc(int(pdfvsc),2*mc_iter))
          !pmc=0.
          cnt4=0
          do j=1,mc_iter                                                          ! Model age distribution; This loop runs mc_iter times (usually ~10000)

            if (j.eq.1 .or. mod(j,500).eq.0) call log_message('Iteration ' + j)

            !call init_random_seed()
            jf=real(j)

! Randomly grab olc (n) grains from model distribution
            do k=1,olc
              call random_number(randflt)                                         ! Generate random number [0,1)
              rint=int(randflt*(cnt))+1                                           ! Get random integer value within range of size of scaled age dist.
              !rint=int(randflt*(cnt2))+1                                           ! Get random integer value within range of size of scaled age dist.
              pagemc(k)=pagesc(rint)                                              ! Add random age to monte carlo age array
              pageumc(k)=pageusc(rint)                                            ! Add associated uncertainty to monte carlo uncertainty array
              !pagemc(k)=opdfv(rint)
              !pageumc(k)=opdfv(rint)*pagemu/100
            enddo

! Generate MC model PDF using olc grains
            paminmc = nint(minval(pagemc)-10*maxval(pageumc))                        ! Find range of ages + uncertainties
            pamaxmc = nint(maxval(pagemc)+10*maxval(pageumc))
            pnummc = nint((pamaxmc-paminmc)/dx)                                           ! Find number of values in MC model PDF arrays
!             do k=1,olc                                                            ! Scale uncertainties using optimal scaling factor of 0.6
!               pageumc(k)=pageumc(k)*0.6                                           ! See Brandon, M., Probability Density Plot for Fission-Track Grain-Age Samples,
!             enddo                                                                 ! Radiation Measurements, Vol. 26, No. 5, pp. 663-676, 1996
            allocate(pnmc(pnummc+1),ppsummc(pnummc+1),ppdfmc(pnummc+1))
            allocate(ppmc(pnummc+1))
            do k=1,pnummc+1
              pnmc(k)=paminmc+(k-1)*dx                                            ! Fill age range array
              !pmc(k,2*jf-1)=paminmc+(k-1)*dx                                       ! Fill age range array
            enddo
            ppsummc=0.
            do k=1,olc
              do l=1,pnummc+1
                ppmc(l)=(1./(pageumc(k)*sqrt(2.*pi)))*exp(-0.5*((pnmc(l)-&        ! Fill probability array
                        pagemc(k))/(pageumc(k)))**2.)
                ppsummc(l)=ppsummc(l)+ppmc(l)                                     ! Fill sum array to check area under array curve
              enddo
            enddo
            psummc=0.
            do k=1,pnummc+1
              !ppdfmc(k)=(ppsummc(k)/olc)*dx                                       ! Scale PDF array to normalize area under PDF curve
              ppdfmc(k)=(ppsummc(k)/olc)                                          ! Scale PDF array to normalize area under PDF curve
              !pmc(k,2*jf)=(ppsummc(k)/olc)*dx                                      ! Fill age range array
              psummc=psummc+ppdfmc(k)                                             ! Calculate area under curve
            enddo

! Generate MC model PDF vector
            cnt3=0
            do k=1,pnummc+1
              hm2=nint(pdfvsc*ppdfmc(k))                                           ! Set number of occurances of given age at current probability
              do l=1,hm2
                cnt3=cnt3+1                                                       ! Count total number of ages in PDF vector array for allocation below
              enddo
            enddo
            allocate(ppdfvmc(cnt3))                                                ! Allocate PDF vector array
            cnt3=0
            ppdfvmc=0.
            do k=1,pnummc+1
              hm2=nint(pdfvsc*ppdfmc(k))
              do l=1,hm2
                cnt3=cnt3+1
                ppdfvmc(cnt3)=pnmc(k)                                              ! Fill PDF vector array scaling number of ages by the probability they occur at given age
              enddo
            enddo

! Run Kuiper test to get misfit between data and model
            call kptwo(opdfv,cnt2,ppdfvmc,cnt3,olc,d,prob,h)
            kuiper_res(j)=h                                                       ! Store kuiper test result (0=pass;1=fail) for this iteration in kuiper results array
            if (h.eq.0) cnt4=cnt4+1                                               ! Increment counter for number of models that pass Kuiper test

! Write out data PDF
!             if (j.eq.1) then
!               open(22,file=run//'/mc_pdf_output/data_age_PDF_'//obasin,&
!                   status='unknown')
!               do k=1,onum+1
!                 write(22,*) on(k),opdf(k)
!               enddo
!               close(22)
!             endif

! Write out 500 PDFs that pass the Kuiper test
            if (h.eq.0 .and. cnt4.le.500) then
              write(hc,'(i5)') j
              if (j.lt.10) hc(1:4)='0000'
              if (j.lt.100) hc(1:3)='000'
              if (j.lt.1000) hc(1:2)='00'
              if (j.lt.10000) hc(1:1)='0'
              open(23,file=run(1:nrun)//'/mc_pdf_output/Basin_X'//x_basin_char(num_xchars:20)//'_Y'//y_basin_char(num_ychars:20)//'/pass_mc_age_PDF/pass_mc_age_PDF_'//hc//'.dat',&
                  status='unknown')
              do k=1,pnummc+1
                write(23,*) pnmc(k),ppdfmc(k)
              enddo
              close(23)
            endif

! Write out first 500 monte carlo PDFs
            if (j.le.500) then
              write(jc,'(i3)') j
              if (j.lt.10) jc(1:2)='00'
              if (j.lt.100) jc(1:1)='0'
              open(25,file=run(1:nrun)//'/mc_pdf_output/Basin_X'//x_basin_char(num_xchars:20)//'_Y'//y_basin_char(num_ychars:20)//'/mc_age_PDF/mc_age_PDF_'//jc//'.dat',&
                  status='unknown')
              do k=1,pnummc+1
                write(25,*) pnmc(k),ppdfmc(k)
              enddo
              close(25)
            endif

! Deallocate arrays
            deallocate(pnmc,ppsummc,ppdfmc,ppmc,ppdfvmc)                           ! Deallocate arrays reallocated during monte carlo sim
          enddo

          mc_iterf=real(mc_iter)
          kpct=(1-(sum(kuiper_res)/mc_iterf))*100.                             ! Store percent of models that passed kuiper test for given basin

! Write output files
          open(20,file=run(1:nrun)//'/mc_pdf_output/Basin_X'//x_basin_char(num_xchars:20)//'_Y'//y_basin_char(num_ychars:20)//'/kuiper_mc_results.dat',status='unknown')
          do j=1,mc_iter
            write(20,*) kuiper_res(j)                                             ! Write out individual subset kuiper test result values to file for this basin
          enddo
          close(20)

          write(21,'(a25,f10.5)') obasin,kpct                                   ! Write out the summary percent of models that passed the Kuiper test

          !open(24,file='age_pdf_output/all_PDF_results_'//obasin,status='unknown')
          !write(24,*) pmc                                                         ! Write out individual subset kuiper test result values to file for this basin
          !close(24)

! Deallocate arrays
          deallocate(on,opdf,op,opsum)
          deallocate(pagesc,pageusc)
          deallocate(pagemc,pageumc,kuiper_res)
          deallocate(opdfv)!,pmc)

! End of main loop
!         endif
!       enddo

! Close open files
      close(21)

! Exit
      end subroutine detrital_mc
end module m_detrital_mc

