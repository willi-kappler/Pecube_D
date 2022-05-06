module m_catchments_output
    use m_compiler

    contains

subroutine catchments_output (age_info,nstep,nsurf,run,xp,yp,det_calc,edot,nrun,&
                              cascadedir,node_thresh,age_flags,header_info,fnme,num_topo_files,&
                              num_of_age_flags)

  ! This subroutine is used when the input model is from Cascade
  ! The Cascade output topo_tec_****.dat files are needed for this subroutine to function
  ! *************************************************************************
  !
  ! These files MUST be in a directory called 'cascade_topo' within the Pecube working directory
  !
  ! *************************************************************************
  ! The x, y, z coordinates are all read from the Cascade topography file along with
  ! each points' catchment number. If it is not the first timestep (i) the subroutine
  ! will attempt to keep the catchment numbers for the current timestep the same as
  ! the previous timestep by checking the points in the catchment and comparing them
  ! to past catchments. Next, ages for the Cascade mesh need to be interpolated from
  ! the Pecube mesh. Finally, if a PDF is wanted the pdfmaker_for_pecube.f90 subroutine
  ! is called with the proper ages and erosion rates (for scaling) and this subroutine
  ! outputs the Cascade mesh with ages

    use m_bivar
    use m_pdfmaker_for_cascade
    use m_data_structures

  implicit none

  integer i,j,k,m,n,p,head,nrun
  integer max_catch,max_catch_temp
  integer md,ncp,nstep,nsurf,nnode,cstore
  integer,dimension(:),allocatable :: catch,catch_count
  real ( kind = 8 ),dimension(:),allocatable :: x,y,z,eros_rate,new_edot
  real ( kind = 8 ),dimension(:,:),allocatable :: new_ages
  real ( kind = 8 ),dimension(:),allocatable :: pdf_ages,pdf_edot
  real ( kind = 8 ),dimension(:,:),allocatable :: pdf_ages_all,pdf_edot_all
  real(8), intent(in) :: edot(nstep, nsurf)
  real ( kind = 8 ) xp(nstep,nsurf),yp(nstep,nsurf),age(nsurf)
  real ( kind = 8 ) junk1,junk2,junk3,junk4,junk5,junk6
  character count4*4,c3*4,junkc*2,run*100,det_calc*300
  logical catch_set
  integer,dimension(:),allocatable :: catch_prev
  logical,dimension(:),allocatable :: lowest_prev
  integer ncascadedir,node_thresh
  character cascadedir*100,as*5,am*5,al*5,eUl*5,eUm*5,eUh*5
  integer(4) :: num_of_age_flags
  integer age_flags(num_of_age_flags),num_topo_files
  real*8 header_info(6)
  character fnme(num_topo_files)*300


  type(age_info_t) age_info(nstep)

  ! may be uninitialized
  max_catch_temp = 0

  ! Writes the various grain sizes and radiation levels to
  ! character arrays printed to ages output Tecplot file
  write (as,'(f5.1)') header_info(1)
  write (am,'(f5.1)') header_info(2)
  write (al,'(f5.1)') header_info(3)
  write (eUl,'(f5.1)') header_info(4)
  write (eUm,'(f5.1)') header_info(5)
  write (eUh,'(f5.1)') header_info(6)

  md=1                          ! See bivar.f90 for explanation of variable
  ncp=4                         ! See bivar.f90 for explanation of variable

  ncascadedir = 1

  do i=1,100
    if (cascadedir(i:i).ne.' ') ncascadedir=i
  enddo

  open (88,file=cascadedir(1:ncascadedir)//'/topo_tec_'//fnme(1)(13:16)//'.dat',status='unknown')

  read (88,*)
  read (88,*)
  read (88,*)
  read (88,*) junkc,nnode

  allocate (catch_prev(nnode),lowest_prev(nnode))

  rewind (88)
  close (88)

  do i=1,nstep
    head=i

    write (count4,'(i4)') head                  ! 'head' is the timestep

    if (head.lt.1000) count4(1:1)='0'
    if (head.lt.100) count4(1:2)='00'
    if (head.lt.10) count4(1:3)='000'

    if (i.le.num_topo_files) then
      open(88,file=cascadedir(1:ncascadedir)//'/topo_tec_'//fnme(i)(13:16)//'.dat',status='unknown')
    else
      open(88,file=cascadedir(1:ncascadedir)//'/topo_tec_'//fnme(num_topo_files)(13:16)//'.dat',status='unknown')
    endif

    read (88,*)
    read (88,*)
    read (88,*)
    read (88,*) junkc,nnode

    allocate (x(nnode),y(nnode),z(nnode))
    allocate (eros_rate(nnode),catch(nnode))
    allocate (new_ages(num_of_age_flags,nnode),new_edot(nnode))

    ! Reads in the coordinates and catchment numbers
    ! Also finds the maximum catchment number for the current timestep
    max_catch = 0
    do j=1,nnode
      read (88,*) x(j),y(j),z(j),junk1,junk2,junk3,junk4,junk5,eros_rate(j),junk6,catch(j)
      if (catch(j).gt.max_catch) max_catch = catch(j)
    enddo

    ! This loop will iterate through all the catchments in the current timestep and compare them to previous catchments
    ! If thea current and previous catchment are found to have the same outlet node but different catchment numbers
    ! the number of the current catchment is changed to the number of the previous catchment so same catchments do not
    ! change numbers over timesteps. If the current catchment's outlet is the not the same as any of the previous catchments
    ! then a new catchment (number) is created for the current timestep and stored for future comparisons
    if (i .eq. 1) max_catch_temp = max_catch                    ! 'max_catch_temp' holds the actual max catchment number from all the previous timesteps
    do j=1,max_catch                                            ! Iterates through all the catchments in the current timestep
      cstore = 0                                                ! Holds the catchment number to be assigned to the current point
      catch_set = .false.                                       ! Tells if there are any values currently in the catchment
      do n=1,nnode                                              ! Iterates through all Cascade nodes
        if (i .eq. 1) then                                      ! If it is the first timestep
          if (catch(n).eq.j) then                               ! If the catchment number of the node is equal to the current catchment number
            if (.not.catch_set) then                            ! If the catchment does not contain any nodes yet
              lowest_prev(n) = .true.                           ! Sets the current node to the outlet of the current catchment
              catch_set=.true.                                  ! Sets the catchment as not empty
            endif
          endif
        else                                                    ! If it is not the first timestep
          if (catch(n).eq.j) then                               ! If the catchment number of the node is equal to the current catchment number
            if (.not.catch_set) then                            ! If the current catchment is empty
              if (lowest_prev(n)) then                          ! If the node was the outlet of catchment in the previous timestep
                catch_set = .true.                              ! Set the catchment as not empty
                cstore = catch_prev(n)                          ! Set 'cstore' to the previous catchment number of the node
                catch(n) = cstore                               ! Set the current catchment number of the node to the previous number
              else                                              ! If the current node was not the outlet of the catchment in the previous timestep
                catch_set = .true.                              ! Set the catchment is not empty
                max_catch_temp = max_catch_temp + 1             ! Increments the actual amount of catchments through all timesteps (this is when an actual new catchment is created, not just a labeling error)
                cstore = max_catch_temp                         ! Set 'cstore' to the new catchment number
                catch(n) = cstore                               ! Set the node's catchment number to the new catchment number
                lowest_prev(n) = .true.                         ! Set the node as the outlet of the new catchment
              endif
            else                                                ! If the catchment is not empty
              catch(n) = cstore                                 ! Set the node's catchment number to the catchment number in 'cstore'
            endif
          endif
        endif
      enddo
    enddo

    max_catch = max_catch_temp                                  ! Stores the actual maximum catchment number into 'max_catch'

    do j=1,nnode
      catch_prev(j)=catch(j)                                    ! Sets the current catchment numbers as the previous numbers for use in the next timestep
    enddo

    allocate (catch_count(max_catch))
    catch_count = 0
    do n=1,nnode
      catch_count(catch(n)) = catch_count(catch(n)) + 1         ! 'catch_count' holds the number of nodes in the specified catchment
    enddo

    ! The 'age' array is used to cast the real*4 ages to real*8 before passing into the interpolation routine
    ! All calls to idbvip are to interpolate a set of ages from the Pecube mesh onto the Cascade mesh
    do p=1,num_of_age_flags
      age = age_info(i)%all_ages(:,p)
      call idbvip (nsurf,xp(i,:),yp(i,:),age,nnode,x,y,new_ages(p,:))
    enddo

    ! This interpolation is for the erosion rates. Used in scaling the PDF later.
    age = edot(i,:)
    call idbvip (nsurf,xp(i,:),yp(i,:),age,nnode,x,y,new_edot)

    call sys_command('if [ ! -d '//run(1:nrun)//'/pdf_cascade ]; then mkdir '//run(1:nrun)//'/pdf_cascade; fi')     ! Checks if the pdf_pecube folder exists in the associated run output directory; if not, then it creates it

    open (101,file=run(1:nrun)//'/pdf_cascade/Timestep_'//count4//'_Catchments_data.dat',status="unknown")
    write (101,'(A40)') 'Catchment Number       Number of nodes'

    do k=1,max_catch                                                    ! Iterates through all catchments
      if (catch_count(k).ge.node_thresh) then                           ! If the catchment is not empty

        if (det_calc .eq. '1') then                                     ! If a PDF is wanted
          allocate (pdf_ages(catch_count(k)),pdf_edot(catch_count(k)))
          allocate (pdf_ages_all(num_of_age_flags,catch_count(k)),pdf_edot_all(num_of_age_flags,catch_count(k)))
          pdf_ages_all=0.
          pdf_edot_all=0.
          do p=1,num_of_age_flags                                                     ! Iterates through all types of ages
            if (age_flags(p) .eq. 1) then
              n=1
              do j=1,nnode                                                ! Iterates through all Cascade nodes
                if (catch(j).eq.k) then                                   ! If the node is in the current catchment
                  pdf_ages(n) = new_ages(p,j)                             ! Compile array of ages of nodes in the catchment
                  pdf_edot(n) = new_edot(j)                               ! Compile array of erosion rates
                  n = n + 1                                               ! Counter variable
                endif
              enddo
              ! Subroutine to calculate and output PDF for specified age type (p) for current catchment (k) at current timestep (i)
!               call pdfmaker_for_pecube (pdf_ages,catch_count(k),p,i,k,run,pdf_edot,'','',-1,-1,nrun)

              pdf_ages_all(p,:)=pdf_ages
              pdf_edot_all(p,:)=pdf_edot

            endif
          enddo

          call pdfmaker_for_cascade (pdf_ages_all,catch_count(k),i,k,run,pdf_edot_all,&
                                     nrun,age_flags,header_info,num_of_age_flags)

          write (101,'(i10,i20)') k,catch_count(k)

        endif

        write (c3,'(i4)') k

        if (k.lt.1000) c3(1:1)='0'
        if (k.lt.100) c3(1:2)='00'
        if (k.lt.10) c3(1:3)='000'

        ! Checks if the 'catchments' folder exists in the output directory and creates it if it doesn't exist
        call sys_command('if [ ! -d '//run(1:nrun)//'/catchments ]; then mkdir '//run(1:nrun)//'/catchments; fi')

        open(100,file=run(1:nrun)//'/catchments/Timestep_'//count4//'_Catchment_'//c3//'.dat',status='unknown')

        write (100,*) 'TITLE = "Pecube Catchment Output"'
        write (100,'(A59)',ADVANCE="no") 'VARIABLES = "X (km)" "Y (km)" "Z (km)" "Total Erosion Rate"'
        if (age_flags(1).eq.1) write (100,'(A)',ADVANCE="no") ' "AHe Age - Farley, 2000 (Ma)"'
        if (age_flags(2).eq.1) write (100,'(A)',ADVANCE="no") ' "AHe Age - a='//as//' um (Ma)"'
        if (age_flags(3).eq.1) write (100,'(A)',ADVANCE="no") ' "AHe Age - a='//am//' um (Ma)"'
        if (age_flags(4).eq.1) write (100,'(A)',ADVANCE="no") ' "AHe Age - a='//al//' um (Ma)"'
        if (age_flags(5).eq.1) write (100,'(A)',ADVANCE="no") ' "AHe Age - eU='//eUl//' ppm (Ma)"'
        if (age_flags(6).eq.1) write (100,'(A)',ADVANCE="no") ' "AHe Age - eU='//eUm//' ppm (Ma)"'
        if (age_flags(7).eq.1) write (100,'(A)',ADVANCE="no") ' "AHe Age - eU='//eUh//' ppm (Ma)"'
        if (age_flags(9).eq.1) write (100,'(A)',ADVANCE="no") ' "AFT Age (Ma)"'
        if (age_flags(12).eq.1) write (100,'(A)',ADVANCE="no") ' "KFeldAr Age (Ma)"'
        if (age_flags(8).eq.1) write (100,'(A)',ADVANCE="no") ' "ZHe Age (Ma)"'
        if (age_flags(10).eq.1) write (100,'(A)',ADVANCE="no") ' "ZFT Age (Ma)"'
        if (age_flags(14).eq.1) write (100,'(A)',ADVANCE="no") ' "MuscAr Age (Ma)"'
        if (age_flags(13).eq.1) write (100,'(A)',ADVANCE="no") ' "BioAr Age (Ma)"'
        if (age_flags(11).eq.1) write (100,'(A)',ADVANCE="no") ' "MAr Age (Ma)"'
        if (age_flags(15).eq.1) write (100,'(A)',ADVANCE="no") ' "HornAr Age (Ma)"'
        if (age_flags(16).eq.1) write (100,'(A)',ADVANCE="no") ' "APb (Ma)"'
        if (age_flags(17).eq.1) write (100,'(A)',ADVANCE="no") ' "Biotite (Ma)"'
        if (age_flags(18).eq.1) write (100,'(A)',ADVANCE="no") ' "RUPb (Ma)"'
        if (age_flags(19).eq.1) write (100,'(A)',ADVANCE="no") ' "TUPb (Ma)"'
        if (age_flags(20).eq.1) write (100,'(A)',ADVANCE="no") ' "ZUPb (Ma)"'
        if (age_flags(21).eq.1) write (100,'(A)',ADVANCE="no") ' "TUTh/He (Ma)"'

        write (100,*) 'ZONE T = "Pecube"'
        write (100,*) 'I=',catch_count(k),', J=1, K=1, ZONETYPE=Ordered'
        write (100,*) 'DATAPACKING=POINT'
        write (100,'(A)',ADVANCE="no") ' DT=(DOUBLE DOUBLE DOUBLE DOUBLE'
        do j=1,num_of_age_flags
          if(age_flags(j) .eq. 1) write (100,'(A)',ADVANCE="no") ' DOUBLE'
        enddo
        write (100,'(A)') ')'

        do m=1,nnode
          if (catch(m).eq.k) then
            write (100,'(4f12.4)',ADVANCE="no") x(m),y(m),z(m),new_edot(m)

            do j=1,7
              if (age_flags(j).eq.1) then
                write (100,'(f12.4)',ADVANCE="no") new_ages(j,m)
              endif
            enddo
            if (age_flags(9).eq.1) write (100,'(f12.4)',ADVANCE="no") new_ages(9,m)
            if (age_flags(12).eq.1) write (100,'(f12.4)',ADVANCE="no") new_ages(12,m)
            if (age_flags(8).eq.1) write (100,'(f12.4)',ADVANCE="no") new_ages(8,m)
            if (age_flags(10).eq.1) write (100,'(f12.4)',ADVANCE="no") new_ages(10,m)
            if (age_flags(14).eq.1) write (100,'(f12.4)',ADVANCE="no") new_ages(14,m)
            if (age_flags(13).eq.1) write (100,'(f12.4)',ADVANCE="no") new_ages(13,m)
            if (age_flags(11).eq.1) write (100,'(f12.4)',ADVANCE="no") new_ages(11,m)
            if (age_flags(15).eq.1) write (100,'(f12.4)',ADVANCE="no") new_ages(15,m)
            write(100,*)

          endif
        enddo
        close (100)
        deallocate (pdf_ages,pdf_edot)
        deallocate (pdf_ages_all,pdf_edot_all)
      endif
    enddo

    close (101)
    close (88)
    deallocate (x,y,z,eros_rate,catch,catch_count)
    deallocate (new_ages,new_edot)

  enddo

  return
end subroutine catchments_output

end module m_catchments_output
