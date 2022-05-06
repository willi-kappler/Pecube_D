module m_pdfmaker_for_cascade
contains

subroutine pdfmaker_for_cascade (age,number,t,c,run,edot,&
                                nrun,age_flags,header_info,&
                                num_of_age_flags)

  ! % April 26, 2005  Greg Stock
  ! % Modified by Dave Whipp - 08/07
  ! % Converted to FORTRAN by cspath 03/08
  ! %
  ! % Calculates synoptic probability density function (SPDF) to create basin-wide PDFs 
  ! % from individual detrital apatite He ages; based on Equation 2 from Ruhl and Hodges (2005).
  ! % Input is a two-column Excel (xls) file of measured He ages and one sigma uncertainties.
  ! % Assumes Gaussian PDF distribution for each grain age, then sums all grain PDF's to create  
  ! % basin-wide PDF.
  !
  ! Input ages are passed as argument to this subroutine
  ! The ages are either found using the output catchment file from Cascade (if model is from Cascade)
  ! or they are found using the TAPES_G algorithm for a user specified basin outlet
  ! The subroutines that call this subroutine are catchments_output.f90 and pdfmaker_for_data2.f, respectively
  ! The error is a set percentage depending on the type of age being used

  implicit none

  real*8 dx,pi,i,start,finish,avg_edot
  integer number,counter,k,t,c,j
  real*8,dimension(:,:),allocatable:: P,Psum,normPDF
  real*8,dimension(:),allocatable:: n
  integer(4) :: num_of_age_flags
  real*8 age(num_of_age_flags,number),error(num_of_age_flags,number),per_error(num_of_age_flags),edot(num_of_age_flags,number)
  character timestep*4,catchment*4,run*100
  integer nrun,age_flags(num_of_age_flags)
  character as*5,am*5,al*5,eUl*5,eUm*5,eUh*5
  real*8 header_info(6)

  per_error(1)=0.05             !% AHe 5% error
  per_error(2)=0.05             !% AHe 5% error
  per_error(3)=0.05             !% AHe 5% error
  per_error(4)=0.05             !% AHe 5% error
  per_error(5)=0.05             !% AHe 5% error
  per_error(6)=0.05             !% AHe 5% error
  per_error(7)=0.05             !% AHe 5% error
  per_error(9)=0.1              !% AFT 10% error
  per_error(8)=0.05             !% ZHe 5% error
  per_error(10)=0.1             !% ZFT 10% error
  per_error(11)=0.04            !% MAr 4% error

  ! TODO: add error for other ages
  per_error(12)=0.1             !
  per_error(13)=0.1             !
  per_error(14)=0.1             !
  per_error(15)=0.1             !

  dx = 0.01
  pi = 3.14159265

  ! Calculates absolute error of the age data based on age type
  do k=1,num_of_age_flags
    do j=1,number
      if (age_flags(k) .eq. 1) then
        error(k,j) = age(k,j) * per_error(k)
      else
        error(k,j)=0.
      endif
    enddo
  enddo

  start = minval(age)-10*maxval(error)-mod(minval(age)-10*maxval(error),dx)           ! Finds minimum possible value
  finish = maxval(age)+10*maxval(error)-mod(maxval(age)+10*maxval(error),dx)        ! Finds maximum possible value

!   start = nint(minval(age(:))-2*maxval(error(:)))
!   finish = nint(maxval(age(:))+10*maxval(error(:)))

  !% Range of n values to consider, rounded to whole numbers
  counter = 0
  i = start
  do while (i .lt. finish)
    counter = counter + 1
    i = i + dx
  enddo

  allocate (n(counter),P(num_of_age_flags,counter),Psum(num_of_age_flags,counter),&
            normPDF(num_of_age_flags,counter))

  counter = 0
  i = start
  do while (i .lt. finish)
    counter = counter + 1
    n(counter) = i
    i = i + dx
  enddo

  do j=1,num_of_age_flags
    Psum(j,:) = 0.
    do k=1,number
      !% Calculate probability for range of n values from Eq. 2a of Ruhl and
      !% Hodges (2005)
      P(j,:)=(1./(error(j,k)*sqrt(2*pi)))*exp(-0.5*((n-age(j,k))/(error(j,k)))**2)*edot(j,k)
      !% Sum the resulting PDF with previous PDF's from Eq. 2b of Ruhl and
      !% Hodges (2005)
      Psum(j,:)=Psum(j,:)+P(j,:)
    enddo

  ! Calculates the average erosion rate
  ! Used to normalize area under curve to 1
  avg_edot=0.
  do k=1,number
    avg_edot = avg_edot + edot(j,k)
  enddo
  avg_edot = avg_edot/number

  !% Normalize area under PDF curve to 1
  normPDF(j,:)=(Psum(j,:)/number)/avg_edot       !% Normalize PDF so area under curve = 1

  enddo

      write (timestep,'(i4)') t
      if (t.lt.1000) timestep(1:1)='0'
      if (t.lt.100) timestep(1:2)='00'
      if (t.lt.10) timestep(1:3)='000'

    ! c will be -1 when this function is called from pdfmaker_for_data2.f because c refers to a catchment number
    ! and there is only catchment numbers when using Cascade models
    ! No significance to using c to determine what to output but simpler than adding another flag to be passed in

      write (catchment,'(i4)') c
      if (c.lt.1000) catchment(1:1)='0'
      if (c.lt.100) catchment(1:2)='00'
      if (c.lt.10) catchment(1:3)='000'

      write (as,'(f5.1)') header_info(1)
      write (am,'(f5.1)') header_info(2)
      write (al,'(f5.1)') header_info(3)
      write (eUl,'(f5.1)') header_info(4)
      write (eUm,'(f5.1)') header_info(5)
      write (eUh,'(f5.1)') header_info(6)

!       write (age_flag,'(i2)') flag
! 
!       if (flag.lt.10) age_flag(1:1)='0'

      ! Creates new file with associated time step, catchment number, and age type (1-8)
      ! Outputs range of values and their probabilities
      open (106,file=run(1:nrun)//'/pdf_cascade/Timestep_'//timestep//'_Catchment_'//catchment//'_pdf.dat',status='unknown')

        write (106,*) 'TITLE = "Cascade PDF"'
        write (106,'(A)',ADVANCE="no") 'VARIABLES = "Age Range (Ma)"'
        if (age_flags(1).eq.1) write (106,'(A)',ADVANCE="no") ' "AHe - Farley, 2000 Frequency"'
        if (age_flags(2).eq.1) write (106,'(A)',ADVANCE="no") ' "AHe - a='//as//' um Frequency"'
        if (age_flags(3).eq.1) write (106,'(A)',ADVANCE="no") ' "AHe - a='//am//' um Frequency"'
        if (age_flags(4).eq.1) write (106,'(A)',ADVANCE="no") ' "AHe - a='//al//' um Frequency"'
        if (age_flags(5).eq.1) write (106,'(A)',ADVANCE="no") ' "AHe - eU='//eUl//' ppm Frequency"'
        if (age_flags(6).eq.1) write (106,'(A)',ADVANCE="no") ' "AHe - eU='//eUm//' ppm Frequency"'
        if (age_flags(7).eq.1) write (106,'(A)',ADVANCE="no") ' "AHe - eU='//eUh//' ppm Frequency"'
        if (age_flags(9).eq.1) write (106,'(A)',ADVANCE="no") ' "AFT Frequency"'
        if (age_flags(8).eq.1) write (106,'(A)',ADVANCE="no") ' "ZHe Frequency"'
        if (age_flags(10).eq.1) write (106,'(A)',ADVANCE="no") ' "ZFT Frequency"'
        if (age_flags(11).eq.1) write (106,'(A)',ADVANCE="no") ' "MAr Frequency"'
        if (age_flags(12).eq.1) write (106,'(A)',ADVANCE="no") ' "KFeldAr Frequency"'
        if (age_flags(13).eq.1) write (106,'(A)',ADVANCE="no") ' "BioAr Frequency"'
        if (age_flags(14).eq.1) write (106,'(A)',ADVANCE="no") ' "MusAr Frequency"'
        if (age_flags(15).eq.1) write (106,'(A)',ADVANCE="no") ' "HornAr Frequency"'

        write (106,*) 'ZONE T = "Cascade PDF"'
        write (106,*) 'I=',counter,', J=1, K=1, ZONETYPE=Ordered'
        write (106,*) 'DATAPACKING=POINT'
        write (106,'(A)',ADVANCE="no") ' DT=(DOUBLE'
        do j=1,num_of_age_flags
          if(age_flags(j) .eq. 1) write (106,'(A)',ADVANCE="no") ' DOUBLE'
        enddo
        write (106,'(A)') ')'



      do k=1,counter
        write (106,'(f10.3)',advance="no") n(k)
        do j=1,7
          if (age_flags(j) .eq. 1) write (106,'(f15.5)',advance="no") normPDF(j,k)
        enddo
        if (age_flags(9) .eq. 1) write (106,'(f15.5)',advance="no") normPDF(9,k)
        if (age_flags(8) .eq. 1) write (106,'(f15.5)',advance="no") normPDF(8,k)
        if (age_flags(10) .eq. 1) write (106,'(f15.5)',advance="no") normPDF(10,k)
        if (age_flags(11) .eq. 1) write (106,'(f15.5)',advance="no") normPDF(11,k)
        if (age_flags(12) .eq. 1) write (106,'(f15.5)',advance="no") normPDF(12,k)
        if (age_flags(13) .eq. 1) write (106,'(f15.5)',advance="no") normPDF(13,k)
        if (age_flags(14) .eq. 1) write (106,'(f15.5)',advance="no") normPDF(14,k)
        if (age_flags(15) .eq. 1) write (106,'(f15.5)',advance="no") normPDF(15,k)
        write (106,*)
      enddo
      close (106)

    deallocate (n,P,Psum,normPDF)

  return

end subroutine pdfmaker_for_cascade

end module m_pdfmaker_for_cascade

