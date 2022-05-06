module m_pdfmaker_for_pecube
use m_compiler
use m_logger

contains

subroutine pdfmaker_for_pecube (age,number,flag,t,c,run,edot,x_basin_char,y_basin_char,&
                                num_xchars,num_ychars,nrun)

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
  integer number,counter,flag,k,t,c,j
  real*8,dimension(:),allocatable::n,P,Psum,normPDF
  real*8 age(number),error(number),per_error(11),edot(number)
  character timestep*4,catchment*4,run*100,age_flag*2
  character x_basin_char*100,y_basin_char*100
  integer num_xchars,num_ychars,nrun

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

  dx = 0.01
  pi = 3.14159265

  ! Calculates absolute error of the age data based on age type
  do j=1,number
    error(j) = age(j) * per_error(flag)
  enddo

  start = minval(age(:))-10*maxval(error(:))-mod(minval(age(:))-10*maxval(error(:)),dx)           ! Finds minimum possible value
  finish = maxval(age(:))+10*maxval(error(:))-mod(maxval(age(:))+10*maxval(error(:)),dx)        ! Finds maximum possible value

  call log_message("pdfmaker_for_pecube.f90: start: " + start + ", finish: " + finish)

!   start = nint(minval(age(:))-2*maxval(error(:)))
!   finish = nint(maxval(age(:))+10*maxval(error(:)))

  !% Range of n values to consider, rounded to whole numbers
  counter = 0
  i = start
  do while (i .lt. finish)
    counter = counter + 1
    i = i + dx
  enddo

  allocate (n(counter),P(counter),Psum(counter),normPDF(counter))

  counter = 0
  i = start
  do while (i .lt. finish)
    counter = counter + 1
    n(counter) = i
    i = i + dx
  enddo

  Psum(:) = 0.
  do k=1,number
    !% Calculate probability for range of n values from Eq. 2a of Ruhl and
    !% Hodges (2005)
    P=(1./(error(k)*sqrt(2*pi)))*exp(-0.5*((n-age(k))/(error(k)))**2)*edot(k)
    !% Sum the resulting PDF with previous PDF's from Eq. 2b of Ruhl and
    !% Hodges (2005)
    Psum=Psum+P
  enddo

  ! Calculates the average erosion rate
  ! Used to normalize area under curve to 1
  avg_edot=0.
  do k=1,number
    avg_edot = avg_edot + edot(k)
  enddo
  avg_edot = avg_edot/number

  !% Normalize area under PDF curve to 1
  normPDF=(Psum/number)/avg_edot       !% Normalize PDF so area under curve = 1


      write (timestep,'(i4)') t
      if (t.lt.1000) timestep(1:1)='0'
      if (t.lt.100) timestep(1:2)='00'
      if (t.lt.10) timestep(1:3)='000'

    ! c will be -1 when this function is called from pdfmaker_for_data2.f because c refers to a catchment number
    ! and there is only catchment numbers when using Cascade models
    ! No significance to using c to determine what to output but simpler than adding another flag to be passed in
    if (c .ne. -1) then

      write (catchment,'(i4)') c
      if (c.lt.1000) catchment(1:1)='0'
      if (c.lt.100) catchment(1:2)='00'
      if (c.lt.10) catchment(1:3)='000'

      write (age_flag,'(i2)') flag

      if (flag.lt.10) age_flag(1:1)='0'

      call sys_command('if [ ! -d '//run(1:nrun)//'/pdf_cascade ]; then mkdir '//run(1:nrun)//'/pdf_cascade; fi')     ! Checks if the pdf_pecube folder exists in the associated run output directory; if not, then it creates it

      ! Creates new file with associated time step, catchment number, and age type (1-8)
      ! Outputs range of values and their probabilities
      open (106,file=run(1:nrun)//'/pdf_cascade/Timestep_'//timestep//'_Catchment_'//catchment//'_Agetype_'//age_flag//'_pdf.dat',status='unknown')
      do k=1,counter
        write (106,'(2f10.5)') n(k),normPDF(k)
      enddo
      close (106)

    else

      call sys_command('if [ ! -d '//run(1:nrun)//'/pdf_pecube ]; then mkdir '//run(1:nrun)//'/pdf_pecube; fi')        ! Checks if the pdf_data folder exists and creates it if it doesn't

      write (age_flag,'(i2)') flag

      if (flag.lt.10) age_flag(1:1)='0'

      ! Creates file with associated basin outlet coordiantes as specified by the user and the cooresponding age type
      ! Writes the range of age values and their probabilities
      open (204,file=run(1:nrun)//'/pdf_pecube/Timestep_'//timestep//'_Basin_X_'//x_basin_char(num_xchars:100)//'_Y_'//y_basin_char(num_ychars:100)//'_Agetype_'//age_flag//'_pdf.dat',status='unknown')
      do k=1,counter
        write (204,'(2f10.5)') n(k),normPDF(k)
      enddo
      close (204)

    endif

    deallocate (n,P,Psum,normPDF)

  return

end subroutine

end module m_pdfmaker_for_pecube
