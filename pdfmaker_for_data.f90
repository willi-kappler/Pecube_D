module m_pdfmaker_for_data
use m_compiler

contains

subroutine pdfmaker_for_data (age,error,number,age_type,run,xoutlet,youtlet,nrun)

  !% April 26, 2005  Greg Stock
  !% Calculates synoptic probability density function (SPDF) to create basin-wide PDFs
  !% from individual detrital apatite He ages; based on Equation 2 from Ruhl and Hodges (2005).
  !% Input is a two-column Excel (xls) file of measured He ages and one sigma uncertainties.
  !% Assumes Gaussian PDF distribution for each grain age, then sums all grain PDF's to create
  !% basin-wide PDF.
  !
  ! Converted to FORTRAN by cspath 03/08
  ! Measured ages and errors are read in by Pecube.f90 and passed into this subroutine
  ! where A PDF for the specified basin is written

  implicit none

  integer number,k,age_type,j,num_xchars,num_ychars
  real*8 dx,i,start,finish,pi,xoutlet,youtlet
  real*8 age(number)
  real*8 error(number)
  integer counter,nrun
  real*8,dimension(:),allocatable::n,P,Psum,normPDF
  character run*100,x_basin_char*100,y_basin_char*100,age_char*2

  num_xchars = 0
  num_ychars = 0
  dx=0.0001
  pi=3.14159

  start = minval(age(:))-10*maxval(error(:))-mod(minval(age(:))-10*maxval(error(:)),dx)           ! Calculates minimal age value
  finish = maxval(age(:))+10*maxval(error(:))-mod(maxval(age(:))+10*maxval(error(:)),dx)        ! Calculates maximum age value

  ! The 'counter' variable holds how many age values will have an associated probability value
  counter=0
  i=start
  do while (i .lt. finish)
    counter=counter+1
    i=i+dx
  enddo

  allocate (n(counter),P(counter),Psum(counter),normPDF(counter))

  counter=0
  i=start
  do while (i .lt. finish)
    counter=counter+1
    n(counter)=i                   ! The 'n' array holds the actual age values for which a probability will be calculated
    i=i+dx
  enddo

  Psum(:)=0.
  do k=1,number
    !% Calculate probability for range of n values from Eq. 2a of Ruhl and
    !% Hodges (2005)
    P=(1./(error(k)*sqrt(2*pi)))*exp(-0.5*((n-age(k))/(error(k)))**2)
    !% Sum the resulting PDF with previous PDF's from Eq. 2b of Ruhl and
    !% Hodges (2005)
    Psum=Psum+P
  enddo

  !% Normalize area under PDF curve to 1
  normPDF=(Psum/number)          !% Normalize PDF so area under curve = 1
  !check=sum(normPDF)            !% Check that total area under new PDF curve = 1.

  call sys_command('if [ ! -d '//run(1:nrun)//'/pdf_data ]; then mkdir '//run(1:nrun)//'/pdf_data; fi')  ! Checks if the pdf_data folder exists in the run output directory and creates it if it does not exist

  do j=1,100
    x_basin_char(j:j)=' '
    y_basin_char(j:j)=' '
  enddo
  write (x_basin_char,'(f100.2)') xoutlet
  write (y_basin_char,'(f100.2)') youtlet
  do j=1,100
    if (x_basin_char(j:j).ne.'0' .and. x_basin_char(j:j).ne.'1' .and. x_basin_char(j:j).ne.'2' &
        .and. x_basin_char(j:j).ne.'3' .and. x_basin_char(j:j).ne.'4' .and. x_basin_char(j:j).ne.'5' &
        .and. x_basin_char(j:j).ne.'6' .and. x_basin_char(j:j).ne.'7' .and. x_basin_char(j:j).ne.'8' &
        .and. x_basin_char(j:j).ne.'9' .and. x_basin_char(j:j).ne.'.') num_xchars=j+1
    if (y_basin_char(j:j).ne.'0' .and. y_basin_char(j:j).ne.'1' .and. y_basin_char(j:j).ne.'2' &
        .and. y_basin_char(j:j).ne.'3' .and. y_basin_char(j:j).ne.'4' .and. y_basin_char(j:j).ne.'5' &
        .and. y_basin_char(j:j).ne.'6' .and. y_basin_char(j:j).ne.'7' .and. y_basin_char(j:j).ne.'8' &
        .and. y_basin_char(j:j).ne.'9' .and. y_basin_char(j:j).ne.'.') num_ychars=j+1
  enddo

  write (age_char,'(i2)') age_type
  if (age_type.lt.10) age_char(1:1)='0'

  ! Creates a new file for the specified basin outlet and age type
  ! Writes the age values and their associated probabilities
  open (106,file=run(1:nrun)//'/pdf_data/Basin_X_'//x_basin_char(num_xchars:100)//'_Y_'//y_basin_char(num_ychars:100)//'_AgeType_'//age_char//'_pdf.dat',status='unknown')
  do k=1,counter
    write (106,'(2f10.5)') n(k), normPDF(k)
  enddo
  close (106)

  return

end subroutine

end module m_pdfmaker_for_data

