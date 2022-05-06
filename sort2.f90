module m_sort2
contains

! Subroutine sort2 from Numerical Recipes in Fortran 90
! Slightly modified to include passing-in of array sizes
!
! dwhipp - 04/08

      SUBROUTINE sort2(arr,n1,slave,n2)
      USE nrutil, ONLY : assert_eq
      use m_indexx
      IMPLICIT NONE
      REAL*8,DIMENSION(n1),INTENT(INOUT) :: arr
      REAL*8,DIMENSION(n2),INTENT(INOUT) :: slave
      INTEGER :: ndum,n1,n2
      INTEGER,DIMENSION(n1) :: index
      ndum=assert_eq(size(arr),size(slave),'sort2')
      call indexx(arr,n1,index)
      arr=arr(index)
      slave=slave(index)
      END SUBROUTINE sort2

end module m_sort2

