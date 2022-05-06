module m_deplss
    use m_logger
contains

! All of these subroutines have been taken from the TAPES_G
! program and specifically the tapes_g.f file. Below is source reference:
! '******************************************************'
! '*                                                    *'
! '*                    TAPES_GRD                       *'
! '*                                                    *'
! '* TERRAIN ANALYSIS PROGRAMS FOR THE ENVIRONMENTAL    *'
! '*            SCIENCES - GRID VERSION                 *'
! '*                                                    *'
! '*               VERSION 6.3   July 1997              *'
! '*                                                    *'
! '******************************************************'
! '*         Developed by:                              *'
! '*                                                    *'
! '*                   IAN D. MOORE                     *'
! '*    Centre for Resource & Environmental Studies     *'
! '*        The Australian National University          *'
! '*          Canberra, ACT 0200, Australia             *'
! '*                                                    *'
! '*      and is now maintained by John Gallant         *'
! '*             of the same address                    *'
! '******************************************************'
! Note: Most of the subroutines have been modified to handle the new UPSTREAMJ and
! UPSTREAMI arrays which hold the x and y indices, respectively, of the points which
! flow into the associated point

!C =================================================================

      SUBROUTINE DEPLSS(NR,NC,LL,CAREA,DI,DIS,ZDD,ZD,CHAN,MFP, &
      Z_temp,UPSTREAMI,UPSTREAMJ)

      IMPLICIT INTEGER*4 (A-Y)
      PARAMETER(LJ=5000,LI=5000)
      INTEGER*4 JOUT(20),IOUT(20),KTEST
      INTEGER*4 IDUM,LL,NR,NC,MM,JN,IN,MFP,IDUM9
      ! Z(-2:LJ+3,-2:LI+3)
      REAL*8 WT, FD(8), SLOPE,DUM1,DUM2,DUM3,DUM4, &
        DUM5,DUM7,DUM8,ASPECT

      ! COUNT(-2:LJ+3,-2:LI+3)
      ! ASP(LJ,LI)
      ! FPATH(-2:LJ+3,-2:LI+3)
      ! SUM(-2:LJ+3,-2:LI+3)

      LOGICAL ACTIVITY,GOAGAIN,DRD8
      ! ACTIVE(-2:LI+3)
      DIMENSION  SELECT(256)
      ! DIR(-2:LJ+3,-2:LI+3)
      ! DIRA(-2:LJ+3,-2:LI+3)
      ! SDIR(-2:LJ+3,3)
      ! DDIR(LJ)
      REAL*8 CAREA, DI, DIS, ZDD, ZD, CHAN, Z_temp(NC,NR)
      INTEGER*4 UPSTREAMI(NC,NR,8),UPSTREAMJ(NC,NR,8),i,j,LASTL
      INTEGER*4 PASS, FIRSTL, IM1, IM2, DUM6

      integer(4), dimension(:,:), allocatable :: Z, DIR, SDIR, DIRA
      integer(4), dimension(:), allocatable :: DDIR
      real(8), dimension(:,:), allocatable :: COUNT, ASP, FPATH, SUM
      logical, dimension(:), allocatable :: ACTIVE


!      EXTERNAL LINK
!       COMMON ACTIVE
!       COMMON/CON/COUNT
!       COMMON/PATH/FPATH
!       COMMON/DIREC/DIR
!       COMMON/DIREE/DIRA
!       COMMON/OUTLET/JOUT,IOUT,NOUT,MCELL
!       COMMON/ASPEC/ASP
      DATA SELECT/  0,  1,  2,  2,  4,  1,  2,  2,  8,  1, &
       8,  2,  8,  4,  4,  2, 16, 16, 16,  2, 16,  4,  4, &
       2,  8,  8,  8,  8,  8,  8,  8,  4, 32,  1,  2,  2, &
       4,  4,  2,  2, 32,  8,  8,  2,  8,  8,  4,  4, 32, &
      32, 32, 32, 16, 32,  4,  2, 16, 16, 16, 16,  8, 16, &
       8,  8, 64, 64, 64,  1, 64,  1,  2,  2, 64, 64,  8, &
       2,  8,  8,  4,  2, 16, 64, 64,  2, 16, 64,  2,  2, &
      16,  8,  8,  8,  8,  8,  8,  4, 32, 64, 32,  1, 32, &
      32, 32,  2, 32, 32, 32,  2, 32,  8,  4,  4, 32, 32, &
      32, 32, 32, 32, 32, 32, 32, 32, 16, 16, 16, 16,  8, &
       8,128,128,128,  1,  4,  1,  2,  2,128,128,  2,  1, &
       8,  4,  4,  2, 16,128,  2,  1,  4,128,  2,  1,  8, &
     128,  8,  1,  8,  8,  4,  2, 32,128,  1,  1,128,128, &
       2,  1, 32,128, 32,  1,  8,128,  4,  2, 32, 32, 32, &
       1, 32,128, 32,  1, 16, 16, 16,  1, 16, 16,  8,  4, &
     128,128,128,128,128,128,  2,  1,128,128,128,  1,128, &
     128,  4,  2, 64,128,128,  1,128,128,128,  1,  8,128, &
       8,  1,  8,  8,  8,  2, 64,128, 64,128, 64,128, 64, &
     128, 32, 64, 64,128, 64, 64, 64,  1, 32, 64, 64,128, &
      64, 64, 64,128, 32, 32, 32, 64, 32, 32, 16,128/

!    ************************************************************
!        This program fills depressions in a DEM, computes flow
!     directions and cumulative cell counts.  The subroutine is
!     based on the methods described in the following references
!     reference.  [See also JENSON & DOMINGUE (1988)]
!C
!        Creating Depressionless DEMs:
!C
!     JENSON & TRAUTWEIN, 1987. Methods and applications in
!     surface depression analysis. Proc. Auto-Carto 8,
!     This is the D8 component of this algorithm.
!C
!        Users have the option of using the Rho8 algorithm for
!     calculating flow directions instead of the D8 algorithm.
!     It uses a random component for flow paths in NE, SE, SW or
!     NW directions.  It is a modification of the method described
!     in the following reference. This method overcomes some of
!     the problems in being limited to only 8 flow directions.
!C
!     FAIRFIELD & LEYMARIE, 1991. Drainage networks from grid
!     digital elevation models. Water Resour. Res. 27(5): 709-
!     717.
!C
!        Both the D8 and Rho8 algorithms can be used with a
!     multiple drainage path algorithm.  The multiple flow path
!     algorithm is used for upslope contributing areas less than
!     a defined critical area and either the D8 or Rho8 elsewhere.
!     The algorithm is based in part on the following references.
!C
!     FREEMAN, 1991.  Calculating catchment area with divergent
!     flow based on a regular grid. Computers & Geosci. 17(3):
!     413-422.
!C
!     QUINN et al., 1991.  The prediction of hillslope flow paths
!     for distributed hydrological modelling using digital terrain
!     models. Hydrological Processes 5(1): 59-79.
!C
!        An alternative approach to the above is to compute drain-
!     age areas using a grid-based stream-tube analogy.  The method
!     has been adapted to that described in;
!C
!     COSTA-CABRAL & BURGESS, 1993.  DEMON (Digital Elevation
!     Model Networks): a model of flow over hillslopes for comput-
!     ation of specific contributing and dispersal areas.  Water
!     Resour. Res. (submitted).
!     ************************************************************

      allocate(Z(-2:LJ+3,-2:LI+3))
      allocate(DIR(-2:LJ+3,-2:LI+3))
      allocate(SDIR(-2:LJ+3,3))
      allocate(DIRA(-2:LJ+3,-2:LI+3))
      allocate(DDIR(LJ))

      allocate(COUNT(-2:LJ+3,-2:LI+3))
      allocate(ASP(LJ,LI))
      allocate(FPATH(-2:LJ+3,-2:LI+3))
      allocate(SUM(-2:LJ+3,-2:LI+3))

      allocate(ACTIVE(-2:LI+3))

      COUNT = 0.
      FPATH = 0.
      DIR = 0
      JOUT = 0
      IOUT = 0
      NOUT = 0
      MCELL = 2
      DRD8 = .FALSE.

      call log_message("subroutine deplss")

      do j=1,NR
        do i=1,NC
          Z(i,j) = NINT(Z_temp(i,j) * 100.0) + 1
        enddo
      enddo
!     --------------------------------------------------------------

      WRITE(6,850)

!     Begin flow direction computations
!C
!     DIR codes as ... 7  8  1
!                      6     2
!                      5  4  3
!C

      DO J=0,NC+1
        Z(J,0)=0
        Z(J,NR+1)=0
        DIR(J,0)=0
        DIR(J,NR+1)=0
        SUM(J,0)=0.
        SUM(J,NR+1)=0.
        COUNT(J,0)=-1.
        COUNT(J,NR+1)=-1.
      ENDDO
      DO I=0,NR+1
        COUNT(0,I)=-1.
        COUNT(NC+1,I)=-1.
        SUM(0,I)=0.
        SUM(NC+1,I)=0.
        DIR(0,I)=0
        DIR(NC+1,I)=0
        Z(0,I)=0
        Z(NC+1,I)=0
      ENDDO
      IDUM=-1

      DO I=1,NR
      DO J=1,NC
         IF(Z(J,I).EQ.0) GO TO 65
!C
!        NCELL=1: Cells on boundary will drain outside boundary
!        NCELL=2: Cells on boundary will drain to a cell inside the
!                 catchment boundary
!C
         NCELL=MCELL
!          IF(NCELL.EQ.2) THEN
!             DO 62 K=1,NOUT
! 62          IF(J.EQ.JOUT(K).AND.I.EQ.IOUT(K)) NCELL=1
!          ENDIF
         DIR(J,I)=THEDIR(NCELL,Z(J,I),Z(J+1,I-1),Z(J+1,I),&
              Z(J+1,I+1),Z(J,I+1),Z(J-1,I+1),Z(J-1,I),Z(J-1,I-1),&
              Z(J,I-1),IDUM,DRD8)
!     -------------------------------------------------------------
65       CONTINUE
         IF(DIR(J,I).LT.0) GOTO 70
         DIR(J,I)=SELECT(DIR(J,I)+1)
 70   CONTINUE
      ENDDO
      ENDDO
!C
!     Now make a pass resolving non-flats with more than one down
!     link.  Iterate linking in the flats
!C
      DO I=1,NR
      ACTIVE(I)=.TRUE.
      ENDDO
      ACTIVE(0)=.FALSE.
      ACTIVE(NR+1)=.FALSE.
      I1=1
      I2=2
      I3=3
      FIRSTL=1
      LASTL=NR
      PASS=0
!C
!     Process the downward pass
!C
 80   ACTIVITY=.FALSE.
      DO J=0,NC+1
      SDIR(J,I1)=DIR(J,FIRSTL-1)
      SDIR(J,I2)=DIR(J,FIRSTL)
      SDIR(J,I3)=DIR(J,FIRSTL+1)
      ENDDO
      PASS=PASS+1
      WRITE(6,855) PASS,FIRSTL,LASTL
      I=FIRSTL
 90   ACTIVE(I)=.FALSE.
 95   GOAGAIN = .FALSE.
      DO J=1,NC
      DDIR(J)=SDIR(J,I2)
      IF(SDIR(J,I2).LT.0) DDIR(J)=LINK2(SDIR(J,I2),ACTIVE(I),&
       SDIR(J+1,I1),SDIR(J+1,I2),SDIR(J+1,I3),SDIR(J,I3),&
       SDIR(J-1,I3),SDIR(J-1,I2),SDIR(J-1,I1),SDIR(J,I1),SELECT,&
       ACTIVITY,GOAGAIN)
      ENDDO
      IF (GOAGAIN) GOTO 95
      DO J=1,NC
      DIR(J,I)=DDIR(J)
      ENDDO
!C
! --  ROTATE TO THE NEXT LINE
!C
110   I=I+1
      IF(I.GE.LASTL+1) GOTO 115
      ITEMP=I1
      I1=I2
      I2=I3
      I3=ITEMP
      IF(.NOT.ACTIVE(I).AND..NOT.ACTIVE(I+1).AND..NOT.&
                  ACTIVE(I+2)) GOTO 110
      DO J=1,NC
      SDIR(J,I3)=DIR(J,I+1)
      ENDDO
      IF(.NOT.ACTIVE(I)) GOTO 110
      GOTO 90
!C
!     Done with this iteration, update FIRSTL and LASTL and go again
!C
115   DO I=FIRSTL,LASTL
      IF(ACTIVE(I)) GOTO 125
      CONTINUE
      ENDDO
!C
!     All done
!C
      GOTO 300
125   IF (ACTIVITY)GOTO 130
      WRITE(6,860)
      GOTO 300
130   FIRSTL=I
      DO I=LASTL,FIRSTL,-1
      IF(ACTIVE(I)) GOTO 140
      CONTINUE
      ENDDO
140   LASTL=I
!C
!     Process the upward pass
!C
      DO J=1,NC
      SDIR(J,I3)=DIR(J,LASTL+1)
      SDIR(J,I2)=DIR(J,LASTL)
      SDIR(J,I1)=DIR(J,LASTL-1)
      ENDDO
      PASS=PASS+1
      WRITE(6,855) PASS,FIRSTL,LASTL
      I=LASTL
160   ACTIVE(I)=.FALSE.
165   GOAGAIN = .FALSE.
      DO J=1,NC
      DDIR(J)=SDIR(J,I2)
      IF(SDIR(J,I2).LT.0) DDIR(J)=LINK2(SDIR(J,I2),ACTIVE(I),&
        SDIR(J+1,I1),SDIR(J+1,I2),SDIR(J+1,I3),SDIR(J,I3),&
        SDIR(J-1,I3),SDIR(J-1,I2),SDIR(J-1,I1),SDIR(J,I1),&
        SELECT,ACTIVITY,GOAGAIN)
      ENDDO
      IF (GOAGAIN) GOTO 165
      DO J=1,NC
      DIR(J,I)=DDIR(J)
      ENDDO
!C
!     Rotate to the next line
!C
180   I=I-1
      IF(I.LE.FIRSTL-1) GOTO 185
      ITEMP=I3
      I3=I2
      I2=I1
      I1=ITEMP
      IM1=MAX(I-1,1)
      IM2=MAX(I-2,1)
      IF(.NOT.ACTIVE(I).AND..NOT.ACTIVE(IM1).AND..NOT.&
                   ACTIVE(IM2)) GOTO 180
      DO J=1,NC
      SDIR(J,I1)=DIR(J,I-1)
      ENDDO
      IF(.NOT.ACTIVE(I)) GOTO 180
      GOTO 160
!C
!     Done with this iteration, update FIRSTL and LASTL and go again
!C
185   DO I=LASTL,FIRSTL,-1
      IF(ACTIVE(I)) GOTO 200
      CONTINUE
      ENDDO
!C
!     All done
!C
      GOTO 300
200   IF (ACTIVITY) GOTO 205
      WRITE(6,860)
      GOTO 300
205   LASTL=I
      DO I=FIRSTL,LASTL
      IF(ACTIVE(I)) GOTO 215
      CONTINUE
      ENDDO
215   FIRSTL=I
!C
!     End of upward pass
!C
      GOTO 80
300   CONTINUE
      DO I=1,NR
      DO J=1,NC
         IF(DIRA(J,I).GT.0) DIR(J,I)=DIRA(J,I)
      CONTINUE
      ENDDO
      ENDDO
!     --------------------------------------------------------------
!      WRITE(6,810)
!      READ(5,*) WT
      WT=1.0
      !WRITE(6,815)
      !READ(5,*) MFP
      GOTO (520,510,530) MFP
!     ---------------------------------------------------------------
510   continue
      !WRITE(6,818)
      !READ(5,*) CAREA
      CAREA=CAREA/DI/DI
      IF(CAREA.EQ.0.0) THEN
         MFP=1
         GOTO 520
      ENDIF
      DO I=1,NR
      DO J=1,NC
         IF(Z(J,I).EQ.0) THEN
            SUM(J,I)=0.0
         ELSE
            NCELL=MCELL
!             IF(NCELL.EQ.2) THEN
!                DO 301 K=1,NOUT
! 301            IF(J.EQ.JOUT(K).AND.I.EQ.IOUT(K)) NCELL=1
!             ENDIF
            SUM(J,I)=THEFLW(NCELL,Z(J,I),Z(J+1,I-1),Z(J+1,I),&
             Z(J+1,I+1),Z(J,I+1),Z(J-1,I+1),Z(J-1,I),Z(J-1,I-1),&
             Z(J,I-1))
         ENDIF
      CONTINUE
      ENDDO
      ENDDO
      GOTO 521
!     -----------------------------------------------------------
!     Begin counting routine - for MFP = 1 and 2
!C
!     For the first pass, find all the mask cells and all the cells
!     that nothing points to and set their counts to a value of
!     'WT'.  When COUNT = -1 it hasn't been solved for yet.  A
!     DIR = 0 indicates a mask.
!C
520   continue
      !IF (DRD8) THEN
      !   WRITE(RECORD, *) 'OPTION Flow accumulation Rho8'
      !ELSE
      !   WRITE(RECORD, *) 'OPTION Flow accumulation D8'
      !ENDIF
      !METAOK = METAWRITERECORD(METAUN, RECORD, METAFTYPE)
521   continue
      WRITE(6,865)
      MM=0
      DO I=1,NR
      DO J=1,NC
      IF(MFP.EQ.2) THEN
         CALL FLOW(Z(J,I),Z(J+1,I-1),Z(J+1,I),Z(J+1,I+1),Z(J,I+1),&
          Z(J-1,I+1),Z(J-1,I),Z(J-1,I-1),Z(J,I-1),FD,SUM(J+1,I-1),&
          SUM(J+1,I),SUM(J+1,I+1),SUM(J,I+1),SUM(J-1,I+1),&
          SUM(J-1,I),SUM(J-1,I-1),SUM(J,I-1),DIR(J+1,I-1),&
          DIR(J+1,I),DIR(J+1,I+1),DIR(J,I+1),DIR(J-1,I+1),&
          DIR(J-1,I),DIR(J-1,I-1),DIR(J,I-1))
      ENDIF
      COUNT(J,I)=THEZEROS(DIR(J,I),DIR(J+1,I-1),DIR(J+1,I),&
        DIR(J+1,I+1),DIR(J,I+1),DIR(J-1,I+1),DIR(J-1,I),&
        DIR(J-1,I-1),DIR(J,I-1),WT,MFP,FD)
      IF(COUNT(J,I).GE.0.0) MM=MM+1
      CONTINUE
      ENDDO
      ENDDO
      WRITE(6,870) MM,(LL-MM)
!C
!     Now iterating on the CFILE, computing COUNTs
!C
      PASS=0
      KTEST=1
310   ACTIVITY=.FALSE.
      PASS=PASS+1
      WRITE(*,820) PASS
      DO I=1,NR
315   GOAGAIN=.FALSE.
      DO J=1,NC
      IF(COUNT(J,I).LT.0) THEN
         IF(MFP.EQ.2.AND.KTEST.EQ.1) THEN
            CALL FLOW(Z(J,I),Z(J+1,I-1),Z(J+1,I),Z(J+1,I+1),Z(J,I+1),&
             Z(J-1,I+1),Z(J-1,I),Z(J-1,I-1),Z(J,I-1),FD,SUM(J+1,I-1),&
             SUM(J+1,I),SUM(J+1,I+1),SUM(J,I+1),SUM(J-1,I+1),&
             SUM(J-1,I),SUM(J-1,I-1),SUM(J,I-1),DIR(J+1,I-1),&
             DIR(J+1,I),DIR(J+1,I+1),DIR(J,I+1),DIR(J-1,I+1),&
             DIR(J-1,I),DIR(J-1,I-1),DIR(J,I-1))
         ENDIF
         COUNT(J,I)=LADO(COUNT(J,I),COUNT(J+1,I-1),COUNT(J+1,I),&
          COUNT(J+1,I+1),COUNT(J,I+1),COUNT(J-1,I+1),&
          COUNT(J-1,I),COUNT(J-1,I-1),COUNT(J,I-1),DIR(J+1,I-1),&
          DIR(J+1,I),DIR(J+1,I+1),DIR(J,I+1),DIR(J-1,I+1),&
          DIR(J-1,I),DIR(J-1,I-1),DIR(J,I-1),SUM(J,I),CAREA,&
          ACTIVITY,GOAGAIN,FPATH(J,I),FPATH(J+1,I-1),FPATH(J+1,I),&
          FPATH(J+1,I+1),FPATH(J,I+1),FPATH(J-1,I+1),FPATH(J-1,I),&
          FPATH(J-1,I-1),FPATH(J,I-1),WT,MFP,FD,J,I,UPSTREAMI,&
          UPSTREAMJ,NC,NR)
         IF(COUNT(J,I).GE.0.0) MM=MM+1
      ENDIF
      CONTINUE
      ENDDO
      IF (GOAGAIN) GOTO 315
      CONTINUE
      ENDDO
      WRITE(6,870) MM,(LL-MM)
!C
!     Done with this iteration
!C
      IF (.NOT.ACTIVITY) GOTO 400
      ACTIVITY=.FALSE.
      PASS=PASS+1
      WRITE(6,830) PASS
      DO I=NR,1,-1
330   GOAGAIN=.FALSE.
      DO J=1,NC
      IF(COUNT(J,I).LT.0) THEN
         IF(MFP.EQ.2.AND.KTEST.EQ.1) THEN
            CALL FLOW(Z(J,I),Z(J+1,I-1),Z(J+1,I),Z(J+1,I+1),Z(J,I+1),&
             Z(J-1,I+1),Z(J-1,I),Z(J-1,I-1),Z(J,I-1),FD,SUM(J+1,I-1),&
             SUM(J+1,I),SUM(J+1,I+1),SUM(J,I+1),SUM(J-1,I+1),&
             SUM(J-1,I),SUM(J-1,I-1),SUM(J,I-1),DIR(J+1,I-1),&
             DIR(J+1,I),DIR(J+1,I+1),DIR(J,I+1),DIR(J-1,I+1),&
             DIR(J-1,I),DIR(J-1,I-1),DIR(J,I-1))
         ENDIF
         COUNT(J,I)=LADO(COUNT(J,I),COUNT(J+1,I-1),COUNT(J+1,I),&
          COUNT(J+1,I+1),COUNT(J,I+1),COUNT(J-1,I+1),COUNT(J-1,I),&
          COUNT(J-1,I-1),COUNT(J,I-1),DIR(J+1,I-1),DIR(J+1,I),&
          DIR(J+1,I+1),DIR(J,I+1),DIR(J-1,I+1),DIR(J-1,I),&
          DIR(J-1,I-1),DIR(J,I-1),SUM(J,I),CAREA,ACTIVITY,GOAGAIN,&
          FPATH(J,I),FPATH(J+1,I-1),FPATH(J+1,I),FPATH(J+1,I+1),&
          FPATH(J,I+1),FPATH(J-1,I+1),FPATH(J-1,I),FPATH(J-1,I-1),&
          FPATH(J,I-1),WT,MFP,FD,J,I,UPSTREAMI,UPSTREAMJ,NC,NR)
         IF(COUNT(J,I).GE.0.0) MM=MM+1
      ENDIF
      CONTINUE
      ENDDO
      IF (GOAGAIN) GOTO 330
      CONTINUE
      ENDDO
      WRITE(6,870) MM,(LL-MM)
!C
!     Done with this iteration
!C
400   IF (ACTIVITY) GOTO 310
      IF(MM.NE.LL.AND.KTEST.EQ.1) THEN
         WRITE(6,880) (LL-MM)
         KTEST=2
         DO I=1,NR
         DO J=1,NC
            SUM(J,I)=0.0
         CONTINUE
         ENDDO
         ENDDO
         GOTO 310
      ELSE IF(MM.NE.LL.AND.KTEST.EQ.2) THEN
         WRITE(6,885)
      ENDIF

      GOTO 700

!      deallocate(Z)
!      deallocate(DIR)
!      deallocate(SDIR)
!      deallocate(DIRA)
!      deallocate(DDIR)

!      deallocate(COUNT)
!      deallocate(ASP)
!      deallocate(FPATH)
!      deallocate(SUM)

!      deallocate(ACTIVE)

      RETURN
!     --------------------------------------------------------------
!     Compute catchment areas and specific catchment areas using
!     a stream-tube approach.  The flow direction (ASP) is measured
!     in degrees clockwise from north: E-90, S-180, W-270, N-0,360
!     and is the aspect of a cell with nodal point (J,I) at its
!     centroid.
!C
530   CONTINUE
      !WRITE(RECORD, *) 'OPTION Flow accumulation DEMON'
      !METAOK = METAWRITERECORD(METAUN, RECORD, METAFTYPE)
!              For all pits set DIR=0, Boundary of DEM has DIR=0
      DO I=1,NR
      DO J=1,NC
         IF(DIR(J,I).LT.0) DIR(J,I)=0
      CONTINUE
      ENDDO
      ENDDO
      DO I=1,NR
      DO J=1,NC
         JN=J
         IN=I
         CALL ATRIB(JN,IN,SLOPE,DUM1,ASPECT,DUM2,DUM3,DUM4,ZDD,ZD,&
                     DIS,DUM5,CHAN,DUM6,DUM7,DUM8,IDUM9,2,Z,&
                     COUNT,DIR)
         ASP(J,I)=ASPECT
         IF(SLOPE.LE.0.0025.AND.DIR(J,I).NE.0) ASP(J,I)=45.*(1.+&
                     ALOG(REAL(DIR(J,I)))/ALOG(2.))
         IF(DIRA(J,I).GT.0) ASP(J,I)=45.*(1.+ALOG(REAL(DIRA(J,I)))/&
                 ALOG(2.))
      CONTINUE
      ENDDO
      ENDDO
      WRITE(6,890)
      !CALL DEMON(NR,NC,LL,DI)
      WRITE(6,895)
!     --------------------------------------------------------------
820   FORMAT(3X,'Downward pass ',I4)
830   FORMAT(3X,'Upward pass ',I4)
850   FORMAT(/'COMPUTING PRIMARY FLOW DIRECTIONS'/)
855   FORMAT(3X,'Downward pass',I3,' FIRSTL =',I6,' LASTL =',I6)
860   FORMAT(/'COULD NOT SOLVE FOR ALL CELLS'/)
865   FORMAT(/'COMPUTING COUNTS (CATCHMENT AREAS)'/)
870   FORMAT(3X,'Resolved =',I7,' Unresolved =',I7)
880   FORMAT(/'UNABLE TO RESOLVE',I4,' CELLS - changing to D8'/&
       10X,'flow direction algorithm for these cells')
885   FORMAT(/'*******************************************'/&
             '** WARNING - UNABLE TO RESOLVE ALL CELLS **'/&
             '**  DEPRESSIONS STILL EXIST IN THE DEM   **'/&
             '*******************************************')
890   FORMAT(/'Begin DEMON calculations')
895   FORMAT('DEMON completed')
!     --------------------------------------------------------------

700   deallocate(Z)
      deallocate(DIR)
      deallocate(SDIR)
      deallocate(DIRA)
      deallocate(DDIR)

      deallocate(COUNT)
      deallocate(ASP)
      deallocate(FPATH)
      deallocate(SUM)

      deallocate(ACTIVE)

      RETURN
      END subroutine


! ================================================================
      FUNCTION THEDIR(NCELL,MID,N1,N2,N3,N4,N5,N6,N7,N8,IDUM,DRD8)
      IMPLICIT INTEGER*4 (A-Z)
      INTEGER*4 M(8),HEDIR,THEDIR
      INTEGER*4 NCELL
      REAL N(8),MAXDROP,RNUM,RHO8
      LOGICAL DRD8
      M(1)=N1
      M(2)=N2
      M(3)=N3
      M(4)=N4
      M(5)=N5
      M(6)=N6
      M(7)=N7
      M(8)=N8
!C
!     Return a < 0 mask if the paths are flat
!C
      HEDIR=0
      DO I=1,7,2
        IF(DRD8) THEN
           CALL RAN_tapesg(IDUM,RNUM)
           RHO8=1.0/(2.0-RNUM)
        ELSE
           RHO8=1.0/SQRT(2.0)
        ENDIF
      N(I)=RHO8*(MID-M(I))
      ENDDO
      DO I=2,8,2
      N(I)= MID-M(I)
      ENDDO
      MAXDROP=-6000.
      DO I=1,8
        IF(NCELL.EQ.2.AND.M(I).EQ.0) GOTO 8
        IF(M(I).EQ.0) THEN
           J=I+4
           IF(J.GT.8) J=J-8
           IF(N(J).GT.0.0) GOTO 8
        ENDIF
        IF(MAXDROP.LT.N(I).AND.M(I).NE.0) MAXDROP=N(I)
  8   CONTINUE
      ENDDO
      DO I=1,8
        IF(NCELL.EQ.2.AND.M(I).EQ.0) GOTO 10
        IF (N(I).EQ.MAXDROP) HEDIR=HEDIR+2**(I-1)
 10   CONTINUE
      ENDDO
      IF(MAXDROP.EQ.0) HEDIR=-HEDIR
!C
!     A pit will be a -300
!C
      IF(MAXDROP.LT.0) HEDIR=-300
      THEDIR=HEDIR
      RETURN
      END function



! ==============================================================
      FUNCTION THEFLW(NCELL,MID,M1,M2,M3,M4,M5,M6,M7,M8)
      IMPLICIT INTEGER*4 (A-Z)
      INTEGER*4 M(8)
      INTEGER*4 NCELL
      REAL THEFLW,N(8),SUM
!C
      M(1)=M1
      M(2)=M2
      M(3)=M3
      M(4)=M4
      M(5)=M5
      M(6)=M6
      M(7)=M7
      M(8)=M8
      SUM=0.0
      DO I=1,7,2
      N(I)=(MID-M(I))/1.414
      ENDDO
      DO I=2,8,2
      N(I)=MID-M(I)
      ENDDO
      DO I=1,8
        IF(NCELL.EQ.2.AND.M(I).EQ.0) GOTO 8
        IF(N(I).GT.0.0) THEN
           IF(M(I).GT.0.0) THEN
              SUM=SUM+(N(I)**1.1)
           ELSE
              SUM=0.
              GOTO 10
           ENDIF
        ENDIF
  8   CONTINUE
      ENDDO
 10   CONTINUE
      THEFLW=SUM
      RETURN
      END function



! ===============================================================
      FUNCTION LINK2(CENTER,ACTIVE,D1,D2,D3,D4,D5,D6,D7,D8,SELECT,&
        ACTIVITY,GOAGAIN)
      IMPLICIT INTEGER*4(A-Z)
      INTEGER*4 VERT(8),HORI(8), OUTF
      LOGICAL ACTIVE,C(8),ACTIVITY,GOAGAIN
      DIMENSION SELECT(256),BITMASK(8)
      DATA BITMASK/1,2,4,8,16,32,64,128/
      LINK2=CENTER

      HORI(1) = 1
      HORI(2) = 1
      HORI(3) = 1
      HORI(4) = 0
      HORI(5) = -1
      HORI(6) = -1
      HORI(7) = -1
      HORI(8) = 0

      VERT(1) = -1
      VERT(2) = 0
      VERT(3) = 1
      VERT(4) = 1
      VERT(5) = 1
      VERT(6) = 0
      VERT(7) = -1
      VERT(8) = -1

!C
!     Check if it is a pit
!C
      IF(LINK2.EQ.-300) GOTO 100
      CWORK=-CENTER
      DO I=8,1,-1
        C(I)=.FALSE.
        IF(CWORK-BITMASK(I).LT.0) GOTO 5
        CWORK=CWORK-BITMASK(I)
        C(I)=.TRUE.
 5    CONTINUE
      ENDDO
!C
!     Check for downstream linkks
!C
      OUTF=0
      IF(D1.NE.16.AND.D1.GT.0.AND.C(1)) then
        OUTF=OUTF+1
!         UPSTREAMJ(Jval+1,Ival-1,5) = Jval+1+HORI(5)
!         UPSTREAMI(Jval+1,Ival-1,5) = Ival-1+VERT(5)
      ENDIF
      IF(D2.NE.32.AND.D2.GT.0.AND.C(2)) then
        OUTF=OUTF+2
!         UPSTREAMJ(Jval+1,Ival,6) = Jval+1+HORI(6)
!         UPSTREAMI(Jval+1,Ival,6) = Ival+VERT(6)
      ENDIF
      IF(D3.NE.64.AND.D3.GT.0.AND.C(3)) then
        OUTF=OUTF+4
!         UPSTREAMJ(Jval+1,Ival+1,7) = Jval+1+HORI(7)
!         UPSTREAMI(Jval+1,Ival+1,7) = Ival+1+VERT(7)
      ENDIF
      IF(D4.NE.128.AND.D4.GT.0.AND.C(4)) then
        OUTF=OUTF+8
!         UPSTREAMJ(Jval,Ival+1,8) = Jval+HORI(8)
!         UPSTREAMI(Jval,Ival+1,8) = Ival+1+VERT(8)
      ENDIF
      IF(D5.NE.1.AND.D5.GT.0.AND.C(5)) then
        OUTF=OUTF+16
!         UPSTREAMJ(Jval-1,Ival+1,1) = Jval-1+HORI(1)
!         UPSTREAMI(Jval-1,Ival+1,1) = Ival+1+VERT(1)
      ENDIF
      IF(D6.NE.2.AND.D6.GT.0.AND.C(6)) then
        OUTF=OUTF+32
!         UPSTREAMJ(Jval-1,Ival,2) = Jval-1+HORI(2)
!         UPSTREAMI(Jval-1,Ival,2) = Ival+VERT(2)
      ENDIF
      IF(D7.NE.4.AND.D7.GT.0.AND.C(7)) then
        OUTF=OUTF+64
!         UPSTREAMJ(Jval-1,Ival-1,3) = Jval-1+HORI(3)
!         UPSTREAMI(Jval-1,Ival-1,3) = Ival-1+VERT(3)
      ENDIF
      IF(D8.NE.8.AND.D8.GT.0.AND.C(8)) then
        OUTF=OUTF+128
!         UPSTREAMJ(Jval,Ival-1,4) = Jval+HORI(4)
!         UPSTREAMI(Jval,Ival-1,4) = Ival-1+VERT(4)
      ENDIF
      IF(OUTF.EQ.0) GOTO 10
      CENTER=SELECT(OUTF+1)
      LINK2=CENTER
      ACTIVITY=.TRUE.
      GOAGAIN=.TRUE.
      GOTO 100
 10   ACTIVE=.TRUE.
 100  CONTINUE
      RETURN
      END function


! =================================================================
      SUBROUTINE RAN_tapesg(IDUM,RNUM)
      DIMENSION R(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=1./M1)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=1./M2)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
      DATA IFF /0/
!     *************************************************************
!        Subroutine returns a uniform random deviate between 0.0
!     & 1.0.  Set IDUM to any negative value to initialize or
!     reinitialize the sequence.  This routine is derived from
!     an algorithm described in the following reference.
!C
!     PRESS W.H., FLANNERY B.P., TEUKOLSKY S.A. & VETTERLING W.T.,
!     1989. Numerical Recipes (Fortran version), Cambridge Uni-
!     versity Press, Sydney, pp 196-197.
!     *************************************************************
      ! may be uninitialized
      ix3 = 0
      IF(IDUM.LT.0.OR.IFF.EQ.0) THEN
     IFF=1
     IX1=MOD(IC1-IDUM,M1)
     IX1=MOD(IA1*IX1+IC1,M1)
     IX2=MOD(IX1,M2)
     IX1=MOD(IA1*IX1+IC1,M1)
     IX3=MOD(IX1,M3)
     DO J=1,97
       IX1=MOD(IA1*IX1+IC1,M1)
       IX2=MOD(IA2*IX2+IC2,M2)
       R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
       CONTINUE
     ENDDO
     IDUM=1
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
!      IF(J.GT.97.OR.J.LT.1) PAUSE
      RNUM=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
      END subroutine


! ============================================================
      FUNCTION THEZEROS(MID,N1,N2,N3,N4,N5,N6,N7,N8,WT,MFP,F)
      IMPLICIT INTEGER*4 (A-Z)
      INTEGER*4 MFP
      REAL*8 THEZEROS,WT,F(8)
!     ********************************************************
!     THE N'S ARE THE DIRS ARE ARRANGED AS N7 N8 N1
!                                          N6    N2
!                                          N5 N4 N3
!
!                         AND ENCODED AS   64 128 1
!                                          32     2
!                                          16   8 4
!     ********************************************************
      THEZEROS=-1.
      IF(MID.NE.0) GOTO 5
      THEZEROS=WT
      GOTO 100
!
!     If nothing points to me make me zero
!
 5    IF(N1.EQ.16.OR.N2.EQ.32.OR.N3.EQ.64.OR.N4.EQ.128.OR.&
        N5.EQ.1.OR.N6.EQ.2.OR.N7.EQ.4.OR.N8.EQ.8) GOTO 100
      IF(MFP.EQ.2) THEN
         DO I=1,8
           IF(F(I).GT.0.0) GOTO 100
         CONTINUE
         ENDDO
      ENDIF
      THEZEROS=WT
100   CONTINUE
      RETURN
      END function
! ===============================================================
      SUBROUTINE FLOW(MID,M1,M2,M3,M4,M5,M6,M7,M8,F,S1,S2,S3,S4,&
        S5,S6,S7,S8,DD1,DD2,DD3,DD4,DD5,DD6,DD7,DD8)
      real(8) :: S1, S2, S3, S4, S5, S6, S7, S8
      INTEGER*4 M(8)
      INTEGER*4 D(8),DD(8),DD1,DD2,DD3,DD4,DD5,DD6,DD7,DD8
      REAL*8 F(8),SUM(8),N(8)
      DATA D/16,32,64,128,1,2,4,8/
!
      M(1)=M1
      M(2)=M2
      M(3)=M3
      M(4)=M4
      M(5)=M5
      M(6)=M6
      M(7)=M7
      M(8)=M8
      DD(1)=DD1
      DD(2)=DD2
      DD(3)=DD3
      DD(4)=DD4
      DD(5)=DD5
      DD(6)=DD6
      DD(7)=DD7
      DD(8)=DD8
      SUM(1)=S1
      SUM(2)=S2
      SUM(3)=S3
      SUM(4)=S4
      SUM(5)=S5
      SUM(6)=S6
      SUM(7)=S7
      SUM(8)=S8
      DO I=1,7,2
      N(I)=(M(I)-MID)/1.414
      ENDDO
      DO I=2,8,2
      N(I)=M(I)-MID
      ENDDO
      DO I=1,8
      IF(SUM(I).EQ.0.0) THEN
         IF(DD(I).EQ.D(I)) THEN
            F(I)=1.0
         ELSE
            F(I)=0.0
         ENDIF
      ELSEIF(N(I).LE.0.0) THEN
         F(I)=0.0
      ELSE
         F(I)=(N(I)**1.1)/SUM(I)
      ENDIF
      CONTINUE
      ENDDO
      RETURN
      END subroutine
! ===============================================================
      FUNCTION LADO(CENTER,C1,C2,C3,C4,C5,C6,C7,C8,D1,D2,D3,D4,&
       D5,D6,D7,D8,SUM,CAREA,ACTIVITY,GOAGAIN,PP,P1,P2,P3,P4,&
       P5,P6,P7,P8,WT,MFP,F,Jval,Ival,UPSTREAMI,UPSTREAMJ,&
       NC,NR)
      IMPLICIT INTEGER*4(A-Z)
      ! PARAMETER(LJ=5000,LI=5000)
      INTEGER*4 MFP,counter,HORI(8),VERT(8), NC, NR
      INTEGER*4 Jval, Ival
      REAL*8 CAREA
      REAL*8 LADO,CENTER,C1,C2,C3,C4,C5,C6,C7,C8,C(8),TOTAL,PP,P1,&
       P2,P3,P4,P5,P6,P7,P8,P(8),PT(8),WT,A,F(8),SUM
      DIMENSION N(8),D(8)
      LOGICAL ACTIVITY,GOAGAIN
      INTEGER*4 UPSTREAMI(NC,NR,8),UPSTREAMJ(NC,NR,8)
      DATA N/16,32,64,128,1,2,4,8/
!
      C(1)=C1
      C(2)=C2
      C(3)=C3
      C(4)=C4
      C(5)=C5
      C(6)=C6
      C(7)=C7
      C(8)=C8
      D(1)=D1
      D(2)=D2
      D(3)=D3
      D(4)=D4
      D(5)=D5
      D(6)=D6
      D(7)=D7
      D(8)=D8
      P(1)=P1
      P(2)=P2
      P(3)=P3
      P(4)=P4
      P(5)=P5
      P(6)=P6
      P(7)=P7
      P(8)=P8

      HORI(1) = 1
      HORI(2) = 1
      HORI(3) = 1
      HORI(4) = 0
      HORI(5) = -1
      HORI(6) = -1
      HORI(7) = -1
      HORI(8) = 0

      VERT(1) = -1
      VERT(2) = 0
      VERT(3) = 1
      VERT(4) = 1
      VERT(5) = 1
      VERT(6) = 0
      VERT(7) = -1
      VERT(8) = -1

!
!     C's are counts, D's are directions, P's are flow paths
!
      LADO=CENTER
      TOTAL=0.0
      DO I=1,8
      PT(I)=0.0
      ENDDO
!     ------------------------------------------------
      IF(MFP.EQ.2) THEN
         DO I=1,8
           IF(C(I).GT.CAREA.OR.SUM.EQ.0.0) THEN
              IF(D(I).NE.N(I)) THEN
                 F(I)=0.0
              ELSE
                 F(I)=1.0
              ENDIF
           ELSE
              IF(C(I).LT.0.0.AND.F(I).GT.0.0) GOTO 100
           ENDIF
         CONTINUE
         ENDDO
      ENDIF
!     ----------------------------------------------
      counter=1
      DO I=1,8
      A=1.0
      IF(I.EQ.1.OR.I.EQ.3.OR.I.EQ.5.OR.I.EQ.7) A=SQRT(2.0)
      IF(MFP.EQ.2) THEN
         IF(F(I).EQ.0.0) GOTO 20
         IF(C(I).LT.0.0) GOTO 100
         TOTAL=TOTAL+F(I)*C(I)
         UPSTREAMJ(Jval,Ival,I) = Jval+HORI(I)
         UPSTREAMI(Jval,Ival,I) = Ival+VERT(I)
         IF(D(I).EQ.N(I)) PT(I)=P(I)+A
      ELSE
         IF(D(I).NE.N(I)) GOTO 20
         IF(C(I).LT.0.0) GOTO 100
         TOTAL=TOTAL+C(I)
         UPSTREAMJ(Jval,Ival,I) = Jval+HORI(I)
         UPSTREAMI(Jval,Ival,I) = Ival+VERT(I)
         PT(I)=P(I)+A
      ENDIF

20    CONTINUE
      ENDDO
      LADO=TOTAL+WT
      PP=AMAX1(PT(1),PT(2),PT(3),PT(4),PT(5),PT(6),PT(7),PT(8))
      ACTIVITY=.TRUE.
      GOAGAIN=.TRUE.
100   CONTINUE
      RETURN
      END function

! =================================================================
      SUBROUTINE ATRIB(J,I,SLOPE,SLOPED8,ASPECT,PROFC,PLANC,&
       GRADC,ZDD,ZD,DIS,DIS2,CHAN,FLOWD,DASDS,WIDTH,MFP,INDX,Z,&
       COUNT,DIR)
      PARAMETER(LJ=5000,LI=5000)
      INTEGER*4 Z(-2:LJ+3,-2:LI+3),HEDIR,FLOWD
      REAL*8 ZDD, ZD, DIS, CHAN, SLOPE, SLOPED8, ASPECT
      real(8) :: PROFC, PLANC, GRADC, DIS2, WIDTH
      real(8) :: SXX, SYY, SXY, DASDS, COUNT, CNT
      DIMENSION S(8),II(8),JJ(8),COUNT(-2:LJ+3,-2:LI+3)
      INTEGER*4 DIR(-2:LJ+3,-2:LI+3),D(8),N(8)
!       COMMON/CON/COUNT
!       COMMON/DIREC/DIR
      DATA II/-1,0,1,1,1,0,-1,-1/
      DATA JJ/1,1,1,0,-1,-1,-1,0/
      DATA D/1,2,4,8,16,32,64,128/
      DATA N/16,32,64,128,1,2,4,8/
!     *************************************************************
!        Compute SLOPE, ASPECT, PLANE CURVATURE, PROFILE CURVATURE
!     & TANGENTAIAL CURVATURE using a central finite difference
!     approximation and the rate of change of specific catchment
!     area in the flow direction using a forward difference
!     approximation.  Method produces identical results to method of
!     Zevenbergen and Thorne (1987), except for plan curvature,
!     who fit a 9 term polynomial to a 3x3 moving grid and
!     determine the slope, aspect, etc. of the point in the centre
!     of the grid.
!     *************************************************************
!     JDUM     Variable determining location in DEM or on edge of
!              DEM field
!     SX,SY    First partial derivatives in x and y directions
!     SXX,SYY  Second partial derivatives in x and y directions
!     SXY      d2z/(dxdy)
!     P        SX**2 + SY**2
!     Q        P+1
!     Z        Elevation (x100)
!     *************************************************************
!
!     Compute slopes  SLOPE= SQRT(P} - in m/m
!                     Slope angle = artan{SLOPE) - in degrees
!

      ! may be uninitialized
      sxy = 0.0

      JDUM=1
      IDUM=1
      IF(Z(J-1,I).EQ.ZD) JDUM=2
      IF(Z(J+1,I).EQ.ZD) JDUM=3
      IF(Z(J+1,I).EQ.ZD.AND.Z(J-1,I).EQ.ZD) JDUM=4
      IF(Z(J,I-1).EQ.ZD) IDUM=2
      IF(Z(J,I+1).EQ.ZD) IDUM=3
      IF(Z(J,I-1).EQ.ZD.AND.Z(J,I+1).EQ.ZD) IDUM=4
      IF(JDUM.EQ.1) THEN
         SX=real(0.5*FLOAT(Z(J+1,I)-Z(J-1,I))/DIS)
      ELSEIF (JDUM.EQ.2) THEN
         SX=real(FLOAT(Z(J+1,I)-Z(J,I))/DIS)
      ELSEIF (JDUM.EQ.3) THEN
         SX=real(FLOAT(Z(J,I)-Z(J-1,I))/DIS)
      ELSEIF (JDUM.EQ.4) THEN
         SX=0.0
      ENDIF
      IF(IDUM.EQ.1) THEN
         SY=real(0.5*FLOAT(Z(J,I+1)-Z(J,I-1))/DIS)
      ELSEIF (IDUM.EQ.2) THEN
         SY=real(FLOAT(Z(J,I+1)-Z(J,I))/DIS)
      ELSEIF (IDUM.EQ.3) THEN
         SY=real(FLOAT(Z(J,I)-Z(J,I-1))/DIS)
      ELSE
         SY=0.0
      ENDIF
      SX=real(SX*CHAN)
      SY=real(SY*CHAN)
      P=SX*SX+SY*SY
      Q=P+1.0
      SLOPE=SQRT(ABS(P))
!     ------------------------------------------------------------
!     Compute aspect in units of degrees clockwise from north
!              ASPECT = arctan(Sy/Sx)
!     where ASPECT=0 in westerly direction
!
      IF(IDUM.EQ.4.OR.JDUM.EQ.4) THEN
         ASPECT=real(ZDD)
      ELSE
         IF(SX.NE.0.0) THEN
            ASPECT=180.0*(1.0-0.5*ATAN(-SY/SX)/ACOS(0.0))+90.0*&
                            SX/ABS(SX)
         ELSE
            ASPECT=0.0
            IF(SY.LT.0.0) ASPECT=180.0
            IF(P.EQ.0.0) ASPECT=real(ZDD)
         ENDIF
      ENDIF
      IF(INDX.EQ.2) RETURN
!     -----------------------------------------------------------
!     Compute slope based on actual maximum slope to eight nearest
!     neighbours (SLOPED8) and flow directions (FLOWD)
!                K's are orgainized as      7 8 1
!                                           6   2
!                                           5 4 3
!
      DO K=1,8
         IF(Z(J+JJ(K),I+II(K)).NE.ZD) THEN
            S(K)=real(FLOAT((Z(J,I)-Z(J+JJ(K),I+II(K))))/DIS)
         ELSE
            S(K)=0.0
         ENDIF
      CONTINUE
      ENDDO
      DO K=1,7,2
      S(K)=S(K)/SQRT(2.0)
      ENDDO
      SMAX=-60000.
      SLOPED8=0.0
      DO K=1,8
        IF(S(K).GT.SMAX) SMAX=S(K)
!        IF(ABS(S(K)).GT.SLOPED8) SLOPED8=ABS(S(K))
        IF(S(K).GT.SLOPED8) SLOPED8=S(K)
      CONTINUE
      ENDDO
      HEDIR=0
      DO K=1,8
        IF(S(K).GE.SMAX.AND.S(K).GE.0.0) HEDIR=HEDIR+2**(K-1)
      CONTINUE
      ENDDO
      SLOPED8=real(SLOPED8*CHAN)
      FLOWD=HEDIR
!     ------------------------------------------------------------
!     Calculation of conponents of the second derivatives of
!     the elevation surface - for use in estimating the curvature
!     parameters
!
      IF(JDUM.EQ.1) THEN
         SXX=(Z(J+1,I)+Z(J-1,I)-2.0*Z(J,I))/DIS2
      ELSEIF (JDUM.EQ.2) THEN
         SXX=(Z(J+2,I)+Z(J,I)-2.0*Z(J+1,I))/DIS2
      ELSEIF (JDUM.EQ.3) THEN
         SXX=(Z(J,I)+Z(J-2,I)-2.0*Z(J-1,I))/DIS2
      ELSEIF (JDUM.EQ.4) THEN
         SXX=0.0
      ENDIF
      IF(IDUM.EQ.1) THEN
         SYY=(Z(J,I+1)+Z(J,I-1)-2.0*Z(J,I))/DIS2
      ELSEIF (IDUM.EQ.2) THEN
         SYY=(Z(J,I+2)+Z(J,I)-2.0*Z(J,I+1))/DIS2
      ELSEIF (IDUM.EQ.3) THEN
         SYY=(Z(J,I)+Z(J,I-2)-2.0*Z(J,I-1))/DIS2
      ELSE
         SYY=0.0
      ENDIF
      IF(JDUM.EQ.1.AND.IDUM.EQ.1) THEN
         SXY=0.25*(Z(J+1,I+1)-Z(J-1,I+1)-Z(J+1,I-1)+Z(J-1,I-1))/DIS2
      ELSEIF(JDUM.EQ.1.AND.IDUM.NE.1) THEN
         IF(IDUM.EQ.2) THEN
            SXY=0.5*(Z(J+1,I+1)-Z(J-1,I+1)-Z(J+1,I)+Z(J-1,I))/DIS2
         ELSEIF(IDUM.EQ.3) THEN
            SXY=0.5*(Z(J+1,I)-Z(J-1,I)-Z(J+1,I-1)+Z(J-1,I-1))/DIS2
         ENDIF
      ELSEIF(JDUM.EQ.2) THEN
         IF(IDUM.EQ.1) THEN
            SXY=0.5*(Z(J+1,I+1)-Z(J,I+1)-Z(J+1,I-1)+Z(J,I-1))/DIS2
         ELSEIF(IDUM.EQ.2) THEN
            SXY=0.5*(Z(J+1,I+1)-Z(J,I+1)-Z(J+1,I)+Z(J,I))/DIS2
         ELSEIF(IDUM.EQ.3) THEN
            SXY=0.5*(Z(J+1,I)-Z(J,I)-Z(J+1,I-1)+Z(J,I-1))/DIS2
         ENDIF
      ELSEIF(JDUM.EQ.3) THEN
         IF(IDUM.EQ.1) THEN
            SXY=0.5*(Z(J,I+1)-Z(J-1,I+1)-Z(J,I-1)+Z(J-1,I-1))/DIS2
         ELSEIF(IDUM.EQ.2) THEN
            SXY=0.5*(Z(J,I+1)-Z(J-1,I+1)-Z(J,I)+Z(J-1,I))/DIS2
         ELSEIF(IDUM.EQ.3) THEN
            SXY=0.5*(Z(J,I)-Z(J-1,I)-Z(J,I-1)+Z(J-1,I-1))/DIS2
         ENDIF
      ENDIF
      IF(JDUM.EQ.4.OR.IDUM.EQ.4) SXY=0.0
      SXX=real(SXX*CHAN)
      SYY=real(SYY*CHAN)
      SXY=real(SXY*CHAN)
!     ------------------------------------------------------------
!     Compute plan, profile and tangent curvature
!     PROFC = [Sxx(Sx**2) + 2SxySxSy + Syy(Sy**2)]/[P(Q**1.5)]
!     PLANC = [Sxx(Sy**2) - 2SxySxSy + Syy(Sx**2)]/(P**1.5)
!     GRADC = [Sxx(Sy**2) - 2SxySxSy + Syy(Sx**2)]/(P*Q**0.5)
!
      IF(P.EQ.0.0) THEN
         PROFC=0.0
         PLANC=0.0
         GRADC=0.0
      ELSEIF(JDUM.EQ.4.OR.IDUM.EQ.4) THEN
         PROFC=real(ZDD)
         PLANC=real(ZDD)
         GRADC=real(ZDD)
      ELSE
         PROFC=(SXX*SX*SX+2.0*SXY*SX*SY+SYY*SY*SY)/(P*Q**1.5)
         PLANC=(SXX*SY*SY-2.0*SXY*SX*SY+SYY*SX*SX)/(P**1.5)
         GRADC=(SXX*SY*SY-2.0*SXY*SX*SY+SYY*SX*SX)/(P*Q**0.5)
!     GRADC is a modified plan curvature- curvature of normal
!           plane section in a direction perpendicular to gradient
!           (direction of tangent to the contour line) - tangent
!           curvature
      ENDIF
!     -------------------------------------------------------------
!     Calculate rate of change of specific catchment area in the
!     primary flow direction and flow widths
!
      KD=0
      DS1=1.0
      DS2=1.0
      CNT=0
      DO K=1,8
         JK=J+JJ(K)
         IK=I+II(K)
         IF(DIR(JK,IK).EQ.N(K)) THEN
            IF(COUNT(JK,IK).GT.CNT) THEN
               CNT=COUNT(JK,IK)
               IF(K.EQ.1.OR.K.EQ.3.OR.K.EQ.5.OR.K.EQ.7) DS1=SQRT(2.0)
            ENDIF
         ENDIF
         IF(DIR(J,I).EQ.D(K)) THEN
            IF(K.EQ.1.OR.K.EQ.3.OR.K.EQ.5.OR.K.EQ.7) DS2=SQRT(2.0)
         ENDIF
      CONTINUE
      ENDDO
      CALL FWIDTH(J,I,MFP,WIDTH,ASPECT,Z,DIR)
      DS=0.5*(DS1+DS2)
      DASDS=(COUNT(J,I)-CNT)/DS/WIDTH
      RETURN
      END subroutine

! ===============================================================
      SUBROUTINE FWIDTH(J,I,MFP,WIDTH,ASPECT,Z,DIR)
      real(8) :: ASPECT, ASP1, alpha, beta, WIDTH
      PARAMETER(LJ=5000,LI=5000)
      PARAMETER (PI=3.14159265358979323846)
      PARAMETER (DTOR=PI/180.0)
      INTEGER*4 DIR(-2:LJ+3,-2:LI+3)
      INTEGER*4 Z(-2:LJ+3,-2:LI+3)
!       COMMON/DIREC/DIR
!     ************************************************************
!     This subroutine computes the flow width for the element.
!     ************************************************************
      diag=SQRT(2.)
      ! may be uninitialized
      alpha = 0.0
      IF(MFP.EQ.1) THEN
         ASP1=45.*(1.+ALOG(REAL(DIR(J,I)))/ALOG(2.))
      ELSEIF(MFP.EQ.2) THEN
         WIDTH = FD8WIDTH(J,I,Z)
         RETURN
      ELSE
     ASP1=ASPECT
      ENDIF
      IF(ASP1.LE.90.) alpha=90.-ASP1
      IF(ASP1.GT.270..AND.ASP1.LE.360.) alpha=ASP1-270.
      IF(ASP1.GT.180..AND.ASP1.LE.270.) alpha=ASP1-180.
      IF(ASP1.GT.90..AND.ASP1.LE.180.) alpha=ASP1-90.
      beta=ABS(alpha-45.)
      WIDTH=diag*COS(beta*DTOR)
      RETURN
      END subroutine
! ==============================================================
      FUNCTION FD8WIDTH(JJ,II,Z)
      IMPLICIT INTEGER*4 (A-Z)
      PARAMETER(LJ=5000,LI=5000)
      INTEGER*4 Z(-2:LJ+3,-2:LI+3)
      INTEGER*4 M(8),MID
      REAL FD8WIDTH,N(8),SUM,CSUM,A
!
!   This subroutine calculates a modified flow width for the
!   FD8/FRho8 multiple flow direction algorithm. It is based
!   on the assumption that the desired specific catchment area
!   is the average specific catchment area leaving the cell.
!   Flow width must therefore have a minimum value of 1 for
!   convergent and planar flows, but increase to a maximum
!   of 4 when flow in all directions is equal, as at a
!   hilltop.
!
!   Added John Gallant August 1995
!
      I = II
      J = JJ
      MID = Z(J,I)
      M(1)=Z(J+1,I-1)
      M(2)=Z(J+1,I)
      M(3)=Z(J+1,I+1)
      M(4)=Z(J,I+1)
      M(5)=Z(J-1,I+1)
      M(6)=Z(J-1,I)
      M(7)=Z(J-1,I-1)
      M(8)=Z(J,I-1)
      SUM=0.0
      DO I=1,7,2
      N(I)=(MID-M(I))/1.414
      ENDDO
      DO I=2,8,2
      N(I)=MID-M(I)
      ENDDO
      DO I=1,8
        IF(N(I).GT.0.0.AND.M(I).GT.0.0) THEN
           N(I) = N(I)**1.1
           SUM=SUM+N(I)
        ENDIF
      CONTINUE
      ENDDO
      CSUM = 0.0
      DO I=1,8
        IF(N(I).GT.0.0.AND.M(I).GT.0.0) THEN
           A = N(I) / SUM
           CSUM = CSUM + A*A
        ENDIF
      CONTINUE
      ENDDO
      CONTINUE
      IF (CSUM.EQ.0) THEN
         A = 1.0
      ELSE
         A = 1.0 / CSUM
      ENDIF
      IF (A.LT.3) THEN
         FD8WIDTH = 1.0
      ELSE
         FD8WIDTH = 1.0 + 3.0*((A-3.0)/5.0)
      ENDIF
      RETURN
      END function

end module m_deplss

