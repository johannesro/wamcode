      PROGRAM readbound
! Read ECWAM spectral boundary file and print to screen.
! Stolen from bouint.F
! 
! CCC must have real*8 on Njord
! xlf -O2 readbound_ecwam.f -o readbound_ecwam
!
! Linux: gfortran readbound_ecwam.f -o readbound_ecwam
! 
! USAGE:
! readbound_ecwam /work/carrasco/wamrun/B0120070726000000
!
! 
! 2011-12-19, Oyvind.Breivik@met.no
! 
      !USE YOWTEXT  , ONLY : PATH

      IMPLICIT NONE

      ! Functions called
      INTEGER iargc

      INTEGER :: narg
      INTEGER :: K, M, IJ 
      INTEGER :: NANG, NFRE, MBMAX
      INTEGER :: KL, ML, NBOUNC, IDELPRC
      INTEGER :: IU01
      INTEGER :: NTIME

      ! Strings
      CHARACTER*256 fin
      CHARACTER*14 CDATE1
      !CHARACTER CDATE1(14)

      ! CCC *8
      REAL*8 :: XANG, XFRE, TH0, FR1, CO, XBOU, XDELC, XDELF
      REAL*8, ALLOCATABLE, DIMENSION(:) :: XLAT, XLON
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: FMEAN, EMEAN, THQ 
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: F


!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *USERID*    CHARACTER USER IDENTIFIER.
!      *RUNDI*     CHARACTER RUN  IDENTIFIER OF INPUT DATA.
!      *FILEDI*    CHARACTER FILE IDENTIFIER OF INPUT DATA.
!      *PATHI*     CHARACTER DIRECTORY OF INPUT DATA.
!      *RUNDO*     CHARACTER RUN  IDENTIFIER OF OUTPUT DATA.
!      *FILEDO*    CHARACTER FILE IDENTIFIER OF OUTPUT DATA.
!      *PATHO*     CHARACTER DIRECTORY OF OUTPUT DATA.
! ----------------------------------------------------------------------

!*    1. INITIALISATION.
!        ---------------
      NTIME = 1
 
!*    1.1 UNITS.
!         ------
 
      IU01 = 1
      NANG=-1
      NFRE=-1
      MBMAX=-1

      ! Command line
      narg=iargc()
      if (narg < 1) then
         print *,"USAGE: readbound_ecwam wamboundary.dat"
         print *,"Print ECMWF WAM boundary file to screen."
         stop
      endif

      ! ECWAM boundary file
      call getarg(1,fin)
      open (IU01,file=fin,form="unformatted",err=1001)
      write (*,*) "CCC opened IU01"

!*    2.2 READ BOUNDARY FILE HEADER.
!         --------------------------
 
      READ (IU01, ERR=1001, END=1001)
     &     XANG, XFRE, TH0, FR1, CO, XBOU, XDELC
      KL = NINT(XANG)
      ML = NINT(XFRE)
      NBOUNC  = NINT(XBOU)
      IDELPRC = NINT(XDELC)
      WRITE(*,*) ' '
      WRITE(*,*) ' INPUT FILE HEADER:'
      WRITE(*,*) ' NO. OF DIRECTIONS IS      KL     = ', KL
      WRITE(*,*) ' NO. OF FREQUENCIES IS     ML     = ', ML
      WRITE(*,*) ' FIRST DIRECTION IS        TH0    = ', TH0
      WRITE(*,*) ' FIRST FREQUENCY IS        FR(1)  = ', FR1
      WRITE(*,*) ' FREQUENCY RATIO IS        CO     = ', CO
      WRITE(*,*) ' NO. OF BOUNDRAY POINTS IS NBOUNC = ', NBOUNC
      WRITE(*,*) ' TIME STEP OF DATA IS      IDELPRC= ', IDELPRC

!*    2.3 CHECK DIMENSIONS.
!         -----------------
 
      IF (NANG .EQ. -1 .AND. NFRE .EQ. -1 .AND. MBMAX.EQ.-1) THEN
        NANG=KL
        NFRE=ML
        MBMAX=NBOUNC
        ALLOCATE(F(MBMAX,NANG,NFRE,2))
        ALLOCATE(FMEAN(MBMAX,2))
        ALLOCATE(EMEAN(MBMAX,2))
        ALLOCATE(THQ(MBMAX,2))
        ALLOCATE(XLAT(MBMAX))
        ALLOCATE(XLON(MBMAX))
      ELSE IF (KL.GT.NANG .OR. ML.GT.NFRE .OR. NBOUNC.GT.MBMAX) THEN
         WRITE(*,*) '*******************************************'
         WRITE(*,*) '*                                         *'
         WRITE(*,*) '*    FATAL ERROR PROGRAM READBOUND_ECWAM  *'
         WRITE(*,*) '*    ============================         *'
         WRITE(*,*) '* ONE OR MORE DIMENSIONS ARE TO SMALL.    *'
         WRITE(*,*) '* NO. OF DIRECTIONS IS      KL     = ', KL
         WRITE(*,*) '*    DIMENSION IS           NANG   = ', NANG
         WRITE(*,*) '* NO. OF FREQUENCIES IS     ML     = ', ML
         WRITE(*,*) '*             DIMENSION IS  NFRE   = ', NFRE
         WRITE(*,*) '* NO. OF BOUNDRAY POINTS IS NBOUNC = ', NBOUNC
         WRITE(*,*) '*              DIMENSION IS MBMAX  = ', MBMAX
         WRITE(*,*) '*                                         *'
         WRITE(*,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
         WRITE(*,*) '*                                         *'
         WRITE(*,*) '*******************************************'
         CALL ABORT
      ENDIF

!*    3.1 READ BOUNDARY VALUES.
!         ---------------------
 
999   CONTINUE
      DO IJ=1,NBOUNC
         READ (IU01, ERR=1001, END=1001)
     &     XLON(IJ), XLAT(IJ), CDATE1, EMEAN(IJ,1),
     &     THQ(IJ,1), FMEAN(IJ,1)
         READ (IU01) ((F(IJ,K,M,1),K=1,KL),M=1,ML)
      ENDDO
      WRITE(*,*) ' '
      WRITE(*,*) ' NTIME = ', NTIME
      WRITE(*,*) ' CDATE1 = ', CDATE1
      WRITE(*,*) ' NBOUNC= ', NBOUNC
      WRITE(*,*) ' KL, ML = ', KL, ML
      WRITE(*,*) ' XANG ...',XANG, XFRE, TH0, FR1, CO, XBOU, XDELC
      WRITE(*,*) ' EMEAN(1,1), FMEAN(1,1)', EMEAN(1,1), FMEAN(1,1)
      WRITE(*,*) ' F(1,1,1,1)', F(1,1,1,1)
      ntime = ntime+1
      GOTO 999
      
      ! Summary
1000  CONTINUE
      GOTO 1030

      ! Error handling
1001  continue
      write (*,*) "ERROR: readbound_ecwam: Open error unit ", IU01
      goto 1030

      ! Close files
1030  continue
      close (IU01)

      END ! PROGRAM readbound_ecwam
