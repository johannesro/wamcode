PROGRAM make_spec_netcdf_metno
! ---------------------------------------------------------------------------- !
!                                                                              !
!     ROTATES THE SPECTRA TO THE TRUE NORTH AND PRINT THEM IN A NETCDF         !
!                                                                              !
!      H. GUNTHER     GKSS/ECMWF  DECEMBER 1989                                !
!                     HZG         DECEMBER 2010      RE-ORGANISED              !
!      A. CARRASCO    METNO       SEPTEMBER 2014                               !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!        POSTPROCESSING OF WAM MODEL SPECTRA OUTPUT.                           !
!                                                                              !
!                                                                              !
!     INTERFACE.                                                               !
!     ----------                                                               !
!                                                                              !
!          IU01    INPUT UNIT OF SPECTRA FILE.                                !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THIS PROGRAM READS THE WAM MODEL SPECTRA OUTPUTS AND EXTRACTS          !
!       SPECTRA AT SPECIFIED LOCATIONS AND TIMES. THEN IT ROTATES THE SPECRA   !
!       TO  THE TRUE NORTH AND WRITES THEM IN A NETCDF FILE.                   !
!       THIS PROGRAM ONLY WORKS FOR THE SPECIFIC ROTATION OF WAM10 AT METNO.   !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!      EXTERNALS.                                                              !
!     -----------                                                              !
USE WAM_COORDINATE_MODULE


USE WAM_GENERAL_MODULE, ONLY:  &
&       INCDATE,               &  !! UPDATES A DATE/TIME GROUP.
&       OPEN_FILE                 !! OPEN A FILE.
USE WAM_PRINT_MODULE,   ONLY:  &
&       PRINT_SPECTRA_USER        !! PRINT A PROTOCOLL OF USER SETTINGS.

use wam_spec_netcdf_metno_module

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE,  ONLY: IU05, FILE05, IU06, FILE06, ITEST
USE WAM_PRINT_MODULE, ONLY: CDATEA, CDATEE, IDELDO,                            &
&                           IU01, FILE01,  CDTFILE, IDFILE,                    &
&                           NOUT_S, PFLAG_S, TITL_S, NOUTT, COUTT,             &
&                           NOUTP, OUTLONG, OUTLAT, NAME, CFLAG_S,             &
&                           KL, ML, CO, FR, THETA,                             &
&                           SPEC_LAT, SPEC_LON, SPEC_DATE,                     &
&                           SPEC,                                              &
&                           SPEC_SEA,                                          &
&                           SPEC_SWELL


IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER     :: IFAIL,NFAIL     !! OPEN ERROR
LOGICAL     :: IEOF            !! END OF FILE ENCOUNTED IN SUB. READ_SPECTRUM
LOGICAL,SAVE :: FRSTIME = .TRUE.
INTEGER            :: I, tstep
character (len=16) :: xfile
CHARACTER (LEN=40) :: HEADER

REAL, PARAMETER :: PI = 3.1415927    !! PI.
REAL, PARAMETER :: ZPI = 2.*PI       !! 2.* PI.
REAL, PARAMETER :: DEG = 180./PI     !! COVERTION FROM RADIANS TO DEGREE
REAL, PARAMETER :: RAD = PI/180.     !! COVERTION FROM DEGREE TO RADIANS
REAL, PARAMETER    :: xcen  =  -40.000000  !! Longitude of center point in rotated sph. grid for WAM10 at METNO
REAL, PARAMETER    :: ycen  =   68.000000  !! Latitude of center point in rotated sph. grid for WAM10 at METNO
real, allocatable, dimension (:)     :: TRUELONG
real, allocatable, dimension (:)     :: TRUELAT

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITALISATION.                                                        !
!        --------------                                                        !

!     1.1 SET USER INPUT AND PROTOCOLL FILE NAMES.                             !

FILE05 = 'Spectra_User'
FILE06 = 'Spectra_Prot'

!     1.2  OPEN USER FILE AND READ USER INPUT.                                 !

OPEN (UNIT=IU06, FILE=FILE06, FORM='FORMATTED', STATUS="UNKNOWN")
CALL READ_SPECTRA_USER
CALL PRINT_SPECTRA_USER

if (.not.allocated(TRUELONG))  allocate(TRUELONG(1:NOUTP))
if (.not.allocated(TRUELAT))  allocate(TRUELAT(1:NOUTP))
CALL ROTATE_COORDINATES(TRUELONG,TRUELAT)


!     1.3 FIRST AND LAST OUTPUT DATE.                                         !!

IF (NOUTT.GT.0) THEN
   CDATEE = COUTT(1)
   CDATEA = COUTT(1)
   DO I = 1,NOUTT
      IF (COUTT(I).LT.CDATEA) CDATEA = COUTT(I)
      IF (COUTT(I).GT.CDATEE) CDATEE = COUTT(I)
   END DO
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. LOOP OVER INPUT FILES.                                                !
!        ----------------------                                                !
NFAIL=0
tstep = 0
FILES: DO

!     2.1 FETCH FILE.                                                          !

   CALL OPEN_FILE (IU06, IU01, FILE01, CDTFILE, 'OLD', IFAIL)
!   IF (IFAIL.NE.0) STOP
   
! tbruns 26.01.2012
   IF (IFAIL.EQ.0) THEN

!     2.2  LOOP OVER OUTPUT TIMES.                                             !

   TIMES: DO

!     2.2.1 READ IN ONE OUTPUT TIME.                                           !

      CALL READ_SPECTRA_FILE (IU01, IEOF)

!     2.2.2 END OF FILE ENCOUNTED?                                             !

      IF (IEOF) EXIT TIMES   !! END OF FILE ENCOUNTERED
      IF (ITEST.GT.0) THEN
         WRITE (IU06,*) 'SUB. READ_SPECTRA_FILE DONE'
         WRITE (IU06,*) 'NEXT OUTPUT DATE, SPEC_DATE, SPEC_LAT, SPEC_LON: ',   &
&                        CDATEA, SPEC_DATE, SPEC_LAT, SPEC_LON
      END IF

!     2.2.3 OUTPUT TIME FOUND?                                                 !

      IF (SPEC_DATE.LT.CDATEA) CYCLE TIMES
      DO WHILE (SPEC_DATE.GT.CDATEA)
         CALL NEXT_OUTPUT_TIME
         tstep = tstep + 1
         IF (CDATEA.GT.CDATEE) EXIT FILES
         IF (SPEC_DATE.LT.CDATEA) CYCLE TIMES
      END DO

!     2.2.4 OUTPUT LOCATION?                                                   !

      WRITE (IU06,*) ' '
      IF (FRSTIME) THEN

!*       open NetCDF file
!        ----------------

         xfile( 1: 3) = 'SPC'
         xfile( 4:13) = cdatea(1:10)
         xfile(14:16) = '.nc'
         write (iu06,*) ' +++'
         write (iu06,*) ' +++ NetCDF file has been opened - name is: ', trim(xfile)
      
         

         call wkncospc(xfile, cdatea, CFLAG_S,TRUELONG,TRUELAT)
         FRSTIME = .FALSE.
         WRITE (IU06,*) ' ' 
      END IF

       
      LOCATION: DO I = 1,NOUTP
         IF (MOD(OUTLONG(I)-SPEC_LON+2*M_S_PER,+M_S_PER).EQ.0 .AND.            &
&           OUTLAT(I).EQ.SPEC_LAT) THEN
            IF (PFLAG_S(1) .AND. CFLAG_S(1)) THEN
            call ROTATE_SPECTRUM(IU06, SPEC_DATE, SPEC_LON, SPEC_LAT,   &
&                  HEADER, FR, THETA, SPEC,I,tstep,1)
            END IF
            IF (PFLAG_S(2) .AND. CFLAG_S(2)) THEN
            call ROTATE_SPECTRUM(IU06, SPEC_DATE, SPEC_LON, SPEC_LAT,   &
&                  HEADER, FR, THETA, SPEC_SEA,I,tstep,2)
            END IF
            IF (PFLAG_S(3) .AND. CFLAG_S(3)) THEN
            call ROTATE_SPECTRUM(IU06, SPEC_DATE, SPEC_LON, SPEC_LAT,   &
&                  HEADER, FR, THETA, SPEC_SWELL,I,tstep,3)
            END IF

           
            CYCLE TIMES
         END IF
       END DO LOCATION

   END DO TIMES
  
   CLOSE (UNIT=IU01, STATUS='KEEP')   !! CLOSE OLD FILE
   ELSE
     NFAIL=NFAIL+1
     IF(NFAIL.gt.100) EXIT FILES
   ENDIF

   IF (CDATEA.EQ.CDATEE) EXIT FILES   !! ALL DONE?
   IF (IDFILE.GT.0) THEN
      CALL INCDATE (CDTFILE, IDFILE)  !! INCREMENT DATE FOR THE NEXT FILE.
      WRITE(IU01,*)'CDTFILE',CDTFILE,'CDATEA',CDATEA,'IDFILE',IDFILE
   ELSE
      EXIT FILES
   END IF
END DO FILES

call wknccspc                              !! close NetCDF file
999 stop

WRITE (*,*) ' PROGRAM PRINT_SPECTRA: ALL DONE'

! ---------------------------------------------------------------------------- !

CONTAINS
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   SUBROUTINE NEXT_OUTPUT_TIME

   CHARACTER (LEN=14) :: IHH

      IF (NOUTT.EQ.0) THEN
         CALL INCDATE (CDATEA,IDELDO)
      ELSE
         IHH = '99999999999999'
         DO I=1,NOUTT
            IF (COUTT(I).GT.CDATEA .AND. COUTT(I).LT.IHH) IHH = COUTT(I)
         END DO
         CDATEA = IHH
      END IF

   END SUBROUTINE NEXT_OUTPUT_TIME

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE ROTATE_COORDINATES(TRUELONG,TRUELAT)

REAL               ::oldlatdeg,latrad
REAL               ::oldlongdeg,longrad
REAL               :: tlong
REAL               :: tlat
REAL               :: TRUELONG(:)
REAL               ::TRUELAT(:)
INTEGER            :: I

DO I = 1,NOUTP
   !LAT and LONG come in seconds*100
   oldlatdeg = (REAL(OUTLAT(I))/(60*60*100)) !LAT converted to degrees
   latrad = oldlatdeg*RAD
   oldlongdeg = REAL(OUTLONG(I))/(60*60*100) !LON converted to degrees
   longrad = oldlongdeg*RAD
   ! Compute coordinates in the new geografic grid (or non rotated grid)
   call SPHROT(longrad,latrad,tlong,tlat,xcen*RAD,ycen*RAD,-1)
   TRUELONG(I)=tlong*DEG
   TRUELAT(I)=tlat*DEG
ENDDO

END SUBROUTINE  ROTATE_COORDINATES


! ---------------------------------------------------------------------------- !

!
 SUBROUTINE ROTATE_SPECTRUM(IUNIT, CDATE, LONG, LAT, TITL, FR, THETA, SPEC, &
&                            IPOS,tstep,ivariable)


!     INTERFACE VARIABLES.
!     --------------------



INTEGER,            INTENT(IN) :: IUNIT     !! OUTPUT UNIT.
CHARACTER (LEN=14), INTENT(IN) :: CDATE     !! DATE OF SPECTRUM (YYYYMMDDHHMMSS).
INTEGER,            INTENT(IN) :: LONG      !! LONGITUDE OF SPECTRUM (DEGREE).
INTEGER,            INTENT(IN) :: LAT       !! LATITUDE OF SPECTRUM (DEGREE).
CHARACTER (LEN=40), INTENT(IN) :: TITL      !! TITLE.
REAL,               INTENT(IN) :: SPEC(:,:) !! SPECTRUM.
REAL,               INTENT(IN) :: FR(:)     !! FREQUENCY ARRAY IN HERTZ.
REAL,               INTENT(IN) :: THETA(:)  !! DIRECTION ARRAY IN RAD.


! ---------------------------------------------------------------------------- !
! 
!     LOCAL VARIABLES.
!     ----------------

INTEGER, PARAMETER :: IPDIR = 12   !! NUMBER OF DIRECTIONS PRINTED PER LINE.

INTEGER            :: KL              !! NUMBER OF DIRECTIONS.
INTEGER            :: ML              !! NUMBER OF FREQUENCIES.
INTEGER            :: IPE, IP, M, LEN
REAL               :: DELTH
REAL               :: SPECR(SIZE(THETA),SIZE(FR)) !!SPECT. ROTATED TO TRUE NORTH
REAL               :: ODSPEC(SIZE(FR))  !! 1-D SPECTRUM.  (M*M/HERTZ)

REAL               :: ANG(SIZE(THETA))     !! DIRECTIONS IN DEGREE.
REAL               ::frac
INTEGER            ::k2
INTEGER            ::k
INTEGER            ::i,IPOS,tstep
INTEGER            ::kplus
INTEGER            ::dk,ivariable
REAL               ::longrad
REAL               ::oldlongdeg
REAL               ::latrad
REAL               ::oldlatdeg
REAL               ::truelong
REAL               ::truelat
REAL               ::angle
REAL               ::delang, delang2
REAL               ::delta

! ---------------------------------------------------------------------------- !
! 
!     1. INITIALISE DIRECTIONS.
!        ----------------------

KL = SIZE(THETA)
ML = SIZE(FR)
DELTH = ZPI/REAL(KL)

!     1.1 ROTATE THE SPECTRA TO THE TRUE NORTH

!LAT and LONG come in seconds*100
oldlatdeg = (REAL(LAT)/(60*60*100)) !LAT converted to degrees
latrad = oldlatdeg*RAD
oldlongdeg = REAL(LONG)/(60*60*100) !LON converted to degrees
longrad = oldlongdeg*RAD
! Local distortion angle in old grid
call ANCORR(delta,oldlongdeg,oldlatdeg,xcen,ycen)
delang  = delta  

write(IUNIT,*)'COORDINATES IN ROTATED GRID  LONG = ',oldlongdeg,'LAT =', oldlatdeg
write(IUNIT,*)'ROTATED ANGLE  =     ',delta


! Compute coordinates in the new geografic grid (or non rotated grid)
call SPHROT(longrad,latrad,truelong,truelat,xcen*RAD,ycen*RAD,-1)
truelong=truelong*DEG
truelat=truelat*DEG
write(IUNIT,*)'TRUE COORDINATES   LONG = ',truelong,'LAT =', truelat
write(IUNIT,*)'ROTATED ANGLE  =     ',delta


! Local distortion angle in new grid geographic grid
call ANCORR(delta,truelong,truelat,0.0,0.0)
delang2  = delta 
 

! Compute angle between old (rotated) and new (no-rotated) grid
! Angle correction from rotated meridian to true meridian
angle = ang360(delang2-delang)
WRITE(IUNIT,*)'NEW ANGLE              = ',angle


! Compute gap measured in directional bins
dk = int(angle/DELTH)
frac = 1.0-mod(angle,DELTH)/DELTH

! Rotate spectra to true north
do k2 = 1, KL
   k = mod(k2+dk-1, KL) + 1              
   kplus = mod(k, KL) + 1
   do m = 1, ML
      SPECR(k2,m) = frac*SPEC(k,m)+(1.0-frac)*SPEC(kplus,m)
   enddo
enddo

!writes as text (ascii)
!WRITE(IUNIT,*)' SPECR',SPECR

!writes in a netcdf file
call wkncwspc (ivariable,SPECR,cdatea,tstep,IPOS)

! ---------------------------------------------------------------------------- !
!
!     2. COMPUTE 1-D SPECTRUM.
!        ---------------------
!do k2 = 1, KL
!   ANG(k2) = ang360(THETA(k2)*DEG - angle)
!enddo

!ANG = THETA*DEG

!ODSPEC = SUM(SPECR, DIM=1)
!ODSPEC = ODSPEC*DELTH





END SUBROUTINE ROTATE_SPECTRUM

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
SUBROUTINE ANCORR(delta,rlon,rlat,xcen,ycen)
!
!Direction angle correction, from rotated Wam grid to geographic
!
!CODED BY: J.E.HAUGEN, A.FOSS   - DNMI/R&D  APRIL 1998
!
integer IU
real delta,rlon,rlat,xcen,ycen
real :: zpir18,xca,yca,zsyca,zcyca
real :: x2,y2,x3,y3
real*8 :: zsxsph,zcxsph,zsysph,zcysph
real*8 :: zxmxc, zsxmxc, zcxmxc,zsxrot
real*8 :: zcxrot,zsyrot,zcyrot
reaL*8 :: za1,za2,za3,za4
real dda,ua,va,u,v,dd

      zpir18 = 2.0*asin(1.0)/180.
!


      xca = xcen*zpir18
      yca = ycen*zpir18
      x3  = rlon*zpir18
      y3  = rlat*zpir18
!     CAREFUL the arguments are different than the routine of ØB
      call sphrot(x3,y3,x2,y2,xca,yca,-1)
!
      zsyca = sin(yca)
      zcyca = cos(yca)
!
      zsxsph = sin(x2)
      zcxsph = cos(x2)
      zsysph = sin(y2)
      zcysph = cos(y2)
      zxmxc  = x2 - xca
      zsxmxc = sin(zxmxc)
      zcxmxc = cos(zxmxc)
      zsxrot = sin(x3)
      zcxrot = cos(x3)
      zsyrot = sin(y3)
      zcyrot = cos(y3)
      za1 = zcxmxc*zcxrot + zcyca*zsxmxc*zsxrot
      za2 = zcyca*zsxmxc*zcxrot*zsyrot + zsyca*zsxmxc*zcyrot
      za2 = za2 - zcxmxc*zsxrot*zsyrot
      za3 =-zsyca*zsxrot/zcysph
      za4 = (zcyca*zcyrot - zsyca*zcxrot*zsyrot)/zcysph
!
      dda = 45.
      ua  = -1.
      va  = -1.
!
      u = za1*ua + za2*va
      v = za3*ua + za4*va
!
      dd=270.-atan2(v,u)/zpir18
      if(dd.gt.+180.) dd=dd-360.
      if(dd.lt.-180.) dd=dd+360.
!
      delta=dd-dda
      if(delta.gt.+180.) delta=delta-360.
      if(delta.lt.-180.) delta=delta+360.
!
END SUBROUTINE ANCORR

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE SPHROT(xin,yin,xout,yout,xcen,ycen,icall)
! conversion between spherical (xsph,ysph) and spherical rotated
!  (xrot,yrot) coordinates. (xcen,ycen) is the position of the
!  rotated equator/greenwich in terms of (longitude,latitude).
!  all values are given in radians. Øyvind Breivik 

integer icall
real xcen,ycen
real xout,yout,xin,yin
real xsph,ysph,xrot,yrot
real  zsycen,zcycen,zxmxc,zsxmxc,zcxmxc
real zsysph,zcysph,zsyrot,zcyrot
real zcxrot,zsxrot


zsycen = sin(ycen)
zcycen = cos(ycen)

if (icall.eq.1) then
!
!  compute spherical rotated coordinates as function of
!  spherical coordinates
!
      xsph = xin;   !xin=log
      ysph = yin;   !yin=lat

      zxmxc  = xsph - xcen
      zsxmxc = sin(zxmxc)
      zcxmxc = cos(zxmxc)
      zsysph = sin(ysph)
      zcysph = cos(ysph)
      zsyrot = zcycen*zsysph - zsycen*zcysph*zcxmxc
      zsyrot = max(zsyrot,-1.0)
      zsyrot = min(zsyrot,+1.0)
      yrot   = asin(zsyrot)
      zcyrot = cos(yrot)
      zcxrot = (zcycen*zcysph*zcxmxc + zsycen*zsysph)/zcyrot
      zcxrot = max(zcxrot,-1.0)
      zcxrot = min(zcxrot,+1.0)
      zsxrot = (zcysph*zsxmxc)/zcyrot
      xrot   = acos(zcxrot)

      if (zsxrot.lt.0.0) xrot = -xrot
      xout = xrot
      yout = yrot

elseif (icall.eq.-1) then
!
!  compute spherical coordinates as function of
!  spherical rotated coordinates
  
      xrot = xin  !xin=xrot
      yrot = yin  !yin=yrot     
      zsxrot = sin(xrot)
      zcxrot = cos(xrot)
      zsyrot = sin(yrot)
      zcyrot = cos(yrot)
      zsysph = zcycen*zsyrot + zsycen*zcyrot*zcxrot
      zsysph = max(zsysph,-1.0)
      zsysph = min(zsysph,+1.0)
      ysph   = asin(zsysph)
      zcysph = cos(ysph)
      zcxmxc = (zcycen*zcyrot*zcxrot-zsycen*zsyrot)/zcysph     
      zcxmxc = max(zcxmxc,-1.0)
      zcxmxc = min(zcxmxc,+1.0)
      zsxmxc = (zcyrot*zsxrot)/zcysph        
      zxmxc  = acos(zcxmxc)
     
      if (zsxmxc.lt.0.0) zxmxc = -zxmxc
      xsph = zxmxc + xcen

      xout = xsph
      yout = ysph

else
   
   write(*,'(1x,''invalid icall in sphrot'')')
   stop

endif

END SUBROUTINE SPHROT
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
REAL FUNCTION ang360 (ANG)
IMPLICIT NONE

      REAL ang    ! [deg]
      ang360 = mod(ang, 360.0) - (sign(1.0,ang)-1.0)*180.0
   
END FUNCTION ang360
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !




END PROGRAM make_spec_netcdf_metno
