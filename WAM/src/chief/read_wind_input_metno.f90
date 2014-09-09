SUBROUTINE READ_WIND_INPUT

! ---------------------------------------------------------------------------- !
!                                                                              !
!   READ_WIND_INPUT - ROUTINE TO READ WINDFIELDS.                              !
!                                                                              !
!     HEINZ GUNTHER    GKSS    JANUARY 2001                                    !
!     ANA CARRASCO             SEPTEMBER 2013                                  !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO READ A WIND FIELD AND TRANSFER IT TO THE WAM_WIND_MODULE.           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!        FORMATTED READ FROM UNIT IU01, FILE01.                                !
!                                                                              !
!       FILE01, WHICH IS DEFINED IN THE USER INPUT, IS ASSIGNED TO IU01.       !
!       THE FILE MUST STORE:                                                   !
!           1. RECORD: THE WIND DATA HEADER.                                   !
!           2. RECORD: THE WIND DATA TIME.                                     !
!           FOLLOWING RECORDS: THE WIND DATA MATRIX.                           !
!           (RECORD 2 AND FOLLOWING RECORDS ARE REPEATED FOR NEXT WIND FIELD.) !
!                                                                              !
!       THE FOLLOWING INFORMATION HAS TO BE TRANSFERRED BY                     !
!       SUB. SET_WIND_HEADER TO THE WAM_WIND_MODULE:                           !
!          INTEGER :: N_LON      !! NUMBER OF LONGITUDES IN GRID.              !
!          INTEGER :: N_LAT      !! NUMBER OF LATITUDES IN GRID.               !
!          REAL*8  :: D_LAT      !! LATITUDE INCREMENT OF GRID [DEG].          !
!          REAL*8  :: D_LON      !! LONGITUDE INCREMENT OF GRID [DEG].         !
!          REAL*8  :: SOUTH      !! SOUTH LATITUDE OF GRID [DEG].              !
!          REAL*8  :: NORTH      !! NORTH LATITUDE OF GRID [DEG].              !
!          REAL*8  :: WEST       !! WEST LONGITUDE OF GRID [DEG].              !
!          REAL*8  :: EAST       !! EAST LONGITUDE OF GRID [DEG].              !
!          INTEGER :: ICODE      !! WIND CODE: 1= USTAR; 2= USTRESS; 3= U10    !
!                                                                              !
!       THE FOLLOWING INFORMATION HAS TO BE TRANSFERRED BY                     !
!       SUB. SET_WIND_FIELD TO THE WAM_WIND_MODULE:                            !
!          CHARACTER (LEN=14) :: CDTWIR     !! DATE/TIME OF WIND FIELD.        !
!          REAL               :: U_MAP(:,:) !! U COMPONENT OF WIND MAP [M/S].  !
!          REAL               :: V_MAP(:,:) !! V COMPONENT OF WIND MAP [M/S].  !
!                                                                              !
!       THE WINDS MUST BE ON A REGULAR LATITUDE-LONGITUDE GRID ARRANGED        !
!       FROM  WEST TO EAST AND FROM SOUTH TO NORTH, WHICH IS:                  !
!       IN THE ARRAYS "U_MAP(I,K)" AND  "V_MAP(I,K)" THE CORNER POINTS ARE:    !
!                 (    1,    1 ) <==> SOUTH WEST                               !
!                 (N_LON,    1 ) <==> SOUTH EAST                               !
!                 (    1, N_LAT) <==> NORTH WEST                               !
!                 (N_LON, N_LAT) <==> NORTH EAST                               !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     EXTERNALS.                                                               !
!     ----------                                                               !



USE WAM_GENERAL_MODULE, ONLY:  &
&       INCDATE                  !! UPDATES A DATE/TIME GROUP.


USE WAM_WIND_MODULE,       ONLY: &
&       SET_WIND_HEADER,         & !! SETS WIND HEADER 
&       SET_WIND_FIELD,          & !! SETS WIND FIELD 
&       PRINT_WIND_STATUS          !! PRINTS WIND MODULE STATUS

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE, ONLY: IU06, ITEST, IU01, FILE01

USE NETCDF

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !


INTEGER, PARAMETER :: KIND_D = 8

INTEGER, SAVE         :: ICODE = 3  !! WIND CODE: 1= USTAR; 2= USTRESS; 3= U10
INTEGER, SAVE         :: N_LON      !! NUMBER OF LONGITUDES IN GRID.
INTEGER, SAVE         :: N_LAT      !! NUMBER OF LATITUDES IN GRID.
REAL (KIND=KIND_D)    :: D_LAT      !! LATITUDE INCREMENT OF GRID [DEG].
REAL (KIND=KIND_D)    :: D_LON      !! LONGITUDE INCREMENT OF GRID [DEG].
INTEGER, SAVE         :: D_TIME     !! TIME INCREMENT [SEC].
REAL (KIND=KIND_D)    :: SOUTH      !! SOUTH LATITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: NORTH      !! NORTH LATITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: WEST       !! WEST LONGITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: EAST       !! EAST LONGITUDE OF GRID [DEG].
REAL,    ALLOCATABLE  :: U_MAP(:,:) !! 1. COMPONENT OF WIND MAP [M/S].
REAL,    ALLOCATABLE  :: V_MAP(:,:) !! 2. COMPONENT OF WIND MAP [M/S].
REAL*8,  ALLOCATABLE  :: time(:)   !! TIME
CHARACTER (LEN=14)    :: CDTWIR     !! DATE/TIME OF WIND FIELD
CHARACTER (LEN=14)     :: firstdate

LOGICAL, SAVE  :: FRSTIME = .TRUE.
INTEGER        :: LENF

INTEGER, SAVE  :: ITIMES
INTEGER, SAVE  :: U10id, V10id, NTIMES 
INTEGER, SAVE  :: ncid,timeid
! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. FOR FIRST CALL: OPEN FILE AND READ HEADER.                            !
!        ------------------------------------------                            !




IF (FRSTIME) THEN
   LENF = LEN_TRIM(FILE01)
   ITIMES = 0
   
  
  !Open existing  netcdf file
   WRITE(IU06,*)'open wind file', FILE01(1:LENF)                                             !

   CALL Pf(nf90_open(path = FILE01(1:LENF), mode = nf90_nowrite, ncid = ncid))
   

   !Read the headers
 
   CALL readheaderWINDSnc(ncid,N_LON,N_LAT,D_LON, &
&                         D_LAT,  D_TIME, NORTH, SOUTH,EAST, WEST,U10id, &
&                         V10id,timeid, NTIMES,firstdate) 
   
   CDTWIR=firstdate
   CALL SET_WIND_HEADER (WEST=WEST,   SOUTH=SOUTH,    &
&                        EAST=EAST,   NORTH=NORTH,    &
&                        D_LON=D_LON, D_LAT=D_LAT,    &
&                        N_LON=N_LON, N_LAT=N_LAT,    &
&                        CODE=ICODE)
   IF (ITEST.GT.0) CALL PRINT_WIND_STATUS
   
   FRSTIME = .FALSE.
   
END IF
!

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. ALLOCATE WIND INPUT ARRAYS.                                           !
!        ---------------------------                                           !

IF (.NOT.ALLOCATED(U_MAP) ) ALLOCATE(U_MAP(N_LON,N_LAT))
IF (.NOT.ALLOCATED(V_MAP) ) ALLOCATE(V_MAP(N_LON,N_LAT))
IF (.NOT.ALLOCATED(time) )  ALLOCATE(time(NTIMES))     
!  

! ---------------------------------------------------------------------------- !
!                                                                              !
!    3. READ WIND FIELD.                                                       !
!       -------------------
  
   ITIMES = 1 + ITIMES
!   CALL Pf(nf90_get_var(ncid, timeid, time(:), start = (/ ITIMES-1 /), &
!&                          count = (/ ITIMES /))
   WRITE(IU06,*)'N_LON, N_LAT =',N_LON, N_LAT                                               !
 
   CALL Pf(nf90_get_var(ncid, V10id, V_MAP(:, :), start = (/ 1, 1, ITIMES /),&
&                       count = (/N_LON , N_LAT, 1/))) 
   CALL Pf(nf90_get_var(ncid, U10id, U_MAP(:, :), start = (/ 1, 1, ITIMES /),&
&                        count = (/N_LON , N_LAT, 1/)))
!    3. GET date    
  
   WRITE(IU06,*)'D_TIME =',D_TIME                                               !
   CALL INCDATE(CDTWIR, D_TIME)  
   

   CALL SET_WIND_FIELD (CDTWIR, U_MAP, V_MAP)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    4. WRITE TEST OUTPUT AND DEALLOCATE ARRAYS.                               !
!       ----------------------------------------                               !

IF (ITEST.GT.1) THEN
   WRITE(IU06,*) ' READ_WIND_INPUT -  WIND FIELD FOR THE CDTWIR = ', CDTWIR
   WRITE(IU06,'(1X,24F5.2)') U_MAP(1:MIN(24,N_LON),1:MIN(5,N_LAT))
   WRITE(IU06,*) ' '
   WRITE(IU06,'(1X,24F5.2)') V_MAP(1:MIN(24,N_LON),1:MIN(5,N_LAT))
END IF

  
DEALLOCATE(U_MAP)
DEALLOCATE(V_MAP)
DEALLOCATE(time)

RETURN

CONTAINS

SUBROUTINE  readheaderWINDSnc(ncid,N_LON,N_LAT,D_LON, &
&                         D_LAT,  D_TIME, NORTH, SOUTH,EAST, WEST,U10id, &
&                         V10id,timeid, NTIMES,firstdate) 
USE NETCDF

implicit none


CHARACTER (LEN=12)     :: x_varname, y_varname 
CHARACTER (LEN=12)     :: lon_name, lat_name     !! NAME OF THE SECOND-VARIABLE IN NETCDF (AC)
CHARACTER (LEN=4)      :: yyyy
CHARACTER (LEN=2)      :: mm
CHARACTER (LEN=2)      :: dd
CHARACTER (LEN=2)      :: hh
CHARACTER (LEN=2)      :: mi
CHARACTER (LEN=2)      :: ss
CHARACTER (LEN=14)     :: firstdate
CHARACTER (LEN=40)     :: time_units
integer :: ncid                    !!ID of the netcdf file opened (in)
integer :: U10id, V10id, timeid,lonid, latid  !!ID of the variables to be picked
integer :: numDims, numAtts
integer :: NTIMES
real*8,    dimension(:),allocatable ::xlon
real*8,    dimension(:),allocatable ::ylat
integer*8,    dimension(:),allocatable ::time
real*8  :: DT
integer ::N_LON,  N_LAT
real*8  ::D_LON, D_LAT, NORTH, SOUTH,EAST, WEST
integer ::D_TIME
integer, dimension(nf90_max_var_dims) :: dimIDs
integer    ::fact
integer*8  ::sectodate

!Name of the wind components 
  x_varname='x_wind'
  y_varname='y_wind'
 
!Check if there is a variable named  x_wind  and get the ID   
      CALL Pf(nf90_inq_varid(ncid, x_varname, U10id))
!Check if there is a variable named  y_wind  and get the ID   
      CALL Pf(nf90_inq_varid(ncid, y_varname, V10id))    
!Check if there is a variable named time and get the ID   
      CALL Pf(nf90_inq_varid(ncid,"time",timeid))
!Get the Number of dimensions of wind variables
      CALL Pf(nf90_inquire_variable(ncid, V10id, ndims = numDims, natts = numAtts))
      
      WRITE(IU06,*)'numDims=',numDims
!Get the Dimensions
      CALL Pf(nf90_inquire_variable(ncid, V10id, dimids = dimIDs(:numDims)))

!Assume first dimension is longitude,
      CALL Pf(nf90_inquire_dimension(ncid, dimIDs(1), name = lon_name, len = N_LON))
      CALL Pf(nf90_inq_varid(ncid,lon_name,lonid))
      WRITE(IU06,*)'N_LON=',N_LON,' lon_name=',lon_name

!Assume second dimension is latitude    
      CALL Pf(nf90_inquire_dimension(ncid, dimIDs(2), name = lat_name, len = N_LAT))
      CALL Pf(nf90_inq_varid(ncid,lat_name,latid))
      WRITE(IU06,*)'N_LAT=',N_LAT,' lat_name=',lat_name

      CALL Pf(nf90_inquire_dimension(ncid, dimIDs(numDims), len = NTIMES))
      WRITE(IU06,*)'NTIMES=',NTIMES

!WRITE(IU06,*)'N_LON N_LAT NTIMES',N_LON,N_LAT,NTIMES
!Get the values of the lons 
      allocate(xlon(N_LON))
      CALL Pf(nf90_get_var(ncid,lonid,xlon))
      D_LON=xlon(2)-xlon(1)
      EAST=MAXVAL(xlon)
      WEST=MINVAL(xlon)


!Get the values of the lats
      allocate(ylat(N_LAT))
      CALL Pf(nf90_get_var(ncid,latid,ylat))       
      D_LAT=ylat(2)-ylat(1)
      NORTH=MAXVAL(ylat)
!remove XXX values in  0.000000XXX 
      NORTH=1.0* INT(NORTH*100000)
      NORTH=NORTH/100000
      SOUTH=MINVAL(ylat)
     

!Get the values of the times
      allocate(time(NTIMES))
      CALL Pf(nf90_get_var(ncid, timeid, time))   
      write(IU06,*)'time(1),time(2)',time(1),time(2)
      DT=(time(2)-time(1))
WRITE(IU06,*)'EAST,WEST,NORTH,SOUTH, D_TIME',EAST,WEST,NORTH,SOUTH, DT
!Get the units of the time
      time_units=''
      CALL Pf(nf90_get_att(ncid, timeid, 'units', time_units))     
      write(IU06,*)'time_units: ', time_units 
      if (time_units(1:4).eq.'days') then

         yyyy=time_units(12:15)
         mm=time_units(17:18)
         dd=time_units(20:21)
         hh=time_units(23:24)
         mi=time_units(26:27)
         ss=time_units(29:30)
         fact=24*60*60 !because time must be in sec
      elseif (time_units(1:4).eq.'hour') then
         yyyy=time_units(13:16)
         mm=time_units(18:19)
         dd=time_units(21:22)
         hh=time_units(24:25)
         mi=time_units(27:28)
         !ss=time_units(30:32)
         ss='00' !!!OJO 
         fact=60*60 !because time must be in sec
    
      elseif (time_units(1:4).eq.'seco') then
         yyyy=time_units(15:18)
         mm=time_units(20:21)
         dd=time_units(23:24)
         hh=time_units(26:27)
         mi=time_units(29:30)
         ss=time_units(32:33)
         fact= 1
      endif
      WRITE(IU06,*)'fact',fact
      WRITE(IU06,*)'DT*fact',DT*fact
      D_TIME=NINT(DT*fact)
      WRITE(IU06,*)'D_TIME=',D_TIME
      firstdate=yyyy//mm//dd//hh//mi//ss
      WRITE(IU06,*)'firstdate=', firstdate
      
      sectodate= time(1)*fact-D_TIME 
     
      WRITE(IU06,*)'sectodate',sectodate
      CALL INCDATELONGAGO (firstdate, sectodate)
      WRITE(IU06,*)'firstdate=', firstdate

         !days since 2000-01-01 00:00:00"
         !12345678901234567890
         !seconds since 2013-09-09 00:00:00"
         !123456789012345678901234567890123
         !hours since 1900-01-01 00:00:0.0" 
         !12345678901234567890123456789012



end SUBROUTINE readheaderWINDSnc
!################################################
!# 	pf      write a NetCDF-error message    #
!#	en	NetCDF error number             #
!################################################

SUBROUTINE pf(en)
use netcdf
integer :: en
if (en/=0) write (iu06,*) 'IN READ_WIND_INPUT err: ', NF90_STRERROR(en)
end SUBROUTINE pf

SUBROUTINE INCDATELONGAGO (CDATE, ISHIFT)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   INCDATE - TO UPDATE DATE TIME GROUP                                        !
!                                                                              !
!     L. BERTOTTI, P.JANSSEN.                                                  !
!                                                                              !
!     H. GUNTHER   ECMWF  NOVEMBER 1989    NEGATIVE INCREMENTS.                !
!     H. GUNTHER   GKSS   NOVEMBER 2001    FT90 AND CENTURY AND SECONDS.       !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       UPDATING DATE TIME GROUP.                                              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
! 
!     INTERFACE VARIABLES.
!     --------------------

CHARACTER (LEN=14), INTENT(INOUT) :: CDATE  !! DATE TIME GROUP (YYYYMMDDHHMMSS).
INTEGER*8,          INTENT(IN)    :: ISHIFT !! TIME INCREMENT IN SECONDS.

! ---------------------------------------------------------------------------- !
! 
!     LOCAL VARIABLES.
!     ----------------

INTEGER, SAVE ::  MON(12) =(/31,28,31,30,31,30,31,31,30,31,30,31/)

INTEGER*8 ::  YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
INTEGER*8 ::  DELT, MDAY

! ---------------------------------------------------------------------------- !
! 
!   1.0 RETURN IF TIME INCREMENT IS ZERO.
!       ---------------------------------

DELT = ISHIFT
IF (ABS(DELT).EQ.0) RETURN

! ---------------------------------------------------------------------------- !
!                                                                              !
!*    2.0 SPLITE DATE TIME GROUP INTO SECONDS, MINUTE, HOUR, DAY, MONTH, YEAR. !
!         -------------------------------------------------------------------- !

READ (CDATE,'(I4,5I2)') YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2.0 ADD AND CHECK SECONDS.                                               !
!         ----------------------                                               !
!                                                                              !
!   2.1 IF SECONDS ARE BETWEEN 0 AND 60 RETURN
!       --------------------------------------

SECOND = SECOND + DELT
IF (SECOND.GE.0. .AND. SECOND.LT.60.) THEN
   WRITE(CDATE,'(I4.4,5I2.2)') YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
   RETURN
END IF

!*  2.2 NEW MIMUTES AND SECONDS.
!       ------------------------

DELT = MODULO(SECOND,60)
MINUTE = MINUTE +(SECOND-DELT)/60
SECOND = DELT

! ---------------------------------------------------------------------------- !
! 
!   3.0 CHECK MINUTES.
!       --------------

IF (MINUTE.GE.60) THEN
   HOUR = HOUR + MINUTE/60        !! MINUTES > 59 ==> NEW HOURS.
ELSE IF (MINUTE.LT.0) THEN
   HOUR = HOUR + (MINUTE-59)/60   !! MINUTES < 0  ==> NEW HOURS.
ELSE
   WRITE (CDATE,'(I4.4,5I2.2)') YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
   RETURN                         !! ALL DONE  ==>  RETURN
END IF 
MINUTE = MODULO(MINUTE,60)        !! NEW MINUTES.

! ---------------------------------------------------------------------------- !
! 
!   4.0 CHECK HOURS.
!       ------------

IF (HOUR.GE.24) THEN
   DAY =  DAY + HOUR/24           !! HOURS > 23 ==> NEW DAYS.
ELSE IF (HOUR.LT.0) THEN
   DAY =  DAY + (HOUR-23)/24      !! HOURS < 0  ==> NEW DAYS.
ELSE
   WRITE (CDATE,'(I4.4,5I2.2)') YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
   RETURN                         !! ALL DONE  ==>  RETURN
END IF 
HOUR = MODULO(HOUR,24)            !! NEW HOURS.

! ---------------------------------------------------------------------------- !
! 
!   5.0 CHECK DAYS.
!       ------------
! 
!   5.1 IF DAYS ARE GREATER THAN DAYS OF MONTH. NEW DAY AND MONTH AND YEAR.
!       -------------------------------------------------------------------

MDAY = MON(MONTH)
IF (MONTH.EQ.2 .AND. MOD(YEAR,4).EQ.0) THEN
   IF (MOD(YEAR,400).EQ.0 .OR. MOD(YEAR,100).NE.0) MDAY = 29
END IF

DO WHILE (DAY > MDAY)
   DAY =  DAY - MDAY
   MONTH = MONTH+1
   IF (MONTH.GE.13) THEN
      YEAR = YEAR+1
      MONTH = MONTH-12
   END IF
   MDAY = MON(MONTH)
   IF (MONTH.EQ.2 .AND. MOD(YEAR,4).EQ.0) THEN
      IF (MOD(YEAR,400).EQ.0 .OR. MOD(YEAR,100).NE.0) MDAY = 29
   END IF
END DO

!   5.2 IF DAYS ARE LESS THAN 1. NEW DAY AND MONTH AND YEAR.
!       ----------------------------------------------------

DO WHILE ( DAY < 1)
   MONTH = MONTH-1
   IF (MONTH.EQ.0) THEN
      MONTH = 12
      YEAR = YEAR-1
   END IF
   MDAY = MON(MONTH)
   IF (MONTH.EQ.2 .AND. MOD(YEAR,4).EQ.0) THEN
      IF (MOD(YEAR,400).EQ.0 .OR. MOD(YEAR,100).NE.0) MDAY = 29
   END IF
   DAY = DAY + MDAY
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     6.0 COMPOSE NEW DATE TIME GROUP.                                         !
!         ----------------------------                                         !

WRITE (CDATE,'(I4.4,5I2.2)') YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

END SUBROUTINE INCDATELONGAGO

END SUBROUTINE READ_WIND_INPUT
