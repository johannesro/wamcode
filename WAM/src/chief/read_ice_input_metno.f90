SUBROUTINE READ_ICE_INPUT

! ---------------------------------------------------------------------------- !
!                                                                              !
!   READ_ICE_INPUT - READ AN ICE MAP.                                          !
!                                                                              !
!     HEINZ GUNTHER    GKSS    JANUARY 1995                                    !
!     ERIK MYKLEBUST           NOVEMBER 2004                                   !
!                                                                              !
!     PURPOSE                                                                  !
!     -------                                                                  !
!                                                                              !
!       TO READ AN ICE MAP.                                                    !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!        FORMATTED READ FROM UNIT IU03, FILE03.                                !
!                                                                              !
!       FILE03, WHICH IS DEFINED IN THE USER INPUT, IS ASSIGNED TO IU03.       !
!       THE FILE MUST STORE:                                                   !
!           1. RECORD: THE ICE DATA HEADER.                                    !
!           2. RECORD: THE ICE DATA TIME.                                      !
!           FOLLOWING RECORDS: THE ICE DATA MATRIX.                            !
!           (RECORD 2 AND FOLLOWING RECORDS ARE REPEATED FOR NEXT ICE FIELD.)  !
!                                                                              !
!       THE FOLLOWING INFORMATION HAS TO BE TRANSFERRED BY                     !
!       SUB. SET_ICE_HEADER TO THE WAM_WIND_MODULE:                            !
!          INTEGER :: N_LON      !! NUMBER OF LONGITUDES IN GRID.              !
!          INTEGER :: N_LAT      !! NUMBER OF LATITUDES IN GRID.               !
!          REAL*8  :: D_LAT      !! LATITUDE INCREMENT OF GRID [DEG].          !
!          REAL*8  :: D_LON      !! LONGITUDE INCREMENT OF GRID [DEG].         !
!          REAL*8  :: SOUTH      !! SOUTH LATITUDE OF GRID [DEG].              !
!          REAL*8  :: NORTH      !! NORTH LATITUDE OF GRID [DEG].              !
!          REAL*8  :: WEST       !! WEST LONGITUDE OF GRID [DEG].              !
!          REAL*8  :: EAST       !! EAST LONGITUDE OF GRID [DEG].              !
!                                                                              !
!       THE FOLLOWING INFORMATION HAS TO BE TRANSFERRED BY                     !
!       SUB. SET_ICE TO THE WAM_ICE_MODULE:                                    !
!          CHARACTER (LEN=14) :: CDTICE     !! DATE/TIME OFICE FIELD.           !
!          INTEGER            :: ICE_GRID(:,:) !! ICE MAP.                     !
!                                                                              !
!       THE ICE MAP MUST BE ON A REGULAR LATITUDE-LONGITUDE GRID ARRANGED      !
!       FROM  WEST TO EAST AND FROM SOUTH TO NORTH, WHICH IS:                  !
!                       (1,  1 ) <==> SOUTH WEST                               !
!                       (NX, 1 ) <==> SOUTH EAST                               !
!                       (1,  NY) <==> NORTH WEST                               !
!                       (NX, NY) <==> NORTH EAST                               !
!                                                                              !
!        ICE POINTS MUST BE EQUAL TO 1 IN ICE_GRID.                            !
!        THE SAME SUB. SET_ICE CAN BE USED, IF ICE_GRID IS DEFINED AS          !
!        REAL OR LOGICAL. IN THIS CASE ICE POINTS ARE WHERE                    !
!           NINT(ICE_GRID) = 1  OR ICE_GRID =.TRUE. , RESPECTIVELY.            !
!                                                                              !
!     REFERENCES                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE                                                                   !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     EXTERNALS.                                                               !
!     ----------                                                               !

USE WAM_GENERAL_MODULE,   ONLY:  &
&       INCDATE                  !! UPDATES A DATE/TIME GROUP.

USE WAM_ICE_MODULE,       ONLY:  &
&       SET_ICE,                 & !! ICE INPUT INTO MODULE.
&       SET_ICE_HEADER             !! ICE INPUT HEADER INTO MODULE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FILE_MODULE,  ONLY: IU06, ITEST, IU03, FILE03

USE NETCDF

IMPLICIT NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER, PARAMETER :: KIND_D = 8

INTEGER, SAVE         :: NX_ICE        !! ICE MAP DIMENSIONS
INTEGER, SAVE         :: NY_ICE        !! ICE MAP DIMENSIONS
REAL (KIND=KIND_D)    :: D_LAT         !! LATITUDE INCREMENT OF GRID [DEG].
REAL (KIND=KIND_D)    :: D_LON         !! LONGITUDE INCREMENT OF GRID [DEG].
REAL (KIND=KIND_D)    :: SOUTH         !! SOUTH LATITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: NORTH         !! NORTH LATITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: WEST          !! WEST LONGITUDE OF GRID [DEG].
REAL (KIND=KIND_D)    :: EAST          !! EAST LONGITUDE OF GRID [DEG].
INTEGER, ALLOCATABLE  :: ICE_GRID(:,:) !! ICE MAP (ice coverage as 0 or 1 mask
REAL, ALLOCATABLE     :: ICE_conc(:,:) !! Ice concentration between 0 and 1 
CHARACTER (LEN=14), SAVE :: CDTICE         !! ICE DATE
CHARACTER (LEN=14)    :: firstdate      !! 
LOGICAL, SAVE         :: FRSTIME = .TRUE.

integer, SAVE  :: D_TIME
INTEGER        :: LENF
INTEGER, SAVE  :: ITIMES
INTEGER, SAVE  :: NTIMES 
INTEGER, SAVE  :: ICEid,ncid,timeid
! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. OPEN ICE DATA FILE.                                                   !
!        -------------------                                                   !

IF (FRSTIME) THEN
   LENF = LEN_TRIM(FILE03)
 !Open existing  netcdf file
   CALL Pf(nf90_open(path = FILE03(1:LENF), mode = nf90_nowrite, ncid = ncid))
   
! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. READ ICE DATA HEADER.                                                 !
!        ---------------------                                                 !

 
   CALL readheaderICEnc(ncid,NX_ICE,NY_ICE,D_LON, &
&                         D_LAT,D_TIME,NORTH,SOUTH,EAST,WEST,ICEid, &
&                         timeid, NTIMES,firstdate) 
   
   CDTICE=firstdate
   ITIMES=0
   CALL SET_ICE_HEADER (WEST, SOUTH, EAST, NORTH, D_LON, D_LAT)

   FRSTIME = .FALSE.
END IF


! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. READ ICE DATA.                                                        !
!        --------------                                                        !

IF (.NOT.ALLOCATED(ICE_GRID)) ALLOCATE (ICE_GRID(1:NX_ICE,1:NY_ICE))
IF (.NOT.ALLOCATED(ICE_conc)) ALLOCATE (ICE_conc(1:NX_ICE,1:NY_ICE))

IF (NTIMES.gt.1) THEN
   ITIMES = 1 + ITIMES
   CALL Pf(nf90_get_var(ncid, ICEid, ICE_conc(:, :),& 
&          start = (/ 1, 1, ITIMES /),&                      
&          count = (/NX_ICE , NY_ICE, 1/))) 
   CALL INCDATE(CDTICE, D_TIME)  
ELSE
   CALL Pf(nf90_get_var(ncid, ICEid, ICE_conc(:, :),&
&          start = (/ 1, 1/),&     
&          count = (/NX_ICE , NY_ICE/))) 
ENDIF

! convert ice concentration to ice cover
ICE_GRID = 0
where (ICE_conc .gt. 0.5) ICE_GRID = 1 


! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. TRANSFER ICE DATA TO MODULE.                                          !
!        ----------------------------                                          !

CALL SET_ICE (CDTICE, ICE_GRID)

! ---------------------------------------------------------------------------- !
!                                                                              !
!      5. WRITE TEST OUTPUT AND DEALLOCATE ARRAY.                              !
!         ---------------------------------------                              !

IF (ITEST.GE.1) THEN
   WRITE(IU06,*) ' '
   WRITE(IU06,*) '     SUB. READ_ICE_INPUT: FIELD FOR CDTICE = ',CDTICE
END IF   

IF (ALLOCATED(ICE_GRID)) DEALLOCATE (ICE_GRID)
IF (ALLOCATED(ICE_conc)) DEALLOCATE (ICE_conc)

RETURN

CONTAINS

! ----------------
SUBROUTINE  readheaderICEnc(ncid,N_LON,N_LAT,D_LON, &
&                         D_LAT,  D_TIME, NORTH, SOUTH,EAST, WEST,ICEid, &
&                         timeid, NTIMES,firstdate) 
USE NETCDF

implicit none


CHARACTER (LEN=21)      :: xvar     !! NAME OF ICE CONCENTRATION VARIABLE
CHARACTER (LEN=12)     :: lon_name, lat_name     !! NAME OF THE coordinate-VARIABLE IN NETCDF
CHARACTER (LEN=4)      :: yyyy
CHARACTER (LEN=2)      :: mm
CHARACTER (LEN=2)      :: dd
CHARACTER (LEN=2)      :: hh
CHARACTER (LEN=2)      :: mi
CHARACTER (LEN=2)      :: ss
CHARACTER (LEN=14)     :: firstdate
CHARACTER (LEN=40)     :: time_units
integer :: ncid                    !!ID of the netcdf file opened (in)
integer :: ICEid, rtimeid, timeid,lonid, latid  !!ID of the variables
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
integer*8  ::rtime, sectodate

!Name of ice varible
  xvar='sea_ice_concentration'
 
!Check if there is a variable named sea_ice_concentration   and get the ID   
      CALL Pf(nf90_inq_varid(ncid, xvar, ICEid))
 !Check if there is a variable named time and get the ID   
      CALL Pf(nf90_inq_varid(ncid,"time",timeid))
      write(IU06,*)'timeid: ',timeid
!Get the Number of dimensions  
      CALL Pf(nf90_inquire_variable(ncid, ICEid, ndims = numDims, natts = numAtts))

!Get the Dimensions
      CALL Pf(nf90_inquire_variable(ncid, ICEid, dimids = dimIDs(:numDims)))
      WRITE(IU06,*)'numDims for ICE =',numDims

!Assume first dimension is longitude,
      CALL Pf(nf90_inquire_dimension(ncid, dimIDs(1), name = lon_name, len = N_LON))
      CALL Pf(nf90_inq_varid(ncid,lon_name,lonid))
      WRITE(IU06,*)'N_LON=',N_LON,' lon_name=',lon_name

!Assume second dimension is latitude    
      CALL Pf(nf90_inquire_dimension(ncid, dimIDs(2), name = lat_name, len = N_LAT))
      CALL Pf(nf90_inq_varid(ncid,lat_name,latid))
      WRITE(IU06,*)'N_LAT=',N_LAT,' lat_name=',lat_name

!Get the number of time steps 
      IF (numDims.gt.2) THEN
         CALL Pf(nf90_inquire_dimension(ncid, dimIDs(numDims), len = NTIMES))
      ELSE
         NTIMES=1
      ENDIF

      WRITE(IU06,*)'N_LON N_LAT NTIMES',N_LON,N_LAT,NTIMES


      allocate(xlon(N_LON))
      write(IU06,*)'get longitudes:'
      CALL Pf(nf90_get_var(ncid,lonid,xlon))
      D_LON=xlon(2)-xlon(1)
      EAST=MAXVAL(xlon)
      WEST=MINVAL(xlon)
      write(IU06,*)'EAST WEST: ',EAST,WEST


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
      

      IF (NTIMES.gt.1) THEN   
         write(IU06,*)'timeid:',timeid

         write(IU06,*)'get time:'
         CALL Pf(nf90_get_var(ncid, timeid, time))
         write(IU06,*)'time: ',time

         write(IU06,*)'time(1),time(2)',time(1),time(2)
         DT=INT(time(2)-time(1))
         write(IU06,*)'D_TIME=',D_TIME
      ELSEIF (NTIMES.eq.1) THEN 
         !get the   forecast_reference_time 
         CALL Pf(nf90_inq_varid(ncid,"forecast_reference_time",rtimeid))
       
         CALL Pf(nf90_get_var(ncid, rtimeid, rtime))
         write(IU06,*)'In ICE file forecast_reference_time = ',rtime
         DT=0
      ENDIF

      WRITE(IU06,*)'EAST,WEST,NORTH,SOUTH, D_TIME',&
&      EAST,WEST,NORTH,SOUTH, D_TIME

!Get the units of the time

      time_units=''
      CALL Pf(nf90_get_att(ncid, timeid, 'units', time_units))     
      write(IU06,*)'time_units: ',time_units 
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

      IF (NTIMES.gt.1) THEN 
         sectodate= time(1)*fact-D_TIME 
      ELSEIF (NTIMES.eq.1) THEN
          sectodate= rtime*fact
      ENDIF

      WRITE(IU06,*)'sectodate',sectodate
      CALL INCDATELONGAGO (firstdate, sectodate)
      WRITE(IU06,*)'firstdate=', firstdate
      WRITE(IU06,*)'D_TIME=', D_TIME

         !days since 2000-01-01 00:00:00"
         !12345678901234567890
         !seconds since 2013-09-09 00:00:00"
         !123456789012345678901234567890123
         !hours since 1900-01-01 00:00:0.0" 
         !12345678901234567890123456789012



end SUBROUTINE readheaderICEnc
!################################################
!# 	pf      write a NetCDF-error message    #
!#	en	NetCDF error number             #
!################################################

SUBROUTINE pf(en)
use netcdf
integer :: en
if (en/=0) write (iu06,*) 'IN READ_ICE_INPUT err: ', NF90_STRERROR(en)
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


END SUBROUTINE READ_ICE_INPUT
