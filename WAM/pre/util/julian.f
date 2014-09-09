!!!!! FUNCTION julday !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Returns the Julian day number that begins at noon of the calendar date
! specified by month mm, day id, and year iyyy, all integers. Negative year
! signifies B.C.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION julday(iyyy,mm,id)

      IMPLICIT NONE

      ! Interface
      INTEGER iyyy, mm, id ! INTENT (IN)
      INTEGER julday       ! INTENT (OUT)

      ! Constant
      INTEGER IGREG
      PARAMETER ( IGREG=15+31*(10+12*1582) )

      ! Locals
      INTEGER ja, jm, jy

      ! Main
      jy = iyyy

      if (jy .EQ. 0) then
        print *, "JULDAY: There is no year zero"
        STOP
      endif

      if (jy .LT. 0) then
        jy = jy+1
      endif

      if (mm .GT. 2) then
        jm = mm+1
      else
        jy = jy-1
        jm = mm+13
      endif

      julday = 365*jy+int(0.25*jy+2000.0)+int(30.6001*jm)+id+1718995

      if (id+31*(mm+12*iyyy) .GE. IGREG) then
        ja = int(0.01*jy)
        julday = julday+2-ja+int(0.25*ja)
      endif

      END ! FUNCTION julday

!!!!!! SUBROUTINE caldat !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inverse of function julday.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE caldat(julian,iyyy,mm,id)

      IMPLICIT NONE

      ! Interface
      INTEGER julian       ! INTENT (IN)
      INTEGER iyyy, mm, id ! INTENT (OUT)

      ! Constant
      INTEGER IGREG
      PARAMETER ( IGREG=2299161 )

      ! Locals
      INTEGER ja, jalpha, jb, jc, jd, je

      ! Main
      if (julian .GE. IGREG) then
        jalpha = int(((julian-1867216)-0.25)/36524.25)
        ja = julian + 1 + jalpha - int(0.25*jalpha)
      else
        ja = julian
      endif

      jb = ja+1524
      jc = int(6680.0+((jb-2439870)-122.1)/365.25)
      jd = 365*jc + int(0.25*jc)
      je = int((jb-jd)/30.6001)
      id = jb - jd - int(30.6001*je)
      mm = je-1

      if (mm .GT. 12) then
        mm = mm-12
      endif

      iyyy = jc-4715

      if (mm .GT. 2) then
        iyyy = iyyy-1
      endif

      if (iyyy .LE. 0) then
        iyyy = iyyy-1
      endif

      END ! SUBROUTINE caldat
