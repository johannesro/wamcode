!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE getwamspec(iunit, ident, ldata, fdata, ierr)
!     Read sequential binary 2D WAM spectrum
!     (a) field identification (one record)
!     (b) field data           (one record)
!
!     Input:
!        iunit: file unit no.
!        ldata: length of fdata (max field size)
!
!     Output:
!        ident(20): field identification
!        fdata(1:ldata): field (unscaled, according to identification)
!        ierr = 0:  read OK
!               1:  read error
!               2:  end of file
!
      IMPLICIT NONE

      ! Interface
      INTEGER   iunit, ldata ! intent in
      INTEGER   ident(20)    ! intent out
      REAL      fdata(ldata) ! intent out
      INTEGER   ierr         ! intent out

      ! Locals
      INTEGER MFSIZE
      PARAMETER (MFSIZE=24*25*4)

      REAL sc
      INTEGER i, nang, nfre, nword, iscale, ios
      INTEGER*2 idata2(MFSIZE), ident2(20)

      ierr=0

      ! Read ident(20) header
      read (iunit,iostat=ios,err=910,end=920) (ident2(i), i=1,20)
      do i = 1, 20
         ident(i) = ident2(i)
      enddo

      nang=ident(10)
      nfre=ident(11)
      nword=nang*nfre

      if ( (nword > ldata) .OR. (nword > mfsize) ) then
         write (6,*) ' **getwamspec** field length too big'
         write (6,*) ' **          ldata = ',ldata
         write (6,*) ' ** ident: ',(ident(i),i=1,11)
         write (6,*) ' **        ',(ident(i),i=12,20)
         write (6,*) ' ** nang,nfre,nang*nfre: ',nang,nfre,nword
         ierr=1
         goto 990
      endif

      ! Read data
      read (iunit,iostat=ios,err=910,end=920) (idata2(i), i=1,nword)

      ! Scale data
      iscale=ident(20)
      sc=10.**iscale
      do i = 1, nword
         fdata(i) = sc*idata2(i)
      enddo

      ierr=0
      goto 990

  910 ierr=1
      write(6,*) ' **getwamspec** read error. file, iostat: ',iunit,ios
      goto 990

  920 ierr=2
      goto 990

  990 CONTINUE
      END ! SUBROUTINE getwamspec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE putwamspec(iunit, ident, ldata, fdata, ierr)
!     Write sequential binary 2D WAM spectrum
!     (a) field identification (one record)
!     (b) field data           (one record)
!
!     Input:
!        iunit: file unit no.
!        ident(20): field identification
!        ldata: length of fdata (max field size)
!        fdata(1:ldata): field (unscaled, according to identification)
!
!     Output:
!        ierr = 0: write OK
!               1: write error
!
      IMPLICIT NONE

      ! Interface
      INTEGER iunit, ldata  ! intent in
      INTEGER ident(20)     ! intent in
      REAL    fdata(ldata)  ! intent in
      INTEGER ierr          ! intent out

      ! Locals
      INTEGER MFSIZE
      REAL UNDEF, UDEF
      PARAMETER (UNDEF=+1.0E+35)
      PARAMETER (UDEF=UNDEF*0.9)
      PARAMETER (MFSIZE=24*25*4)

      REAL sc, fmax
      INTEGER i, ifmax, nang, nfre, nword, iscale, ios
      INTEGER*2 idata2(MFSIZE), ident2(20)
      INTEGER idento(20)

      ierr=0

      do i = 1, 20
         idento(i) = ident(i)
      enddo

      nang=ident(10)
      nfre=ident(11)
      nword=nang*nfre

      if ( (nword > ldata) .OR. (nword > mfsize) ) then
         write (6,*) ' **putwamspec** field length too big'
         write (6,*) ' **          ldata = ',ldata
         write (6,*) ' ** ident: ',(ident(i),i=1,11)
         write (6,*) ' **        ',(ident(i),i=12,20)
         write (6,*) ' ** nang,nfre,nang*nfre: ',nang,nfre,nword
         ierr=1
         goto 990
      endif

      ! Autoscale?
      iscale=ident(20)
      if (iscale == -32767) then
         fmax = 0.0
         do i = 1, nword
            if (fdata(i) < UDEF) fmax=max(fmax,abs(fdata(i)))
         enddo
         if (fmax > 0.0) then
            iscale = log10(fmax)-4.0
            ifmax = nint(fmax*10.0**(-iscale))
            if (ifmax < 3278) then
               iscale=iscale-1
               ifmax=nint(fmax*10.**(-iscale))
            end if
            if (ifmax > 32766) iscale=iscale+1
            iscale=max(iscale,-30)
         else
            iscale=0
         endif
         idento(20) = iscale
      endif ! autoscale

      ! Scale data, handling undefs
      sc=10.0**(-iscale)
      do i = 1, nword
         if (fdata(i) < UDEF) then
            idata2(i) = nint(sc*fdata(i))
         else
            idata2(i) = -32767
         endif
      enddo
      
      do i = 1, 20
         ident2(i) = idento(i)
      enddo

      ! Write field identification and field data
      write (iunit, iostat=ios, err=900) (ident2(i), i=1,20)
      write (iunit, iostat=ios, err=900) (idata2(i), i=1,nword)

      ierr=0
      goto 990

  900 ierr=1
      write(6,*) ' **putwamspec** write error. file, iostat: ',iunit,ios
      goto 990

  990 CONTINUE
      END ! SUBROUTINE putwamspec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!! SUBROUTINE hiresdeg !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Returns degrees*100+minutes and a remainder*10000 for more accurate
! positioning
!
! Input:
!   r: decimal degrees
!
! Output:
!   idegmin: degrees*100+minutes
!   iremain: remainder*10000
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE hiresdeg(r, idegmin, iremain)

      IMPLICIT NONE

      ! Interface
      REAL    r                ! intent in
      INTEGER idegmin, iremain ! intent out

      ! Locals
      REAL x, d

      x = r*100.0
      idegmin = nint(x)
      d = x-real(idegmin)
      iremain = nint(d*10000.0)

      END ! SUBROUTINE hiresdeg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!! FUNCTION degrees !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Returns decimal degrees. Input is degrees*100+minutes and a remainder*1000000
! for more accurate positioning.
!
! Input:
!   idegmin: degrees*100+minutes
!   iremain: remainder*10000
!
! Output:
!   degrees: decimal degrees
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION degrees(idegmin, iremain)

      IMPLICIT NONE

      INTEGER idegmin, iremain ! intent in

      REAL degrees             ! intent out

      degrees = real(idegmin)/100.0 + real(iremain)/1000000.0

      END ! FUNCTION degrees
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
