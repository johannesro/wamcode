!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE getspec(iunit, ident, maxang, maxfre, fdata, ierr)
!     Read sequential binary 2D spectrum, DNMI style
!     This routine should replace the older getwamspec.
!
!     (a) field identification (one record)
!     (b) field data           (one record)
!
!     Input:
!        iunit: file unit no.
!        maxang, maxfre: max dimensions of fdata
!
!     Output:
!        ident(20): field identification
!        fdata: field (unscaled, according to identification)
!        ierr = 0:  read OK
!               1:  read error
!               2:  end of file
!
! 2006-11-16, Oyvind.Breivik@met.no
!
      IMPLICIT NONE

      ! Interface
      INTEGER   iunit, maxang, maxfre ! intent in
      INTEGER   ident(20)             ! intent out
      REAL      fdata(maxang, maxfre) ! intent out
      INTEGER   ierr                  ! intent out

      ! Vars
      REAL sc
      INTEGER i, j, k, nang, nfre, nword, iscale, ios
      INTEGER*2 idata2(maxang*maxfre), ident2(20)

      ierr=0

      ! Read ident(20) header
      read (iunit,iostat=ios,err=910,end=920) (ident2(i), i=1,20)
      do i = 1, 20
         ident(i) = ident2(i)
      enddo

      nang=ident(10)
      nfre=ident(11)
      nword=nang*nfre

      if (nword > maxang*maxfre) then
         write (6,*) ' **getspec** field length too big'
         write (6,*) ' **          maxang*maxfre = ',maxang*maxfre
         write (6,*) ' ** ident: ',(ident(i),i=1,11)
         write (6,*) ' **        ',(ident(i),i=12,20)
         write (6,*) ' ** nang,nfre,nang*nfre: ',nang,nfre,nword
         ierr=1
         goto 990
      endif

      ! Read data
      read (iunit,iostat=ios,err=910,end=920) (idata2(k), k=1,nword)

      ! Scale data
      iscale=ident(20)
      sc=10.0**iscale
      k=0
      do j=1,nfre
         do i=1,nang
            k=k+1
            fdata(i,j) = sc*idata2(k)
         enddo
      enddo

      ierr=0
      goto 990

  910 ierr=1
      write(6,*) ' **getspec** read error. file, iostat: ',iunit,ios
      goto 990

  920 ierr=2
      goto 990

  990 CONTINUE
      END ! SUBROUTINE getspec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE putspec(iunit, ident, maxang, maxfre, fdata, ierr)
!     Write sequential binary 2D spectrum, DNMI style
!     (a) field identification (one record)
!     (b) field data           (one record)
!
!     Input:
!        iunit: file unit no.
!        ident(20): field identification
!        maxang, maxfre: max dimensions of fdata
!        fdata: field (unscaled, according to identification)
!
!     Output:
!        ierr = 0: write OK
!               1: write error
!
      IMPLICIT NONE

      ! Interface
      INTEGER iunit, maxang, maxfre ! intent in
      INTEGER ident(20)             ! intent in
      REAL    fdata(maxang,maxfre)  ! intent in
      INTEGER ierr                  ! intent out

      ! Locals
      INTEGER MFSIZE
      REAL UNDEF, UDEF
      PARAMETER (UNDEF=+1.0E+35)
      PARAMETER (UDEF=UNDEF*0.9)

      REAL sc, fmax
      INTEGER i, j, k, ifmax, nang, nfre, nword, iscale, ios
      INTEGER*2 idata2(maxang*maxfre), ident2(20)
      INTEGER idento(20)

      ierr=0

      do k = 1, 20
         idento(k) = ident(k)
      enddo

      nang=ident(10)
      nfre=ident(11)
      nword=nang*nfre

      if (nword > maxang*maxfre) then
         write (6,*) ' **putspec** field length too big'
         write (6,*) ' **          maxang*maxfre = ',maxang*maxfre
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
         !do i = 1, nword
         k=0
         do j=1,nfre
            do i=1,nang
               k=k+1
               if (fdata(i,j) < UDEF) fmax=max(fmax,abs(fdata(i,j)))
            enddo
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
      k=0
      do j=1,nfre
         do i=1,nang
            k=k+1
            if (fdata(i,j) < UDEF) then
               idata2(k) = nint(sc*fdata(i,j))
            else
               idata2(k) = -32767
            endif
         enddo
      enddo
      
      do k = 1, 20
         ident2(k) = idento(k)
      enddo

      ! Write field identification and field data
      write (iunit, iostat=ios, err=900) (ident2(k), k=1,20)
      write (iunit, iostat=ios, err=900) (idata2(k), k=1,nword)

      ierr=0
      goto 990

  900 ierr=1
      write(6,*) ' **putspec** write error. file, iostat: ',iunit,ios
      goto 990

  990 CONTINUE
      END ! SUBROUTINE putspec
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
