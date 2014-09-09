      PROGRAM ecwambound

! Create ECMWF-WAM spectral boundary file from sequential binary spectra 
! (DNMI format).
!
! Note: only rotated and unrotated lat-lon projections accepted.
!
! USAGE: ecwambound inspec.seq wambound.inp outspec.wam 
!
! Input  file  - inspec.seq:  sequential binary 2D spectra - DNMI format
! Input file   - wambound.inp: boundary and resolution parameters
! Output file  - outspec.wam: ECMWF-WAM 2D spectral binary
! Output file2 - outspec.wam: ECMWF-WAM 2D spectral binary
! Example:
! 
!
! 2012-01-27 - oyvind.breivik@met.no
! 
! guest@overcast:/metno/routine/metop/wave/cray/wam/inp/preproc.inp_8km 
! contains the name of the fine-run input boundary file:
! ../runwam/CBOUT_50km.dat
! Temporal resolution of boundary spectra is 1 h.
!
! wambound.inp: - serves three purposes:
! (1) boundary parameters,
!     "coarse grid" resolution parameters,
!     projection parameters (xcen, ycen)
!     to wamboundinp.f
! (2) boundary parameters,
!     "coarse grid" resolution parameters
!     to wambound.f
! (3) dtresh and maxneighbours passed through to invdist.f
!
! Requires rotspher.f, specio.f, julian.f
!
! Compile:
! Linux: gfortran -O2 ecwambound.f rotspher.f specio.f julian.f -o ecwambound
! Njord: xlf -O2 ecwambound.f rotspher.f specio.f julian.f -o ecwambound
! ECMWF Linux: pgf90 -O2 ecwambound.f rotspher.f specio.f julian.f -o ecwambound
! CCC test: pgf90 -O2 ecwambound_test_nested_grids_info.f rotspher.f specio.f julian.f -o ecwambound_test

      IMPLICIT NONE

      ! Functions called
      INTEGER iargc
      INTEGER julday
      REAL degrees

      ! Constants
      INTEGER MAXANG, MAXFRE, MAXARR
      INTEGER IU02, IU03, IU21, IU25, IDUM
      REAL RDUM

      PARAMETER ( MAXANG=50, MAXFRE=51, MAXARR=10000 )
      PARAMETER ( IU02=2, IU03=3, IU21=21, IU25=25, IDUM = 999 )
      PARAMETER ( RDUM = 999.9 ) 

      ! Vars
      INTEGER fchh, fchh0, fchhfirst, hh, hh0, jd0, jd
      INTEGER immdd0, immdd, ihhmi0, ihhmi
      INTEGER idelprc, idtpro, yyyy, mm, dd
      INTEGER ntime, nskip, nang, nfre, narg
      INTEGER i, k, ij, m, nx, ny, idgrid, idummy
      INTEGER ierr
      INTEGER nbounc
      REAL xang, xfre, xbou, xdelc
      REAL th0, co, co1, delt25, fd, rdummy
      REAL dlamac, dphiac, amosoc, amonoc, amoeac, amowec
      REAL pi, rad, deg
      REAL lonrad, latrad, rlon, rlat
      REAL xcen, ycen
      REAL emean, thq, fmean
      REAL delth
      REAL temp
      REAL si, ci, sih, cih

      LOGICAL lreal, lswan

      ! Arrays
      INTEGER ident(20)
      REAL    blngc(MAXARR), blatc(MAXARR)
      REAL    slon(MAXARR), slat(MAXARR)
      REAL    sp(MAXANG,MAXFRE)
      REAL    th(MAXANG), dfim(MAXFRE), fr(MAXFRE)

      ! Strings
      CHARACTER*256 fin, fout, foutnest, flist, str
      CHARACTER*14 cdate

      ! Initialize
      lreal = .FALSE.
      lswan = .FALSE. ! CCC for the time being

      nbounc=0
      ierr=0
      ntime=0
      nskip=0
      pi=4.0*atan(1.0)
      rad=pi/180.0
      deg=180.0/pi

      ! Command line
      narg=iargc()
      if (narg < 3) then
         write (*,*) 
     2   "USAGE: ecwambound inspec.seq wambound.inp outspec.wam
     2wam_nested_grids_info"
         write (*,*)
     2   " Create WAM spectral boundary file from sequential binary"
         write (*,*) 
     2   " spectra (DNMI format)"
         stop
      endif

      ! Open binary sequential input file with 2D spectra
      call getarg(1,fin)
      open (IU25, file=fin, access="sequential", status="old",
     2      form="unformatted", err=1025)

      ! Open ASCII input file
      call getarg(2,flist)
      open (IU21, file=flist, form="formatted", status="old",
     2      err=1021)

      ! Open WAM boundary binary output file
      call getarg(3,fout)
      open (IU02, file=fout, form="unformatted", err=1009)

      ! Open WAM nested grids binary output file
      call getarg(4,foutnest)
      open (IU03, file=foutnest, form="unformatted", err=1009)

      !!! Read input list !!!

      ! Read grid.sph identificator line for coarse grid boundary
      read (IU21,"(a)",err=1021) str
      read (str(10:100),*) 
     2 idgrid, nx, ny, amowec, amosoc, dlamac, dphiac, xcen, ycen
      read (IU21,*,err=1021) rdummy, idummy, fr(1), fchhfirst

      ! Compute ROTATED eastern and northern boundaries [deg]
      amoeac = amowec+(nx-1)*dlamac
      amonoc = amosoc+(ny-1)*dphiac


!!! Initialize spectral parameters from first spectrum !!!

      ! Read first spectrum
      call getspec(IU25,ident,MAXANG,MAXFRE,sp,ierr)
      if (ierr /= 0) goto 1010

      nang = ident(10)
      nfre = ident(11)

      ! Convert analysis date to Julian days
      jd0 = julday(ident(1),ident(2)/100,mod(ident(2),100))
      hh0 = ident(3)/100

      delth = 360.0/real(nang)

      ! Frequency array
      co = real(ident(13))/100.0
      do m = 2, nfre
         fr(m) = fr(m-1)*co
      enddo
      delt25 = fr(nfre)/4.0*delth*rad

      ! Delta-frequency array
      co1 = 0.5*(co-1.0)*delth*rad
      dfim(1)= co1*fr(1)
      do m = 2, nfre-1
         dfim(m) = co1*(fr(m)+fr(m-1))
      enddo
      dfim(nfre) = co1*fr(nfre-1)

      ! Direction array
      th0 = real(ident(19))/100.0
      do k = 1, nang
         th(k)= rad*(real(k-1)*delth+th0)
      enddo

!!! Count spectral locations nbounc (number of wet boundary points) !!!

      fchh0 = ident(4) ! first forecast time
      fchh = fchh0
      ihhmi0 = ident(3)
      ihhmi = ihhmi0
      immdd0 = ident(2)
      immdd = immdd0
      
      do while ( ierr==0 .AND. fchh==fchh0 .AND. ihhmi==ihhmi0 .AND. 
     +           immdd==immdd0 )
         nbounc = nbounc+1
         slat(nbounc) = degrees(ident(5),ident(9))
         slon(nbounc) = degrees(ident(6),ident(17))
         latrad = slat(nbounc)*rad
         lonrad = slon(nbounc)*rad ! CCC slon, slat not needed

         ! Compute rotated spherical coordinates
         call sphrot(lonrad,latrad,rlon,rlat,1,xcen*rad,ycen*rad,1)
         ! Radians to degrees
         blngc(nbounc)=rlon*deg
         blatc(nbounc)=rlat*deg

         call getspec(IU25,ident,MAXANG,MAXFRE,sp,ierr)
         fchh = ident(4)
         ihhmi = ident(3)
         immdd = ident(2)

      enddo ! while

      idelprc = (fchh-fchh0)*3600 ! timestep [s]

      ! Rewind spectral input file
      ierr = 0
      rewind IU25


!!! Write nested_grids_info (inspired by mbounc.F) !!!
      
      ! Information about number of nested grids in "coarse" model
      write (IU03) 1       ! gbounc, number of nested grids in "coarse" model
      write (IU03) nbounc  ! IPOGBN(I),I=1,GBOUNC - the number of boundary points in each nested subdomain

      ! Information about our specific fine grid appearance in the "coarse" grid
      write (IU03) nbounc ! was NBOUNC-INI+1
      write (IU03) (1,i=1,nbounc) ! the block number of boundary points - always 1 these days
      write (IU03) (i,i=1,nbounc) ! fake indices
      write (IU03) dlamac, dphiac, amosoc, amonoc, amoeac, amowec,
     &     (blngc(i),i=1,nbounc), (blatc(i),i=1,nbounc)

!!! Write header !!!

      xang = float(nang)
      xfre = float(nfre)
      xbou = float(nbounc)
      xdelc = float(idelprc)
      if (lreal .or. lswan) then
         write (IU02) real(xang), real(xfre), real(th0), 
     &   real(fr(1)), real(co), real(xbou), real(xdelc)
      else
         write (IU02) xang, xfre, th0, fr(1), co, xbou, xdelc
      endif

!!! Time loop !!!

      do while ( ierr==0 )

         ntime = ntime+1

       !!! Loop over spectral locations !!!

         do ij = 1, nbounc

          !!! Read spectrum !!!

            call getspec(IU25,ident,MAXANG,MAXFRE,sp,ierr)
            if (ierr /= 0) GOTO 1000

          !!! Skip records before forecast time fchhfirst

            fchh = ident(4)
            if (fchh < fchhfirst) then 
               if  (ij == 1) then
                  nskip = nskip+1
               endif
               GOTO 110
            endif

          !!! Compute date-string YYYYMMDDHHMISS !!!

            if (ij == 1) then
               jd = jd0+(hh0+fchh)/24
               hh = mod(hh0+fchh,24)
               ! Negative lead time?
               if (hh < 0) then
                  hh = hh+24
                  jd = jd-1
               endif
               ! Convert Julian day to civil date
               call caldat(jd,yyyy,mm,dd)
               write (cdate,'(i4.4,5i2.2)') yyyy, mm, dd, hh, 0, 0
            endif

            !rlon = blngc(ij)
            !rlat = blatc(ij)

          !!! Compute emean - see SEMEAN !!!

            emean=0.0
            ! Loop through all spectral bins
            do m = 1,nfre
               temp = 0.0
               do k = 1,nang
                  temp = temp+sp(k,m)
               enddo
               emean = emean+temp*dfim(m)
            end do     
            ! Compute tail energy
            emean = emean+delt25*temp

          !!! Compute mean direction - see STHQ !!!

            si = 0.0
            ci = 0.0
            ! Integrate over frequencies and directions.
            do k = 1,NANG
               cih = cos(th(k))
               sih = sin(th(k))
               temp = 0.0
               do m = 1,nfre
                  temp = temp+sp(k,m)*dfim(m)
               enddo
               si = si+sih*temp
               ci = ci+cih*temp
            enddo

            if (ci == 0.0) ci = 0.1E-33
            thq = atan2(si,ci)
            if (thq < 0.0) thq = thq + 2*pi


          !!! Compute mean frequency - in deep water - see FEMEAN !!!

            fmean = 0.0
            do m = 1, nfre
               fd = dfim(m)/fr(m)
               temp = 0.0
               do k = 1,nang
                  temp = temp+sp(k,m)
               enddo     
               fmean = fmean+fd*temp
            enddo

            ! Add tail energy.
            fmean = fmean+delt25*temp
            fmean = emean/fmean

            if (lreal) then
               write (IU02) real(blngc(ij)), real(blatc(ij)), cdate,
     &               emean, thq, fmean
               write (IU02) (((sp(k,m)), k=1,nang), m=1,nfre)
            elseif (lswan) then
               write (IU02) real(blngc(ij)), real(blatc(ij)),
     &               cdate(3:14), emean, thq, fmean
               write (IU02) (((sp(k,m)), k=1,nang), m=1,nfre)
            else
               write (IU02) blngc(ij), blatc(ij), cdate,
     &               emean, thq, fmean
               write (IU02) (((sp(k,m)), k=1,nang), m=1,nfre)
            endif

110         CONTINUE

         enddo ! ij

      enddo ! while ntime

!!! End time loop !!!


      ! Summary
1000  CONTINUE
      ! Adjust time counter
      ntime = ntime-1
      if (nbounc==1) then
         write (*,*) "wambound: wrote one spectrum at", ntime-nskip,
     2               " times"
      elseif (nbounc>1) then
         write (*,*) "wambound: wrote", nbounc," spectra at",
     2               ntime-nskip," times"
      endif
      write (*,*) "wambound: skipped", nskip,
     2" times before", fchhfirst,"H"

      goto 1030

      ! Error handling
1009  continue
      write (*,*) "ERROR: wambound: Open error unit ", IU02
      goto 1030
1010  continue
      write (*,*) "ERROR: wambound: Error in getspec, ierr == ", ierr
      goto 1030
1021  continue
      write (*,*) "ERROR: wambound: Open error unit ", IU21
1025  continue
      write (*,*) "ERROR: wambound: Open error unit ", IU25
      goto 1030
      goto 1030

      ! Close files
1030  continue
      close (IU25)
      close (IU02)
      close (IU03)

      END ! PROGRAM wambound
