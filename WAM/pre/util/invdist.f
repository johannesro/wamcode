      PROGRAM invdist
!
! Interpolate binary sequential WAM spectra (DNMI format) using inverse 
! distance weighting.
!
! Note: handles only (rotated) spherical co-ordinates.
!
! USAGE: invdist inspec.seq invdist.inp outspec.seq [xcen ycen -silent]
!
! Input file 1 - inspec.seq: sequential binary 2D spectra - DNMI format
! Input file 2 - invdist.inp: ASCII list of desired spectral locations
!                Note: first line of invdist.inp contains information on the
!                influence radius of neighbouring spectra, max number of 
!                neighbours to be used, and possibly spectral information used by 
!                other programs
! 
! Output file  - outspec.seq: same format as input file 1, spectra interpolated
!                to desired locations
! [xcen ycen]  - lon and lat of center pos in rotated spherical grid,
!                default to xcen=0.0, ycen=0.0 (unrotated grid) if omitted,
!                xcen=0.0, ycen=65.0 applies for the WAM 50 km grid.
! [-silent]    - do not report missing spectra
! 
! Requires:
!    sphere.f
!    rotspher.f
!    specio.f
!
! Compile IBM: 
! f90 -qfixed -O2 -C invdist.f sphere.f rotspher.f specio.f -o invdist
! Compile Linux: 
! pgf90 -O2 -C invdist.f sphere.f rotspher.f specio.f -o invdist
! or
! gfortran -O2 -C invdist.f sphere.f rotspher.f specio.f -o invdist
! make -f invdist.make
! 
!
! 2006m01d23, Oyvind.Breivik@met.no
! 2007-03-09, allows spectra to be set to zero if outside range
! 2008-11-03, silent option
! 2009-05-27, handles multiple analysis times
! 2009-05-28, bug fix: handles single time
!
! Test:
! invdist spk00.dat spk00int.seq invdist.inp 0.0 65.0 -silent

      IMPLICIT NONE

      ! Functions called
      REAL degrees
      REAL spheredist
      REAL ang360
      INTEGER iargc

      ! Constants
      INTEGER MAXANG, MAXFRE, MAXARR, MAXNEIGHBOURS
      INTEGER IU25, IU26, IU27
      REAL TOL, POW, DISTMAX, EPS, XCEN0, YCEN0

      PARAMETER ( MAXANG=50, MAXFRE=51, MAXARR=105000,MAXNEIGHBOURS=15 )
      PARAMETER ( IU25=25, IU26=26, IU27=27 )
      PARAMETER ( TOL=0.1, POW=2.0, DISTMAX=40000.0 )
      PARAMETER ( XCEN0=0.0, YCEN0=0.0 )

      ! Vars
      INTEGER i, j, k, l, m, ntime, nang, nfre
      INTEGER fchh, fchh0
      INTEGER immdd0, immdd, ihhmi0, ihhmi
      INTEGER nwish, nspec, narg, miss, neighbours
      INTEGER ierr
      REAL pi, rad, deg
      REAL th0, co, co1, delth, tmp, si, ci
      REAL dtresh, delta
      REAL lonrad, latrad, rlon, rlat, rlonrad, rlatrad
      REAL xcen, ycen, xcrad, ycrad
      REAL th, thq, thw
      REAL u10, ustar, u10eastw, u10northw, useastw, usnorthw
      CHARACTER c
      LOGICAL nullspec, silent

      ! Arrays
      INTEGER ident(20), identnew(20), gpdummy(6)
      INTEGER idx(MAXARR,MAXNEIGHBOURS)
      REAL dfim(MAXFRE), fr(MAXFRE)
      REAL slon(MAXARR), slat(MAXARR), lon(MAXARR), lat(MAXARR)
      REAL dist(MAXARR), sumw(MAXARR), delang(MAXARR)
      REAL w(MAXARR,MAXNEIGHBOURS), distmin(MAXARR,MAXNEIGHBOURS)
      REAL sp(MAXANG,MAXFRE), spw(MAXANG,MAXFRE)
      REAL sparr(MAXARR,MAXANG,MAXFRE)
      REAL u10east(MAXARR), u10north(MAXARR)
      REAL useast(MAXARR), usnorth(MAXARR)

      ! Strings
      CHARACTER*256 fin, fout, fwish, str

      ! Initialize
      silent = .FALSE.
      xcen = XCEN0
      ycen = YCEN0
      miss=0
      nwish = 0
      ierr = 0
      ntime=0
      nspec=0
      pi=4.0*atan(1.0)
      rad=pi/180.0
      deg=180.0/pi
      nullspec = .FALSE.

      ! Command line
      narg=iargc()
      if (narg < 3) then
         write (*,*) 
     2"USAGE: invdist spec.seq outspec.seq invdist.inp [xcen ycen slnt]"
         write (*,*) 
     2   " Interpolate 2D spectra using inverse distance weighting"
         write (*,*) 
     2   " spec.seq: binary sequential input file with 2D spectra"
         write (*,*) 
     2   " outspec.seq: ditto output file for selected locations"
         write (*,*) 
     2   " invdist.inp line 1: dtresh [km], maxneighbours"
         write (*,*) 
     2   "  dtresh is the max allowed distance to neighbouring spectra"
         write (*,*) 
     2   "  if dtresh<0, spectra outside range are set to null"
         write (*,*) 
     2   "  maxneighbours is max number of neighbours used for interp"
         write (*,*) 
     2   " subsequent lines: list of wanted locations (lat, lon)"
         write (*,*) 
     2   " xcen ycen: Lon and lat of center pos in rot. spherical grid,"
         write (*,*) 
     2   "  default to xcen=0.0, ycen=0.0 (unrotated grid) if omitted."
         write (*,*) 
     2   "  xcen=0.0, ycen=65.0 applies to the WAM 50km grid."
         write (*,*) 
     2   "  -slnt: do not report missing spectra."
         GOTO 1001
      endif

      ! Open binary sequential file with 2D spectra
      call getarg(1,fin)
      open (IU25, file=fin, access="sequential", status="old",
     2      form="unformatted", err=1025)

      ! Output file
      call getarg(2,fout)
      open (IU26, file=fout, access="sequential",
     2      form="unformatted", err=1026)

      ! Open ASCII wish list (desired spectral locations)
      call getarg(3,fwish)
      open (IU27, file=fwish, form="formatted", status="old",
     2      err=1027)

      ! Read grid parameters xcen, ycen
      if (narg >= 5) then
         call getarg(4,str)
         read (str,*) xcen 
         call getarg(5,str)
         read (str,*) ycen 
      endif
      xcrad = xcen*rad
      ycrad = ycen*rad

      if (narg > 5) then
         silent = .TRUE.
      endif

!!! Read wish list !!!

      ! Loop until end of wish list
      do while ( .TRUE. )
         read (IU27,"(a)",end=100) str
         c = str(1:1)
         ! List may contain comments, first character must be c, *, ! or #
         if ( c/='c' .AND. c/='*' .AND. c/='!' .AND. c/='#') then
            ! First line contains 
            ! dtresh - treshold [km] distance to neighbours and
            ! neighbours - max number of neighbours used for interpolation
            if (nwish == 0) then
               read (str,*) dtresh, neighbours
               if (dtresh < 0.0) nullspec = .TRUE.
               dtresh = abs(dtresh)
               if (neighbours > MAXNEIGHBOURS) then
                  write (*,*) "invdist: Neighbours cannot exceed ",
     2                        MAXNEIGHBOURS
                  GOTO 1030
               endif
            else
            ! Following lines contain lat, lon [deg] of desired locations
               read (str,*) lat(nwish), lon(nwish)
               latrad = lat(nwish)*rad
               lonrad = lon(nwish)*rad
               ! Compute rotated coordinates
               call
     2          sphrot(lonrad,latrad,rlonrad,rlatrad,1,xcrad,ycrad,1)
               rlon=rlonrad*deg
               rlat=rlatrad*deg
               ! Compute local distortion angle [deg]
               call ancorr(delta,rlon,rlat,xcen,ycen)
               delang(nwish) = delta
               ! Count desired locations
            endif
            nwish = nwish+1
         endif
      enddo
100   CONTINUE
      nwish = nwish-1


!!! Count spectral grid points !!!

      ! Read first spectrum
      call getspec(IU25, ident, MAXANG, MAXFRE, sp, ierr)
      if (ierr /= 0) GOTO 1010

      nang = ident(10)
      nfre = ident(11)
      delth = 360.0/real(nang)

      fr(1) = real(ident(12))/1000.0
      th0 = real(ident(19))/100.0

      ! Frequency array
      co = real(ident(13))/100.0
      do m = 2,nfre
         fr(m) = fr(m-1)*co
      enddo

      ! Delta-frequency array
      co1 = 0.5*(co-1.0)*delth*rad
      dfim(1)= co1*fr(1)
      do m = 2,nfre-1
         dfim(m) = co1*(fr(m)+fr(m-1))
      enddo
      dfim(nfre) = co1*fr(nfre-1)

      fchh0 = ident(4) ! first forecast time
      fchh = fchh0
      ihhmi0 = ident(3)
      ihhmi = ihhmi0
      immdd0 = ident(2)
      immdd = immdd0
      
      do while ( ierr==0 .AND. fchh==fchh0 .AND. ihhmi==ihhmi0 .AND. 
     +           immdd==immdd0 )
         nspec = nspec+1
         slat(nspec) = degrees(ident(5),ident(9))
         slon(nspec) = degrees(ident(6),ident(17))

         call getspec(IU25, ident, MAXANG, MAXFRE, sp, ierr)
         !if (ierr /= 0) goto 1010
         fchh = ident(4)
         ihhmi = ident(3)
         immdd = ident(2)
      enddo ! while


!!! Find nearest neighbours !!!

      ! Loop over wish list
      do j = 1, nwish

         ! Compute distance to all spectral locations [km]
         do l = 1, nspec
            dist(l) = spheredist(slon(l),slat(l),lon(j),lat(j))/1000.0
            if (dist(l) < TOL) then
               dist(l) = TOL
            endif
         enddo

         ! Find nearest neighbours
         do i = 1, neighbours
            distmin(j,i) = DISTMAX

            ! Loop over spectral locations
            do l = 1, nspec
               if (dist(l) < distmin(j,i)) then
                  distmin(j,i) = dist(l)
                  idx(j,i) = l
               endif
            enddo ! l

            ! Remove nearest spectral location from next search
            dist(idx(j,i)) = DISTMAX
         enddo ! i

         ! Compute weights
         i = 1
         sumw(j)=0.0
         do while ( (i <= neighbours) .AND. (distmin(j,i) <= dtresh) )
            l = idx(j,i)
            w(j,i) = distmin(j,i)**(-POW)
            sumw(j) = sumw(j) + w(j,i)
            i = i+1
         enddo ! while

         ! Flag loners
         if ( (distmin(j,1) > dtresh) .AND. .NOT. silent) then
            write (*,"(a,i5,a,f11.6,a,f12.6)")
     2       "WARNING: invdist: No spectra within range for pos", j,
     2       " at lat",lat(j),", lon",lon(j)
            miss = miss+1
         endif
      enddo ! j

      ! Rewind spectral input file
      ierr = 0
      rewind IU25


!!! Time loop !!!

      do while ( ierr==0 )

         ntime = ntime+1
         ! Loop over spectral locations
         do l = 1, nspec

            ! Read spectrum
            call getspec(IU25, ident, MAXANG, MAXFRE, sp, ierr)
            if (ierr /= 0) goto 1000
            rlat=degrees(ident(5),ident(9))
            rlon=degrees(ident(6),ident(17))
            if (spheredist(slon(l),slat(l),rlon,rlat) > 10.0) then
               write (*,*)
     2         "ERROR: invdist: input sequence mismatch at time ",
     2         ntime, nspec
               goto 1030
            endif

            ! Array of spectra
            do k=1,nang
               do m=1,nfre
                  sparr(l,k,m) = sp(k,m)
               enddo
            enddo

            ! Wind speed [m/s]
            u10 = real(ident(7))/100.0
            ! Wind dir [rad]
            ! clockw. from rotated N, going to (oceanographic convention) [rad]
            thw = rad*ang360(real(ident(8)/10.0))
            ! Rot east component wind speed [m/s]
            u10east(l)  = u10*sin(thw)
            ! Rot north component wind speed [m/s]
            u10north(l) = u10*cos(thw)
            ! Friction velocity [m/s]
            ustar = real(ident(16))/1000.0
            ! Rot east component friction velocity [m/s]
            useast(l) = ustar*sin(thw)
            ! Rot north component friction velocity [m/s]
            usnorth(l) = ustar*cos(thw)

         enddo ! l

         ! Loop over wish list
         do j = 1, nwish

            ! Initialize to null
            u10eastw =  0.0
            u10northw = 0.0
            useastw =   0.0
            usnorthw =  0.0
            u10 = 0.0
            ustar = 0.0
            thw = 0.0
            thq = 0.0

            ! Prepare new header
            do i = 1, 20
               identnew(i) = ident(i)
            enddo

            ! Spectrum set to null
            do k=1,nang
               do m=1,nfre
                  spw(k,m) = 0.0
               enddo 
            enddo

            ! Loop over list of neighbours
            i = 1
            do while ( i <= neighbours .AND. distmin(j,i) <= dtresh )

               ! Spectral location index
               l = idx(j,i)

               ! Spectral component weights
               do k=1,nang
                  do m=1,nfre
                     spw(k,m) = spw(k,m) + sparr(l,k,m)*w(j,i)
                  enddo
               enddo

               ! Scalar weights
               u10eastw = u10eastw+u10east(l)*w(j,i)
               u10northw = u10northw+u10north(l)*w(j,i)
               useastw = useastw+useast(l)*w(j,i)
               usnorthw = usnorthw+usnorth(l)*w(j,i)

               i = i+1

            enddo ! while i

            ! Normalize spectrum and compute scalars
            if (distmin(j,1) <= dtresh) then

               ! Divide spectrum by sum of weights
               do k=1,nang
                  do m=1,nfre
                     spw(k,m) = spw(k,m)/sumw(j)
                  enddo
               enddo
               u10eastw = u10eastw/sumw(j)
               u10northw = u10northw/sumw(j)
               useastw = useastw/sumw(j)
               usnorthw = usnorthw/sumw(j)

               ! Wind speed [m/s]
               u10 = sqrt(u10eastw**2.0+u10northw**2.0)
               ! Wind dir [deg], cw from N, going to (oceanog. conv.)
               thw = atan2(u10eastw,u10northw) ! [rad]
               thw = ang360(thw*deg) ! adjust angle and convert to deg
               ! Friction velocity [m/s]
               ustar = sqrt(useastw**2.0+usnorthw**2.0)

               ! Wave mean direction thq [rad], see STHQ.f
               si=0.0
               ci=0.0
               do k=1,nang
                  tmp = 0.0
                  do m=1,nfre
                     tmp = tmp+spw(k,m)*dfim(m)
                  enddo
                  th = real(k-1)*delth + 0.5*delth
                  si = si+sin(th*rad)*tmp
                  ci = ci+cos(th*rad)*tmp
               enddo
               if (ci==0.0) ci = 0.1e-33
               thq = atan2(si,ci)*deg
               thq = ang360(thq)

            endif ! distmin <= dtresh

            ! Write spectrum if location is within range or nullspec is true
            if ( (distmin(j,1) <= dtresh) .OR. nullspec) then
               ! Update ident
               identnew(7) = nint(u10*100.0)
               identnew(8) = nint(thw*10.0)    ! [deg/10]
               identnew(14) = nint(delang(j)*100.0)
               identnew(15) = 0 ! WAM index undefined, set to zero
               identnew(16) = nint(ustar*1000.0)
               identnew(18) = nint(thq*deg*100.0) ! [deg/100]
               ! Longitude and latitude [deg]
               call hiresdeg(lat(j), identnew(5), identnew(9))
               call hiresdeg(lon(j), identnew(6), identnew(17))
               ! Automatic scaling
               identnew(20) = -32767

               ! Store ident and spectrum
               call putspec(IU26, identnew, MAXANG, MAXFRE, spw, ierr)
            endif

         enddo ! j = 1, nwish

      enddo ! while l

!!! End time loop !!!

      ! Summary
1000  CONTINUE
      ! Adjust time counter
      ntime = ntime-1
      if (miss==1) then
         write (*,*) "WARNING: invdist: one location was out of bounds"
      else if (miss>1) then
         write (*,*) "WARNING: invdist:", miss, 
     2               " locations were out of bounds"
      endif

      if (nullspec) miss=0

      if (nwish-miss==0) then
         write (*,*) "WARNING: invdist: no spectra written"
      endif
      if (nwish-miss==1) then
         write (*,*) "invdist: wrote", nwish-miss, " spectrum at",ntime,
     2               " times"
      else if (nwish-miss>1) then
         write (*,*) "invdist: wrote", nwish-miss, " spectra at",ntime,
     2               " times"
      endif
         
      goto 1030

      ! Error handling
1010  CONTINUE
      write (*,*) "ERROR: invdist: Error in getspec, ierr == ", ierr
      GOTO 1030
1025  CONTINUE
      write (*,*) "ERROR: invdist: Open error unit ", IU25
      GOTO 1030
1026  CONTINUE
      write (*,*) "ERROR: invdist: Open error unit ", IU26
      GOTO 1030
1027  CONTINUE
      write (*,*) "ERROR: invdist: Open error unit ", IU27
      GOTO 1030

      ! Close files
1030  CONTINUE
      close (IU27)
      close (IU25)
      close (IU26)

1001  CONTINUE
      END
