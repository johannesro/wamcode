      PROGRAM boundwam

! Create spectral boundary file from sequential binary spectra 
! (DNMI format).
!
! Note: only rotated and unrotated lat-lon projections accepted.
!
! USAGE: boundwam inspec.seq wambound_ecwam10.inp outspec.wam
!
! Input  file - inspec.seq:  sequential binary 2D spectra - DNMI format
! Output file - outspec.wam: ECMWF-WAM 2D spectral binary
! Input file  - wambound_ecwam10.inp: boundary and resolution parameters
! Example:
! 
!
! 2012-02-08 - oyvind.breivik@met.no
! 2013-11-01 - anac@met.no .. Modified so ir works with mywave WAM version
!
! wambound_ecwam10.inp: - serves three purposes:
! (1) boundary parameters,
!     "coarse grid" resolution parameters,
!     projection parameters (xcen, ycen)
!     to wamboundinp.f
! (2) boundary parameters,
!     "coarse grid" resolution parameters
!     to boundwam.f
! (3) dtresh and maxneighbours passed through to invdist.f
!
! Requires rotspher.f, specio.f, julian.f
!
! Compile:
! Vilje: 
! module load intelcomp/13.0.1 mpt/2.06 netcdf/4.3.0
! ifort  -O3 -real_size 64 -traceback -ip -align all   -xAVX -ftz -fno-alias -no-prec-div -no-prec-sqrt -o boundwam  boundwam.f rotspher.f specio.f julian.f
!
! Linux: gfortran -O2 boundwam.f rotspher.f specio.f julian.f -o boundwam
!
!
!
!
! CCC Must be real*8 

      IMPLICIT NONE

      ! Functions called
      INTEGER iargc
      INTEGER julday
      REAL degrees

      ! Constants
      INTEGER MAXANG, MAXFRE, MAXARR
      INTEGER IU02, IU04, IU21,IU22, IU25, IDUM
      REAL RDUM

      PARAMETER ( MAXANG=50, MAXFRE=51, MAXARR=10000 )
      PARAMETER ( IU02=2, IU04=4, IU21=21, IU22=22, IU25=25)
      PARAMETER ( IDUM = 999)
      PARAMETER ( RDUM = 999.9 ) 

      ! Vars
      INTEGER*8 fchh, fchh0, fchhfirst, hh, hh0, jd0,jd
      INTEGER*8 immdd0, immdd, ihhmi0, ihhmi
      INTEGER*8 idelprc, idtpro
      INTEGER   yyyy, mm, dd
      INTEGER*8 ntime, nskip, nang, nfre, narg
      INTEGER*8 i, j, k, ij, m, nx, ny, idgrid, idummy
      INTEGER*8 ierr
      INTEGER*8 nbounc
      INTEGER*8 bdepth
      INTEGER*8 depth(MAXARR)
      REAL*8  xang, xfre, xbou, xdelc, xdelf
      REAL*8 th0, co, co1, delt25, fd, rdummy
      REAL*8 dlamac, dphiac, amosoc, amonoc, amoeac, amowec
      REAL*8 pi, rad, deg
      REAL*8 lonrad, latrad, rlon, rlat
      REAL*8 xcen, ycen
      REAL*8 emean, thq, fmean
      REAL*8 delth
      REAL*8 temp
      REAL*8 si, ci, sih, cih
      REAL*8 brotlat, brotlon

     

      ! Arrays
      INTEGER ident(20)
      REAL*8    blngc(MAXARR), blatc(MAXARR)
      REAL*8    slon(MAXARR), slat(MAXARR)
      REAL*8    sp(MAXANG,MAXFRE)
      REAL     th(MAXANG), dfim(MAXFRE), fr(MAXFRE)

      ! Strings
      CHARACTER*256 fin, fout, flist, str
      CHARACTER*14 cdate

      ! Initialize
     

      nbounc=0
      ierr=0
      ntime=0
      nskip=0
      pi=4.0*atan(1.0)
      rad=pi/180.0
      deg=180.0/pi
      xdelf = 3 * 60 * 60 !time to save files in sec 
      ! Command line
      narg=iargc()
      if (narg < 3) then
         write (*,*) 
     2   "USAGE: boundwam inspec.seq boundwam.inp boundwam.out"
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

      ! Open ASCII boundary depth input file
      call getarg(3,flist)
      open (IU22, file=flist, form="formatted", status="old",
     2      err=1021)
      write(*,*)'boundary depth input file',flist
      

      ! Open WAM binary output file
      call getarg(4,fout)
      open (IU02, file=fout, form="unformatted", status="unknown", 
     2      err=1009)
      ! Open asscii output file
      open (IU04, file='TMP', form="formatted", err=1009)

      !!! Read input list !!!

      ! Read grid.sph identificator line for coarse grid boundary
      ! CCC not needed
      read (IU21,"(a)",err=1021) str
      read (str(10:100),*) 
     2 idgrid, nx, ny, amowec, amosoc, dlamac, dphiac, xcen, ycen
      read (IU21,*,err=1021) rdummy, idummy, fr(1), fchhfirst

      ! Compute ROTATED eastern and northern boundaries [deg] !CCC
      !amoeac = amowec+(nx-1)*dlamac
      !amonoc = amosoc+(ny-1)*dphiac


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

      ! Read the depths at the boundary
      do ij=1,nbounc
          read (IU22,125)  brotlat, brotlon, bdepth
          if ( (abs( brotlat - blatc(ij) ).le.0.00001).and. 
     &         (abs( brotlon - blngc(ij) ).le.0.00001)) then
                depth(ij) = bdepth         
          endif
      enddo
      write(*,*)'depth(1:10)',depth(1:10)
125   format(2F12.6,2x,I10)
       
!!! Write header !!! 
      xang = float(nang)
      xfre = float(nfre)
      xbou = float(nbounc)
     
      xdelc=xdelf !boundary input data are every 3h. The data are 
                  !stored in files which have a time step of 3h 
                  !too, which is that there is one field per file.
     !!! OJO I repeat two times xdlef 
         write (IU02) real(xang,4), real(xfre,4), real(th0*rad,4), 
     &   real(fr(1),4), real(co,4), real(xbou,4), real(xdelc,4), 
     &   real(xdelf,4)  
        
         write (IU04,*)real(xang,4), real(xfre,4), real(th0*rad,4), 
     &   real(fr(1),4), real(co,4), real(xbou,4), real(xdelc,4), 
     &   real(xdelf,4)  

         write (IU02) ( real(depth(m),4), m=1,nbounc )
 
         write (IU04,*)( real(depth(m),4), m=1,nbounc )


 200     format(f12.5,3x, f12.5, 3x, e11.4, 2x, e11.4, f11.6,f16.4,
     &   f15.2,f15.2)
 201     format(5(3x,f8.3))

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

          
               write (IU02) real(blngc(ij),4), real(blatc(ij),4),
     &               cdate,  real(emean,4), real(thq,4), real(fmean,4)

               write (IU04,*)real(blngc(ij),4), real(blatc(ij),4),
     &               cdate,  real(emean,4), real(thq,4), real(fmean,4)

 203           format(2f16.5,2x, a14,3f13.7)

               write (IU02) ( (real(sp(k,m),4), k=1,nang), m=1,nfre)
               write (IU04,*) ( (real(sp(k,m),4), k=1,nang), m=1,nfre)

           

110         CONTINUE

         enddo ! ij

      enddo ! while ntime

!!! End time loop !!!


      ! Summary
1000  CONTINUE

     
 
      ! Adjust time counter
      ntime = ntime-1
      if (nbounc==1) then
         write (*,*) "boundwam: wrote one spectrum at", ntime-nskip,
     2               " times"
      else if (nbounc>1) then
         write (*,*) "boundwam: wrote", nbounc," spectra at",
     2               ntime-nskip," times"
      endif
      write (*,*) "boundwam: skipped", nskip,
     2" times before", fchhfirst,"H"

      goto 1030

      ! Error handling
1009  continue
      write (*,*) "ERROR: boundwam: Open error unit ", IU02
      goto 1030
1010  continue
      write (*,*) "ERROR: boundwam: Error in getspec, ierr == ", ierr
      goto 1030
1021  continue
      write (*,*) "ERROR: boundwam: Open error unit ", IU21
1025  continue
      write (*,*) "ERROR: boundwam: Open error unit ", IU25
      goto 1030
      goto 1030

      ! Close files
1030  continue
      close (IU25)
      close (IU02)
      close (IU22)
      close (IU04)
      END ! PROGRAM boundwam

!##compiling boundwam
!#!module load cmkl intelcomp mpt
!#!ifort  -O3 -real_size 64 -traceback -ip -align all\
!#!             -xAVX -ftz -fno-alias -no-prec-div -no-prec-sqrt -o boundwam  boundwam.f rotspher.f specio.f julian.f


