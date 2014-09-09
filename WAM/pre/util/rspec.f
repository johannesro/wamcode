      PROGRAM rspec
! Read sequential binary 2D spectrum file, DNMI style.
! Print ident(20).
!
! See also rseq, statseq, se_sekv.
!
! Compile: gfortran -O2 rspec.f specio.f -o rspec
!
! 2006-11-16, oyvind.breivik@met.no
! 2007-03-12, handling empty spectra, OB
! 2009-05-27, handles multiple analysis times
! 2009-05-28, handles single time step

      IMPLICIT NONE

      ! Constants
      REAL EPS
      INTEGER MAXANG, MAXFRE
      INTEGER IFIN

      PARAMETER ( EPS=0.000001 )
      PARAMETER ( MAXANG=50, MAXFRE=51 )
      PARAMETER ( IFIN=25 )

      ! Functions called
      INTEGER iargc

      ! Vars
      REAL pi, deg, rad
      REAL th0, co, co1, delt25, delang, fd
      REAL hm0, etot, thpeak, fmean, tmean, tpeak, emax
      REAL delth
      REAL temp, xx, denom
      INTEGER narg, n, nspec, ierr, fchh, fchh0
      INTEGER immdd0, immdd, ihhmi0, ihhmi
      INTEGER nang, nfre, k, l, m, ntime, mmax, kmax
      
      ! Arrays
      INTEGER ident(20)
      REAL    sp(MAXANG,MAXFRE)
      REAL    th(MAXANG), dfim(MAXFRE), fr(MAXFRE)
      REAL    esumth(MAXFRE), esumfr(0:MAXANG+1), theta(0:MAXANG+1)
      CHARACTER*132 fname

      ! Initialize
      n=0
      nspec=0
      pi=4.0*atan(1.0)
      rad=pi/180.0
      deg=180.0/pi

      ! Command line arguments
      narg=iargc()
      call getarg(1,fname)
      if (narg==0) then
         print *,"USAGE: rspec <wam.seq>"
         print *," Lists ident(1:20) of DNMI sequential binary"
         print *," containing spectral information."
         print *," Computes total energy, mean frequency and mean dir."
         print *," See also rseq, statseq, se_sekv."
         GOTO 1001
      endif

      open(IFIN, file = fname, status = 'old', access = 'sequential',
     >     form = 'unformatted')

      ! Write header line.
      write (6,101)


!!! Initialize spectral parameters from first spectrum !!!

      ! Read first spectrum
      call getspec(IFIN,ident,MAXANG,MAXFRE,sp,ierr)
      if (ierr /= 0) goto 999

      nang = ident(10)
      nfre = ident(11)
      delth = 360.0/real(nang)
      fr(1) = real(ident(12))/1000.0
      th0 = real(ident(19))/100.0
      do k=1,nang
         theta(k)= real(k-1)*delth !+ th0
      enddo
      theta(0)=theta(1)-delth
      theta(nang+1)=theta(nang)+delth  

      ! Frequency array
      co = real(ident(13))/100.0
      do m = 2,nfre
         fr(m) = fr(m-1)*co
      enddo
      delt25 = fr(nfre)/4.0*delth*rad

      ! Delta-frequency array
      co1 = 0.5*(co-1.0)*delth*rad
      dfim(1)= co1*fr(1)
      do m = 2,nfre-1
         dfim(m) = co1*(fr(m)+fr(m-1))
      enddo
      dfim(nfre) = co1*fr(nfre-1)

!!! Count spectral locations (number of wet boundary points) !!!

      fchh0 = ident(4) ! first forecast time
      fchh = fchh0
      ihhmi0 = ident(3)
      ihhmi = ihhmi0
      immdd0 = ident(2)
      immdd = immdd0
      
      do while ( ierr==0 .AND. fchh==fchh0 .AND. ihhmi==ihhmi0 .AND. 
     +           immdd==immdd0 )
         nspec = nspec+1

         call getspec(IFIN,ident,MAXANG,MAXFRE,sp,ierr)
         !if (ierr /= 0) goto 999
         fchh = ident(4)
         ihhmi = ident(3)
         immdd = ident(2)

      enddo ! while

      ! Rewind spectral input file
      rewind IFIN
      ierr = 0


!!! Time loop !!!

      do while ( ierr==0 )

         ntime = ntime+1


       !!! Loop over spectral locations !!!

         do l = 1, nspec


          !!! Read spectrum !!!

            call getspec(IFIN,ident,MAXANG,MAXFRE,sp,ierr)
            if (ierr /= 0) GOTO 999

            ! Set position-specific parameters
            delang = real(ident(14))/100.0

          !!! Compute etot !!!

            etot=0.0
            ! Loop through all spectral bins
            do m = 1,nfre
               temp = 0.0
               do k = 1,nang
                  temp = temp+sp(k,m)
               enddo
               etot = etot+temp*dfim(m)
            end do     
            ! Compute tail energy.
            etot = etot+delt25*temp
            hm0 = 4.0*sqrt(etot)


          !!! Compute mean frequency in deep water !!!

            fmean = 0.0
            do m = 1,nfre
               fd = dfim(m)/fr(m)
               temp = 0.0
               do k= 1,nang
                  temp = temp+sp(k,m)
               enddo     
               fmean = fmean+fd*temp
            enddo

            ! Add tail energy
            fmean = fmean+delt25*temp
            fmean = max(EPS,fmean)
            fmean = etot/fmean

            ! Check for division by zero
            if (fmean > EPS) then
               tmean = 1.0/fmean
            else
               tmean = 0.0
            endif

          !!! Find peak period !!!

            emax=0.0
            mmax=0
            do m = 1,nfre
               esumth(m)=0.0
               do k = 1,nang
                  esumth(m)=esumth(m)+sp(k,m)
               enddo
               if (esumth(m) > emax) then 
                  emax=esumth(m)
                  mmax=m
               endif
            enddo

            ! More accurate computation of peak period
            tpeak=0.0
            if (mmax==0) then
               tpeak=0.0 ! Null spectrum
            else if (mmax==1) then
               tpeak=1.0/fr(1)
            else if (mmax == nfre) then
               tpeak=1.0/fr(nfre)
               if (hm0<tpeak**2/25.0) tpeak = 5.0*sqrt(hm0)
            else
               xx = (esumth(mmax+1)-esumth(mmax-1))/
     $              (esumth(mmax+1)-esumth(mmax)-esumth(mmax)+
     $               esumth(mmax-1))
               tpeak = (1.0+xx*(.066+2.04315e-3*(1.0-xx)))/fr(mmax)
            endif

         !!! Compute peak direction !!!

            emax=0.0
            do k=1,nang
               esumfr(k)=0.0
               do m=1,nfre
                  esumfr(k)=esumfr(k)+sp(k,m)
               enddo
               if (esumfr(k) > emax) then
                  emax=esumfr(k)
                  kmax=k
               endif
            enddo         

            ! More accurate computation of peak direction
            thpeak=0.0
            if (esumfr(kmax) > 0.0) then
               esumfr(0) = esumfr(nang)
               esumfr(nang+1) = esumfr(1)
               denom = esumfr(kmax+1)-esumfr(kmax)*2.0+esumfr(kmax-1)
               thpeak = theta(kmax)-0.5*delth*(esumfr(kmax+1)
     2                  -esumfr(kmax-1))/denom
               thpeak = thpeak + th0 + delang + 180.0
               if (thpeak < 0.0) thpeak = thpeak+360.0
               if (thpeak >= 360.0) thpeak = thpeak-360.0
            endif          

            n=n+1
            write
     +       (6,"(i7,1x,i4.4,1x,i4.4,i5,i6,i7,i8,3i6,i5,i4,i5,
     +        i6,i7,3i8,i7,i7,i7,f6.2,2f6.2,f6.1)")
     +        n, (ident(k),k=1,20), hm0, tmean, tpeak, thpeak

         enddo ! l

      enddo ! while ntime

 101  FORMAT(
     +"   NREC YYYY MMDD HHMM  FCHH    LAT     LON   W10  WDIR   LAX ",
     +"NANG NFRE FR(1)  CO   DANG     IDX   USTAR     LOX    THQ",
     +"    TH0  SCALE   HM0 TMEAN    TP   DDP")
 999  close (IFIN)
      write (6, 101)
      print *, "Read", n, " records from", nspec, " locations at", 
     +         ntime-1, " times"

1001  CONTINUE

      END ! PROGRAM rspec
