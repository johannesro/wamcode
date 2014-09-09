      PROGRAM CHANGE_IDENT
C
C        ***** Changes field identification on met.no SEQUENTIAL *****
C        ***** feltfiles according to input list                 *****
C----------------------------------------------------------------------------
C    Syntax for compiling: <F77-compiler> -o change_ident change_ident.f
C----------------------------------------------------------------------------
C        DNMI/FoU   December 2007  Harald Engedahl: First version
C----------------------------------------------------------------------------
C
      implicit none
c
c  mfsize : maximum field size, input and output
c
      integer   mfsize,ldata
      parameter (mfsize=2000000,ldata=mfsize)
c
      character*132 navn1, navn2
c
      integer   iunit1,iunit2,iunit3,ierr
      integer   identf1(20),identf2(20),identtest1(20),identtest2(20)
      integer   i,ios,ihours1,ihours2,iprog,ierr1,ierr2
      integer   nxin,nyin,nword,nw,ierror
      integer*2 ident(20),idata(mfsize)
      integer   lgeom
c
c===========================================================================
c

c***** Start read input list file

c  Read input and output file names:
      read(5,*) navn1
      write(6,*)'Input SEQUENTIAL  feltfile : ',navn1
      read(5,*) navn2
      write(6,*)'Output SEQUENTIAL feltfile : ',navn2
      write(6,*)
c
c  Read list of 20 words of identification.
c  All input idents gt. -32760 is changed on output file.
      do i=1,20
        read(5,*) identtest1(i),identtest2(i)
        write(6,*) identtest1(i),identtest2(i)
      end do

c***** End read input list file *******
c
      iunit1=14 ! input file
      iunit2=15 ! output file
c
      open(iunit1, file = navn1, status = 'old', access = 'sequential',
     >     form = 'unformatted')
      open(iunit2, file = navn2, status = 'new', access = 'sequential',
     >     form = 'unformatted')
c
c------------------------------------------------------------------------------------------
c
      write(6,*)' Start scanning file.....'
      write(6,*)
c
 11   continue
c
c Reads ident... 
      read(iunit1,iostat=ios,err=911,end=921) (ident(i),i=1,20)
      do i=1,20
        identf1(i)=ident(i)
      end do
c
c Check if extra ident:
c
      lgeom=0
      if(ident(9).ge.1000)
     ,   lgeom=ident(9)-(ident(9)/1000)*1000
c
      nxin=ident(10)
      nyin=ident(11)
      nword=nxin*nyin
      nw=nword+lgeom
cc      write(6,*)
cc      write(6,*)' nx, ny, lgeom : ',nxin,nyin,lgeom
cc      write(6,*)
c
      if(nword.gt.ldata) then
        write(6,*) ' ** field length too big',
     +                        ' (input ldata too small)'
        write(6,*) ' **          ldata = ',ldata
        write(6,*) ' ** ident: ',(ident(i),i=1,11)
        write(6,*) ' **        ',(ident(i),i=12,20)
        write(6,*) ' ** nx,ny,nx*ny: ',nxin,nyin,nword
        ierr=1
        goto 990
      end if
c
      if(nw.gt.mfsize) then
        write(6,*) ' ** field length + extra ident too big',
     +                        ' (input ldata too small)'
        write(6,*) ' **         mfsize = ',mfsize
        write(6,*) ' ** ident: ',(ident(i),i=1,11)
        write(6,*) ' **        ',(ident(i),i=12,20)
        write(6,*) ' ** nx,ny,nx*ny,lgeom: ',nxin,nyin,nword,lgeom
        ierr=1
        goto 990
      end if
c
c Reads the field + extra ident...
      read(iunit1,iostat=ios,err=911,end=921) (idata(i),i=1,nw)	
c
c Changes ident according to input list...
      do i=1,20
        if(identtest1(i).eq.identf1(i)) then
          identf2(i)=identtest2(i)
        else
          identf2(i)=identf1(i)
        endif
        ident(i)=identf2(i)
      enddo
c
      write(iunit2,iostat=ios,err=900) (ident(i),i=1,20)
      write(iunit2,iostat=ios,err=900) (idata(i),i=1,nw)

cc    write(6,100) (ident(i),i=1,20)
cc100 format(20i6)
c
c Next field...
      go to 11
c
c------------------------------------------------------------------------------------------
c
      goto 990
c
  900 ierr=1
      write(6,*) ' ** write error. file,iostat: ',iunit2,ios
c
  911 ierr=1
      write(6,*) ' ** read error. file,iostat: ',iunit1,ios
      goto 990
  921 ierr=2
      write(6,*) ' ** end_of_file. file,iostat: ',iunit1,ios
      goto 990
c
  990 continue
C
      END
