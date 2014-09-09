C *****************************************************************
C
      subroutine getflt(mode,iunit,itime,ident,ldata,fdata,
     +                  igtype,gparam,nfskip,ierr)
c
c  16 bit input and unpack (not handling 'undefined')
c
c  input:
c     mode:   0 = read field
c             1 = read field, skip fields until time > itime
c             2 = read field if field time = itime
c                 (otherwise next read starts at the same field)
c           100 = read field identification
c                 (next read starts at the same field)
c           101 = read field identification, skip fields
c                 until time > itime
c                 (next read starts at the same (last) field)
c           102 = read field identification, skip fields
c                 until time >= itime
c                 (next read starts at the same (last) field)
c           200 = scan rest of the file and read field with
c                 matching identification, specified identification
c                 input in ident(1:20) where -32767 means any value
c           201 = scan the whole file and read field with
c                 matching identification, specified identification
c                 input in ident(1:20) where -32767 means any value
c            -1 = clean up after a file is closed, and the same
c                 file unit no. is used for another file.
c     iunit: file unit no.
c     itime: used if mode is 1, 2, 101 or 102.
c            itime(1) - year
c                 (2) - month
c                 (3) - day
c                 (4) - hour (utc)
c                 (5) - forecast time in hours
c                       (added to previous date and time)
c     ident:  used if mode is 200 or 201
c     ldata:  length of fdata (max field size)
c
c  output:
c     ident(20): field identification
c     fdata(..): field (unscaled, according to identification)
c     igtype:    grid type (1=polar. 2=geo. 3+rot.sph. ...)
c     gparam(6): grid parameters
c     nfskip:    no. of fields skipped (for mode=2,101,102)
c     ierr = 0:  read o.k.
c            1:  read error
c            2:  read error, end_of_file
c            3:  no field returned due to mode=2 and itime spec.
c
c  note: mode=100,101,102 and possibly mode=2 will only read
c        a field's identification and not the field (data).
c        to avoid (cray) backspace problems, the identification
c        is stored until the next read.
c
c  warning: using file unit no. (not file name) to identify
c           files when storing field identification.
c           if more than one file is opened with the same
c           unit, use 'call getflt(-1,...)' after closing
c           a file to avoid errors.
c
c  no computer dependant i/o methodes
c                  ------
c                  Does also work with CRAY's from PrgEnv 3.0.2.1
c                  and above, provided an assign -F f77 -N ieee u:iunit
c                  is applied somewhere.
c                  NTNU/ITEA 1998-06-23 Jorn Amundsen
c                  ------
c
c  DNMI/FoU  12.04.1992  Anstein Foss
c  DNMI/FoU  19.08.1993  Anstein Foss
c  DNMI/FoU  16.04.1998  Anstein Foss ... cray.t3e
c  DNMI/FoU  02.07.1998  Anstein Foss ... integer*2 also for Cray
c
      include 'inout.inc'
c
c  input/output
      integer   mode,iunit,ldata,igtype,nfskip,ierr
      integer   itime(5),ident(20)
      real      fdata(ldata),gparam(6)
c
c  local
      integer*2 ident2(20),idata2(mfsize)
c
      integer   itimef(5),idents(20)
c
      parameter (maxsav=5,maxs20=maxsav*20)
      integer   iusave(maxsav),idsave(20,maxsav)
c
      parameter (mlgeom=50)
      integer   lgtype
      integer*2 idgeom(20+1+mlgeom)
      real      gplast(6)
c
      data iusave/maxsav*0/
      data idsave/maxs20*0/
c
      data lgtype/0/
      data idgeom/71*0/
      data gplast/6*999999./
c
c
      if(iunit.lt.1  .or. iunit.gt.99 .or.
     +  (mode.ne.0   .and. mode.ne.1   .and.  mode.ne.2   .and.
     +   mode.ne.100 .and. mode.ne.101 .and.  mode.ne.102 .and.
     +   mode.ne.200 .and. mode.ne.201 .and.  mode.ne.-1)) then
        write(6,*) ' **getflt** error ** mode,iunit: ',mode,iunit
        stop
      end if
c
      if(mode.eq.200 .or. mode.eq.201) then
        do i=1,20
          idents(i)=ident(i)
        end do
      end if
c
      if(mode.eq.201) rewind(iunit)
c
c        a) read field identification (one record)
c        b) read field data           (one record)
c
      ierr=0
      nfskip=0
c
      nsave=0
      isave=0
      do i=1,maxsav
        if(iusave(i).gt.0)     nsave=nsave+1
        if(iusave(i).eq.iunit) isave=i
      end do
c
      if(mode.eq.-1) then
        if(isave.gt.0) iusave(isave)=0
        goto 990
      end if
c
      if(mode.eq.201 .and. isave.gt.0) then
        iusave(isave)=0
        isave=0
        nsave=nsave-1
      end if
c
      if(isave.gt.0) then
        do i=1,20
          ident(i)=idsave(i,isave)
        end do
        iusave(isave)=0
        nsave=nsave-1
        goto 30
      end if
c
   20 continue
c
      read(iunit,iostat=ios,err=910,end=920) (ident2(i),i=1,20)
      do i=1,20
        ident(i)=ident2(i)
      end do
c
   30 continue
c
c..iexit set to 1 if the next field is not read
c..iskip set to 1 if the next field is skipped
      iexit=0
      iskip=0
c
      if(mode.eq.100) then
        iexit=1
      elseif(mode.eq.200 .or. mode.eq.201) then
        n=0
        do i=1,20
          if(idents(i).eq.-32767 .or. idents(i).eq.ident(i)) n=n+1
        end do
        if(n.ne.20) iskip=1
      elseif(mode.ne.0) then
        itimef(1)=ident(12)
        itimef(2)=ident(13)/100
        itimef(3)=ident(13)-(ident(13)/100)*100
        itimef(4)=ident(14)/100
        itimef(5)=ident(4)
        call hrdiff(0,0,itime,itimef,ihours,ierr1,ierr2)
        if(mode.eq.1   .and. ihours.le.0) iskip=1
        if(mode.eq.2   .and. ihours.ne.0) iexit=1
        if(mode.eq.101 .and. ihours.le.0) iskip=1
        if(mode.eq.101 .and. ihours.gt.0) iexit=1
        if(mode.eq.102 .and. ihours.lt.0) iskip=1
        if(mode.eq.102 .and. ihours.ge.0) iexit=1
      end if
c
      if(iexit.eq.1) then
c  return without reading the field (data). save identification.
        if(nsave.eq.maxsav) iusave(1)=0
        isave=0
        do n=1,maxsav
          if(iusave(n).gt.0) then
            isave=isave+1
            iusave(isave)=iusave(n)
            do i=1,20
              idsave(i,isave)=idsave(i,n)
            end do
          end if
        end do
        isave=maxsav
        iusave(isave)=iunit
        do i=1,20
          idsave(i,isave)=ident(i)
        end do
      end if
c
      if(iexit.eq.1) then
        ierr=0
        if(mode.eq.2) ierr=3
        goto 990
      end if
c
      nxin=ident(10)
      nyin=ident(11)
      nword=nxin*nyin
c
      lgeom=0
      if(ident(9).ge.1000) lgeom=ident(9)-(ident(9)/1000)*1000
c
      if(nword.gt.ldata .and. iskip.eq.0) then
        write(6,*) ' **getflt** field length too big',
     +                        ' (input ldata too small)'
        write(6,*) ' **          ldata = ',ldata
        write(6,*) ' ** ident: ',(ident(i),i=1,11)
        write(6,*) ' **        ',(ident(i),i=12,20)
        write(6,*) ' ** nx,ny,nx*ny: ',nxin,nyin,nword
        ierr=1
        goto 990
      end if
c
      if(nword+lgeom.gt.mfsize) then
        write(6,*) ' **getflt** field length too big',
     +                        ' (buffer too small)'
        write(6,*) ' **         mfsize = ',mfsize
        write(6,*) ' ** ident: ',(ident(i),i=1,11)
        write(6,*) ' **        ',(ident(i),i=12,20)
        write(6,*) ' ** nx,ny,nx*ny,lgeom: ',nxin,nyin,nword,lgeom
        ierr=1
        goto 990
      end if
c
      nw=nword+lgeom
      read(iunit,iostat=ios,err=910,end=920) (idata2(i),i=1,nw)	
c
      if(iskip.eq.0) then
        iscale=ident(20)
        scale=10.**iscale
c  following (fast) loop does not handle 'undefined' values
        do i=1,nword
          fdata(i)=scale*idata2(i)
        end do
ccc.............................................................
ccc the next (slow) loop handles 'undefined' values
ccc     do i=1,nword
ccc       if(idata2(i).ne.-32767) then
ccc         fdata(i)=scale*idata2(i)
ccc       else
ccc         fdata(i)=undef
ccc       end if
ccc     end do
ccc.............................................................
      end if
c
      if(iskip.ne.0) then
        nfskip=nfskip+1
        goto 20
      end if
c
      if(ident(9).gt.0) then
        newgrid=0
        if(ident(9).ne.idgeom(9)) then
	  newgrid=1
	else
	  do i=15,18
	    if(ident(i).ne.idgeom(i)) newgrid=1
	  end do
	  do i=1,lgeom
	    if(idata2(nword+i).ne.idgeom(21+i)) newgrid=1
	  end do
	end if
	if(newgrid.eq.1) then
	  do i=1,20
	    idgeom(i)=ident(i)
	  end do
	  idgeom(10)=1
	  idgeom(11)=1
	  idgeom(21)=0
	  do i=1,lgeom
	    idgeom(21+i)=idata2(nword+i)
	  end do
          call gridpars(+1,20+1+mlgeom,idgeom,
     +                  lgtype,ix,iy,gplast,ierror)
          if(ierror.ne.0) then
            write(6,*) 'GETFLT: GRIDPARS ERROR ',ierror
	    idgeom(9)=-9999
            ierr=2
            goto 990
          end if
	end if
	igtype=lgtype
	do i=1,6
	  gparam(i)=gplast(i)
	end do
      else
	igtype   =ident(9)
	gparam(1)=ident(15)
	gparam(2)=ident(16)
	gparam(3)=ident(17)
	gparam(4)=ident(18)
	gparam(5)=0.
	gparam(6)=0.
      end if
c
      ierr=0
      goto 990
c
  910 ierr=1
      write(6,*) ' **getflt** read error. file,iostat: ',iunit,ios
      if(mode.eq.200 .or. mode.eq.201) goto 940
      goto 990
  920 ierr=2
      write(6,*) ' **getflt** end_of_file. file,iostat: ',iunit,ios
      if(mode.eq.200 .or. mode.eq.201) goto 940
      goto 990
c
  940 do i=1,20
        ident(i)=idents(i)
      end do
c
  990 continue
      return
      end
C
C *****************************************************************
C
      subroutine putflt(iunit,ident,igtype,gparam,
     +                  ldata,fdata,scale,ierr)
c
c  16 bit pack and output (handling 'undefined')
c
c  input:
c     iunit: file unit
c     ident: field identification,
c	     (set ident(20)=-32767 for automatic scaling)
c     igtype: grid type
c     gparam: grid specification parameters
c     ldata: length of the field
c            (field size given in identification used when
c             writing the field)
c     fdata: the field, unscaled
c     scale: additional scaling of the field before output,
c            usually 1. (this is in addition to the scaling given
c            in the identification, for changing the basic unit)
c
c  output:
c     ierr = 0:  write o.k.
c            1:  write error
c
c  no computer dependant i/o methodes.
c                  ------
c                  Does also work with CRAY's from PrgEnv 3.0.2.1
c                  and above, provided an assign -F f77 -N ieee u:iunit
c                  is applied somewhere.
c                  NTNU/ITEA 1998-06-23 Jorn Amundsen
c                  ------
c
c  DNMI/FoU  12.04.1992  Anstein Foss
c  DNMI/FoU  19.08.1993  Anstein Foss
c  DNMI/FoU  10.10.1997  Anstein Foss ... added extra geometry ident.
c  DNMI/FoU  16.04.1998  Anstein Foss ... cray.t3e
c  DNMI/FoU  02.07.1998  Anstein Foss ... integer*2 also for Cray
c
      include 'inout.inc'
c
      parameter (mldata=mfsize+50)
c
c  input/output
      integer   iunit,igtype,ident(20),ldata,ierr
      real      gparam(6),fdata(ldata),scale
c
c  local
      integer   idento(20)
      integer*2 ident2(20),idata2(mldata)
c
      parameter (mlgeom=50)
      integer   lgtype,llgeom
      integer*2 idgeom(20+1+mlgeom)
      real      gplast(6)
c
      data lgtype,llgeom/-9999,0/
      data idgeom/71*0/
      data gplast/6*999999./
c
      do i=1,20
        idento(i)=ident(i)
      end do
c
      if(igtype.gt.0) then
c
        newgrid=0
        if(igtype.ne.lgtype) newgrid=1
        do i=1,6
          if(gparam(i).ne.gplast(i)) newgrid=1
        end do
        if(newgrid.eq.1) then
          lgtype=igtype
          do i=1,6
            gplast(i)=gparam(i)
          end do
          ix=1
          iy=1
          call gridpars(-1,20+1+mlgeom,idgeom,
     +                  igtype,ix,iy,gparam,ierror)
          if(ierror.ne.0) then
            write(6,*) 'PUTFLT: GRIDPARS ERROR ',ierror
	    lgtype=-999
            ierr=2
            goto 990
          end if
          llgeom=0
          if(idgeom(9).ge.1000) llgeom=idgeom(9)-(idgeom(9)/1000)*1000
        end if
c
        idento( 9)=idgeom( 9)
        idento(15)=idgeom(15)
        idento(16)=idgeom(16)
        idento(17)=idgeom(17)
        idento(18)=idgeom(18)
        lgeom=llgeom
c
      else
c
c  in case of "something" stored as if it were a field, "geometry"
c  identification is left unchanged (extended "geometry" not possible)
c
        lgeom=0
c
      end if
c
      nxout=ident(10)
      nyout=ident(11)
      nword=nxout*nyout
c
      ierr=0
c
      if(nword+lgeom.gt.mldata) then
        write(6,*) ' **putflt** field length too big',
     +                        ' (buffer too small)'
        write(6,*) ' **         mfsize = ',mfsize
        write(6,*) ' **         mldata = ',mldata
        write(6,*) ' ** ident: ',(idento(i),i=1,11)
        write(6,*) ' **        ',(idento(i),i=12,20)
        write(6,*) ' ** nx,ny,nx*ny,lgeom: ',nxout,nyout,nword,lgeom
        ierr=1
        goto 990
      end if
c
      sc=scale
      if(sc.eq.0.) sc=1.
      iscale=ident(20)
c
      if(iscale.eq.-32767) then
c..automatic scaling
	fmax=0.
ccc.............................................................
ccc  fast code not handling undefines values
ccc     do i=1,nword
ccc       fmax=max(fmax,abs(fdata(i)))
ccc     end do
ccc.............................................................
c  slow code handling undefines values
        do i=1,nword
          if(fdata(i).lt.udef) fmax=max(fmax,abs(fdata(i)))
        end do
	if(fmax.gt.0.) then
          fmax=fmax*abs(sc)
	  iscale=log10(fmax)-4.
	  ifmax=nint(fmax*10.**(-iscale))
	  if(ifmax.lt.3278) then
	    iscale=iscale-1
	    ifmax=nint(fmax*10.**(-iscale))
	  end if
	  if(ifmax.gt.32766) iscale=iscale+1
	  iscale=max(iscale,-30)
	else
	  iscale=0
	end if
	idento(20)=iscale
      end if
c
c  scale and 'pack' data
c
      sc=sc*10.**(-iscale)
c
ccc.............................................................
ccc  fast code not handling undefines values
ccc   do i=1,nword
ccc     idata2(i)=nint(sc*fdata(i))
ccc   end do
ccc.............................................................
c  slow code handling undefines values
      do i=1,nword
        if(fdata(i).lt.udef) then
          idata2(i)=nint(sc*fdata(i))
        else
          idata2(i)=-32767
        end if
      end do
ccc.............................................................
ccc   if(ident(20).ne.-32767) then
ccc  fixed scaling, check integer*2 value range
ccc	do i=1,nword
ccc	  idata2(i)=max(min(idata2(i),+32767),-32768)
ccc	end do
ccc   end if
ccc.............................................................
      if(lgeom.gt.0) then
        do i=1,lgeom
          idata2(nword+i)=idgeom(21+i)
        end do
        nword=nword+lgeom
      end if
c
      do i=1,20
        ident2(i)=idento(i)
      end do
c
c  write field identification and field data
c
      write(iunit,iostat=ios,err=900) (ident2(i),i=1,20)
      write(iunit,iostat=ios,err=900) (idata2(i),i=1,nword)
c
      goto 990
c
  900 ierr=1
      write(6,*) ' **putflt** write error. file,iostat: ',iunit,ios
c
  990 continue
c
      return
      end
C
C ********************************************************************
C
      subroutine hrdiff(iup1,iup2,itime1,itime2,ihours,ierr1,ierr2)
c
c  calculate interval in hours between 'itime1' and 'itime2',
c  ihours = 'itime2' - 'itime1'  (positive,zero or negative)
c
c  itime.: year,month,day,time in hours, forecast length in hours
c                                      (positive,zero or negative)
c
c  if iup1=1: 'itime1' is updated to give valid date,time
c  if iup2=1: 'itime2' is updated to give valid date,time
c                                         (forecast length = 0)
c
c  ierr1,ierr2: 0 = o.k. itime1,itime2
c               1 = not o.k. itime1,itime2
c
c  DNMI/FoU  xx.xx.1990  Anstein Foss
c
      integer itime1(5),itime2(5)
c
      integer mdays(12),it(5,2)
      data mdays/31,28,31,30,31,30,31,31,30,31,30,31/
c
      do 10 i=1,5
        it(i,1)=itime1(i)
        it(i,2)=itime2(i)
   10 continue
c
c  put  prog.time  into  year,month,day,time
c
      call vtime(it(1,1),ierr1)
      call vtime(it(1,2),ierr2)
c
c  'it' now gives "veryfing" time of 'itime1' and 'itime2'
      if(iup1.eq.1 .and. ierr1.eq.0) then
        do 40 i=1,5
   40   itime1(i)=it(i,1)
      endif
      if(iup2.eq.1 .and. ierr2.eq.0) then
        do 42 i=1,5
   42   itime2(i)=it(i,2)
      endif
c
      if(ierr1.ne.0 .or. ierr2.ne.0) then
        ihours=-32767
        return
      endif
c
      do 50 i=1,4
        if(it(i,1).ne.it(i,2)) goto 55
   50 continue
c  no time difference
      ihours=0
      return
c
   55 nt1=1
      nt2=2
      if(it(i,1).gt.it(i,2)) then
        nt1=2
        nt2=1
      endif
      nhh=0
c
      if(it(1,nt1).eq.it(1,nt2)) goto 70
      iy=it(1,nt1)
c  remaining hours first year:
      do 60 im=it(2,nt1),12
        md=mdays(im)
        if(im.eq.2) then
          if(iy/4*4.eq.iy) md=29
          if(iy/100*100.eq.iy .and. iy/400*400.ne.iy) md=28
        endif
        nhh=nhh+md*24
   60 continue
      nhh=nhh-(it(3,nt1)-1)*24-it(4,nt1)
c  one year steps
      do 65 iy=it(1,nt1)+1,it(1,nt2)-1
        nd=365
        if(iy/4*4.eq.iy) nd=366
        if(iy/100*100.eq.iy .and. iy/400*400.ne.iy) nd=365
        nhh=nhh+nd*24
   65 continue
      it(1,nt1)=it(1,nt2)
      it(2,nt1)=1
      it(3,nt1)=1
      it(4,nt1)=0
c
   70 if(it(2,nt1).eq.it(2,nt2)) goto 80
c  remaining hours first month
      iy=it(1,nt1)
      im=it(2,nt1)
      md=mdays(im)
      if(im.eq.2) then
        if(iy/4*4.eq.iy) md=29
        if(iy/100*100.eq.iy .and. iy/400*400.ne.iy) md=28
      endif
      nhh=nhh+(md-it(3,nt1)+1)*24-it(4,nt1)
c  one month steps
      do 75 im=it(2,nt1)+1,it(2,nt2)-1
        md=mdays(im)
        if(im.eq.2) then
          if(iy/4*4.eq.iy) md=29
          if(iy/100*100.eq.iy .and. iy/400*400.ne.iy) md=28
        endif
        nhh=nhh+md*24
   75 continue
      it(2,nt1)=it(2,nt2)
      it(3,nt1)=1
      it(4,nt1)=0
c
   80 if(it(3,nt1).ne.it(3,nt2)) then
        nhh=nhh+(it(3,nt2)-it(3,nt1))*24-it(4,nt1)
        it(3,nt1)=it(3,nt2)
        it(4,nt1)=0
      endif
c  hours last day
      nhh=nhh+it(4,nt2)-it(4,nt1)
      it(4,nt1)=it(4,nt2)
c
      if(nt1.eq.1) then
        ihours=nhh
      else
        ihours=-nhh
      endif
c
      return
      end
C
C ********************************************************************
C
      subroutine vtime(itime,ierror)
c
c  'itime' is updated to give "valid" date,time (with prog.time = 0)
c
c  input:  itime(5) - itime(1): year
c                     itime(2): month (1-12)
c                     itime(3): day (1-28/29/30/31)
c                     itime(4): time in hours (00-23)
c                     itime(5): time in hours of prognosis
c                                     (negative, zero or positive)
c  output: itime(5) -  as above, itime(5)=0
c          ierror   -  0 = o.k. input date/time
c                      1 = not o.k. input date/time
c                          ('itime' not changed)
c
c  DNMI/FoU  xx.xx.1992  Anstein Foss
c
      integer itime(5)
c
      integer mdays(12)
      data mdays/31,28,31,30,31,30,31,31,30,31,30,31/
c
      iy=itime(1)
      im=itime(2)
      id=itime(3)
      ih=itime(4)
c
c  test input time
      ierror=0
      if(im.lt.1 .or. im.gt.12) then
        ierror=1
      else
        md=mdays(im)
        if(im.eq.2) then
          if(iy/4*4.eq.iy) md=29
          if(iy/100*100.eq.iy .and. iy/400*400.ne.iy) md=28
        endif
        if(id.lt.1 .or. id.gt.md) ierror=1
      endif
      if(ih.lt.00 .or. ih.gt.23) ierror=1
c
      if(ierror.ne.0) return
c
      ih=ih+itime(5)
      if(ih.ge.0 .and. ih.le.23) goto 50
      if(ih.lt.0) goto 30
c
      nd=ih/24
      ih=ih-24*nd
      do 20 n=1,nd
        id=id+1
        if(id.gt.md) then
          im=im+1
          if(im.gt.12) then
            iy=iy+1
            im=1
            md=mdays(im)
          elseif(im.eq.2) then
            md=mdays(im)
            if(iy/4*4.eq.iy) md=29
            if(iy/100*100.eq.iy .and. iy/400*400.ne.iy) md=28
          else
            md=mdays(im)
          endif
          id=1
        endif
   20 continue
      goto 50
c
   30 nd=(-ih+23)/24
      ih=ih+24*nd
      do 40 n=1,nd
        id=id-1
        if(id.lt.1) then
          im=im-1
          if(im.lt.1) then
            iy=iy-1
            im=12
            md=mdays(im)
          elseif(im.eq.2) then
            md=mdays(im)
            if(iy/4*4.eq.iy) md=29
            if(iy/100*100.eq.iy .and. iy/400*400.ne.iy) md=28
          else
            md=mdays(im)
          endif
          id=md
        endif
   40 continue
c
   50 itime(1)=iy
      itime(2)=im
      itime(3)=id
      itime(4)=ih
      itime(5)=0
c
      return
      end
c
c
c
      subroutine ancorr(delta,rlon,rlat,xcen,ycen)
c
c     Direction angle correction, from rotated Wam grid to geographic.
c
c     CODED BY: J.E.HAUGEN, A.FOSS   - DNMI/R&D  APRIL 1998
c
      real delta,rlon,rlat,xcen,ycen
c
      zpir18 = 2.0*asin(1.0)/180.
c
      xca = xcen*zpir18
      yca = ycen*zpir18
      x3  = rlon*zpir18
      y3  = rlat*zpir18
c
      call sphrot(x2,y2,x3,y3,1,xca,yca,-1)
c
      zsyca = sin(yca)
      zcyca = cos(yca)
c
      zsxsph = sin(x2)
      zcxsph = cos(x2)
      zsysph = sin(y2)
      zcysph = cos(y2)
      zxmxc  = x2 - xca
      zsxmxc = sin(zxmxc)
      zcxmxc = cos(zxmxc)
      zsxrot = sin(x3)
      zcxrot = cos(x3)
      zsyrot = sin(y3)
      zcyrot = cos(y3)
      za1 = zcxmxc*zcxrot + zcyca*zsxmxc*zsxrot
      za2 = zcyca*zsxmxc*zcxrot*zsyrot + zsyca*zsxmxc*zcyrot -
     +      zcxmxc*zsxrot*zsyrot
      za3 =-zsyca*zsxrot/zcysph
      za4 = (zcyca*zcyrot - zsyca*zcxrot*zsyrot)/zcysph
c
      dda = 45.
      ua  = -1.
      va  = -1.
c
      u = za1*ua + za2*va
      v = za3*ua + za4*va
c
      dd=270.-atan2(v,u)/zpir18
      if(dd.gt.+180.) dd=dd-360.
      if(dd.lt.-180.) dd=dd+360.
c
      delta=dd-dda
      if(delta.gt.+180.) delta=delta-360.
      if(delta.lt.-180.) delta=delta+360.
c
      return
      end
c
c
c
      subroutine sphrot(xsph,ysph,xrot,yrot,n,
     +                  xcen,ycen,icall)
c
c  conversion between spherical (xsph,ysph) and spherical rotated
c  (xrot,yrot) coordinates. (xcen,ycen) is the position of the
c  rotated equator/greenwich in terms of (longitude,latitude).
c  all values are given in radians.
c
      real xsph(n), ysph(n), xrot(n), yrot(n)
c
      zsycen = sin(ycen)
      zcycen = cos(ycen)
c
      if (icall.eq.1) then
c
c  compute spherical rotated coordinates as function of
c  spherical coordinates
c
      do j = 1,n
         zxmxc  = xsph(j) - xcen
         zsxmxc = sin(zxmxc)
         zcxmxc = cos(zxmxc)
         zsysph = sin(ysph(j))
         zcysph = cos(ysph(j))
         zsyrot = zcycen*zsysph - zsycen*zcysph*zcxmxc
         zsyrot = max(zsyrot,-1.0)
         zsyrot = min(zsyrot,+1.0)
         yrot(j) = asin(zsyrot)
         zcyrot = cos(yrot(j))
         zcxrot = (zcycen*zcysph*zcxmxc +
     +             zsycen*zsysph)/zcyrot
         zcxrot = max(zcxrot,-1.0)
         zcxrot = min(zcxrot,+1.0)
         zsxrot = zcysph*zsxmxc/zcyrot
         xrot(j) = acos(zcxrot)
         if (zsxrot.lt.0.0) xrot(j) = -xrot(j)
      enddo
c
      elseif (icall.eq.-1) then
c
c  compute spherical coordinates as function of
c  spherical rotated coordinates
c
      do j = 1,n
         zsxrot = sin(xrot(j))
         zcxrot = cos(xrot(j))
         zsyrot = sin(yrot(j))
         zcyrot = cos(yrot(j))
         zsysph = zcycen*zsyrot + zsycen*zcyrot*zcxrot
         zsysph = max(zsysph,-1.0)
         zsysph = min(zsysph,+1.0)
         ysph(j) = asin(zsysph)
         zcysph = cos(ysph(j))
         zcxmxc = (zcycen*zcyrot*zcxrot -
     +             zsycen*zsyrot)/zcysph
         zcxmxc = max(zcxmxc,-1.0)
         zcxmxc = min(zcxmxc,+1.0)
         zsxmxc = zcyrot*zsxrot/zcysph
         zxmxc  = acos(zcxmxc)
         if (zsxmxc.lt.0.0) zxmxc = -zxmxc
         xsph(j) = zxmxc + xcen
      enddo
c
      else
c
      write(*,'(1x,''invalid icall in sphrot'')')
      stop
c
      endif
c
      return
      end
c
c
c
      subroutine gridpars(icall,ldata,idata,igtype,nx,ny,grid,ierror)
cccccc............gridpar(icall,ldata,idata,igtype,nx,ny,grid,ierror)
c
c  NAME:
c     gridpar
c
c  PURPOSE:
c     Conversion between integer*2 field identification and
c     variables used in programs. Note that some gridtypes also
c     has extended (geometry) identification behind the field data.
c
c  SYNOPSIS:
c     subroutine gridpar(icall,ldata,idata,igtype,nx,ny,grid,ierror)
c     integer   icall,ldata,igtype,nx,ny,ierror
c     integer*2 idata(ldata)
c     real      grid(6)
c
c  INPUT:
c     icall  - +1 : from idata to igtype,nx,ny,grid
c              -1 : from igtype,nx,ny,grid to idata
c     ldata  - length of the idata array
c
c  INPUT/OUTPUT:
c     idata  - field identification, field data (not touched)
c              and possibly extra geometry specifications
c              (max 20 words in this version)
c     igtype - grid type, 1 = polarstereographic grid (true at 60N)
c                         2 = geographic
c                         3 = spherical rotated grid
c                         4 = polarstereographic grid
c                         5 = mercator grid (unrotated)
c                         * = unknown grid type (-32767 to +32 accepted)
c              (note that 'codes' are used in the field identification
c               when extra geometry identification is stored after
c               the field data, the reason for the +32 limit)
c     nx     - no. of gridpoints in x direction (1 - 32767)
c     ny     - no. of gridpoints in y direction (1 - 32767)
c     grid   - grid parameters
c                igtype=1,4, polarstereographic grid:
c                          grid(1) = x position of north pole (xp)
c                          grid(2) = y position of north pole (yp)
c                          grid(3) = no. of grid units between
c                                    North Pole and Equator
c                          grid(4) = grid rotation angle (degrees),
c                                    positive east, negative west
c                          grid(5) = projection latitude
c                                    (degrees) standard is 60 (60 deg. N)
c                          grid(6) = 0.  (not used)
c                igtype=2,3, geographic or spherical rotated grid:
c                          grid(1) = western boundary (degrees)
c                                    (longitude for x=1)
c                          grid(2) = southern boundary (degrees)
c                                    (latitude for y=1)
c                          grid(3) = longitude increment (degrees)
c                          grid(4) = latitude  increment (degrees)
c                          grid(5) = longitude position of rotated equator
c                                    (degrees)  (0 if geographic grid)
c                          grid(6) = latitude  position of rotated equator
c                                    (degrees)  (0 if geographic grid)
c                igtype=5, mercator (unrotated) grid:
c                          grid(1) = western boundary (degrees)
c                                    (longitude for x=1)
c                          grid(2) = southern boundary (degrees)
c                                    (latitude for y=1)
c                          grid(3) = x (longitude) increment (km)
c                          grid(4) = y (latitude)  increment (km)
c                          grid(5) = reference (construction) latitude
c                                    (degrees)
c                          grid(6) = 0.  (not used)
c                igtype=*, unknown grid type,
c                          only use grid type less than 1 if the
c                          grid parameters have no meaning:
c                          grid(1:6) : unknown grid parameters
c
c  OUTPUT:
c     ierror - error status: 0 = no error
c                            1 = bad value in input identification
c                                or in input grid parameters etc.
c                            2 = ldata too small
c                            3 = unknown icall
c                            4 = ldata too small for needed extra
c                                geometry identification (icall=-1),
c                                but the best possible identification
c                                is done (see NOTES below)
c
c  DOCUMENTATION:
c     Basic document on felt files:
c          FILE STRUKTUR FOR "SANNTIDS" LAGRING AV GRID-DATA
c               Forskningsavdeling DNMI, oktober 1982
c     See also /usr/local/doc/felt.doc (at DNMI)
c
c  NOTES:
c   - This routine maintain compability with old formats
c     (see comments in the source code),
c     new formats added when required for some reason.
c   - Specyfing ldata too short to hold extra geometry identification
c     (for icall=-1) will restrict this routine from making this even
c     when it seems appropriate. In order to be compatible with
c     old unextended (or less extended) formats you then may
c     disregard the returned ierror=4 (in special cases).
c   - Avoid calling this routine more often than necessary,
c     i.e. only after reading the first field or only before output
c     of the first field.
c     In models with several calling routines you may do a 'setup' call
c     and keep (part of) the first 20 words of identification and
c     possibly extra geometry specification locally (see note below).
c   - The value of nx and ny (in idata if icall=+1) does not influence
c     conversion of grid parameters.
c     You may then use this routine to convert between field
c     identification and grid parameters with nx=ny=1 and possibly
c     find extra identification in word 22,23,... in idata.
c     (see DOCUMENTATION to find format description).
c
c-----------------------------------------------------------------------
c  DNMI/FoU  05.05.1995  Anstein Foss
c  DNMI/FoU  09.06.1995  Anstein Foss
c  DNMI/FoU  14.05.1996  Anstein Foss ... mercator (unrotated)
c  DNMI/FoU  02.09.1996  Anstein Foss ... gridtype 2012 -> 2 (bad test)
c  DNMI/FoU  15.10.1996  Anstein Foss ... even better scaling when needed
c  DNMI/FoU  17.02.1997  Anstein Foss ... and again (for 'image' fields)
c-----------------------------------------------------------------------
c
      implicit none
c
c..input/output
      integer   icall,ldata,igtype,nx,ny,ierror
      integer*2 idata(ldata)
      real      grid(6)
c
c..local
      integer   ngw,lgeom,i,igr,ig1,ig2,ld,igscale,kgeom,lgeom1,lgeom2
      integer*2 igeom1(12),igeom2(20)
      real      gscale,glog,gws,dglim,dgmax,dgprv,grx
      real      gw(6),gr(6)
c
      ierror=0
c
      if(icall.eq.+1) then
c
c..idata(ldata) -> igtype,nx,ny,grid(6)
c
        igtype=0
        nx=0
        ny=0
        do i=1,6
          grid(i)=0.
        end do
        if(ldata.lt.20) then
          ierror=2
          return
        end if
        igtype=idata(9)
        lgeom =0
        if(igtype.gt.999) then
          i=igtype
          igtype=igtype/1000
          lgeom =i-igtype*1000
        end if
        nx=idata(10)
        ny=idata(11)
        if(nx.lt.1 .or. ny.lt.1) then
          ierror=1
          return
        end if
c
        if(igtype.eq.1 .or. igtype.eq.4) then
c
c..polarstereographic (type 1 always true at 60 degrees North)
c
          if(idata(17).gt.0) then
c..the standard
            grid(1)=idata(15)*0.01
            grid(2)=idata(16)*0.01
            grid(3)=idata(17)*0.1
            grid(4)=idata(18)
          else
c..an old extension
            grid(1)=idata(15)
            grid(2)=idata(16)
            grid(3)=-idata(17)*0.1
            grid(4)=idata(18)
          end if
          grid(5)=60.
c
        elseif(igtype.eq.2 .or. igtype.eq.3) then
c
c..geographic (2) or spherical rotated grid (3)
c
          grid(1)=idata(16)*0.01
          grid(2)=idata(15)*0.01
          grid(3)=idata(18)*0.01
          grid(4)=idata(17)*0.01
c
        elseif(igtype.eq.5) then
c
c..mercator grid (5)
c
          grid(1)=idata(15)*0.01
          grid(2)=idata(16)*0.01
          grid(3)=idata(17)*0.1
          grid(4)=idata(18)*0.1
c
        else
c
c..unknown/undefined grid type
c
          grid(1)=idata(15)
          grid(2)=idata(16)
          grid(3)=idata(17)
          grid(4)=idata(18)
c
        end if
c
        if(lgeom.gt.0) then
c
          if(igtype.eq.1 .or. igtype.eq.4) then
            gscale=100.
            ngw=5
          elseif(igtype.eq.2 .or. igtype.eq.3) then
            gscale=10000.
            ngw=6
          elseif(igtype.eq.5) then
            gscale=10000.
            ngw=5
          else
            gscale=100.
            ngw=6
          end if
c
          ld=20+nx*ny
c
          if(lgeom.eq.ngw*2 .and. ld+lgeom.le.ldata) then
c
c..first extended method
            do i=1,ngw
              ig1=idata(ld+1)
              ig2=idata(ld+2)
              ld=ld+2
              grid(i)=float(ig1*10000+ig2)/gscale
            end do
c
          elseif(lgeom.eq.2+ngw*3 .and. ld+lgeom.le.ldata) then
c
c..second extended method
            if(idata(ld+1).eq.ngw .and. idata(ld+2).eq.3) then
              ld=ld+2
              do i=1,ngw
                igscale=idata(ld+1)
                ig1    =idata(ld+2)
                ig2    =idata(ld+3)
                ld=ld+3
                gscale=10.**igscale
                grid(i)=float(ig1*10000+ig2)*gscale
              end do
            else
              ierror=2
            end if
c
          else
c
            ierror=2
c
          end if
c
        end if
c
        if(ierror.eq.0) then
c
          if(igtype.eq.1 .or. igtype.eq.4) then
c..the DNMI "standard" (150km grid => an=grid(3)=79.)
            if(grid(3).ne.0.) grid(3)=79.*150./grid(3)
            if(grid(3).eq.0.   .or. grid(5).eq.0. .or.
     +         grid(5).lt.-90. .or. grid(5).gt.+90.) ierror=1
          elseif(igtype.eq.2 .or. igtype.eq.3) then
            if(grid(3).eq.0. .or. grid(4).eq.0.) ierror=1
          elseif(igtype.eq.5) then
            if(grid(3).eq.0. .or. grid(4).eq.0.) ierror=1
          end if
c
        end if
c
      elseif(icall.eq.-1) then
c
c..igtype,nx,ny,grid(6) -> idata(ldata)
c
        if(ldata.lt.20) then
          ierror=2
          return
        end if
c
        if(igtype.gt.32 .or.
     +     nx.lt.1 .or. nx.gt.32767 .or.
     +     ny.lt.1 .or. ny.gt.32767) then
          ierror=1
          return
        end if
c
        idata( 9)=igtype
        idata(10)=nx
        idata(11)=ny
        do i=1,6
          gw(i)=grid(i)
        end do
c
        if(igtype.eq.1 .or. igtype.eq.4) then
c
c..polarstereographic (type 1 always true at 60 degrees North)
c
          if(gw(3).eq.0.) then
            ierror=1
            return
          end if
c
c..the DNMI "standard" (150km grid => an=grid(3)=79.)
          gw(3)=150.*79./gw(3)
c
          if(abs(gw(1)).lt.327.66 .and.
     +       abs(gw(2)).lt.327.66) then
c..the standard
            idata(15)=nint(gw(1)*100.)
            idata(16)=nint(gw(2)*100.)
            idata(17)=nint(gw(3)*10.)
            idata(18)=nint(gw(4))
            gr(1)=float(idata(15))*0.01
            gr(2)=float(idata(16))*0.01
            gr(3)=float(idata(17))*0.1
            gr(4)=float(idata(18))
            gr(5)=60.
          elseif(abs(gw(1)).lt.32766. .and.
     +           abs(gw(2)).lt.32766.) then
c..an old extension
            idata(15)=nint(gw(1))
            idata(16)=nint(gw(2))
            idata(17)=-nint(gw(3)*10.)
            idata(18)=nint(gw(4))
            gr(1)=float(idata(15))
            gr(2)=float(idata(16))
            gr(3)=-float(idata(17))*0.1
            gr(4)=float(idata(18))
            gr(5)=60.
          else
c..old impossible case
            idata(15)=32767
            idata(16)=32767
            idata(17)=32767
            idata(18)=32767
            gr(1)=327.67
            gr(2)=327.67
            gr(3)=3276.7
            gr(4)=32767.
            gr(5)=60.
          end if
c
          gscale=100.
          ngw=5
c
        elseif(igtype.eq.2 .or. igtype.eq.3) then
c
c..geographic (2) or spherical rotated grid (3)
c
          if(gw(3).eq.0. .or. gw(4).eq.0.) then
            ierror=1
            return
          end if
c
          if(gw(1).gt.+180.) gw(1)=gw(1)-360.
          if(gw(1).lt.-180.) gw(1)=gw(1)+360.
c
          idata(15)=nint(gw(2)*100.)
          idata(16)=nint(gw(1)*100.)
          idata(17)=nint(gw(4)*100.)
          idata(18)=nint(gw(3)*100.)
          gr(1)=float(idata(16))*0.01
          gr(2)=float(idata(15))*0.01
          gr(3)=float(idata(18))*0.01
          gr(4)=float(idata(17))*0.01
          gr(5)=0.
          gr(6)=0.
c
          gscale=10000.
          ngw=6
c
        elseif(igtype.eq.5) then
c
c..mercator grid (5)
c
          if(gw(3).eq.0. .or. gw(4).eq.0.) then
            ierror=1
            return
          end if
c
          idata(15)=nint(gw(1)*100.)
          idata(16)=nint(gw(2)*100.)
          idata(17)=nint(gw(3)*10.)
          idata(18)=nint(gw(4)*10.)
          gr(1)=float(idata(15))*0.01
          gr(2)=float(idata(16))*0.01
          gr(3)=float(idata(17))*0.1
          gr(4)=float(idata(18))*0.1
          gr(5)=0.
c
          gscale=10000.
          ngw=5
c
        else
c
c..unknown/undefined grid type
c
          idata(15)=nint(gw(1))
          idata(16)=nint(gw(2))
          idata(17)=nint(gw(3))
          idata(18)=nint(gw(4))
          gr(1)=float(idata(15))
          gr(2)=float(idata(16))
          gr(3)=float(idata(17))
          gr(4)=float(idata(18))
          gr(5)=0.
          gr(6)=0.
c
          gscale=100.
          ngw=6
c
        end if
c
c..check if the standard packing above was good enough
c..or if the first or second extended method should be used,
c..for compability with old formats the least extended is preferred
c
c..a limit to avoid high precision complications
        dglim=1.e-8
c
        dgmax=0.
c
        if(igtype.gt.0) then
          do i=1,ngw
            if(gw(i).ne.0.) dgmax=max(dgmax,abs((gr(i)-gw(i))/gw(i)))
          end do
        end if
c
        if(dgmax.gt.dglim) then
c
          kgeom=0
          ld=20+nx*ny
c
c..first extended method, same scaling for all grid parameters
c..but don't use it unless it gives better results
c..(kept in the code for compability with old formats, possibly
c.. in models etc. not using this routine)
c
          dgprv=dgmax
          dgmax=0.
          lgeom1=0
c
          do i=1,ngw
            gws=gw(i)*gscale
c..check overflow
            if(gws.ne.0. .and. abs(gws).lt.3.e+8) then
              igr=nint(gws)
              ig1=igr/10000
              ig2=igr-ig1*10000
              igeom1(lgeom1+1)=ig1
              igeom1(lgeom1+2)=ig2
              grx=float(ig1*10000+ig2)/gscale
              dgmax=max(dgmax,abs((grx-gw(i))/gw(i)))
            else
c..zero or value not possible for this method
              igeom1(lgeom1+1)=0
              igeom1(lgeom1+2)=0
              if(gws.ne.0.) dgmax=1.
            end if
            lgeom1=lgeom1+2
          end do
c
          if(dgmax.lt.dgprv) then
            if(ldata.ge.ld+lgeom1) then
              kgeom=1
            else
              ierror=4
            end if
          else
            dgmax=dgprv
          end if
c
          if(dgmax.gt.dglim .and. ierror.eq.0) then
c
c..second extended method, good enough for any real*4 precision value
c..but don't use it unless it gives better results
c
            dgprv=dgmax
            dgmax=0.
c
            igeom2(1)=ngw
            igeom2(2)=3
            lgeom2=2
            do i=1,ngw
              if(gw(i).ne.0.) then
c..8 decimals precision (more than enough for real*4)
                glog=log10(abs(gw(i)))-8.
                igscale=int(glog)
                if(glog.gt.0.) igscale=igscale+1
c..keep scaling within real*4 range
                igscale=min(igscale,+25)
                igscale=max(igscale,-25)
                gscale=10.**(-igscale)
                igr=nint(gw(i)*gscale)
                ig1=igr/10000
                ig2=igr-ig1*10000
                igeom2(lgeom2+1)=igscale
                igeom2(lgeom2+2)=ig1
                igeom2(lgeom2+3)=ig2
                grx=float(ig1*10000+ig2)/gscale
                dgmax=max(dgmax,abs((grx-gw(i))/gw(i)))
              else
                igeom2(lgeom2+1)=0
                igeom2(lgeom2+2)=0
                igeom2(lgeom2+3)=0
              end if
              lgeom2=lgeom2+3
            end do
c
            if(dgmax.lt.dgprv) then
              if(ldata.ge.ld+lgeom2) then
                kgeom=2
              else
                ierror=4
              end if
ccc         else
ccc           dgmax=dgprv
            end if
c
          end if
c
          if(kgeom.eq.1) then
            idata(9)=igtype*1000+lgeom1
            do i=1,lgeom1
              idata(ld+i)=igeom1(i)
            end do
          elseif(kgeom.eq.2) then
            idata(9)=igtype*1000+lgeom2
            do i=1,lgeom2
              idata(ld+i)=igeom2(i)
            end do
          end if
c
        end if
c
c
      else
c
c..wrong icall
c
        ierror=3
c
      end if
c
      return
      end
C      subroutine earthr(rearth)
c
c  NAME:
c     earthr
c
c  PURPOSE:
c     Return earth radius in unit meter,
c     i.e. the DNMI standard value.
c
c  SYNOPSIS:
c     subroutine earthr(rearth)
c     real rearth
c
c  OUTPUT:
c     rearth  - earth radius i unit meter
c
c-----------------------------------------------------------------------
c  DNMI/FoU  xx.xx.19xx  ............ ... rearth = 6368.00 km in models
c  DNMI/FoU  25.08.1995  Anstein Foss ... rearth = 6371.22 km
c  DNMI/FoU  21.03.1996  Anstein Foss ... rearth = 6371.00 km
c-----------------------------------------------------------------------
c
C      real rearth
c
C      rearth = 6371000.
c
C      return
c      end
