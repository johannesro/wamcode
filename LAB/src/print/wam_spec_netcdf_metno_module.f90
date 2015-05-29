!
!########################################################
!#							#
!#	NetCDF Outputmodule				#
!#                                                      #
!#      Wolfgang Koch  GKSS    July 2009                #
!#      Arno Behrens   GKSS    October 2010             #
!#	Ana Carrasco   METNO   September 2014           #
!#                                                      #
!#                                                      # 
!#							#         
!#                                                      #
!########################################################
!
module wam_spec_netcdf_metno_module

use wam_file_module,          only: iu06
use wam_print_module,         only: KL, ML, CO, FR, THETA, IDELDO, NOUTP 



    
implicit none
private
integer, parameter :: nf=40             !! maximum number of MAP fields
integer, parameter :: nfs=3             !! maximum number of SPEC fields
integer, save :: nc = -1
logical, save :: fc = .true.		!! erster Aufruf
integer, save :: no = -1		!! Offset beim schreiben
integer, save :: ncid(0:nf+3) = -1
integer, save :: ncidspc(0:nfs+3) = -1
character (len=14), save :: fsd
     
public wkncospc, wkncwspc, wknccspc, pf, FR, THETA 

contains

!
! ---------------------------------------------------------------------------- !
subroutine wkncospc (name, sd, cflag_s, SPECLONG, SPECLAT)

use netcdf
!########################################################
!#      write header of a  NetCDF file                    #
!########################################################
    
implicit none
character (len=16), intent(in) :: name
character (len=14), intent(in) :: sd
logical, intent(in) :: cflag_s(:)
character (len=33) :: tua   
character (len=21) :: tda
character (len=42) :: metno
character (len=100), dimension (4,nfs) :: vl  
integer, dimension (5)   :: diid
integer, dimension (0:4) :: did
integer, dimension (nfs)  :: flg, ty
integer :: n, m, l, i, loop,II
real*8 :: dtor, dt
real*4, dimension (nfs) :: fw, vmin, vmax
REAL,               INTENT(IN) :: SPECLONG(:)  
REAL,               INTENT(IN) :: SPECLAT(:)
integer :: yy
integer, allocatable, dimension (:) :: xx

metno='Norwegian Meteorological Institute, met.no'
allocate (xx(size(SPECLONG)))
xx = (/ (II, II=1,size(SPECLONG)) /)
yy = 1
  
dt = real(IDELDO)
IF (ncidspc(0)<0) THEN
   fsd = sd
   vl  = ''                          !! integrated parameters
   vl(1, 1) = 'SPEC'
   vl(2, 1) = '2D_SPECTRUM'
   vl(3, 1) = '2-D spectrum of total sea'
   vl(4, 1) = 'm**2 s  '
   vmin(1)  =  0.00
   vmax(1)  =  100.0
    
   vl(1, 2) = 'SPEC_SEA'
   vl(2, 2) = '2D_SPECTRUM_OF_WIND_SEA'
   vl(3, 2) = 'Wind Sea 2-D spectrum'
   vl(4, 2) = 'm**2 s'
   vmin(2)  =  0.0
   vmax(2)  =  100.0

   vl(1, 3) = 'SPEC_SWELL'
   vl(2, 3) = '2D_SPECTRUM_OF_WIND_SWELL'
   vl(3, 3) = 'Swell 2-D spectrum'
   vl(4, 3) = 'm**2 s'
   vmin(3)  =  0.0
   vmax(3)  =  100.0


   flg = 1                   !! prepare table for required parameters only
   do loop=1,nfs
      if (.not.cflag_s(loop)) then
         flg([loop]) = 0
      endif
   enddo

   ty = NF90_FLOAT           !! type of parameter
   fw = -99999.              !! dummy (zmiss)
   did = [5,4,3,2,1]
   i = dt
   WRITE(tda,'("0000-00-00 (",2(i2.2,":"),i2.2,")")')i/3600,MOD(i/60,60),MOD(i,60)
   tua="seconds since "//sd(1:4)//"-"//sd(5:6)//"-"//sd(7:8)//" "//sd(9:10)//":"//sd(11:12)//":"//sd(13:14)

   CALL Pf(NF90_CREATE(name,ior(NF90_CLOBBER,NF90_SHARE),ncidspc(0)))
   CALL Pf(NF90_DEF_DIM(ncidspc(0),'time',NF90_UNLIMITED,diid(1))) 
   CALL Pf(NF90_DEF_DIM(ncidspc(0),'x',NOUTP ,diid(3)))
   CALL Pf(NF90_DEF_DIM(ncidspc(0),'y',1 ,diid(2)))
   CALL Pf(NF90_DEF_DIM(ncidspc(0),'freq',ML,diid(4)))
   CALL Pf(NF90_DEF_DIM(ncidspc(0),'direction',KL,diid(5)))
 
   CALL Pf(NF90_DEF_VAR(ncidspc(0),'direction',NF90_FLOAT,[diid(5)],ncidspc(nfs+1)))
   call pf(nf90_put_att(ncidspc(0),ncidspc(nfs+1),'units','degree'))

   CALL Pf(NF90_DEF_VAR(ncidspc(0),'freq',NF90_FLOAT,[diid(4)],ncidspc(nfs+2)))
   call pf(nf90_put_att(ncidspc(0),ncidspc(nfs+2),'units','1/s'))

   CALL Pf(NF90_DEF_VAR(ncidspc(0),'x',NF90_INT,[diid(3)],ncidspc(nfs+3)))
   call pf(nf90_put_att(ncidspc(0),ncidspc(nfs+3),'axis','x'))
   call pf(nf90_put_att(ncidspc(0),ncidspc(nfs+3),'long_name','x-coordinate in Cartesian system'))
   call pf(nf90_put_att(ncidspc(0),ncidspc(nfs+3),'standard_name','projection_x_coordinate')) 
   call pf(nf90_put_att(ncidspc(0),ncidspc(nfs+3),'units','1'))


   CALL Pf(NF90_DEF_VAR(ncidspc(0),'y',NF90_INT,[diid(2)],ncidspc(nfs+4)))
   call pf(nf90_put_att(ncidspc(0),ncidspc(nfs+4),'axis','y'))
   call pf(nf90_put_att(ncidspc(0),ncidspc(nfs+4),'long_name','y-coordinate in Cartesian system'))
   call pf(nf90_put_att(ncidspc(0),ncidspc(nfs+4),'standard_name','projection_y_coordinate')) 
   call pf(nf90_put_att(ncidspc(0),ncidspc(nfs+4),'units','1'))

   CALL Pf(NF90_DEF_VAR(ncidspc(0),'time',NF90_INT,[diid(1)],ncidspc(nfs+5)))
   CALL Pf(NF90_PUT_ATT(ncidspc(0),ncidspc(nfs+5),'delta_t',tda))
   CALL Pf(NF90_PUT_ATT(ncidspc(0),ncidspc(nfs+5),'units',tua))
   CALL Pf(NF90_PUT_ATT(ncidspc(0),ncidspc(nfs+5),'dt',i))

   CALL Pf(NF90_DEF_VAR(ncidspc(0),'latitude',NF90_FLOAT,[diid(3), diid(2)],ncidspc(nfs+6)))  
   call pf(nf90_put_att(ncidspc(0),ncidspc(nfs+6),'long_name','latitude'))
   call pf(nf90_put_att(ncidspc(0),ncidspc(nfs+6),'standard_name','latitude'))
   call pf(nf90_put_att(ncidspc(0),ncidspc(nfs+6),'units','degree_north'))

   CALL Pf(NF90_DEF_VAR(ncidspc(0),'longitude',NF90_FLOAT,[diid(3), diid(2)],ncidspc(nfs+7)))  
   call pf(nf90_put_att(ncidspc(0),ncidspc(nfs+7),'long_name','longitude'))
   call pf(nf90_put_att(ncidspc(0),ncidspc(nfs+7),'standard_name','longitude'))
   call pf(nf90_put_att(ncidspc(0),ncidspc(nfs+7),'units','degree_east'))

   CALL Pf(NF90_PUT_ATT(ncidspc(0),NF90_GLOBAL,'Conventions','CF-1.0'))
   CALL Pf(NF90_PUT_ATT(ncidspc(0),NF90_GLOBAL,'institution',metno))
   
   DO i=1,nfs
      IF (flg(i)>0) THEN
         CALL Pf(NF90_DEF_VAR(ncidspc(0),vl(1,i),ty(i),diid(did),ncidspc(i)))
         CALL Pf(NF90_PUT_ATT(ncidspc(0),ncidspc(i),'_FillValue',fw(i)))
         call pf(nf90_put_att(ncidspc(0),ncidspc(i),'long_name',vl(3,i)))
         call pf(nf90_put_att(ncidspc(0),ncidspc(i),'standard_name',vl(2,i)))
         call pf(nf90_put_att(ncidspc(0),ncidspc(i),'units',vl(4,i)))
!         call pf(nf90_put_att(ncidspc(0),ncidspc(i),'valid_min',vmin(i)))
!         call pf(nf90_put_att(ncidspc(0),ncidspc(i),'valid_max',vmax(i)))
      ENDIF
   ENDDO

   CALL Pf(NF90_ENDDEF(ncidspc(0)))
   CALL Pf(NF90_PUT_VAR(ncidspc(0),ncidspc(nfs+1),THETA*180/3.14159265359))   
   CALL Pf(NF90_PUT_VAR(ncidspc(0),ncidspc(nfs+2),FR))
   CALL Pf(NF90_PUT_VAR(ncidspc(0),ncidspc(nfs+3),xx))
   CALL Pf(NF90_PUT_VAR(ncidspc(0),ncidspc(nfs+4),yy))
   CALL Pf(NF90_PUT_VAR(ncidspc(0),ncidspc(nfs+7),SPECLONG))
   CALL Pf(NF90_PUT_VAR(ncidspc(0),ncidspc(nfs+6),SPECLAT))
  
ELSE
   write (iu06,*) ' +++ NetCDF-file is open already'
   STOP           ' +++ Problem occurs during NetCDF-Output handling '
ENDIF
 
   
END subroutine  wkncospc 
! ---------------------------------------------------------------------------- !

subroutine wkncwspc (ip,grid,ad,t,IPOS)
use netcdf
    
implicit none
real,    intent(in)    :: grid(:,:)   !! grid of parameter
integer, intent(in)    :: ip          !! parameter number
integer, dimension (5) :: start,count
integer                :: IPOS    
character (len=14), intent(in) :: ad
integer :: vid, j, t, tda


if (fc) then
   no = 1-t                           !! t+no=1 : first output
   fc = .false.
endif
tda = IDELDO
vid = ip

write(iu06,*)'IPOS',IPOS,'t*tda,[t+no]',t*tda,[t+no],' MAX =',maxval(grid)
write(iu06,*)'size1',size(grid,1), 'size2',size(grid,2)

if (ncidspc(vid)>0) then
   start = [1,   1, IPOS , 1, t+no]
   count   = [size(grid,1), size(grid,2), 1 , 1, 1]
   write(iu06,*)'start',start
   CALL Pf(NF90_PUT_VAR(ncidspc(0),ncidspc(vid),grid,start,count))

   CALL Pf(NF90_PUT_VAR(ncidspc(0),ncidspc(nfs+5),t*tda,[t+no]))
else
   write (iu06,*) ' +++ parameter ',vid,' is not available in NetCDF-file'
endif

end subroutine wkncwspc

!##
!##############################
!# 	close NetCDF-file     #
!##############################


subroutine wknccspc
use netcdf
if (ncidspc(0)>=0) CALL Pf (NF90_CLOSE(ncidspc(0)))
ncid = -1
end subroutine wknccspc

!################################################
!# 	pf      write a NetCDF-error message    #
!#	en	NetCDF error number             #
!################################################

subroutine pf(en)
use netcdf
integer :: en
if (en/=0) write (iu06,*) ' +++ Error : ', NF90_STRERROR(en)
end subroutine pf
end module wam_spec_netcdf_metno_module
