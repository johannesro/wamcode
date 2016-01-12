!
!########################################################
!#							#
!#	NetCDF Outputmodule for radiation stress,       #
!#             wave force and stokes drift              #
!#                                                      #
!#      Wolfgang Koch  GKSS    July 2009                #
!#      Arno Behrens   GKSS    October 2010             #
!#      Arno Behrens   HZG     June 2014                #
!#                             rad parameters           #
!#      Ana Carrsco    METNo   Jan 2015                 #
!#                             add metno attributes     #       
!#							#
!########################################################
!
module wam_rad_netcdf_metno_module

use wam_file_module,          only: iu06
use wam_print_module,         only: nx, ny, amosop, amowep, xdello, xdella
use wam_output_set_up_module, only: idelint
use wam_general_module,       only: difdate

    
implicit none
private
integer, parameter :: nf=8              !! maximum number of fields
integer, save :: nc = -1
logical, save :: fc = .true.		!! erster Aufruf
integer, save :: no = -1		!! Offset beim schreiben
integer, save :: ncid(0:nf+3) = -1
character (len=14), save :: fsd

real*8, allocatable, dimension (:) :: lat
real*8, allocatable, dimension (:) :: lon
character (len=14)       :: Basedate
real                     :: hradd      
     
public wknco, wkncw, wkncc, pf, nx, ny, hradd

contains

!##################################################
!# 	open NetCDF file                          #
!#                                                #
!#	name	filename                          #
!#	n	number of longitudes              #
!#	m	number of latitudes               #
!#	l	maximum number of depth layers    #
!#	lon	longitudes                        #
!#	lat	latitudes                         #
!#	ncid	NetCDF file ID-number             #
!#	sd	date [yyyymmddhh]                 #
!#	dt	time step [s]                     #
!#	flg	parameter switch                  #
!##################################################

subroutine wknco (name, sd, cflag_p, amowep, amosop, xdella, xdello)
use netcdf
    
implicit none
character (len=16), intent(in) :: name
character (len=14), intent(in) :: sd
logical, intent(in) :: cflag_p(:)
logical, dimension(5) :: ismember
integer, dimension(5) :: metn

real :: amowep, amosop, xdella, xdello

character (len=33) :: tua
character (len=21) :: tda
character (len=42) :: metno
character (len=100), dimension (4,nf) :: vl  
     
integer, dimension (3)   :: diid
integer, dimension (0:2) :: did
integer, dimension (nf)  :: flg, ty
integer :: n, m, l, i, j, loop
    
real*8 :: dtor, dt
real*4, dimension (nf) :: fw, vmin, vmax

metno='Norwegian Meteorological Institute, MET Norway'    
allocate (lat(ny), lon(nx))
dt = real(idelint)
lat(1) = amosop
do loop=2,ny
   lat(loop) = lat(loop-1)+xdella    !! latitudes
enddo
lon(1) = amowep
do loop=2,nx
   lon(loop) = lon(loop-1)+xdello    !! lonitudes
enddo
 
IF (ncid(0)<0) THEN
   fsd = sd
   vl  = ''                          !! integrated parameters
   vl(1, 1) = 'sxx'
   vl(2, 1) = 'radiation_stress_tensor_sxx'
   vl(3, 1) = 'radiation stress tensor sxx'
   vl(4, 1) = 'kg/s*s'
   vmin(1)  = -1.e07
   vmax(1)  = 1.e07
    
   vl(1, 2) = 'syy'
   vl(2, 2) = 'radiation_stress_tensor_syy'
   vl(3, 2) = 'radiation stress tensor syy'
   vl(4, 2) = 'kg/s*s'
   vmin(2)  = -1.e07
   vmax(2)  = 1.e07
    
   vl(1, 3) = 'sxy'
   vl(2, 3) = 'radiation_stress_tensor_sxy'
   vl(3, 3) = 'radiation stress tensor sxy'
   vl(4, 3) = 'kg/s*s'
   vmin(3)  = -1.e07
   vmax(3)  = 1.e07

   vl(1, 5) = 'wfx'
   vl(2, 5) = 'wave_force_x_component'
   vl(3, 5) = 'wave force x-component'
   vl(4, 5) = 'N/m*m'
   vmin(5)  = -1.
   vmax(5)  = 1.

   vl(1, 6) = 'wfy'
   vl(2, 6) = 'wave_force_y_component'
   vl(3, 6) = 'wave force y-component'
   vl(4, 6) = 'N/m*m'
   vmin(6)  = -1
   vmax(6)  = 1.
    
   vl(1,7) = 'sdx'
   vl(2,7) = 'sea_surface_wave_stokes_drift_x_velocity'
   vl(3,7) = 'stokes drift x-component'
   vl(4,7) = 'm/s'
   vmin(7) = -1.
   vmax(7) = 1.
    
   vl(1,8) = 'sdy'
   vl(2,8) = 'sea_surface_wave_stokes_drift_y_velocity'
   vl(3,8) = 'stokes drift y-component'
   vl(4,8) = 'm/s'
   vmin(8) = -1.
   vmax(8) = 1.
   
   flg = 1                   !! prepare table for required parameters only
   do loop=1,nf
      if (.not.cflag_p(loop)) then
         flg([loop]) = 0
      endif
   enddo
      
   ty = NF90_FLOAT           !! type of parameter
   fw = -999999.             !! dummy (zmiss)
   did = [3,2,1]
   i = dt
   WRITE(tda,'("0000-00-00 (",2(i2.2,":"),i2.2,")")')i/3600,MOD(i/60,60),MOD(i,60)
   Basedate = '19700101000000'
   call DIFFDATEHR(Basedate, sd ,hradd)
   
   tua="hours since "//Basedate(1:4)//"-"//Basedate(5:6)//"-"//Basedate(7:8)//" "//Basedate(9:10)//":"//Basedate(11:12)//":"//Basedate(13:14)
    
   CALL Pf(NF90_CREATE(name,ior(NF90_CLOBBER,NF90_SHARE),ncid(0)))
   CALL Pf(NF90_DEF_DIM(ncid(0),'time',NF90_UNLIMITED,diid(1)))
   CALL Pf(NF90_DEF_DIM(ncid(0),'lat',ny,diid(2)))
   CALL Pf(NF90_DEF_DIM(ncid(0),'lon',nx,diid(3)))
   CALL Pf(NF90_DEF_VAR(ncid(0),'lat',NF90_FLOAT,[diid(2)],ncid(nf+1)))
   CALL Pf(NF90_DEF_VAR(ncid(0),'lon',NF90_FLOAT,[diid(3)],ncid(nf+2)))
   CALL Pf(NF90_DEF_VAR(ncid(0),'time',NF90_FLOAT,[diid(1)],ncid(nf+3)))
   CALL Pf(NF90_PUT_ATT(ncid(0),ncid(nf+3),'standard_name','time'))
   CALL Pf(NF90_PUT_ATT(ncid(0),ncid(nf+3),'calendar',"gregorian"))
   CALL Pf(NF90_PUT_ATT(ncid(0),ncid(nf+3),'delta_t',tda))
   CALL Pf(NF90_PUT_ATT(ncid(0),ncid(nf+3),'units',tua))
   CALL Pf(NF90_PUT_ATT(ncid(0),ncid(nf+3),'dt',i))

   CALL Pf(NF90_PUT_ATT(ncid(0),NF90_GLOBAL,'Conventions','CF-1.6'))
   CALL Pf(NF90_PUT_ATT(ncid(0),NF90_GLOBAL,'institution',metno))

     metn=(/ 1, 2, 3, 5, 6/)
 
   DO i=1,nf
      IF (flg(i)>0) THEN
         CALL Pf(NF90_DEF_VAR(ncid(0),vl(1,i),ty(i),diid(did),ncid(i)))
         CALL Pf(NF90_PUT_ATT(ncid(0),ncid(i),'_FillValue',fw(i)))
         call pf(nf90_put_att(ncid(0),ncid(i),'long_name',vl(3,i)))
          do j=1,size(metn)
            if (i.eq.metn(j)) then
               ismember(j)=.true.
            else
               ismember(j)=.false.
            endif
         end do
         IF (ANY(ismember)) THEN 
            call pf(nf90_put_att(ncid(0),ncid(i),'metno_name',vl(2,i)))
         ELSE
            call pf(nf90_put_att(ncid(0),ncid(i),'standard_name',vl(2,i)))
         ENDIF
         call pf(nf90_put_att(ncid(0),ncid(i),'units',vl(4,i)))
         !call pf(nf90_put_att(ncid(0),ncid(i),'valid_min',vmin(i)))
         !call pf(nf90_put_att(ncid(0),ncid(i),'valid_max',vmax(i)))
      ENDIF
   ENDDO
   CALL Pf(NF90_ENDDEF(ncid(0)))
   CALL Pf(NF90_PUT_VAR(ncid(0),ncid(nf+1),lat))
   CALL Pf(NF90_PUT_VAR(ncid(0),ncid(nf+2),lon))
ELSE
   write (iu06,*) ' +++ NetCDF-file is open already'
   STOP           ' +++ Problem occurs during NetCDF-Output handling '
ENDIF
END subroutine wknco

!########################################################
!#      write a field to NetCDF file                    #
!#                                                      #
!#      nid     pointer to NetCDF file ID-number        #
!#      vid     pointer to NetCDF parameter ID-number   #
!########################################################

subroutine wkncw (ip,grid,ad,t)
use netcdf
    
implicit none
real,    intent(in)    :: grid(:,:)   !! grid of parameter
integer, intent(in)    :: ip          !! parameter number
integer, dimension (3) :: sta
    
character (len=14), intent(in) :: ad
integer :: vid, j, t
real    :: tda


if (fc) then
   no = 1-t                           !! t+no=1 : first output
   fc = .false.
endif
tda = real(idelint)/(60.*60.)
vid = ip
if (ncid(vid)>0) then
   sta = [1,1,t+no]
   CALL Pf(NF90_PUT_VAR(ncid(0),ncid(vid),grid,sta))
   CALL Pf(NF90_PUT_VAR(ncid(0),ncid(nf+3),t*tda + hradd,[t+no])) !time
else
   write (iu06,*) ' +++ parameter ',vid,' is not available in NetCDF-file'
endif
end subroutine wkncw

!##############################
!# 	close NetCDF-file     #
!##############################

subroutine wkncc
use netcdf
if (ncid(0)>=0) CALL Pf (NF90_CLOSE(ncid(0)))
ncid = -1
end subroutine wkncc

!################################################
!# 	pf      write a NetCDF-error message    #
!#	en	NetCDF error number             #
!################################################

subroutine pf(en)
use netcdf
integer :: en
if (en/=0) write (iu06,*) ' +++ Error : ', NF90_STRERROR(en)
end subroutine pf

!######################################
!#                                    #
!# Calculates the difference in hrs   #                       
!# between two dates CDATE1, CDATE2   #                         
!#                                    #                      
!#  Ana Carrasco METNO Jan 2015       #                      
!#                                    #                      
!######################################

SUBROUTINE DIFFDATEHR(CDATE1, CDATE2, hradd)
character (len=14) :: CDT1, CDT2, CDT0,CDATE1,CDATE2
integer itime1(6),itime2(6),itime0(6), SEC,i,yn
real    hradd
integer*8 addsec

!  CHANGE DATE TIME TO ENSURE THAT THE SECOND IS LARGER.         
!

IF (CDATE1 .GT. CDATE2) THEN
   CDT1 = CDATE2
   CDT2 = CDATE1
ELSE IF (CDATE2 .GT. CDATE1) THEN
   CDT1 = CDATE1
   CDT2 = CDATE2
END IF


READ(CDT1,'(I4,5I2)')itime1(1),itime1(2),itime1(3),itime1(4),itime1(5),itime1(6)
READ(CDT2,'(I4,5I2)')itime2(1),itime2(2),itime2(3),itime2(4),itime2(5),itime2(6)

itime0(2)=itime1(2)
itime0(3)=itime1(3)
itime0(4)=itime1(4)
itime0(5)=itime1(5)
itime0(6)=itime1(6)
!LOOP OVER YEARS
yn=abs(itime1(1)-itime2(1))
addsec=0

if (yn.gt.1) then
 do i = 1,yn-2
   !add a year to date CDT1
   itime0(1)=itime1(1)+i
   WRITE (CDT0,'(I4.4,5I2.2)')itime0(1),itime0(2),itime0(3),itime0(4),itime0(5),itime0(6)  
   call difdate (CDT1, CDT0, SEC)   !DATES HAVE TO BE IN CONSECUTIVE YEARS   
   addsec=addsec+SEC 
   CDT1=CDT0
 enddo
endif

call difdate (CDT1, CDT2, SEC)    !DATES HAVE TO BE IN CONSECUTIVE YEARS.  
write(iu06,*)'SEC ',SEC

addsec=addsec+SEC 
write(iu06,*)'addsec ',addsec
hradd=real(addsec)/(60.0*60.0)
write(iu06,*)'hradd',hradd
end subroutine DIFFDATEHR

end module wam_rad_netcdf_metno_module
