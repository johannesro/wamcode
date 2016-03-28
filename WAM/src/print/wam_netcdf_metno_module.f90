!
!########################################################
!#							#
!#	NetCDF Outputmodule				#
!#                                                      #
!#      Wolfgang Koch  GKSS    July 2009                #
!#      Arno Behrens   GKSS    October 2010             #
!#      Arno Behrens   HZG     June 2014                #
!#                             additional parameters    #
!#	Ana Carrasco   METno   Jan 2015                 #
!#                             add METno and CF1.6 Conv.#
!#                             attributes               #
!########################################################
!
module wam_netcdf_metno_module

use wam_file_module,          only: iu06
use wam_print_module,         only: nx, ny, amosop, amowep, xdello, xdella
use wam_output_set_up_module, only: idelint
use wam_general_module,       only: difdate
   
implicit none
private
integer, parameter :: nf=40             !! maximum number of fields
integer, save :: nc = -1
logical, save :: fc = .true.		!! erster Aufruf
integer, save :: no = -1		!! Offset beim schreiben
integer, save :: ncid(0:nf+3) = -1
character (len=14), save :: fsd

real*8, allocatable, dimension (:) :: lat
real*8, allocatable, dimension (:) :: lon
character (len=14)       :: Basedate
real*8                     :: hradd     
                 
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
character (len=17), intent(in) :: name
character (len=14), intent(in) :: sd
logical, intent(in) :: cflag_p(:)
logical, dimension(19) :: ismember
integer, dimension(19) :: metn

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

metno='Norwegian Meteorological Institute'    
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
   vl(1, 1) = 'ff'
   vl(2, 1) = 'wind_speed'
   vl(3, 1) = 'Wind speed'
   vl(4, 1) = 'm s-1'
   vmin(1)  = 0.
   vmax(1)  = 40.
    
   vl(1, 2) = 'dd'
   vl(2, 2) = 'wind_to_direction'
   vl(3, 2) = 'Wind direction'
   vl(4, 2) = 'degree'
   vmin(2)  = 0.
   vmax(2)  = 360.
    
   vl(1, 3) = 'FV'
   vl(2, 3) = 'FV'
   vl(3, 3) = 'friction velocity'
   vl(4, 3) = 'm s-1'
   vmin(3)  = 0.
   vmax(3)  = 3.

   vl(1, 4) = 'DC'
   vl(2, 4) = 'DC'
   vl(3, 4) = 'drag coefficient'
   vl(4, 4) = ''
   vmin(4)  = 0.
   vmax(4)  = 0.01

   vl(1, 5) = 'CP'
   vl(2, 5) = 'Charnock_parameter'
   vl(3, 5) = 'Charnock parameter'
   vl(4, 5) = ''
   vmin(5)  = 0.
   vmax(5)  = 0.1

   vl(1, 6) = 'DEPTH'
   vl(2, 6) = 'water depth'
   vl(3, 6) = 'water depth'
   vl(4, 6) = 'm'
   vmin(6)  = 0
   vmax(6)  = 8000

   vl(1, 9) = 'hs'
   vl(2, 9) = 'sea_surface_wave_significant_height'
   vl(3, 9) = 'Total significant wave height'
   vl(4, 9) = 'm'
   vmin(9)  = 0.
   vmax(9)  = 20.
    
   vl(1,10) = 'tp'
   vl(2,10) = 'sea_surface_wave_peak_period_from_variance_spectral_density'
   vl(3,10) = 'Total peak period'
   vl(4,10) = 's'
   vmin(10) = 1.
   vmax(10) = 30.
    
   vl(1,11) = 'tmp'
   vl(2,11) = 'sea_surface_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment'
   vl(3,11) = 'Total mean period'
   vl(4,11) = 's'
   vmin(11) = 1.
   vmax(11) = 20.
   
   vl(1,12) = 'tm1'
   vl(2,12) = 'sea_surface_wave_mean_period_from_variance_spectral_density_first_frequency_moment'
   vl(3,12) = 'Total m1-period'
   vl(4,12) = 's'
   vmin(12) = 1.
   vmax(12) = 20.

   vl(1,13) = 'tm2'
   vl(2,13) = 'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment'
   vl(3,13) = 'Total m2-period'
   vl(4,13) = 's'
   vmin(13) = 1.
   vmax(13) = 20.
   
   vl(1,14) = 'thq'
   vl(2,14) = 'sea_surface_wave_to_direction'
   vl(3,14) = 'Total mean wave direction'
   vl(4,14) = 'degree'
   vmin(14) = 0.
   vmax(14) = 360.
     
   vl(1,15) = 'ds'
   vl(2,15) = 'sea_surface_wave_directional_spread'
   vl(3,15) = 'Total directional spreed'
   vl(4,15) = 'degree'
   vmin(15) = 0.
   vmax(15) = 120.
    
   vl(1,16) = 'NWS'
   vl(2,16) = 'normalised_wave_stress'
   vl(3,16) = 'normalised_wave_stress'
   vl(4,16) = ''
   vmin(16) = 0.
   vmax(16) = 2.
    
   vl(1,17) = 'hs_sea'
   vl(2,17) = 'sea_surface_wind_wave_significant_height'
   vl(3,17) = 'Sea significant wave height'
   vl(4,17) = 'm'
   vmin(17) = 0.
   vmax(17) = 20.
   
   vl(1,18) = 'tp_sea'
   vl(2,18) = 'sea_surface_windsea_wave_peak_period_from_variance_spectral_density'
   vl(3,18) = 'Sea peak period'
   vl(4,18) = 's'
   vmin(18) = 1.
   vmax(18) = 20.

   vl(1,19) = 'tmp_sea'
   vl(2,19) = 'sea_surface_wind_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment'
   vl(3,19) = 'Sea mean period'
   vl(4,19) = 's'
   vmin(19) = 1.
   vmax(19) = 20.
     
   vl(1,20) = 'tm1_sea'
   vl(2,20) = 'sea_surface_wind_wave_mean_period_from_variance_spectral_density_first_frequency_moment'
   vl(3,20) = 'Sea m1-period'
   vl(4,20) = 's'
   vmin(20) = 1.
   vmax(20) = 20.
    
   vl(1,21) = 'tm2_sea'
   vl(2,21) = 'sea_surface_wind_wave_mean_period_from_variance_spectral_density_second_frequency_moment'
   vl(3,21) = 'Sea m2-period'
   vl(4,21) = 's'
   vmin(21) = 1.
   vmax(21) = 20.
    
   vl(1,22) = 'thq_sea'
   vl(2,22) = 'sea_surface_wind_wave_to_direction'
   vl(3,22) = 'Sea mean wave direction'
   vl(4,22) = 'degree'
   vmin(22) = 0.
   vmax(22) = 360.
    
   vl(1,23) = 'ds_sea'
   vl(2,23) = 'sea_surface_wind_wave_directional_spread'
   vl(3,23) = 'Sea directional spreed'
   vl(4,23) = 'degree'
   vmin(23) = 0.
   vmax(23) = 120.
   
   vl(1,25) = 'hs_swell'
   vl(2,25) = 'sea_surface_swell_wave_significant_height'
   vl(3,25) = 'Swell significant wave height'
   vl(4,25) = 'm'
   vmin(25) = 0.
   vmax(25) = 20.
    
   vl(1,26) = 'tp_swell'
   vl(2,26) = 'sea_surface_swell_wave_peak_period_from_variance_spectral_density'
   vl(3,26) = 'Swell peak period'
   vl(4,26) = 's'
   vmin(26) = 1.
   vmax(26) = 30.
    
   vl(1,27) = 'tmp_swell'
   vl(2,27) = 'sea_surface_swell_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment'
   vl(3,27) = 'Swell mean period'
   vl(4,27) = 's'
   vmin(27) = 1.
   vmax(27) = 25.
     
   vl(1,28) = 'tm1_swell'
   vl(2,28) = 'sea_surface_swell_wave_mean_period_from_variance_spectral_density_first_frequency_moment'
   vl(3,28) = 'Swell m1-period'
   vl(4,28) = 's'
   vmin(28) = 1.
   vmax(28) = 25.
   
   vl(1,29) = 'tm2_swell'
   vl(2,29) = 'sea_surface_swell_wave_mean_period_from_variance_spectral_density_second_frequency_moment'
   vl(3,29) = 'Swell tm2-period'
   vl(4,29) = 's'
   vmin(29) = 1.
   vmax(29) = 25.

   vl(1,30) = 'thq_swell'
   vl(2,30) = 'sea_surface_swell_wave_to_direction'
   vl(3,30) = 'Swell mean wave direction'
   vl(4,30) = 'degree'
   vmin(30) = 0.
   vmax(30) = 360.
    
   vl(1,31) = 'ds_swell'
   vl(2,31) = 'sea_surface_swell_wave_directional_spread'
   vl(3,31) = 'Swell directional spreed'
   vl(4,31) = 'degree'
   vmin(31) = 0.
   vmax(31) = 120.

   vl(1,33) = 'goda'
   vl(2,33) = 'goda_peakness_parameter'
   vl(3,33) = 'goda peakness parameter'
   vl(4,33) = ''
   vmin(33) = 0.
   vmax(33) = 15.

   vl(1,34) = 'kurtosis'
   vl(2,34) = 'kurtosis'
   vl(3,34) = 'kurtosis'
   vl(4,34) = ''
   vmin(34) = 0.
   vmax(34) = 1.

   vl(1,35) = 'BFI'
   vl(2,35) = 'Benjamin_Feir_index'
   vl(3,35) = 'Benjamin Feir index'
   vl(4,35) = ''
   vmin(35) = -10.
   vmax(35) = 10.

   vl(1,36) = 'mHs'
   vl(2,36) = 'normalized_maximum_wave_height'
   vl(3,36) = 'normalized maximum wave height'
   vl(4,36) = 'm'
   vmin(36) = 0.
   vmax(36) = 40.

   vl(1,37) = 'mwp'
   vl(2,37) = 'maximum_wave_period'
   vl(3,37) = 'maximum wave period'
   vl(4,37) = 's'
   vmin(37) = 1.
   vmax(37) = 30.

   vl(1,38) = 'TpI'
   vl(2,38) = 'interpolated_peak_frequency'
   vl(3,38) = 'interpolated peak frequency'
   vl(4,38) = 's'
   vmin(38) = 0.
   vmax(38) = 0.8
 
   vl(1,39) = 'Pdir'
   vl(2,39) = 'peak_direction'
   vl(3,39) = 'peak direction'
   vl(4,39) = 'degree'
   vmin(39) = 0.
   vmax(39) = 360.

   vl(1,40) = 'msqs'
   vl(2,40) = 'mean_square_slope'
   vl(3,40) = 'mean square slope'
   vl(4,40) = ''
   vmin(40) = 0.
   vmax(40) = 1.
    
   flg = 1                   !! prepare table for required parameters only
   do loop=1,nf
      if (.not.cflag_p(loop)) then
         flg([loop]) = 0
      endif
   enddo
      
   ty = NF90_FLOAT           !! type of parameter
   fw = -999.                !! dummy (zmiss)
   did = [3,2,1]
   i = dt
   WRITE(tda,'("0000-00-00 (",2(i2.2,":"),i2.2,")")')i/3600,MOD(i/60,60),MOD(i,60)
   Basedate = '20010101000000'
   call DIFFDATEHR(Basedate, sd ,hradd)
   
   tua="hours since "//Basedate(1:4)//"-"//Basedate(5:6)//"-"//Basedate(7:8)//" "//Basedate(9:10)//":"//Basedate(11:12)//":"//Basedate(13:14)
     
   CALL Pf(NF90_CREATE(name,ior(NF90_CLOBBER,NF90_SHARE),ncid(0)))
   CALL Pf(NF90_DEF_DIM(ncid(0),'time',NF90_UNLIMITED,diid(1)))
   CALL Pf(NF90_DEF_DIM(ncid(0),'rlat',ny,diid(2)))
   CALL Pf(NF90_DEF_DIM(ncid(0),'rlon',nx,diid(3)))
   CALL Pf(NF90_DEF_VAR(ncid(0),'rlat',NF90_FLOAT,[diid(2)],ncid(nf+1)))
   CALL Pf(NF90_PUT_ATT(ncid(0),ncid(nf+1),'axis','Y'))  
   CALL Pf(NF90_DEF_VAR(ncid(0),'rlon',NF90_FLOAT,[diid(3)],ncid(nf+2)))
   CALL Pf(NF90_PUT_ATT(ncid(0),ncid(nf+2),'axis','X'))
!   CALL Pf(NF90_DEF_VAR(ncid(0),'time',NF90_FLOAT,[diid(1)],ncid(nf+3)))
   CALL Pf(NF90_DEF_VAR(ncid(0),'time',NF90_DOUBLE,[diid(1)],ncid(nf+3)))
   CALL Pf(NF90_PUT_ATT(ncid(0),ncid(nf+3),'axis','T'))
   CALL Pf(NF90_PUT_ATT(ncid(0),ncid(nf+3),'standard_name','time'))
   CALL Pf(NF90_PUT_ATT(ncid(0),ncid(nf+3),'calendar',"gregorian"))
   CALL Pf(NF90_PUT_ATT(ncid(0),ncid(nf+3),'delta_t',tda))
   CALL Pf(NF90_PUT_ATT(ncid(0),ncid(nf+3),'units',tua))
   CALL Pf(NF90_PUT_ATT(ncid(0),ncid(nf+3),'dt_sec',i))

   CALL Pf(NF90_PUT_ATT(ncid(0),NF90_GLOBAL,'Conventions','CF-1.6'))
   CALL Pf(NF90_PUT_ATT(ncid(0),NF90_GLOBAL,'institution',metno))

 

   metn=(/ 3, 4, 5, 10, 15, 16, 18, 23, 26, 30, 31, 33, 34, 35, &
&                 36, 37, 38, 39, 40 /)
   
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
tda = real(idelint)/(60*60)
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
real*8    hradd
real*8 addsec

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
addsec=0.0

if (yn.gt.1) then
 do i = 1,yn-1
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
!hradd=real(addsec)/(60.0*60.0)
hradd=addsec/(60.0*60.0)
write(iu06,*)'hradd',hradd
end subroutine DIFFDATEHR


end module wam_netcdf_metno_module
