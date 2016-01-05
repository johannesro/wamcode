program prepTOPO 
  implicit none

  integer                 :: nx,ny,i,j,k,kmax,n
  real                    :: dlon,dlat,lats,latn,lonw,lone
  integer                 :: n1,n2

  real,allocatable        :: bath(:,:)
  integer,allocatable     :: idat(:,:)
  character*1,allocatable :: a(:,:)

  real                    ::     undef,udef
  parameter (undef=+1.e+35)
  parameter (udef=undef*0.9)


  nx=21
  ny=21
  lats=0
  lonw=0
  dlat=0.2
  dlon=0.2
  latn = lats + (ny-1)*dlat
  lone = lonw + (nx-1)*dlon

  allocate(bath(nx,ny))
  allocate(idat(nx,ny))
  allocate(a(nx,ny))

  open(unit=11,file='LABTOPO.DAT')

! Create bathymetry

  do j=1,ny
     do i=1,nx
        bath(i,j)=1000.
!        if (i.eq.1. .or. i.eq.nx)bath(i,j)=+1.e+35
        if (i.eq.1.)bath(i,j)=+1.e+35
!        if (j.eq.1. .or. j.eq.ny)bath(i,j)=+1.e+35
     end do
  end do

  kmax=(nx+11)/12
  do k=1,ny
     do i=1,nx
!        write(*,*)data(i,k)
        if(bath(i,k).lt.udef)then
           idat(i,k)=-nint(abs(bath(i,k)))
           a(i,k)='D'
        else
           idat(i,k)=200
           a(i,k)='E'
        end if
     end do
  end do

  write(*,'(8F10.5)')dlat,dlon,lats,latn,lonw,lone
  write(11,'(8F10.5)')dlat,dlon,lats,latn,lonw,lone
  do n=1,ny
     do k=1,kmax
        n1=12*(k-1)+1
        n2 = min(12*k,nx)
        write(*,'(12(I5,A1))')(idat(i,n),a(i,n),i=n1,n2)
        write(11,'(12(I5,A1))')(idat(i,n),a(i,n),i=n1,n2)
     end do
  end do


  close(11)

end program prepTOPO
