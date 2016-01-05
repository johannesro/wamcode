import numpy as np
import scipy as sp
from netCDF4 import Dataset
import pylab as pl

ntime = 72
nx=21
ny=21
lats = 0.
lonw = 0.
dlat = 0.2
dlon = 0.2

u = sp.ones((ntime,nx,ny))*10.
v = sp.ones((ntime,nx,ny))*0.

lat = lats + sp.arange(ny)*dlat
lon = lonw + sp.arange(nx)*dlon 

#'icec'

icec = sp.ones((ntime,nx,ny))*0.

ofile = Dataset('WIND_INPUT.DAT','w')
ofile.createDimension('lon',size=nx)
ofile.createDimension('lat',size=ny)
ofile.createDimension('time',size=ntime)

nc_u = ofile.createVariable('x_wind',sp.float32,dimensions=('time','lat','lon'))
nc_u[:] = u
nc_v = ofile.createVariable('y_wind',sp.float32,dimensions=('time','lat','lon'))
nc_v[:] = v
nc_ice = ofile.createVariable('icec',sp.float32,dimensions=('time','lat','lon'))
nc_ice[:] = icec
nc_lat = ofile.createVariable('lat', sp.float32, dimensions=('lat'))
nc_lat[:] = lat
nc_lon = ofile.createVariable('lon', sp.float32, dimensions=('lon'))
nc_lon[:] = lon

nc_time = ofile.createVariable('time',sp.float64,dimensions=('time'))
nc_time[:] = sp.arange(0,ntime,1)
nc_time.setncattr('units','hours since 2015-01-01 00:00:00')


ofile.close()
