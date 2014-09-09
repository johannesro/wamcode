#!/bin/sh
#25 Feb 2013
# Prepare the boundaries from the spectra of ECMWF
##  Usage: 
# ./make_boundary.sh  today_spectra 20130503180000 
# ./make_boundary.sh  today_spectra 20130301060000
#                                   YYYYMMDDHHMISE 
# where
# today_spectra =  /opdata/ec/bwp/BWP0226000* 
#
#
# Output file is: BXXYYYYMMDDHHMISE   (boundary to be used by MyWave WAM)

spec_anaplusfc=$1
lastdate=$2 #format YYYYMMDDHHMISE

# Convert from GRIB to sequential binary 
util/grib2spec -i ${spec_anaplusfc} -o Boundary/output_spectra.seq -t 1

cd ./Boundary

## If need it read the seq spectra file by using the rspec program
## ../util/rspec output_spectra.seq >outofreadingspectra
##  Interpolate to desired positions
../util/invdist output_spectra.seq invdist_out.seq invdist_ecwam10.inp
## Rotate to desired model projection - WAM10/HIRLAM10
../util/rotspec invdist_out.seq rotspec_out.seq 0.0 0.0 -40.0 68.0
## Convert to WAM  boundary file format 
../util/boundwam  rotspec_out.seq  wambound_ecwam10.inp Rotlonlatzbound.txt  ../'BXX'${lastdate}

#rm rotspec_out.seq
#rm invdist_out.seq 
#rm output_spectra.seq
#rm outofreadingspectra


#Compiling on vilje
### update modules: module load intelcomp/13.0.1 mpt/2.06 netcdf/4.3.0

#   The way to compile the programs used above are:
##compiling grib2spec on vilje
#!module load cmkl intelcomp mpt
#!ifort  -O3 -real_size 64 -traceback -ip -align all\
#!             -xAVX -ftz -fno-alias -no-prec-div -no-prec-sqrt  -o  grib2spec grib2spec.f specio.f sphere.f -L/prod/forecast/lib/emos_00381 -lemosR64


##compiling invdist  
#!module load cmkl intelcomp mpt #OJO OJO
#!ifort  -O3 -real_size 64 -traceback -ip -align all  \
#!             -xAVX -ftz -fno-alias -no-prec-div -no-prec-sqrt -o invdist invdist.f sphere.f inout.f specio.f

##compiling rotspec  
#! module load cmkl intelcomp mpt #OJO OJO
#! ifort  -O3 -real_size 64 -traceback -ip -align all  \
#!              -xAVX -ftz -fno-alias -no-prec-div -no-prec-sqrt -o rotspec rotspec.f sphere.f rotspher.f specio.f

##compiling boundwam
#!module load cmkl intelcomp mpt
#!ifort  -O3 -real_size 64 -traceback -ip -align all\
#!             -xAVX -ftz -fno-alias -no-prec-div -no-prec-sqrt -o boundwam  boundwam.f rotspher.f specio.f julian.f

##compiling rspec (reads seq spectra file)
#!module load cmkl intelcomp mpt
#!ifort  -O3 -real_size 64 -traceback -ip -align all\
#!             -xAVX -ftz -fno-alias -no-prec-div -no-prec-sqrt  -o rspec rspec.f  specio.f 




