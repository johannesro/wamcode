#!/bin/bash

# location of WAM source code folders (chief,mod,print,preproc)
SRC=src

# makefile
MKF=mk/makefile_ubuntu.mk

# location for executables
EXE=abs

# vilje:
module unload intelcomp/12.0.5.220
module load intelcomp/13.0.1 mpt/2.06 netcdf/4.3.0

mkdir ${EXE}
cp ${SRC}/*/*.f90 ${EXE} # copy source code
cp ${MKF} ${EXE}/makefile # copy makefile

# compile everything
cd ${EXE}
make clean
make all


