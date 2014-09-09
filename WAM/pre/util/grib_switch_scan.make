# Makefile for bilingrb.c
#
# This makefile assumes that the library libgribw.a has been compiled and is 
# found in ../Gribw together with the header files listed in $(INC).
# Also, a ./TMP directory is assumed for storing object files to keep your
# workspace neat & tidy.

VPATH = .:TMP:../Gribw

.SUFFIXES:
.SUFFIXES: .o .c .h

CC = cc
LD = $(CC)
LDPATH = ../../Gribw
LDFLAGS = -I$(LDPATH) -L$(LDPATH)
LIBS = -lgribw -lm
TARGET = grib_switch_scan

# Suffix rule
.c.o:
	cd ./TMP ; $(CC) -c -Wall -O2 -o $*.o ../$*.c

INC = gribw.h bds.h gds.h
SRC = 
OBJ = $(SRC:.c=.o)

$(TARGET): $(TARGET).c $(INC) $(OBJ)
	cd ./TMP; $(LD) ../$(TARGET).c -o ../$(TARGET) $(OBJ) $(LDFLAGS) $(LIBS)
