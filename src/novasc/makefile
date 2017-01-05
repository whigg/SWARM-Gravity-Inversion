#   create archive library

MODS = novas novascon nutation eph_manager readeph0 nutation \
       solsys1 solsys2 solsys3 jplint
OUT  = libnovas.a

#   create archive library

DIR   	= ../..
LDIR  	= $(DIR)/lib
IDIR  	= $(DIR)/include
VPATH 	= .:$(LDIR):

FFLAGS  = -Wall -I$(IDIR) -O3
FC      = gfortran
CFLAGS  = -Wall -I$(IDIR) -O3
CC      = gcc
LDFLAGS = -L. -L$(LDIR)
LDLIBS  = -lm -lgfortran
AR      = ar
ARFLAGS = rcs

LIB  = $(LDIR)/$(OUT)
OBJS = $(MODS:%=%.o)

all : $(LIB)

$(LIB) : $(LIB)($(OBJS))

clean:
	rm $(LIB)
.PHONY : all clean
