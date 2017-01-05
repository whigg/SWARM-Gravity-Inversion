#   create archive library

MODS = itgitp matvec stats dftflt

OUT  = libnumrs.a

#   create archive library

DIR   	= ../..
LDIR  	= $(DIR)/lib
IDIR  	= $(DIR)/include
VPATH 	= .:$(LDIR):

FFLAGS  = -Wall -I$(IDIR) -O3
FC      = gfortran
CFLAGS  = -Wall -I$(IDIR) -O3
CC      = gcc
CXXFLAGS= -Wall -I$(IDIR) -O3
CXX     = g++
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
