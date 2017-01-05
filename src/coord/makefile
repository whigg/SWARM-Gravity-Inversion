#   create archive library

MODS = tedot orbxyz xform testeph dtdb 

OUT  = libcoord.a

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
