#	create executable

#SRCS = aotides gnv2gps geteigen editacc editsca fitgps orb2l1c l1c2l2
SRCS = sp32gps
#LIBS = libnovas.a libnumrs.a libcoord.a libgrvts.a
LIBS =  libgrace.a libgrvts.a libnovas.a libcoord.a libnumrs.a
XDIR    = ../../bin
DDIR	= .
#	create executable

#DIR  	= $(HOME)/my_packages
DIR  	= ../../
LDIR 	= $(DIR)/lib
IDIR 	= $(DIR)/include
VPATH	= .:$(LDIR):

FFLAGS  = -Wall -I$(IDIR) -O3
FC      = gfortran
CFLAGS  = -Wall -I$(IDIR) -O3 -fopenmp
CC      = gcc
CXXFLAGS  = -Wall -I$(IDIR) -O3 -fopenmp 
CXX      = g++
LDFLAGS = -L. -L$(LDIR)
LDLIBS  = -lm -lgfortran
AR      = ar
ARFLAGS = rcs

EXES  	= $(SRCS:%=$(XDIR)/%.e)
DBGS	= $(SRCS:%=$(DDIR)/%.d)

all : $(EXES) 

d : $(DBGS)

$(XDIR)/%.e : %.o $(LIBS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)

$(DDIR)/%.d : %.c $(LIBS)
	$(CXX) -g -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)

clean :
	rm  $(EXES) 

cleand :
	rm  $(DBGS)
.PHONY : all d clean cleand
