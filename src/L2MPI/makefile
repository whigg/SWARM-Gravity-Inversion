#	create executable

SRCS = L2precheck L2solverDMAp L2solverNMAp
XDIR    = .
#XDIR    = ../../exe
DDIR	= .
#	create executable

#DIR  	= $(HOME)/my_packages
#LDIR 	= $(DIR)/lib
#IDIR 	= $(DIR)/include
#VPATH	= .:$(LDIR):

#FFLAGS  = -Wall -I$(IDIR) -O3
#FC      = gfortran
#CFLAGS  = -ipo -O2 -xHost -openmp -Wall $(MKL_CFLAGS)
CFLAGS  = -ipo -O2 -xHost -Wall $(MKL_CFLAGS)
CC      = mpicc
#LDFLAGS = -L. -L$(LDIR)
LDLIBS  = -lm $(MKL_LIBS) $(SCALAPACK_LIBS) -lifcore
#LDLIBS  = -lm $(MKL_LIBS) /usr/local/scalapack/2.0.2-intel/lib/libscalapack.a 
#AR      = ar
#ARFLAGS = rcs

EXES  	= $(SRCS:%=$(XDIR)/%.e)
DBGS	= $(SRCS:%=$(DDIR)/%.d)

all : $(EXES) 

d : $(DBGS)

$(XDIR)/%.e : %.o $(LIBS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS)

$(DDIR)/%.d : %.c $(LIBS)
	$(CC) -g -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS)

clean :
	rm  $(EXES) 

cleand :
	rm  $(DBGS)
.PHONY : all d clean cleand
