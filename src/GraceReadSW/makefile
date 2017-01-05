SRCS =  Bin2AsciiLevel1
XDIR    = ../../bin

EXES  	= $(SRCS:%=$(XDIR)/%.e)

all: $(EXES)
#all: ../../exe/Bin2AsciiLevel1.e

$(EXES):
	tar xvf GraceReadSW_L1_2010-03-31.tar
	cp GRACEsyspath.h RELEASE_2010-03-31
	$(MAKE) -C RELEASE_2010-03-31 all
	mv RELEASE_2010-03-31/Bin2AsciiLevel1.e $(XDIR)
	rm -rf RELEASE_2010-03-31

clean :
	rm  $(EXES) 
#	rm  ../../exe/Bin2AsciiLevel1.e

.PHONY : all clean 
