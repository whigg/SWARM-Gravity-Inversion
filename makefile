
#MODULES = novas coord grvts grace fitgps orb2l1c RELEASE_2010-03-31
#MODULES = novas coord grvts numrs grace GRACEL1C RELEASE_2010-03-31
MODULES = novasc coord grvts numrs grace GRACEL1C GraceReadSW 

DIR = src

MDIRS = $(patsubst %,$(DIR)/%,$(MODULES))


all:
	mkdir -p bin lib
	for dir in $(MDIRS); do \
		$(MAKE) --no-print-directory -C $$dir all; \
	done
#	cp -f $(DIR)/RELEASE_2010-03-31/Bin2AsciiLevel1.e exe/

clean:
	for dir in $(MDIRS); do \
		$(MAKE) -C $$dir clean; \
  	done

.PHONY: all clean
