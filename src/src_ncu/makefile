OPT=-O3
NOVASOBJ=novas.o novascon.o nutation.o eph_manager.o solsys1.o readeph0.o
NOVASHED=novas.h novascon.h nutation.h eph_manager.h solarsystem.h

main.orb: main.l1b2orb.o $(NOVASOBJ)
	gcc $(OPT) -o main.orb main.l1b2orb.o $(NOVASOBJ) -lm
main.l1b2orb.o: main.l1b2orb.c main.l1b2orb.h
	gcc $(OPT) -c main.l1b2orb.c main.l1b2orb.h -lm
$(NOVASOBJ):$(NOVASHED)

.PHONY: clean
clean:
	rm *.o *.gch
