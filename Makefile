MAIN=lgscore.o molecule.o pcons.o src/nrutil.c src/jacobi.c src/eigsrt.c src/nets.c
AE=lgscore.o molecule.o pconsAE.o src/nrutil.c src/jacobi.c src/eigsrt.c src/nets.c

#Uncomment to compile Sscore version 
#GOAL=Sscore
#SCORE=$(GOAL)
#SSCORE =-D$(SCORE)

#comment to compile Sscore version 
GOAL=LGscore
SCORE=$(GOAL)
SSCORE =
##########

Arch = $(shell hostname)
#FLAG=-lm -O3 -funroll-loops -Isrc/ 
#LFLAG=-lm -O1 -Isrc/ 
LFLAG=-O3 -Isrc/ -fopenmp -std=c99 -static $(SSCORE)
#LFLAG=-g -Isrc/ -fopenmp -std=c99
#LFLAG=-g -lm -Isrc/ 

AE: $(AE)
	$(CC) $(LFLAG) -o bin/pconsAE.$(SCORE).$(Arch) $(AE) -lm

MAIN: $(MAIN)
	$(CC) $(LFLAG) -o bin/pcons.$(SCORE).$(Arch) $(MAIN) -lm


pcons.o: src/pcons.c
	$(CC) $(LFLAG) $(CCFLAG) -c src/pcons.c -lm

pconsAE.o: src/pconsAE.c
	$(CC) $(LFLAG) $(CCFLAG) -c src/pconsAE.c -lm

molecule.o: src/molecule.c src/molecule.h
	$(CC) $(LFLAG) $(CCFLAG) -c src/molecule.c -lm

lgscore.o: src/lgscore.c src/lgscore.h
	$(CC) $(LFLAG) $(CCFLAG) -c src/lgscore.c -lm


.c.o:	
	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c -lm

clean:
	rm -rf *.o 
