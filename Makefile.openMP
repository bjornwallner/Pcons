MAIN=lgscore.o molecule.o pcons.o src/nrutil.c src/jacobi.c src/eigsrt.c src/nets.c


#FLAG=-lm -O3 -funroll-loops -Isrc/ 
#LFLAG=-lm -O1 -Isrc/ 
LFLAG=-O3 -Isrc/ -fopenmp -std=c99
#LFLAG=-g -Isrc/ -fopenmp -std=c99
#LFLAG=-g -lm -Isrc/ 

MAIN: $(MAIN)
	$(CC) $(LFLAG) -o bin/pcons.openmp $(MAIN) -lm


pcons.o: src/pcons.c
	$(CC) $(LFLAG) $(CCFLAG) -c src/pcons.c -lm

molecule.o: src/molecule.c src/molecule.h
	$(CC) $(LFLAG) $(CCFLAG) -c src/molecule.c -lm

lgscore.o: src/lgscore.c src/lgscore.h
	$(CC) $(LFLAG) $(CCFLAG) -c src/lgscore.c -lm


.c.o:	
	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c -lm
