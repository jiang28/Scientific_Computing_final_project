OPT= -O2

all: main

main: testdrive.o qrhouse.o house.o applyH.o applyQ.o 
	gcc $(OPT) -o testdrive testdrive.o -lm qrhouse.o -lm -lgslcblas house.o -lm -lgslcblas applyH.o -lm -lgslcblas applyQ.o -lm -lgslcblas

main.o: testdrive.c qrhouse.c house.c applyH.c applyQ.c qrhouse.h house.h applyH.h applyQ.h cblas.h 
	gcc $(OPT) -c testdrive.c qrhouse.c house.c applyH.c applyQ.c cblas_LINUX.a
	
clean:
	rm -f *.o testdrive
