CC=gcc
OPT=-O2
OBJECTS=ms.o data_sumstat.o rand1.o streec.o
CCFLAGS=-c -g -Wall
LIBS=-lm

msABC: $(OBJECTS)
	$(CC) $(OPT) $(OBJECTS) $(LIBS) -o msABC -g -lm
ms.o: ms.c
	$(CC) $(OPT) $(CCFLAGS) ms.c
data_sumstat.o: data_sumstat.c
	$(CC) $(OPT) $(CCFLAGS) data_sumstat.c
rand1.o: rand1.c
	$(CC) $(OPT) $(CCFLAGS) rand1.c
streec.o: streec.c
	$(CC) $(OPT) $(CCFLAGS) streec.c
clean:
	rm -f *.o
	rm -f msABC
