CC = g++
CFLAGS = -Wall -g -std=c++11
LDFLAGS =

MC_Focusing: MC_Focussing.o Hybrid_Taus.o Allocate.o
	$(CC) $(CFLAGS) MC_Focusing.o Hybrid_Taus.o Allocate.o ${LDFLAGS} -o MC_Focusing

Hybrid_Taus.o: Hybrid_Taus.cpp Hybrid_Taus.h
	${CC} ${CFLAGS} -c Hybrid_Taus.cpp

Allocate.o: Allocate.cpp Allocate.h
	${CC} ${CFLAGS} -c Allocate.cpp

MC_Focussing.o: MC_Focusing.cpp MC_Focusing.h Hybrid_Taus.h
	${CC} ${CFLAGS} -c MC_Focusing.cpp

clean: 
	$(RM) count *.o *~
