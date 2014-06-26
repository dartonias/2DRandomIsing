CC=g++ -std=c++0x
CCOPT=-O3 -DNDEBUG
CCFLAGS=-g -DDEBUG

all: Ising.out
debug: Ising_debug.out

Ising_debug.out: *.cpp *.hpp
	$(CC) $(CCFLAGS) main.cpp -o $@
Ising.out: *.cpp *.hpp
	$(CC) $(CCOPT) main.cpp -o $@
clean:
	rm -rf *.o Ising*.out obs_* savobs_*
