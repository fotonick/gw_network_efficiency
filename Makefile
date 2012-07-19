.SECONDARY: network_efficiency.c

network_efficiency: network_efficiency.c Makefile
	gcc -g -fopenmp -O3 -ffast-math -mfpmath=387 $(shell pkg-config --cflags lal lalsimulation gsl) -o $@.o -c $<
	gcc -fopenmp $(shell pkg-config --libs lal lalsimulation gsl) -o $@ $@.o
