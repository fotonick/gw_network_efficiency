network_efficiency: network_efficiency.o
	gcc -fopenmp $(shell pkg-config --libs lal lalsimulation gsl) -o $@ $<

%.o: %.c Makefile
	gcc -fopenmp -std=c99 -O3 -ffast-math -mfpmath=387 $(shell pkg-config --cflags lal lalsimulation gsl) -o $@ -c $<

