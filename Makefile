network_efficiency: network_efficiency.o
	gcc -std=gnu99 -fopenmp $(shell pkg-config --libs lal lalsimulation gsl) -o $@ $<

%.o: %.c Makefile
	gcc -fopenmp -std=gnu99 -O3 -ffast-math -mfpmath=387 $(shell pkg-config --cflags lal lalsimulation gsl) -o $@ -c $<

clean:
	rm -f network_efficiency
	rm -f network_efficiency.o