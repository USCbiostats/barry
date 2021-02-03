.PHONY: docs tests clean phylo
docs:
	doxygen Doxyfile
tests:
	cd tests && $(MAKE) 
clean:
	cd tests && $(MAKE) clean
phylo: examples/09-phylo-simulation.cpp
	cd examples&& g++ -std=c++11 -Wall -pedantic 09-phylo-simulation.cpp -o phylosim && \
		valgrind --leak-check=full ./phylosim

