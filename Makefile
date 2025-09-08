.PHONY: docs tests clean phylo test-portable
docs:
	cp -R mathjax docs/. && doxygen Doxyfile
tests:
	cd tests && $(MAKE) 
test-portable:
	cd tests && $(MAKE) test-portable
coverage:
	cd tests && $(MAKE) coverage
clean:
	cd tests && $(MAKE) clean
phylo: examples/09-phylo-simulation.cpp
	cd examples&& g++ -std=c++11 -Wall -pedantic 09-phylo-simulation.cpp -o phylosim && \
		valgrind --leak-check=full ./phylosim
prof:
	cd examples/&& \
		g++ -std=c++11 -Wall -g -O2 -Wpedantic phylo-suff-stats.cpp -o phylo-suff-stats && \
		valgrind --tool=callgrind ./phylo-suff-stats
.PHONY: prof counts

counts:
	cloc include/** tests/*.cpp examples/

barry.hpp: 
	Rscript --verbose --vanilla barry-hpp.R


