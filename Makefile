.PHONY: docs tests clean phylo test-portable help

# Default target
help:
	@echo "Barry Makefile - Available targets:"
	@echo ""
	@echo "  help           - Show this help message (default)"
	@echo "  tests          - Build and run tests with OpenMP"
	@echo "  test-portable  - Build and run tests without OpenMP (cross-platform)"
	@echo "  coverage       - Run tests with coverage analysis"
	@echo "  docs           - Generate documentation with Doxygen"
	@echo "  clean          - Clean build artifacts"
	@echo "  phylo          - Build and test phylogenetic simulation example"
	@echo "  prof           - Profile phylogenetic example with callgrind"
	@echo "  counts         - Count lines of code"
	@echo "  barry.hpp      - Generate single header file"
	@echo ""
	@echo "For more information, see the README.md file."

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
.PHONY: help docs tests clean phylo test-portable prof counts coverage

counts:
	cloc include/** tests/*.cpp examples/

barry.hpp: 
	Rscript --verbose --vanilla barry-hpp.R


