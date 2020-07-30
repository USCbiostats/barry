ergm:
	R --debugger=valgrind --debugger-args="--leak-check=full --track-origins=yes --show-leak-kinds=all --log-file=07-ergm.log" -e 'Rcpp::sourceCpp("07-ergm.cpp", verbose = TRUE, rebuild = TRUE)' 
phylo:
	rm phylosim; g++ -std=c++11 -Wall -pedantic 09-phylo-simulation.cpp -o phylosim && \
		valgrind --leak-check=full ./phylosim
.PHONY: ergm
