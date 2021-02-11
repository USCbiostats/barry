aphylo:
	g++ -std=c++11 -Wall -pedantic -g aphylo2/phylo-suffstats2.cpp -o aphylo2/pss2 && \
		valgrind --leak-check=full ./aphylo2/pss2


phylo:
	m; g++ -std=c++11 -Wall -pedantic -g phylosim.cpp -o phylosim && \
		valgrind --leak-check=full ./phylosim
ergm:
	R --debugger=valgrind --debugger-args="--leak-check=full --track-origins=yes --show-leak-kinds=all --log-file=07-ergm.log" -e 'Rcpp::sourceCpp("07-ergm.cpp", verbose = TRUE, rebuild = TRUE)' 
clean: 
	rm phylosim

.PHONY: ergm clean
