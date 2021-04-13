BFILES=barr*.gcov geese-*.gcov flock-*.gcov \
       counters-*.gcov model-*.gcov rules-*.gcov \
       cell-*.gcov powerset-*.gcov support-*.gcov \
       statscounter*.gcov statsdb*.gcov

check: main
	./main.o && \
	$(MAKE) clean
main: main.cpp
	g++ -std=c++11 -Wall -pedantic -fprofile-arcs -ftest-coverage\
	       	main.cpp -g -o main.o
checkv: main
	$(MAKE) main;\
       	valgrind --leak-check=full --track-origins=yes ./main.o #valgrind --tool=cachegrind ./main
clean:
	rm main.o; rm *.gcov; rm *gcno; rm *gcda
coverage: main
	./main.o && gcov -r main.cpp && \
		rm gcov-files/*.gcov && \
		mv *.hpp.gcov gcov-files/ && $(MAKE) clean
.PHONY: clean check checkv
