BFILES=barr*.gcov geese-*.gcov flock-*.gcov \
       counters-*.gcov model-*.gcov rules-*.gcov \
       cell-*.gcov powerset-*.gcov support-*.gcov \
       statscounter*.gcov statsdb*.gcov

CPPFLAGS= -std=c++11 -Wall -pedantic

check: main.o
	./main.o

main.o: main.cpp
	g++ $(CPPFLAGS) -O3 main.cpp -g -o main.o

checkv: mainv.o
	valgrind --leak-check=full --track-origins=yes ./mainv.o #valgrind --tool=cachegrind ./main

mainv.o: main.cpp
	g++ $(CPPFLAGS) -g -O1 -fprofile-arcs -ftest-coverage main.cpp -o mainv.o

clean:
	rm main*.o; rm *.gcov; rm *gcno; rm *gcda

coverage: mainv.o
	./mainv.o && gcov -r main.cpp && \
		rm catch*.hpp.gcov && \
		mv *.hpp.gcov gcov-files/ && $(MAKE) clean

.PHONY: clean check checkv
