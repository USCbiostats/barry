CPP=g++
VPATH = ../include/barry
BFILES=barr*.gcov geese-*.gcov flock-*.gcov \
       counters-*.gcov model-*.gcov rules-*.gcov \
       cell-*.gcov powerset-*.gcov support-*.gcov \
       statscounter*.gcov statsdb*.gcov

CPPFLAGS= -std=c++17 -fopenmp -Wall -pedantic -march=native -Wbool-operation # only works with g++10

CPPFLAGSV= -std=c++17 -Wall -pedantic

# Cross-platform compatible flags (no OpenMP, no march=native, no GCC-specific warnings)
CPPFLAGS_PORTABLE= -std=c++17 -Wall -pedantic -O2

check: main.o 
	./main.o

# Cross-platform test target
test-portable: main-portable.o
	./main-portable.o

main-portable.o: main.cpp
	$(CPP) $(CPPFLAGS_PORTABLE) main.cpp -g -o main-portable.o

.PHONY: %.cpp
main.o: %.cpp
	$(CPP) $(CPPFLAGS) -O2 main.cpp -g -o main.o 

mainsan.o: %.cpp
	$(CPP) $(CPPFLAGS) -fsanitize=address -lasan -O1 main.cpp -g -o mainsan.o

checkv: mainv.o
	valgrind --leak-check=full --track-origins=yes ./mainv.o #valgrind --tool=cachegrind ./main

main: main.cpp
	g++ $(CPPFLAGSV) -g -O1 -fprofile-arcs -ftest-coverage main.cpp -DBARRY_VALGRIND -o main

clean:
	rm -f main*.o; rm -f *.gcov; rm -f *gcno; rm -f *gcda

coverage: mainv.o
	./mainv.o && gcov -r mainv.o-main.cpp && \
		rm catch*.hpp.gcov && \
		mv *.hpp.gcov gcov-files/ && $(MAKE) clean

.PHONY: callgrind
cachegrind: main.calgrind.out

main.cachegrind.out: main.o
	valgrind --tool=callgrind --callgrind-out-file=main.callgrind.out ./main && \
		kcachegrind main.callgrind.out &

.PHONY: clean check checkv

03.o: 03b-support-ergm-b.cpp
	$(CPP) $(CPPFLAGS) -O2 -g 03b-support-ergm-b.cpp -o 03.o

13.o: 13-likelihood-benchmark.cpp
	$(CPP) $(CPPFLAGS) -O2 -g 13-likelihood-benchmark.cpp -o 13.o

10.o: 10-geese-predict.cpp
	$(CPP) $(CPPFLAGS) -O2 -g 10-geese-predict.cpp -o 10.o

# Target depends on cpp file
%.o: %.cpp
	$(CPP) $(CPPFLAGS) -O2 -g $< -o $@
