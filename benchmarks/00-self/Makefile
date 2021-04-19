CPPFLAGS = -std=c++11 -O3 -Wall -pedantic -g 
# icpx
INTELCOMP = dpcpp
all: 01-omp-pll.a 01-omp-pll.b
	echo --- Intel pll ------- && ./01-omp-pll.b && \
	echo --- Intel serial ---- && ./01-omp.b && \
	echo --- GNU pll --------- && ./01-omp-pll.a && \
	echo --- GNU serial ------ && ./01-omp.a

01-omp-pll.b: 01-omp.cpp
	$(INTELCOMP) $(CPPFLAGS) -qopenmp -DBARRY_USE_OMP 01-omp.cpp -o 01-omp-pll.b -xCORE-AVX2 && \
	$(INTELCOMP) $(CPPFLAGS) 01-omp.cpp -o 01-omp.b
	
01-omp-pll.a: 01-omp.cpp
	g++ $(CPPFLAGS) -fopenmp -DBARRY_USE_OMP 01-omp.cpp -o 01-omp-pll.a -mavx2 && \
	g++ $(CPPFLAGS) -mno-avx2 01-omp.cpp -o 01-omp.a 
clean:
	rm *.o;rm *.a;rm *.b
