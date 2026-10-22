CPPFLAGS = -std=c++17 -O3 -Wall -pedantic

all: 01-defm-derivatives.a
	./01-defm-derivatives.a 1

# Requires a compiler with OpenMP (e.g., GNU g++ on Linux)
pll: 01-defm-derivatives-pll.a
	./01-defm-derivatives-pll.a 1 && ./01-defm-derivatives-pll.a 4

01-defm-derivatives.a: 01-defm-derivatives.cpp
	$(CXX) $(CPPFLAGS) 01-defm-derivatives.cpp -o 01-defm-derivatives.a

01-defm-derivatives-pll.a: 01-defm-derivatives.cpp
	$(CXX) $(CPPFLAGS) -fopenmp 01-defm-derivatives.cpp -o 01-defm-derivatives-pll.a

clean:
	rm -f *.a
