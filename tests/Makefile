check: main
	./main && $(MAKE) clean
main: 
	g++ -std=c++11 -Wall -g -pedantic main.cpp -o main
checkv: 
	$(MAKE) clean; $(MAKE) main; valgrind --tool=cachegrind ./main
clean:
	rm main
.PHONY: clean check checkv
