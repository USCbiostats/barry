check: main
	./main && $(MAKE) clean
main: 
	g++ -std=c++11 -Wall -pedantic main.cpp -o main
checkv: 
	$(MAKE) clean; $(MAKE) main; valgrind --leak-check=full ./main
clean:
	rm main
.PHONY: clean check checkv
