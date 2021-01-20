docs:
	doxygen Doxyfile
tests:
	cd tests && $(MAKE) 
clean:
	cd tests && $(MAKE) clean
.PHONY: docs tests clean

