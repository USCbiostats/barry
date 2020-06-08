docs:
	doxygen Doxyfile
tests:
	cd tests && $(MAKE) 
.PHONY: docs tests
